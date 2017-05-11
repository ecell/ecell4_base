from __future__ import print_function

import os
import logging
import tempfile
import pickle
import inspect
import textwrap
import re
import types
import itertools
import binascii
import multiprocessing
import copy

import ecell4.extra.sge as sge


def run_serial(target, jobs, n=1, **kwargs):
    return [[target(copy.copy(job), i + 1, j + 1) for j in range(n)] for i, job in enumerate(jobs)]

def run_multiprocessing(target, jobs, n=1, **kwargs):
    def target_wrapper(f, end_send):
        def wf(*_args, **_kwargs):
            end_send.send(f(*_args, **_kwargs))
        return wf

    processes = []
    end_recvs = []
    for i, job in enumerate(jobs):
        for j in range(n):
            end_recv, end_send = multiprocessing.Pipe(False)
            end_recvs.append(end_recv)
            p = multiprocessing.Process(
                target=target_wrapper(target, end_send), args=(job, i + 1, j + 1))
            p.start()
            processes.append(p)

    for p in processes:
        p.join()

    retval = [end_recv.recv() for end_recv in end_recvs]
    return [retval[i: i + n] for i in range(0, len(retval), n)]

def run_sge(target, jobs, n=1, path='.', delete=True, wait=True, environ=None, modules=(), **kwargs):
    logging.basicConfig(level=logging.DEBUG)

    if isinstance(target, types.LambdaType) and target.__name__ == "<lambda>":
        raise RuntimeError("A lambda function is not accepted")

    # src = textwrap.dedent(inspect.getsource(singlerun)).replace(r'"', r'\"')
    src = textwrap.dedent(inspect.getsource(target)).replace(r'"', r'\"')
    if re.match('[\s\t]+', src.split('\n')[0]) is not None:
        raise RuntimeError(
            "Wrong indentation was found in the source translated")

    if not os.path.isdir(path):
        os.makedirs(path)  #XXX: MYOB

    if environ is None:
        environ = {}
        keys = ("LD_LIBRARY_PATH", "PYTHONPATH")
        for key in keys:
            if key in os.environ.keys():
                environ[key] = os.environ[key]
        if "PYTHONPATH" in environ.keys() and environ["PYTHONPATH"].strip() != "":
            environ["PYTHONPATH"] = "{}:{}".format(os.getcwd(), environ["PYTHONPATH"])
        else:
            environ["PYTHONPATH"] = os.getcwd()

    cmds = []
    pickleins = []
    pickleouts = []
    for i, job in enumerate(jobs):
        (fd, picklein) = tempfile.mkstemp(suffix='.pickle', prefix='sge-', dir=path)
        with os.fdopen(fd, 'wb') as fout:
            pickle.dump(job, fout)
        pickleins.append(picklein)

        pickleouts.append([])
        for j in range(n):
            fd, pickleout = tempfile.mkstemp(suffix='.pickle', prefix='sge-', dir=path)
            os.close(fd)
            pickleouts[-1].append(pickleout)
        # pickleouts.append(
        #     [tempfile.mkstemp(suffix='.pickle', prefix='sge-', dir=path)[1]
        #      for j in range(n)])

        cmd = '#!/bin/bash\n'
        for key, value in environ.items():
            cmd += 'export {:s}={:s}\n'.format(key, value)
        cmd += 'python3 -c "\n'
        cmd += 'import sys\n'
        cmd += 'import os\n'
        cmd += 'import pickle\n'
        cmd += 'with open(sys.argv[1], \'rb\') as fin:\n'
        cmd += '    job = pickle.load(fin)\n'
        cmd += 'pass\n'
        for m in modules:
            cmd += "from {} import *\n".format(m)
        cmd += src
        cmd += '\ntid = int(os.environ[\'SGE_TASK_ID\'])'
        cmd += '\nretval = {:s}(job, {:d}, tid)'.format(target.__name__, i + 1)
        cmd += '\nfilenames = {:s}'.format(str(pickleouts[-1]))
        cmd += '\npickle.dump(retval, open(filenames[tid - 1], \'wb\'))'
        cmd += '" {:s}\n'.format(picklein)
        cmds.append(cmd)

    if isinstance(wait, bool):
        sync = 0 if not wait else 10
    elif isinstance(wait, int):
        sync = wait
    else:
        raise ValueError("'wait' must be either 'int' or 'bool'.")

    jobids = sge.run(cmds, n=n, path=path, delete=delete, sync=sync, **kwargs)

    if not (sync > 0):
        return None

    for jobid, name in jobids:
        outputs = sge.collect(jobid, name, n=n, path=path, delete=delete)
        for output in outputs:
            print(output, end='')

    retval = [[pickle.load(open(pickleout, 'rb')) for pickleout in tasks]
              for tasks in pickleouts]

    if delete:
        for picklename in itertools.chain(pickleins, *pickleouts):
            os.remove(picklename)

    return retval

def genseeds(n):
    """
    Return a random number generator seed for ensemble_simulations.
    A seed for a single run is given by ``getseed(rngseed, i)``.

    Parameters
    ----------
    n : int
        A size of the seed.

    Returns
    -------
    rndseed : bytes
        A random number seed for multiple runs.

    """
    return binascii.hexlify(os.urandom(4 * n))

def getseed(myseed, i):
    """
    Return a single seed from a long seed given by `genseeds`.

    Parameters
    ----------
    myseed : bytes
        A long seed given by `genseeds(n)`.
    i : int
        An index less than n.

    Returns
    -------
    rndseed : int
        A seed (less than (2 ** 31))

    """
    rndseed = int(myseed[(i - 1) * 8: i * 8], 16)
    rndseed = rndseed % (2 ** 31)  #XXX: trancate the first bit
    return rndseed

#XXX:
#XXX:
#XXX:

def singlerun(job, job_id, task_id):
    import ecell4.util
    import ecell4.extra.ensemble
    rndseed = ecell4.extra.ensemble.getseed(job.pop('myseed'), task_id)
    job.update({'return_type': 'array', 'rndseed': rndseed})
    data = ecell4.util.run_simulation(**job)
    return data

import ecell4.util.decorator
import ecell4.util.simulation
import ecell4.util.viz
import ecell4.ode

## observers=(), progressbar=0
def ensemble_simulations(
    t, y0=None, volume=1.0, model=None, solver='ode',
    is_netfree=False, species_list=None, without_reset=False,
    return_type='matplotlib', opt_args=(), opt_kwargs=None,
    structures=None, rndseed=None,
    n=1, nproc=1, method=None, errorbar=True,
    **kwargs):
    """
    Run simulations multiple times and return its ensemble.
    Arguments are almost same with ``ecell4.util.run_simulation``.
    `observers` and `progressbar` is not available here.

    Parameters
    ----------
    n : int, optional
        A number of runs. Default is 1.
    nproc : int, optional
        A number of processors. Ignored when method='serial'.
        Default is 1.
    method : str, optional
        The way for running multiple jobs.
        Choose one from 'serial', 'sge' and 'multiprocessing'.
        Default is None, which works as 'serial'.
    **kwargs : dict, optional
        Optional keyword arugments are passed through to `run_serial`,
        `run_sge`, or `run_multiprocessing`.
        See each function for more details.

    Returns
    -------
    value : list, DummyObserver, or None
        Return a value suggested by ``return_type``.
        When ``return_type`` is 'array', return a time course data.
        When ``return_type`` is 'observer', return a DummyObserver.
        DummyObserver is a wrapper, which has the almost same interface
        with NumberObservers.
        Return nothing if else.

    See Also
    --------
    ecell4.util.run_simulation
    ecell4.extra.run_serial
    ecell4.extra.run_sge
    ecell4.extra.run_multiprocessing

    """
    y0 = y0 or {}
    opt_kwargs = opt_kwargs or {}
    structures = structures or {}

    for key, value in kwargs.items():
        if key == 'r':
            return_type = value
        elif key == 'v':
            volume = value
        elif key == 's':
            solver = value
        elif key == 'm':
            model = value

    if model is None:
        model = ecell4.util.decorator.get_model(is_netfree, without_reset)

    if isinstance(model, ecell4.ode.ODENetworkModel):
        raise ValueError('A model with ratelaws is not supported yet.')

    if species_list is None:
        species_list = ecell4.util.simulation.list_species(model, y0.keys())

    if rndseed is None:
        myseed = genseeds(n)
    elif (not isinstance(rndseed, bytes) or len(rndseed) != n * 4 * 2):
        raise ValueError(
            "A wrong seed for the random number generation was given. Use 'genseeds'.")

    jobs = [{'t': t, 'y0': y0, 'volume': volume, 'model': model, 'solver': solver, 'species_list': species_list, 'structures': structures, 'myseed': myseed}]

    if method is None or method.lower() == "serial":
        retval = run_serial(singlerun, jobs, n=n, **kwargs)
    elif method.lower() == "sge":
        retval = run_sge(singlerun, jobs, n=n, **kwargs)
    elif method.lower() == "multiprocessing":
        retval = run_multiprocessing(singlerun, jobs, n=n, **kwargs)
    else:
        raise ValueError(
            'Argument "method" must be one of "serial", "multiprocessing" and "sge".')

    if return_type is None:
        return

    assert len(retval) == len(jobs) == 1

    if return_type in ("array", 'a'):
        return retval[0]

    import numpy

    class DummyObserver:

        def __init__(self, inputs, species_list, errorbar=True):
            if len(inputs) == 0:
                raise ValueError("No input was given.")

            t = numpy.array(inputs[0], numpy.float64).T[0]
            mean = sum([numpy.array(data, numpy.float64).T[1: ] for data in inputs])
            mean /= len(inputs)

            self.__data = numpy.vstack([t, mean]).T

            if errorbar:
                std = sum([(numpy.array(data, numpy.float64).T[1: ] - mean) ** 2
                           for data in inputs])
                std /= len(inputs)
                std = numpy.sqrt(std)
                self.__error = numpy.vstack([t, std]).T
            else:
                self.__error = None

            self.__species_list = [ecell4.Species(serial) for serial in species_list]

        def targets(self):
            return self.__species_list

        def data(self):
            return self.__data

        def t(self):
            return self.__data.T[0]

        def error(self):
            return self.__error

    if return_type in ("matplotlib", 'm'):
        if isinstance(opt_args, (list, tuple)):
            ecell4.util.viz.plot_number_observer_with_matplotlib(
                DummyObserver(retval[0], species_list, errorbar), *opt_args, **opt_kwargs)
        elif isinstance(opt_args, dict):
            # opt_kwargs is ignored
            ecell4.util.viz.plot_number_observer_with_matplotlib(
                DummyObserver(retval[0], species_list, errorbar), **opt_args)
        else:
            raise ValueError('opt_args [{}] must be list or dict.'.format(
                repr(opt_args)))
    elif return_type in ("nyaplot", 'n'):
        if isinstance(opt_args, (list, tuple)):
            ecell4.util.viz.plot_number_observer_with_nya(
                DummyObserver(retval[0], species_list, errorbar), *opt_args, **opt_kwargs)
        elif isinstance(opt_args, dict):
            # opt_kwargs is ignored
            ecell4.util.viz.plot_number_observer_with_nya(
                DummyObserver(retval[0], species_list, errorbar), **opt_args)
        else:
            raise ValueError('opt_args [{}] must be list or dict.'.format(
                repr(opt_args)))
    elif return_type in ("observer", 'o'):
        return DummyObserver(retval[0], species_list, errorbar)
    elif return_type in ("dataframe", 'd'):
        import pandas
        return [
            pandas.concat([
                pandas.DataFrame(dict(Time=numpy.array(data).T[0],
                                      Value=numpy.array(data).T[i + 1],
                                      Species=serial))
                for i, serial in enumerate(species_list)])
            for data in retval[0]]


if __name__ == "__main__":
    # def myrun(job, job_id=0, task_id=0):
    #     import ecell4
    #     print("Hi, I'm in local!")
    #     print("My job id is {:d}, and my task id is {:d}.".format(job_id, task_id))
    #     print("My job is {:s}.".format(str(job)))
    #     return job['x'] + job['y']

    # jobs = [{'x': i, 'y': i ** 2} for i in range(1, 4)]
    # print(run_serial(myrun, jobs, n=2))
    # print(run_multiprocessing(myrun, jobs, n=2))
    # # print(run_sge(myrun, jobs, n=2, delete=False))
    # print(run_sge(myrun, jobs, n=2))

    from ecell4 import *
    from ecell4.extra import ensemble

    with reaction_rules():
        A + B == C | (0.01, 0.3)

    ensemble.ensemble_simulations(
        10.0, {'C': 60}, solver='gillespie', return_type='matplotlib',
        n=30, method='multiprocessing')
