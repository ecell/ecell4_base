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


def run_serial(target, jobs, n=1):
    return [[target(copy.copy(job), i + 1, j + 1) for j in range(n)] for i, job in enumerate(jobs)]

def run_multiprocessing(target, jobs, n=1):
    def target_wrapper(f, end_send):
        def wf(*args, **kwargs):
            end_send.send(f(*args, **kwargs))
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

def run_sge(target, jobs, n=1, path='.', delete=True, environ={}):
    logging.basicConfig(level=logging.DEBUG)

    if isinstance(target, types.LambdaType) and target.__name__ == "<lambda>":
        raise RuntimeError("A lambda function is not accepted")

    src = textwrap.dedent(inspect.getsource(singlerun)).replace(r'"', r'\"')
    if re.match('[\s\t]+', src.split('\n')[0]) is not None:
        raise RuntimeError(
            "Wrong indentation was found in the source translated")

    cmds = []
    pickleins = []
    pickleouts = []
    for i, job in enumerate(jobs):
        (fd, picklein) = tempfile.mkstemp(suffix='.pickle', prefix='sge-', dir=path)
        with os.fdopen(fd, 'wb') as fout:
            pickle.dump(job, fout)
        pickleins.append(picklein)

        pickleouts.append(
            [tempfile.mkstemp(suffix='.pickle', prefix='sge-', dir=path)[1]
             for j in range(n)])

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
        cmd += src
        cmd += '\ntid = int(os.environ[\'SGE_TASK_ID\'])'
        cmd += '\nretval = {:s}(job, {:d}, tid)'.format(target.__name__, i + 1)
        cmd += '\nfilenames = {:s}'.format(str(pickleouts[-1]))
        cmd += '\npickle.dump(retval, open(filenames[tid - 1], \'wb\'))'
        cmd += '" {:s}\n'.format(picklein)
        cmds.append(cmd)

    jobids = sge.run(cmds, n=n, path=path, delete=delete)

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

#XXX:
#XXX:
#XXX:

def singlerun(job, job_id, task_id):
    import ecell4.util
    myseed = job.pop('myseed')
    rndseed = int(myseed[(task_id - 1) * 8: task_id * 8], 16)
    rndseed = rndseed % (2 ** 31)  #XXX: trancate the first bit
    myrng = ecell4.GSLRandomNumberGenerator()
    myrng.seed(rndseed)
    job.update({'return_type': 'array', 'rng': myrng})
    data = ecell4.util.run_simulation(**job)
    return data

import ecell4.util.decorator
import ecell4.util.simulation
import ecell4.util.viz
import ecell4.ode


def ensemble_simulations(
    t, y0={}, volume=1.0, model=None, solver='ode', species_list=None, structures={},
    is_netfree=False, without_reset=False,
    return_type='matplotlib', opt_args=(), opt_kwargs={},
    errorbar=True, n=1, nproc=1, method=None, environ={}):
    """
    observers=(), progressbar=0, rndseed=None,
    """
    # if not isinstance(solver, str):
    #     raise ValueError('Argument "solver" must be a string.')

    if model is None:
        model = ecell4.util.decorator.get_model(is_netfree, without_reset)

    if isinstance(model, ecell4.ode.ODENetworkModel):
        raise ValueError('A model with ratelaws is not supported yet.')

    if species_list is None:
        species_list = ecell4.util.simulation.list_species(model, y0.keys())

    myseed = binascii.hexlify(os.urandom(4 * n))

    jobs = [{'t': t, 'y0': y0, 'volume': volume, 'model': model, 'solver': solver, 'species_list': species_list, 'structures': structures, 'myseed': myseed}]

    if method is None or method.lower() == "serial":
        retval = run_serial(singlerun, jobs, n=n)
    elif method.lower() == "sge":
        retval = run_sge(singlerun, jobs, n=n, environ=environ)
    elif method.lower() == "multiprocessing":
        retval = run_multiprocessing(singlerun, jobs, n=n)
    else:
        raise ValueError(
            'Argument "method" must be one of "serial", "multiprocessing" and "sge".')

    if return_type == "array":
        return retval

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

    if return_type == "matplotlib":
        if isinstance(opt_args, (list, tuple)):
            ecell4.util.viz.plot_number_observer_with_matplotlib(
                *itertools.chain([DummyObserver(inputs, species_list, errorbar) for inputs in retval],
                                 opt_args),
                **opt_kwargs)
        elif isinstance(opt_args, dict):
            # opt_kwargs is ignored
            ecell4.util.viz.plot_number_observer_with_matplotlib(
                *[DummyObserver(inputs, species_list, errorbar) for inputs in retval],
                **opt_args)
        else:
            raise ValueError('opt_args [{}] must be list or dict.'.format(
                repr(opt_args)))
    elif return_type == "nyaplot":
        if isinstance(opt_args, (list, tuple)):
            ecell4.util.viz.plot_number_observer_with_nya(
                *itertools.chain([DummyObserver(inputs, species_list, errorbar) for inputs in retval],
                                 opt_args),
                **opt_kwargs)
        elif isinstance(opt_args, dict):
            # opt_kwargs is ignored
            ecell4.util.viz.plot_number_observer_with_nya(
                *[DummyObserver(inputs, species_list, errorbar) for inputs in retval],
                **opt_args)
        else:
            raise ValueError('opt_args [{}] must be list or dict.'.format(
                repr(opt_args)))
    elif return_type == "observer":
        return [DummyObserver(inputs, species_list, errorbar) for inputs in retval]
    elif return_type == "dataframe":
        import pandas
        return [[
            pandas.concat([
                pandas.DataFrame(dict(Time=numpy.array(data).T[0],
                                      Value=numpy.array(data).T[i + 1],
                                      Species=serial))
                for i, serial in enumerate(species_list)])
            for data in inputs] for inputs in retval]
    else:
        raise ValueError(
            'Invald Argument "return_type" was given [{}].'.format(str(return_type)))


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

    # # environ = {'LD_LIBRARY_PATH': '/home/kaizu/lily_kaizu/src/ecell4/local/lib', 'PYTHONPATH': '/home/kaizu/lily_kaizu/src/ecell4/local/lib/python3.4/site-packages'}
    # environ = {}
    # # print(run_sge(myrun, jobs, n=2, delete=False, environ=environ))
    # print(run_sge(myrun, jobs, n=2, environ=environ))

    from ecell4 import *
    from ecell4.extra import ensemble

    with reaction_rules():
        A + B == C | (0.01, 0.3)

    environ = {'LD_LIBRARY_PATH': '/home/kaizu/lily_kaizu/src/ecell4/local/lib',
               'PYTHONPATH': '/home/kaizu/lily_kaizu/src/ecell4/local/lib/python3.4/site-packages'}
    retval = ensemble.ensemble_simulations(
        10.0, {'C': 60}, solver='gillespie', return_type='matplotlib',
        n=5, method='multiprocessing', environ=environ)

    # import numpy

    # def concatenate(results):
    #     return sum([numpy.array(data) for data in results]) / len(results)
    # retval = [concatenate(results) for results in retval]
    # print(retval)

    # numpy.savetxt('ens.dat', retval[0])
