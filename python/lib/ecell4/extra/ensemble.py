import os
import logging
import tempfile
import pickle
import inspect
import textwrap
import re
import types
import itertools

import multiprocessing

import ecell4.extra.sge as sge

def run_serial(target, jobs, n=1):
    return [[target(job) for i in range(n)] for job in jobs]

def run_multiprocessing(target, jobs, n=1):
    def target_wrapper(f, end_send):
        def wf(j):
            end_send.send(f(j))
        return wf

    processes = []
    end_recvs = []
    for job in jobs:
        for i in range(n):
            end_recv, end_send = multiprocessing.Pipe(False)
            end_recvs.append(end_recv)
            p = multiprocessing.Process(
                target=target_wrapper(target, end_send), args=(job, ))
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
    for job in jobs:
        (fd, picklein) = tempfile.mkstemp(suffix='.pickle', prefix='sge-', dir=path)
        with os.fdopen(fd, 'wb') as fout:
            pickle.dump(job, fout)
        pickleins.append(picklein)

        pickleouts.append(
            [tempfile.mkstemp(suffix='.pickle', prefix='sge-', dir=path)[1]
             for i in range(n)])

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
        cmd += '\nretval = {:s}(job)'.format(target.__name__)
        cmd += '\ntid = int(os.environ[\'SGE_TASK_ID\'])'
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


if __name__ == "__main__":
    def singlerun(job):
        import ecell4
        print("Hi, I'm in local!")
        print("job => {}".format(str(job)))
        return job['x'] + job['y']

    jobs = [{'x': i, 'y': i ** 2} for i in range(1, 4)]
    print(run_serial(singlerun, jobs, n=2))
    print(run_multiprocessing(singlerun, jobs, n=2))

    # environ = {'LD_LIBRARY_PATH': '/home/kaizu/lily_kaizu/src/ecell4/local/lib', 'PYTHONPATH': '/home/kaizu/lily_kaizu/src/ecell4/local/lib/python3.4/site-packages'}
    environ = {}
    # print(run_sge(singlerun, jobs, n=2, delete=False, environ=environ))
    print(run_sge(singlerun, jobs, n=2, environ=environ))

