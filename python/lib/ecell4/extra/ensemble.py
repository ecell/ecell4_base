import os
import logging
import tempfile
import pickle
import inspect
import textwrap
import re
import types

import multiprocessing
import sge

def run_serial(target, jobs):
    for job in jobs:
        target(job)

def run_multiprocessing(target, jobs):
    processes = []
    for job in jobs:
        p = multiprocessing.Process(target=target, args=(job, ))
        p.start()
        processes.append(p)

    for p in processes:
        p.join()

def run_sge(target, jobs, n=1, path='.', delete=True):
    logging.basicConfig(level=logging.DEBUG)

    if isinstance(target, types.LambdaType) and target.__name__ == "<lambda>":
        raise RuntimeError("A lambda function is not accepted")

    src = textwrap.dedent(inspect.getsource(singlerun)).replace(r'"', r'\"')
    if re.match('[\s\t]+', src.split('\n')[0]) is not None:
        raise RuntimeError(
            "Wrong indentation was found in the source translated")

    cmds = []
    pickles = []
    for job in jobs:
        (fd, picklename) = tempfile.mkstemp(suffix='.pickle', prefix='sge-', dir=path)
        with os.fdopen(fd, 'wb') as fout:
            pickle.dump(job, fout)
        pickles.append(picklename)

        cmd = '#!/bin/bash\n'
        # cmd += 'export LD_LIBRARY_PATH=/home/kaizu/lily_kaizu/src/ecell4/local/lib\n'
        # cmd += 'export PYTHONPATH=/home/kaizu/lily_kaizu/src/ecell4/local/lib/python3.4/site-packages\n'
        cmd += 'python3 -c "\n'
        cmd += 'import sys\n'
        cmd += 'import pickle\n'
        cmd += 'with open(sys.argv[1], \'rb\') as fin:\n'
        cmd += '    job = pickle.load(fin)\n'
        cmd += 'pass\n'
        cmd += src
        cmd += '\n{:s}(job)" {:s}\n'.format(target.__name__, picklename)
        cmds.append(cmd)

    jobids = sge.run(cmds, n=n, path=path, delete=delete)

    for jobid, name in jobids:
        outputs = sge.collect(jobid, name, n=n, path=path, delete=delete)
        for output in outputs:
            print(output, end='')

    if delete:
        for picklename in pickles:
            os.remove(picklename)


if __name__ == "__main__":
    def singlerun(job):
        print("Hi, I'm in local!")
        print("job => {}".format(str(job)))

    jobs = [{'x': i, 'y': i ** 2} for i in range(1, 4)]
    # run_serial(jobs)
    # run_multiprocessing(singlerun, jobs)
    run_sge(singlerun, jobs, delete=False)
