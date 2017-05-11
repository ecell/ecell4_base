import subprocess
import time
import re
import tempfile
import os
import os.path
import logging
import collections

# PREFIX = '/usr/bin'
# QSUB_CMD = os.path.join(PREFIX, 'qsub')
# QSTAT_CMD = os.path.join(PREFIX, 'qstat')
# QDEL_CMD = os.path.join(PREFIX, 'qdel')

rcParams = {}
rcParams["PREFIX"] = "/usr/bin"
rcParams["QSUB"] = "qsub"
rcParams["QSTAT"] = "qstat"
rcParams["QDEL"] = "qdel"

def get_logger():
    return logging.getLogger('sge')

def run(jobs, n=1, path='.', wc_queue_list=None, sync=10, delete=True, opts=None):
    if not isinstance(jobs, collections.Iterable):
        return singlerun(jobs, n, path, wc_queue_list, sync, delete, opts)

    retval = []
    for job in jobs:
        retval.append(singlerun(job, n, path, wc_queue_list, 0, delete, opts))
    if sync > 0:
        try:
            wait([jobid for jobid, name, filename in retval], sync)
        finally:
            if delete:
                for jobid, name, filename in retval:
                    os.remove(filename)
    return [(jobid, name) for jobid, name, filename in retval]

def singlerun(job, n=1, path='.', wc_queue_list=None, sync=10, delete=True, opts=None):
    (fd, filename) = tempfile.mkstemp(suffix='.job', prefix='sge-', dir=path, text=True)
    with os.fdopen(fd, 'w') as fout:
        fout.write(job)

    (jobid, name) = submit(filename, n, path, path, wc_queue_list, opts)

    if sync > 0:
        try:
            wait(jobid, sync)
        finally:
            if delete:
                os.remove(filename)

    return (jobid, name, filename)

def collect(jobid, name, n=1, path='.', delete=True):
    outputs = []
    for i in range(n):
        err = False
        filename = os.path.join(path, '{}.e{}.{}'.format(name, jobid, i + 1))
        output = open(filename, 'r').read()
        if output != "":
            # err = True
            for line in output.split('\n'):
                get_logger().error(
                    "A standard error stream [{}] displays: {}".format(filename, line))
        if not err and delete:
            os.remove(filename)

        filename = os.path.join(path, '{}.o{}.{}'.format(name, jobid, i + 1))
        output = open(filename, 'r').read()
        outputs.append(output)
        if not err and delete:
            os.remove(filename)
    return outputs

def submit(job, n=1, epath='.', opath='.', wc_queue_list=None, opts=None):
    output = subprocess.check_output(
        [os.path.join(rcParams["PREFIX"], rcParams["QSUB"]), '-cwd']
        + (['-q', wc_queue_list] if wc_queue_list is not None else [])
        + (opts or []) + ['-e', epath, '-o', opath, '-t', '1-{:d}'.format(n), job])
    output = output.decode('utf-8')
    get_logger().debug(output.strip())

    #XXX: Your job-array 21.1-1:1 ("sge-date") has been submitted
    jobarray = output.split()[2]
    jobid = int(jobarray.split('.')[0])
    name = output.split()[3][2: -2]
    return (jobid, name)

def wait(jobids, interval=10):
    if isinstance(jobids, collections.Iterable):
        jobidstrs = [str(jobid) for jobid in jobids]
    else:
        jobidstrs = [str(jobids)]

    dowait = True
    try:
        while dowait:
            output = subprocess.check_output([os.path.join(rcParams["PREFIX"], rcParams["QSTAT"])])
            output = output.decode('utf-8')
            for line in output.split('\n'):
                get_logger().debug(line)

            dowait = False
            for line in output.split('\n'):
                state = line.split()
                if len(state) < 5 or state[0] not in jobidstrs:
                    continue

                #XXX: job-ID prior   name       user         state submit/start at     queue                          slots ja-task-ID
                jobid = int(state[0])
                if re.search(state[4], 'qwrt'):
                    get_logger().info(
                        'Job {:d} must be queued, running or being transferred'.format(jobid))
                    dowait = True
                    break
                elif re.search(state[4], 'acuE'):
                    get_logger().error('Job {:d} in error state'.format(jobid))
                else:
                    get_logger().error('Unknown state {:s}'.format(state[4]))

            if dowait:
                time.sleep(interval)
                get_logger().info(
                    "Waiting for jobids {:s} to finish".format(str(jobids)))
    finally:
        if dowait:
            output = subprocess.check_output([os.path.join(rcParams["PREFIX"], rcParams["QDEL"])] + jobidstrs)
            get_logger().debug(output.strip())


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    n = 3
    (jobid, name, filename) = singlerun("#!/bin/bash\ndate\nsleep 5\npwd\necho \"puke\"", n)
    print(collect(jobid, name, n))
