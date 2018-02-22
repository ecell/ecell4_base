import subprocess
import time
import re
import tempfile
import os
import os.path
import logging
import collections

rcParams = {}
rcParams["PREFIX"] = "/usr/bin"
# rcParams["QSUB"] = "qsub"
# rcParams["QSTAT"] = "qstat"
# rcParams["QDEL"] = "qdel"
rcParams["QSUB"] = "sbatch"
rcParams["QSTAT"] = "squeue"
rcParams["QDEL"] = "scancel"

def get_logger():
    return logging.getLogger('sge')

def run(jobs, n=1, path='.', sync=10, delete=True, extra_args=None, max_running_tasks=None):
    if not isinstance(jobs, collections.Iterable):
        return singlerun(jobs, n, path, sync, delete, extra_args, max_running_tasks)

    retval = []
    for job in jobs:
        retval.append(singlerun(job, n, path, 0, delete, extra_args, max_running_tasks))
    if sync > 0:
        try:
            wait([jobid for jobid, name, filename in retval], sync)
        finally:
            if delete:
                for jobid, name, filename in retval:
                    os.remove(filename)
    return [(jobid, name) for jobid, name, filename in retval]

def singlerun(job, n=1, path='.', sync=10, delete=True, extra_args=None, max_running_tasks=None):
    (fd, filename) = tempfile.mkstemp(suffix='.job', prefix='sge-', dir=path, text=True)
    with os.fdopen(fd, 'w') as fout:
        fout.write(job)

    (jobid, name) = submit(filename, n, path, path, extra_args, max_running_tasks)

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
        # err = False
        # filename = os.path.join(path, '{}.e{}.{}'.format(name, jobid, i + 1))
        # output = open(filename, 'r').read()
        # if output != "":
        #     # err = True
        #     for line in output.split('\n'):
        #         get_logger().error(
        #             "A standard error stream [{}] displays: {}".format(filename, line))
        # if not err and delete:
        #     os.remove(filename)

        # filename = os.path.join(path, '{}.o{}.{}'.format(name, jobid, i + 1))
        # output = open(filename, 'r').read()
        # outputs.append(output)
        # if not err and delete:
        #     os.remove(filename)

        err = False
        filename = os.path.join(path, 'slurm-{}_e{}.out'.format(jobid, i + 1))
        output = open(filename, 'r').read()
        if output != "":
            # err = True
            for line in output.split('\n'):
                get_logger().error(
                    "A standard error stream [{}] displays: {}".format(filename, line))
        if not err and delete:
            os.remove(filename)

        filename = os.path.join(path, 'slurm-{}_o{}.out'.format(jobid, i + 1))
        output = open(filename, 'r').read()
        outputs.append(output)
        if not err and delete:
            os.remove(filename)

    return outputs

def submit(job, n=1, epath='.', opath='.', extra_args=None, max_running_tasks=None):
    # output = subprocess.check_output(
    #     [os.path.join(rcParams["PREFIX"], rcParams["QSUB"]), '-cwd']
    #     + (['-tc', str(max_running_tasks)] if max_running_tasks is not None else []) + (extra_args or [])
    #     + ['-e', epath, '-o', opath, '-t', '1-{:d}'.format(n), job])
    cmd = ([os.path.join(rcParams["PREFIX"], rcParams["QSUB"])]
        + (extra_args or [])
        + ['-e', os.path.join(epath, 'slurm-%A_e%a.out')]
        + ['-o', os.path.join(opath, 'slurm-%A_o%a.out')]
        + ['-a', '1-{:d}{}'.format(n, '%{}'.format(max_running_tasks) if max_running_tasks is not None else ''), job])

    try:
        output = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as err:
        get_logger().error(err.stdout.decode('utf-8').rstrip())
        raise err

    output = output.decode('utf-8')
    get_logger().debug(output.strip())

    #XXX: Submitted batch job 202
    jobid = int(output.strip().split()[-1])
    name = None

    # #XXX: Your job-array 21.1-1:1 ("sge-date") has been submitted
    # jobarray = output.split()[2]
    # jobid = int(jobarray.split('.')[0])
    # name = output.split()[3][2: -2]
    return (jobid, name)

def wait(jobids, interval=10):
    # """
    # (python3.5) kaizu@lupin:~/src/ecell4-develop/build/test$ squeue -u kaizu
    #              JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
    #              207_1      defq sge-cwfu    kaizu CG       0:00      1 ibn01
    #              207_2      defq sge-cwfu    kaizu CG       0:00      1 ibn01
    #              208_1      defq sge-org6    kaizu CG       0:00      1 ibn01
    #              208_2      defq sge-org6    kaizu CG       0:00      1 ibn01
    #              209_1      defq sge-wp2q    kaizu CG       0:00      1 ibn01
    #              209_2      defq sge-wp2q    kaizu CG       0:00      1 ibn01
    # (python3.5) kaizu@lupin:~/src/ecell4-develop/build/test$ squeue -u kaizu
    #              JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
    #          215_[1-2]      defq sge-n4qk    kaizu PD       0:00      1 (Resources)
    #          216_[1-2]      defq sge-aop9    kaizu PD       0:00      1 (Resources)
    #          217_[1-2]      defq sge-cjkq    kaizu PD       0:00      1 (Resources)
    # """

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

            # for line in output.split('\n'):
            #     state = line.split()
            #     if len(state) < 5 or state[0] not in jobidstrs:
            #         continue

            #     #XXX: job-ID prior   name       user         state submit/start at     queue                          slots ja-task-ID
            #     jobid = int(state[0])
            #     if re.search(state[4], 'qwrt'):
            #         get_logger().info(
            #             'Job {:d} must be queued, running or being transferred'.format(jobid))
            #         dowait = True
            #         break
            #     elif re.search(state[4], 'acuE'):
            #         get_logger().error('Job {:d} in error state'.format(jobid))
            #     else:
            #         get_logger().error('Unknown state {:s}'.format(state[4]))

            for line in output.split('\n'):
                state = line.split()
                if len(state) < 8:
                    continue

                jobid = state[0].split('_')[0]
                if jobid not in jobidstrs:
                    continue

                jobid = int(jobid)
                if state[4] in ('PD', 'R', 'CF', 'CG', 'CD'):
                    get_logger().info(
                        'Job {:d} must be queued, running or being transferred'.format(jobid))
                    dowait = True
                    break
                elif state[4] in ('CA', 'F', 'TO', 'NF', 'RV', 'SE'):
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
