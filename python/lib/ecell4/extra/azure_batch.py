# azure_batch.py
# Copyright (c) 2017 Kazunari Kaizu
# Released under the GNU General Public License

# python_tutorial_client.py
# Copyright (c) 2017 Microsoft Corporation
# Released under the MIT license

from __future__ import print_function
import datetime
import os
import sys
import time
import binascii
import pickle
import textwrap
import inspect
import itertools
import re
import io
from logging import getLogger
_log = getLogger(__name__)

import azure.storage.blob as azureblob
import azure.batch.batch_service_client as batch
import azure.batch.batch_auth as batchauth
import azure.batch.models as batchmodels

_STANDARD_OUT_FILE_NAME = 'stdout.txt'
_STANDARD_ERROR_FILE_NAME = 'stderr.txt'
_SAMPLES_CONFIG_FILE_NAME = 'configuration.cfg'

try:
    import configparser
except ImportError:
    import ConfigParser as configparser


def print_batch_exception(batch_exception):
    """Prints the contents of the specified Batch exception.

    :param batch_exception:
    """
    _log.error('-------------------------------------------')
    _log.error('Exception encountered:')
    if batch_exception.error and \
            batch_exception.error.message and \
            batch_exception.error.message.value:
        _log.error(batch_exception.error.message.value)
        if batch_exception.error.values:
            _log.error('')
            for mesg in batch_exception.error.values:
                _log.error('{}:\t{}'.format(mesg.key, mesg.value))
    _log.error('-------------------------------------------')

def upload_file_to_container(block_blob_client, container_name, file_path):
    """Uploads a local file to an Azure Blob storage container.

    :param block_blob_client: A blob service client.
    :type block_blob_client: `azure.storage.blob.BlockBlobService`
    :param str container_name: The name of the Azure Blob storage container.
    :param str file_path: The local path to the file.
    :rtype: `azure.batch.models.ResourceFile`
    :return: A ResourceFile initialized with a SAS URL appropriate for Batch
    tasks.
    """
    blob_name = os.path.basename(file_path)

    _log.info('Uploading file {} to container [{}]...'.format(file_path, container_name))

    block_blob_client.create_blob_from_path(container_name,
                                            blob_name,
                                            file_path)

    sas_token = block_blob_client.generate_blob_shared_access_signature(
        container_name,
        blob_name,
        permission=azureblob.BlobPermissions.READ,
        expiry=datetime.datetime.utcnow() + datetime.timedelta(hours=2))

    sas_url = block_blob_client.make_blob_url(container_name,
                                              blob_name,
                                              sas_token=sas_token)

    return batchmodels.ResourceFile(file_path=blob_name,
                                    blob_source=sas_url)

def get_container_sas_token(block_blob_client,
                            container_name, blob_permissions):
    """Obtains a shared access signature granting the specified permissions to the
    container.

    :param block_blob_client: A blob service client.
    :type block_blob_client: `azure.storage.blob.BlockBlobService`
    :param str container_name: The name of the Azure Blob storage container.
    :param BlobPermissions blob_permissions:
    :rtype: str
    :return: A SAS token granting the specified permissions to the container.
    """
    # Obtain the SAS token for the container, setting the expiry time and
    # permissions. In this case, no start time is specified, so the shared
    # access signature becomes valid immediately.
    container_sas_token = \
        block_blob_client.generate_container_shared_access_signature(
            container_name,
            permission=blob_permissions,
            expiry=datetime.datetime.utcnow() + datetime.timedelta(hours=2))

    return container_sas_token

def select_latest_verified_vm_image_with_node_agent_sku(
        batch_client, publisher, offer, sku_starts_with):
    """Select the latest verified image that Azure Batch supports given
    a publisher, offer and sku (starts with filter).
    Originally in azure-batch-samples.Python.Batch.common.helpers

    :param batch_client: The batch client to use.
    :type batch_client: `batchserviceclient.BatchServiceClient`
    :param str publisher: vm image publisher
    :param str offer: vm image offer
    :param str sku_starts_with: vm sku starts with filter
    :rtype: tuple
    :return: (node agent sku id to use, vm image ref to use)
    """
    # get verified vm image list and node agent sku ids from service
    node_agent_skus = batch_client.account.list_node_agent_skus()
    # pick the latest supported sku
    skus_to_use = [
        (sku, image_ref) for sku in node_agent_skus for image_ref in sorted(
            sku.verified_image_references, key=lambda item: item.sku)
        if image_ref.publisher.lower() == publisher.lower() and
        image_ref.offer.lower() == offer.lower() and
        image_ref.sku.startswith(sku_starts_with)
    ]
    # skus are listed in reverse order, pick first for latest
    sku_to_use, image_ref_to_use = skus_to_use[0]
    return (sku_to_use.id, image_ref_to_use)

def wrap_commands_in_shell(ostype, commands):
    """Wrap commands in a shell
    Originally in azure-batch-samples.Python.Batch.common.helpers

    :param list commands: list of commands to wrap
    :param str ostype: OS type, linux or windows
    :rtype: str
    :return: a shell wrapping commands
    """
    if ostype.lower() == 'linux':
        return '/bin/bash -c \'set -e; set -o pipefail; {}; wait\''.format(
            ';'.join(commands))
    elif ostype.lower() == 'windows':
        return 'cmd.exe /c "{}"'.format('&'.join(commands))
    else:
        raise ValueError('unknown ostype: {}'.format(ostype))

def create_pool(batch_service_client, pool_id,
                resource_files, publisher, offer, sku,
                task_file, vm_size, node_count):
    """Creates a pool of compute nodes with the specified OS settings.

    :param batch_service_client: A Batch service client.
    :type batch_service_client: `azure.batch.BatchServiceClient`
    :param str pool_id: An ID for the new pool.
    :param list resource_files: A collection of resource files for the pool's
    start task.
    :param str publisher: Marketplace image publisher
    :param str offer: Marketplace image offer
    :param str sku: Marketplace image sku
    :param str task_file: A file name of the script
    :param str vm_size: A type of vm
    :param str node_count: The number of nodes
    """
    _log.info('Creating pool [{}]...'.format(pool_id))

    # Create a new pool of Linux compute nodes using an Azure Virtual Machines
    # Marketplace image. For more information about creating pools of Linux
    # nodes, see:
    # https://azure.microsoft.com/documentation/articles/batch-linux-nodes/

    # Specify the commands for the pool's start task. The start task is run
    # on each node as it joins the pool, and when it's rebooted or re-imaged.
    # We use the start task to prep the node for running our task script.
    task_commands = [
        # Copy the python_tutorial_task.py script to the "shared" directory
        # that all tasks that run on the node have access to. Note that
        # we are using the -p flag with cp to preserve the file uid/gid,
        # otherwise since this start task is run as an admin, it would not
        # be accessible by tasks run as a non-admin user.
        'cp -p {} $AZ_BATCH_NODE_SHARED_DIR'.format(os.path.basename(task_file)),
        # Install pip
        'curl -fSsL https://bootstrap.pypa.io/get-pip.py | python',
        # Install the azure-storage module so that the task script can access
        # Azure Blob storage, pre-cryptography version
        'pip install azure-storage==0.32.0',
        # Install E-Cell 4
        'pip install https://1028-6348303-gh.circle-artifacts.com/0/root/circle/wheelhouse/ecell-4.1.2-cp27-cp27mu-manylinux1_x86_64.whl']

    # Get the node agent SKU and image reference for the virtual machine
    # configuration.
    # For more information about the virtual machine configuration, see:
    # https://azure.microsoft.com/documentation/articles/batch-linux-nodes/
    sku_to_use, image_ref_to_use = \
        select_latest_verified_vm_image_with_node_agent_sku(
            batch_service_client, publisher, offer, sku)
    user = batchmodels.AutoUserSpecification(
        scope=batchmodels.AutoUserScope.pool,
        elevation_level=batchmodels.ElevationLevel.admin)
    new_pool = batch.models.PoolAddParameter(
        id=pool_id,
        virtual_machine_configuration=batchmodels.VirtualMachineConfiguration(
            image_reference=image_ref_to_use,
            node_agent_sku_id=sku_to_use),
        vm_size=vm_size,
        target_dedicated_nodes=0,
        target_low_priority_nodes=node_count,
        start_task=batch.models.StartTask(
            command_line=wrap_commands_in_shell('linux', task_commands),
            user_identity=batchmodels.UserIdentity(auto_user=user),
            wait_for_success=True,
            resource_files=resource_files),
    )

    try:
        batch_service_client.pool.add(new_pool)
    except batchmodels.batch_error.BatchErrorException as err:
        print_batch_exception(err)
        raise

def create_job(batch_service_client, job_id, pool_id):
    """Creates a job with the specified ID, associated with the specified pool.

    :param batch_service_client: A Batch service client.
    :type batch_service_client: `azure.batch.BatchServiceClient`
    :param str job_id: The ID for the job.
    :param str pool_id: The ID for the pool.
    """
    _log.info('Creating job [{}]...'.format(job_id))

    job = batch.models.JobAddParameter(
        job_id,
        batch.models.PoolInformation(pool_id=pool_id))

    try:
        batch_service_client.job.add(job)
    except batchmodels.batch_error.BatchErrorException as err:
        print_batch_exception(err)
        raise

def add_tasks(batch_service_client, job_id, loads,
              output_container_name, output_container_sas_token,
              task_file, acount_name):
    """Adds a task for each input file in the collection to the specified job.

    :param batch_service_client: A Batch service client.
    :type batch_service_client: `azure.batch.BatchServiceClient`
    :param str job_id: The ID of the job to which to add the tasks.
    :param list input_files: A collection of input files. One task will be
     created for each input file.
    :param output_container_name: The ID of an Azure Blob storage container to
    which the tasks will upload their results.
    :param output_container_sas_token: A SAS token granting write access to
    the specified Azure Blob storage container.
    :param str task_file: A file name of the script
    :param str account_name: A storage account
    """

    _log.info('Adding {} tasks to job [{}]...'.format(len(loads), job_id))
    # _log.info('Adding {} tasks to job [{}]...'.format(len(input_files), job_id))

    tasks = list()

    for (input_file, output_file, i, j) in loads:
        command = ['python $AZ_BATCH_NODE_SHARED_DIR/{} '
                   '--filepath {} --output {} --storageaccount {} '
                   '--task_id {} --job_id {} '
                   '--storagecontainer {} --sastoken "{}"'.format(
                       os.path.basename(task_file),
                       input_file.file_path,
                       output_file,
                       acount_name,
                       i, j,
                       output_container_name,
                       output_container_sas_token)]
        _log.debug('CMD : "{}"'.format(command[0]))

        tasks.append(batch.models.TaskAddParameter(
                'topNtask{}-{}'.format(i, j),
                wrap_commands_in_shell('linux', command),
                resource_files=[input_file]
                )
        )

    batch_service_client.task.add_collection(job_id, tasks)

    task_ids = [task.id for task in tasks]
    _log.info('{} tasks were added.'.format(len(task_ids)))
    return task_ids

def wait_for_tasks_to_complete(batch_service_client, job_ids, timeout):
    """Returns when all tasks in the specified job reach the Completed state.

    :param batch_service_client: A Batch service client.
    :type batch_service_client: `azure.batch.BatchServiceClient`
    :param str job_id: The id of the job whose tasks should be to monitored.
    :param timedelta timeout: The duration to wait for task completion. If all
    tasks in the specified job do not reach Completed state within this time
    period, an exception will be raised.
    """
    timeout_expiration = datetime.datetime.now() + timeout

    print("Monitoring all tasks for 'Completed' state, timeout in {}...".format(timeout), end='')

    while datetime.datetime.now() < timeout_expiration:
        print('.', end='')
        sys.stdout.flush()
        # tasks = batch_service_client.task.list(job_id)
        # incomplete_tasks = [task for task in tasks if
        #                     task.state != batchmodels.TaskState.completed]
        for (job_id, _) in job_ids:
            tasks = batch_service_client.task.list(job_id)
            incomplete_tasks = [task for task in tasks if
                                task.state != batchmodels.TaskState.completed]
            if incomplete_tasks:
                break
        if not incomplete_tasks:
            print()
            return True
        else:
            time.sleep(1)

    raise RuntimeError("ERROR: Tasks did not reach 'Completed' state within "
                       "timeout period of " + str(timeout))

def download_blobs_from_container(block_blob_client,
                                  container_name, directory_path,
                                  prefix=None):
    """Downloads all blobs from the specified Azure Blob storage container.

    :param block_blob_client: A blob service client.
    :type block_blob_client: `azure.storage.blob.BlockBlobService`
    :param container_name: The Azure Blob storage container from which to
     download files.
    :param directory_path: The local directory to which to download the files.
    :param str prefix: A name prefix to filter blobs. None as its default
    """
    _log.info('Downloading all files from container [{}]...'.format(container_name))

    container_blobs = block_blob_client.list_blobs(container_name, prefix=None)
    _log.info('{} blobs are found [{}]'.format(len(tuple(container_blobs)), ', '.join(blob.name for blob in container_blobs.items)))

    for blob in container_blobs.items:
        destination_file_path = os.path.join(directory_path, blob.name)

        block_blob_client.get_blob_to_path(container_name,
                                           blob.name,
                                           destination_file_path)

        _log.info('  Downloaded blob [{}] from container [{}] to {}'.format(
            blob.name,
            container_name,
            destination_file_path))

    _log.info('  Download complete!')

def _read_stream_as_string(stream, encoding):
    """Read stream as string
    Originally in azure-batch-samples.Python.Batch.common.helpers

    :param stream: input stream generator
    :param str encoding: The encoding of the file. The default is utf-8.
    :return: The file content.
    :rtype: str
    """
    output = io.BytesIO()
    try:
        for data in stream:
            output.write(data)
        if encoding is None:
            encoding = 'utf-8'
        return output.getvalue().decode(encoding)
    finally:
        output.close()
    raise RuntimeError('could not write data to stream or decode bytes')

def read_task_file_as_string(
        batch_client, job_id, task_id, file_name, encoding=None):
    """Reads the specified file as a string.
    Originally in azure-batch-samples.Python.Batch.common.helpers

    :param batch_client: The batch client to use.
    :type batch_client: `batchserviceclient.BatchServiceClient`
    :param str job_id: The id of the job.
    :param str task_id: The id of the task.
    :param str file_name: The name of the file to read.
    :param str encoding: The encoding of the file. The default is utf-8.
    :return: The file content.
    :rtype: str
    """
    stream = batch_client.file.get_from_task(job_id, task_id, file_name)
    return _read_stream_as_string(stream, encoding)

def print_task_output(batch_client, job_id, task_ids, encoding=None):
    """Prints the stdout and stderr for each task specified.
    Originally in azure-batch-samples.Python.Batch.common.helpers

    :param batch_client: The batch client to use.
    :type batch_client: `batchserviceclient.BatchServiceClient`
    :param str job_id: The id of the job to monitor.
    :param task_ids: The collection of tasks to print the output for.
    :type task_ids: `list`
    :param str encoding: The encoding to use when downloading the file.
    """
    for task_id in task_ids:
        file_text = read_task_file_as_string(
            batch_client,
            job_id,
            task_id,
            _STANDARD_OUT_FILE_NAME,
            encoding)
        print("{} content for task {}: ".format(
            _STANDARD_OUT_FILE_NAME,
            task_id))
        print(file_text)

        file_text = read_task_file_as_string(
            batch_client,
            job_id,
            task_id,
            _STANDARD_ERROR_FILE_NAME,
            encoding)
        print("{} content for task {}: ".format(
            _STANDARD_ERROR_FILE_NAME,
            task_id))
        print(file_text)

def run_azure_batch(target, jobs, n=1, path='.', delete=True, config=None):
    """Execute a function for multiple sets of arguments on Microsoft Azure,
    and return the results as a list.

    :param function target: A target function.
    :param list jobs: A list of sets of arguments given to the target.
    :param int n: The number of repeats running the target. 1 as default.
    :param str path: A path to save temp files. The current path as default.
    :param bool delete: Delete temp files after finishing jobs, or not. True as default.
    :param config: str or configparser.ConfigParser. A config file. An example is the following:

    ```
    [azure]
    batch.name = foo
    batch.key = bar
    batch.url = hoge
    storage.name = fuga
    storage.key = spam
    pool.nodecount = 2
    # pool.id = MyPool
    # pool.vmsize = Standard_D11_v2
    # os.publisher = Canonical
    # os.offer = UbuntuServer
    # os.sku = 16
    # job.id = MyJob
    ```

    :return: A list of results corresponding the `jobs` list.
    :rtype: list
    """

    if config is None:
        raise ValueError('Argument \'config\' must be given.')
    elif isinstance(config, str):
        if not os.path.isfile(config):
            raise FileNotFoundError('A file [{}] could not be found.'.format(config))
        config_filename = config
        config = configparser.ConfigParser()
        config.sections()
        config.read(config_filename)
        config.sections()
    elif not isinstance(config, configparser.ConfigParser):
        raise ValueError('\'config\' must be eighter str or ConfigParser. [{}] was given.'.format(repr(config)))

    if 'azure' not in config:
        raise KeyError('Key \'azure\' could not be found in the given config.')

    for key in ('batch.name', 'batch.key', 'batch.url', 'storage.name', 'storage.key', 'pool.nodecount'):
        if key not in config['azure']:
            raise KeyError('Key \'{}\' could not be found in the \'azure\' section.'.format(key))

    # Update the Batch and Storage account credential strings below with the values
    # unique to your accounts. These are used when constructing connection strings
    # for the Batch and Storage client objects.
    _BATCH_ACCOUNT_NAME   = config['azure']['batch.name']
    _BATCH_ACCOUNT_KEY    = config['azure']['batch.key']
    _BATCH_ACCOUNT_URL    = config['azure']['batch.url']
    _STORAGE_ACCOUNT_NAME = config['azure']['storage.name']
    _STORAGE_ACCOUNT_KEY  = config['azure']['storage.key']
    _POOL_NODE_COUNT      = config['azure']['pool.nodecount']
    _POOL_ID              = config['azure'].get('pool.id', 'MyPool')
    _POOL_VM_SIZE         = config['azure'].get('pool.vmsize', 'Standard_D11_v2')
    _NODE_OS_PUBLISHER    = config['azure'].get('os.publisher', 'Canonical')
    _NODE_OS_OFFER        = config['azure'].get('os.offer', 'UbuntuServer')
    _NODE_OS_SKU          = config['azure'].get('os.sku', '16')
    _JOB_ID               = config['azure'].get('job.id', 'MyJob')

    if not _POOL_NODE_COUNT.isdigit():
        raise ValueError('The wrong pool node count was given [{}]. This must be an integer'.format(_POOL_NODE_COUNT))

    proc_per_node = 2  #XXX: Does this depend on pool vm?
    nproc = int(_POOL_NODE_COUNT) * proc_per_node

    code_header = """
from __future__ import print_function
import argparse
import os
import string
import azure.storage.blob as azureblob

parser = argparse.ArgumentParser()
parser.add_argument('--filepath', required=True,
                    help='The path to the text file to process. The path'
                         'may include a compute node\\'s environment'
                         'variables, such as'
                         '$AZ_BATCH_NODE_SHARED_DIR/filename.txt')
parser.add_argument('--output', required=True,
                    help='The path to the output.')
parser.add_argument('--job_id', type=int, required=True)
parser.add_argument('--task_id', type=int, required=True)
parser.add_argument('--storageaccount', required=True,
                    help='The name the Azure Storage account that owns the'
                         'blob storage container to which to upload'
                         'results.')
parser.add_argument('--storagecontainer', required=True,
                    help='The Azure Blob storage container to which to'
                         'upload results.')
parser.add_argument('--sastoken', required=True,
                    help='The SAS token providing write access to the'
                         'Storage container.')
args = parser.parse_args()

input_file = os.path.realpath(args.filepath)
output_file = args.output

import pickle
with open(input_file, mode='rb') as fin:
    inputs = pickle.load(fin)

"""
    code_footer = """

with open(output_file, mode='wb') as fout:
    pickle.dump(res, fout, protocol=2)

# Create the blob client using the container's SAS token.
# This allows us to create a client that provides write
# access only to the container.
blob_client = azureblob.BlockBlobService(account_name=args.storageaccount,
                                         sas_token=args.sastoken)
output_file_path = os.path.realpath(output_file)
blob_client.create_blob_from_path(args.storagecontainer,
                                  output_file,
                                  output_file_path)
"""

    src = textwrap.dedent(inspect.getsource(target)).replace(r'"', r'\"')
    if re.match('[\s\t]+', src.split('\n')[0]) is not None:
        raise RuntimeError(
            "Wrong indentation was found in the source translated")

    code = code_header
    code += src
    code += 'res = {}(inputs, args.task_id, args.job_id)'.format(target.__name__)
    code += code_footer
    target = code

    suffix = binascii.hexlify(os.urandom(4)).decode()

    start_time = datetime.datetime.now().replace(microsecond=0)
    _log.info('Sample start: {}'.format(start_time))

    if not os.path.isdir(path):
        os.mkdir(path)

    # task_file = target
    # task_file = 'task-{}.py'.format(suffix)
    task_file = '{}/task-{}.py'.format(path, suffix)
    with open(task_file, 'w') as fout:
        fout.write(target)

    # Prepare input pickle files
    input_file_names = []
    output_file_names = []
    for i, job in enumerate(jobs):
        filename = '{}/input-{}_{}.pickle'.format(path, suffix, i)
        input_file_names.append(filename)
        for j in range(n):
            output_file_names.append('output-{}_{}.{}.pickle'.format(suffix, i, j + 1))
        with open(filename, mode='wb') as fout:
            pickle.dump(job, fout, protocol=2)

    # Create the blob client, for use in obtaining references to
    # blob storage containers and uploading files to containers.
    blob_client = azureblob.BlockBlobService(
        account_name=_STORAGE_ACCOUNT_NAME,
        account_key=_STORAGE_ACCOUNT_KEY)

    n_jobs = -(-(len(jobs) * n) // nproc)  # ceil for int
    _log.info('{} jobs will be created.'.format(n_jobs))
    res = None

    try:
        # Use the blob client to create the containers in Azure Storage if they
        # don't yet exist.
        app_container_name = 'application'
        input_container_name = 'input'
        output_container_name = 'output'
        blob_client.create_container(app_container_name, fail_on_exist=False)
        blob_client.create_container(input_container_name, fail_on_exist=False)
        blob_client.create_container(output_container_name, fail_on_exist=False)

        # Paths to the task script. This script will be executed by the tasks that
        # run on the compute nodes.
        application_file_paths = [os.path.realpath(task_file)]

        # The collection of data files that are to be processed by the tasks.
        input_file_paths = [os.path.realpath(filename) for filename in input_file_names]

        # Upload the application script to Azure Storage. This is the script that
        # will process the data files, and is executed by each of the tasks on the
        # compute nodes.
        application_files = [
            upload_file_to_container(blob_client, app_container_name, file_path)
            for file_path in application_file_paths]

        # Upload the data files. This is the data that will be processed by each of
        # the tasks executed on the compute nodes in the pool.
        input_files = [
            upload_file_to_container(blob_client, input_container_name, file_path)
            for file_path in input_file_paths]

        # Obtain a shared access signature that provides write access to the output
        # container to which the tasks will upload their output.
        output_container_sas_token = get_container_sas_token(
            blob_client,
            output_container_name,
            azureblob.BlobPermissions.WRITE)

        # Create a Batch service client. We'll now be interacting with the Batch
        # service in addition to Storage
        credentials = batchauth.SharedKeyCredentials(_BATCH_ACCOUNT_NAME,
                                                     _BATCH_ACCOUNT_KEY)

        batch_client = batch.BatchServiceClient(
            credentials,
            base_url=_BATCH_ACCOUNT_URL)

        # Create the pool that will contain the compute nodes that will execute the
        # tasks. The resource files we pass in are used for configuring the pool's
        # start task, which is executed each time a node first joins the pool (or
        # is rebooted or re-imaged).
        create_pool(batch_client,
                    _POOL_ID + '-' + suffix,
                    application_files,
                    _NODE_OS_PUBLISHER,
                    _NODE_OS_OFFER,
                    _NODE_OS_SKU,
                    task_file,
                    _POOL_VM_SIZE, _POOL_NODE_COUNT)

        # Create the job that will run the tasks.
        loads = []
        for i, input_file in enumerate(input_files):
            for j, output_file in enumerate(output_file_names[i * n: (i + 1) * n]):
                loads.append((input_file, output_file, i + 1, j + 1))

        assert n_jobs == -(-len(loads) // nproc)  # ceil for int
        job_names = []
        for i in range(n_jobs):
            job_name = '{}-{}-{}'.format(_JOB_ID, suffix, i + 1)

            create_job(batch_client, job_name, _POOL_ID + '-' + suffix)

            # Add the tasks to the job. We need to supply a container shared access
            # signature (SAS) token for the tasks so that they can upload their output
            # to Azure Storage.
            task_ids = add_tasks(batch_client,
                                 job_name,
                                 loads[i * nproc: (i + 1) * nproc],
                                 output_container_name,
                                 output_container_sas_token,
                                 task_file,
                                 _STORAGE_ACCOUNT_NAME)

            job_names.append((job_name, task_ids))

        # Pause execution until tasks reach Completed state.
        wait_for_tasks_to_complete(batch_client,
                                   job_names,
                                   datetime.timedelta(minutes=20))

        _log.info("  Success! All tasks reached the 'Completed' state within the specified timeout period.")

        # Download the task output files from the output Storage container to a
        # local directory. Note that we could have also downloaded the output
        # files directly from the compute nodes themselves.
        download_blobs_from_container(blob_client,
                                      output_container_name,
                                      os.path.abspath(path))

        for job_id, task_ids in job_names:
            print_task_output(batch_client, job_id, task_ids)

        # Print out some timing info
        end_time = datetime.datetime.now().replace(microsecond=0)
        _log.info('Sample end: {}'.format(end_time))
        _log.info('Elapsed time: {}'.format(end_time - start_time))

        res = []
        for output_file in output_file_names:
            with open(os.path.join(path, output_file), mode='rb') as fin:
                res.append(pickle.load(fin))
        res = [res[i * n: (i + 1) * n] for i in range(len(jobs))]
    finally:
        # Clean up storage resources
        _log.info('Deleting containers...')
        blob_client.delete_container(app_container_name)
        blob_client.delete_container(input_container_name)
        blob_client.delete_container(output_container_name)

        # Clean up Batch resources (if the user so chooses).
        for i in range(n_jobs):
            job_name = '{}-{}-{}'.format(_JOB_ID, suffix, i + 1)
            _log.info('Deleting job [{}] ...'.format(job_name))
            batch_client.job.delete(job_name)

        _log.info('Deleting pool...')
        batch_client.pool.delete(_POOL_ID + '-' + suffix)

        if delete:
            _log.info('Deleting temporary files...')
            for filename in output_file_names:
                filename = os.path.join(path, filename)
                if os.path.isfile(filename):
                    os.remove(filename)
            for filename in itertools.chain(input_file_paths, application_file_paths):
                if os.path.isfile(filename):
                    os.remove(filename)

    return res

def singlerun(job, task_id=0, job_id=0):
    """This task is for an example."""

    import ecell4
    print('ecell4.__version__ = {:s}'.format(ecell4.__version__))
    print('job={}, task_id={}, job_id={}'.format(str(job), task_id, job_id))

    with ecell4.reaction_rules():
        A + B == C | (0.01, 0.3)

    res = ecell4.run_simulation(
        10.0,
        y0={'A': job[0], 'B': job[1], 'C': job[2]},
        rndseed=job_id,
        solver='gillespie',
        return_type='array')

    print('A simulation was successfully done.')
    return res


if __name__ == '__main__':
    from logging import basicConfig, StreamHandler, DEBUG
    # basicConfig(level=DEBUG)
    handler = StreamHandler()
    handler.setLevel(DEBUG)
    getLogger(__name__).setLevel(DEBUG)
    getLogger(__name__).addHandler(handler)

    # jobs = [(n, n, n) for n in range(10, 70, 10)]
    jobs = [(30, 30, 30), (60, 60, 60)]
    res = run_azure_batch(singlerun, jobs, n=2, path='tmp', config='example.ini')
    print(res)

    import numpy
    import matplotlib.pylab as plt

    for i, dataset in enumerate(res):
        for j, data in enumerate(dataset):
            data = numpy.array(data).T
            plt.plot(data[0], data[3], '-', label='task{}-{}'.format(i, j))
    plt.xlabel('Time')
    plt.ylabel('# of Molecules')
    # plt.legend(loc='best')
    # plt.savefig('res.png')
    plt.show()
