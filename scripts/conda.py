#!/usr/bin/env python
'''

    Script to interact with the existing conda installation.

    It assumes that a conda environment is active, and will
    run conda commands to:

    * Get the install folder for conda
    * Get a list of the existing environment
    * Pin all packages in all conda environments

    References:
    * https://conda.io/docs/user-guide/tasks/manage-pkgs.html#preventing-packages-from-updating-pinning
    * https://github.com/conda/conda/issues/7375

'''

import os
import subprocess
import json


def run_command(statement):
    os.environ.update({'BASH_ENV': os.path.join(os.environ['HOME'],'.bashrc')})
    process = subprocess.Popen(statement,
                               cwd=os.getcwd(),
                               shell=True,
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               env=os.environ.copy())

    (out, err) = process.communicate()

    if len(err) > 0:
        issue = (' Problem running command: {}\n'
                 ' Stderr was: {}\n'
                 .format(statement, err))
        raise Exception(issue)

    return out, err


statement = "conda info --json"
(out, err) = run_command(statement)
cmd_out = json.loads(out)
install_folder = cmd_out['conda_prefix']
#print(install_folder)

statement = "conda env list --json"
(out, err) = run_command(statement)
cmd_out = json.loads(out)
envs = cmd_out['envs']
#print(envs)

for env_folder in envs:
    if env_folder.startswith(install_folder) and \
       env_folder != install_folder:
        env_name = os.path.basename(env_folder)
        statement = "source activate {} && ".format(env_name)
        statement+= "conda list --json"
        (out, err) = run_command(statement)
        cmd_out = json.loads(out)
        conda_meta = install_folder + '/envs/' + env_name + '/conda-meta/pinned'
        with open(conda_meta, 'w') as outf:
            for pkg in cmd_out:
                outf.write(pkg['name'] + ' ==' + pkg['version'] + '\n')


