from paramiko import SSHClient
from paramiko import AutoAddPolicy
import os

'''
This script removes all patch dock results from nova
'''

def main(username, password):
    ssh = SSHClient()
    ssh.set_missing_host_key_policy(AutoAddPolicy())
    ssh.load_system_host_keys()
    ssh.connect('nova.cs.tau.ac.il', username=username, password=password)
    all_complex_codes = next(os.walk('../data'))[1]
    for c_id in all_complex_codes:
        ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command(
            'cd ~/PatchDock; rm ../data/{}/patch_dock/{}.patch_dock_output.*'.format(c_id, c_id))
        print(ssh_stdout.read(), ssh_stderr.read())

