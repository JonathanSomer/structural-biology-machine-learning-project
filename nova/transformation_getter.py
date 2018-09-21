import time
from paramiko import SSHClient
from scp import SCPClient
import sys
import argparse
import os
from paramiko import AutoAddPolicy



#############################################################################
# 
# 1. must have a PatchDock directory on nova in order to run this!
# 
# 
# 
# 
#############################################################################

def progress(filename, size, sent):
    sys.stdout.write("%s\'s progress: %.2f%%   \r" % (filename, float(sent)/float(size)*100) )


def main(complex_code, num_transformations, username, password, test=False, force=False, keep_transformations=False):
	ssh = SSHClient()
	ssh.set_missing_host_key_policy(AutoAddPolicy())
	ssh.load_system_host_keys()
	ssh.connect('nova.cs.tau.ac.il', username=username, password=password)
	scp = SCPClient(ssh.get_transport(), progress = progress)

	c = complex_code
	if not os.path.exists('../data/{}/patch_dock/'.format(c)):
		os.makedirs('../data/{}/patch_dock/'.format(c))

	ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command('cd ~/PatchDock; rm ../data/{}/patch_dock/{}.patch_dock_output.*'.format(c,c))
	print(ssh_stdout.read(), ssh_stderr.read())

	ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command(
		'cd ~/PatchDock; ./transOutput.pl ../data/{}/patch_dock/{}.patch_dock_output {} {}'.format(c, c, str(1),

																							   str(num_transformations)))

	print(ssh_stdout.read(), ssh_stderr.read())
	print("sleeps 20 sec ...")
	time.sleep(20)
	print("start fetching results from nova...")
	for i in range(1, num_transformations+1):
		if (i-1)%100 == 0: print(i)
		try:
			scp.get('~/data/{}/patch_dock/{}.patch_dock_output.{}.pdb'.format(c,c, str(i)), '../data/{}/patch_dock/'.format(c))
		except Exception as e:

			print(e.message)
			print("ERROR with: %s" % complex_code)

	scp.close()
	ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command('cd ~/PatchDock; rm ../data/{}/patch_dock/{}.patch_dock_output.*'.format(c,c))
	print(ssh_stdout.read(), ssh_stderr.read())

	if not keep_transformations:
		for i in range(1, num_transformations+1):
			try:
				os.remove('../data/{}/patch_dock/{}.patch_dock_output.{}.pdb'.format(c,c,str(i)))
			except:
				print("error with transformation number %d" % i)
			

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--complex', help='complex code to fetch transformations for')
    parser.add_argument('-n', '--num_transformations', help='number of transformations to fetch')
    parser.add_argument('-u', '--username', help='universiy username')
    parser.add_argument('-p', '--password', help='universiy username\'s password')
    parser.add_argument('-t', '--test', help='perform on single complex to check flow is alright', default=False, action='store_true')
    parser.add_argument('-f', '--force', help='re uploads and runs everything', default=False, action='store_true')
    parser.add_argument('-k', '--keep_transformations', help='do not delete transformations when done', default=False, action='store_true')
    args = parser.parse_args()


