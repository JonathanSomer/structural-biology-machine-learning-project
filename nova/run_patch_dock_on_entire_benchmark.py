from paramiko import SSHClient
from scp import SCPClient
import sys
import argparse
import os
import json
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



def main(username, password, test=False, force=False):
	ssh = SSHClient()
	ssh.set_missing_host_key_policy(AutoAddPolicy())
	ssh.load_system_host_keys()
	ssh.connect('nova.cs.tau.ac.il', username=username, password=password)

	scp = SCPClient(ssh.get_transport(), progress = progress)

	# copy benchmark unbound receptor and ligand to nova:
	all_complex_codes = next(os.walk('../data'))[1]

	if test:
		all_complex_codes = all_complex_codes[:1]

	with open('./copied_successfully.json', 'r') as f:
		copied_successfully = json.load(f)

	if force:
		copied_successfully = []

	for complex_code in all_complex_codes:
		if len(complex_code) != 4:
			continue

		print("Copying complex: %s" % complex_code)

		try:
			ssh.exec_command('mkdir ~/data/' + complex_code + '/')
			ssh.exec_command('mkdir ~/data/' + complex_code + '/benchmark/')
			ssh.exec_command('mkdir ~/data/' + complex_code + '/patch_dock/')
			scp.put('../data/' + complex_code + '/benchmark/' + complex_code + '_l_u.pdb', '~/data/' + complex_code + '/benchmark/')
			scp.put('../data/' + complex_code + '/benchmark/' + complex_code + '_r_u.pdb', '~/data/' + complex_code + '/benchmark/')
			scp.put("../data/{0}/patch_dock/{0}.patch_dock_output".format(complex_code), '~/data/' + complex_code + '/patch_dock/')
			copied_successfully.append(complex_code)

			with open('./copied_successfully.json', 'w') as f:
				json.dump(copied_successfully, f)
		except:
			print("ERROR with: %s" % complex_code)

	# run patch dock on each benchmark receptor and ligand which doesn't have an ouput file yet
	with open('./patch_dock_ran_successfully.json', 'r') as f:
		patch_dock_ran_successfully = json.load(f)

	if force:
		patch_dock_ran_successfully = []

	for complex_code in all_complex_codes:
		if len(complex_code) != 4 or complex_code in patch_dock_ran_successfully:
			continue

		print("Running PatchDock on: %s" % complex_code)

		try:
			ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command('cd ~/PatchDock; sh run_patch_dock_on_one_complex.sh ' + complex_code)
			print(ssh_stdout.read(), ssh_stderr.read())
			if ssh_stderr.read() == b'':			
				patch_dock_ran_successfully.append(complex_code)
				with open('./patch_dock_ran_successfully.json', 'w') as f:
					json.dump(patch_dock_ran_successfully, f)
		except Exception as e:
			print(e.message)
			print("ERROR with: %s" % complex_code)

	print("DONE")
	scp.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', '--username', help='universiy username')
    parser.add_argument('-p', '--password', help='universiy username\'s password')
    parser.add_argument('-t', '--test', help='perform on single complex to check flow is alright', default=False, action='store_true')
    parser.add_argument('-f', '--force', help='re uploads and runs everything', default=False, action='store_true')
    args = parser.parse_args()

