import os
'''
you can run this on any machine
this script adds headers for benchmark files that makes patchdock results separable between ligand and receptor
'''


LIGAND_SEPARATOR_HEADER = "LIGAND SEPARATOR\n"
LIGAND_HEADER = "LIGAND\n"
RECEPTOR_HEADER = "RECEPTOR\n"

data_path = r"../proccessed_data_v1"

all_complex_codes = next(os.walk(data_path))[1]
for c_id in all_complex_codes:

    with open(data_path + "/{0}/benchmark/{0}_l_u.pdb".format(c_id), 'r+') as l_file:
        f_lines = l_file.readlines()
        if f_lines[0] != LIGAND_SEPARATOR_HEADER:
            header = LIGAND_SEPARATOR_HEADER + LIGAND_HEADER
            l_file.seek(0)
            l_file.write(header + "".join(f_lines))

    with open(data_path + "/{0}/benchmark/{0}_r_u.pdb".format(c_id), 'r+') as r_file:
        f_lines = r_file.readlines()
        if f_lines[0] != RECEPTOR_HEADER:
            r_file.seek(0)
            r_file.write(RECEPTOR_HEADER + "".join(f_lines))