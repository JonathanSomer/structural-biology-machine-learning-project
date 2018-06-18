from Constants import FilePaths



def split_patchdock_result_complex(receptor_pdb_id, ligand_pdb_id, rank):
    org_file_path = FilePaths.get_patch_dock_result_path(receptor_pdb_id, ligand_pdb_id, rank)
    receptor_file_path = FilePaths.get_patch_dock_result_receptor_path(receptor_pdb_id, ligand_pdb_id, rank)
    ligand_file_path = FilePaths.get_patch_dock_result_ligand_path(receptor_pdb_id, ligand_pdb_id, rank)
    with open(org_file_path, 'r') as org_file:
        with open(receptor_file_path, 'w+') as receptor_file:
            header = org_file.readline()
            while("HEADER" not in header):
                header = org_file.readline()
            receptor_file.write(header)
            line = org_file.readline()
            while("HEADER" not in line):
                receptor_file.write(line)
                line = org_file.readline()
        with open(ligand_file_path, 'w+') as ligand_file:
            while(line != ''):
                ligand_file.write(line)
                line = org_file.readline()

def split_all_patchdock_results(receptor_pdb_id, ligand_pdb_id, num_of_results):
    for i in range(1, num_of_results+1):
        split_patchdock_result_complex(receptor_pdb_id, ligand_pdb_id,i)

split_patchdock_result_complex("1bra", "1aap", )