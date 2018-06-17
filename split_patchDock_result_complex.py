from Constants import file_paths


def split_patchdock_result_complex(receptor_pdb_id, ligand_pdb_id, rank):
    org_file_path = file_paths.patchDock_result_path_format.format(receptor_pdb_id, ligand_pdb_id, rank)
    receptor_file_path = file_paths.patchDock_result_receptor_path_format.format(receptor_pdb_id, ligand_pdb_id, rank)
    ligand_file_path = file_paths.patchDock_result_ligand_path_format.format(receptor_pdb_id, ligand_pdb_id, rank)
    with open(org_file_path, 'r') as org_file:
        with open(receptor_file_path, 'w+') as receptor_file:
            line = org_file.readline()
            while("HEADER" not in line):
                receptor_file.write(line)
        with open(ligand_file_path, 'w+') as ligand_file:
            while(line != ''):
                ligand_file.write(line)
                line = org_file.readline()


split_patchdock_result_complex("1bra", "1app", 1)