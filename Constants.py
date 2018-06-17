import os.path

project_dir = os.path.dirname(os.path.abspath(__file__))

class FilePaths(object):

    def __init__(self):

        self.patchDock_results_path_format = os.path.join(project_dir, r"patchDock_results\{}-{}")
        self.patchDock_result_path_format = os.path.join(project_dir, r"patchDock_results\{}-{}\docking.res.{}.pdb")
        self.patchDock_result_receptor_path_format = self.patchDock_results_path_format + "\docking.res.{}.receptor.pdb"
        self.patchDock_result_ligand_path_format = self.patchDock_results_path_format + "\docking.res.{}.ligand.pdb"
        self.pdb_file_path_format = r"pdb_files\{}.pdb"



file_paths = FilePaths()

