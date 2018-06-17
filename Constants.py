


class FilePaths(object):

    def __init__(self):
        self.patchDock_results_path_format = r"patchDock_results/{}-{}"
        self.patchDock_result_path_format = r"patchDock_results/{}-{}/docking.res.{}.pdb"
        self.patchDock_result_receptor_path_format = self.patchDock_results_path_format + "/docking.res.{}.receptor.pdb"
        self.patchDock_result_ligand_path_format = self.patchDock_results_path_format + "/docking.res.{}.ligand.pdb"
        self.pdb_file_path_format = r"pdb_files/{}.pdb"


    def __getattr__(self, item):
        return r"./" + item


file_paths = FilePaths()

