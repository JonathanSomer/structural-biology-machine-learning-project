
class FilePaths(object):

    def __init__(self):
        self.patchDock_result_path_format = r"patchDock_results/{}-{}/docking.res.{}.pdb"
        self.pdb_file_path_format = r"pdb_files/{}.pdb"

    def __getattr__(self, item):
        return r"./" + item


file_paths = FilePaths()
