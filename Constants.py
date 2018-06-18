import os.path

project_dir = os.path.dirname(os.path.abspath(__file__))

class FilePaths(object):

    PATCHDOCK_RESULTS_DIR_PATH_FORMAT = r"data\patch_dock_results\{}-{}"
    PATCHDOCK_RESULT_PATH_FORMAT = PATCHDOCK_RESULTS_DIR_PATH_FORMAT + r"\docking.res.{}.pdb"
    PATCHDOCK_RESULT_RECEPTOR_PATH_FORMAT = PATCHDOCK_RESULTS_DIR_PATH_FORMAT + "\docking.res.{}.receptor.pdb"
    PATCHDOCK_RESULT_LIGAND_PATH_FORMAT = PATCHDOCK_RESULTS_DIR_PATH_FORMAT + "\docking.res.{}.ligand.pdb"
    PDB_FILE_PATH_FORMAT = r"data\pdb_files\{}.pdb"
    RAPTORX_RESULTS_FILE_PATH_FORMAT = r"data\raptorx_results\{}-{}"

    @classmethod
    def get_patch_dock_result_path(cls, receptor_pdb_id, ligand_pdb_id, rank):
        return cls.PATCHDOCK_RESULT_PATH_FORMAT.format(receptor_pdb_id, ligand_pdb_id, rank)

    @classmethod
    def get_patch_dock_result_receptor_path(cls, receptor_pdb_id, ligand_pdb_id, rank):
        return cls.PATCHDOCK_RESULT_RECEPTOR_PATH_FORMAT.format(receptor_pdb_id, ligand_pdb_id, rank)

    @classmethod
    def get_patch_dock_result_ligand_path(cls, receptor_pdb_id, ligand_pdb_id, rank):
        return cls.PATCHDOCK_RESULT_LIGAND_PATH_FORMAT.format(receptor_pdb_id, ligand_pdb_id, rank)

    @classmethod
    def get_pdb_file_path(cls, pdb_id):
        return cls.PDB_FILE_PATH_FORMAT.format(pdb_id)

    @classmethod
    def get_raptorx_results_file_path(cls, receptor_id, ligand_id):
        return cls.RAPTORX_RESULTS_FILE_PATH_FORMAT.format(receptor_id, ligand_id)