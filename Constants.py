import os.path

project_dir = os.path.dirname(os.path.abspath(__file__))

class FilePaths(object):

    BASE_DATA_PATH = os.path.join(project_dir, r"data", r"{}-{}")
    PATCHDOCK_RESULTS_DIR_PATH_FORMAT = os.path.join(BASE_DATA_PATH, r"patch_dock_results")
    PATCHDOCK_RAW_RESULT_PATH_FORMAT = os.path.join(PATCHDOCK_RESULTS_DIR_PATH_FORMAT, r"raw_results", r"docking.res.{}.pdb")
    PATCHDOCK_RESULT_RECEPTOR_PATH_FORMAT = os.path.join(PATCHDOCK_RESULTS_DIR_PATH_FORMAT, r"splitted_results", r"docking.res.{}.receptor.pdb")
    PATCHDOCK_RESULT_LIGAND_PATH_FORMAT = os.path.join(PATCHDOCK_RESULTS_DIR_PATH_FORMAT, r"splitted_results", r"docking.res.{}.ligand.pdb")
    PDB_FILES_DIR_PATH_FORMAT = os.path.join(BASE_DATA_PATH, r"pdb_files")
    RECEPTOR_PDB_FILE_PATH = os.path.join(PDB_FILES_DIR_PATH_FORMAT, r"receptor.pdb")
    LIGAND_PDB_FILE_PATH = os.path.join(PDB_FILES_DIR_PATH_FORMAT, r"ligand.pdb")
    KNOWN_COMPLEX_PDB_FILE_PATH = os.path.join(PDB_FILES_DIR_PATH_FORMAT, r"complex.pdb")
    RAPTORX_MATRIX_FILE_PATH = os.path.join(BASE_DATA_PATH, r"raptorx_results", r"raptorx_matrix")

    @classmethod
    def get_patch_dock_raw_result_path(cls, receptor_pdb_id, ligand_pdb_id, rank):
        return cls.PATCHDOCK_RAW_RESULT_PATH_FORMAT.format(receptor_pdb_id, ligand_pdb_id, rank)

    @classmethod
    def get_patch_dock_result_receptor_path(cls, receptor_pdb_id, ligand_pdb_id, rank):
        return cls.PATCHDOCK_RESULT_RECEPTOR_PATH_FORMAT.format(receptor_pdb_id, ligand_pdb_id, rank)

    @classmethod
    def get_patch_dock_result_ligand_path(cls, receptor_pdb_id, ligand_pdb_id, rank):
        return cls.PATCHDOCK_RESULT_LIGAND_PATH_FORMAT.format(receptor_pdb_id, ligand_pdb_id, rank)

    @classmethod
    def get_receptor_pdb_path(cls, receptor_pdb_id, ligand_pdb_id):
        return cls.RECEPTOR_PDB_FILE_PATH.format(receptor_pdb_id, ligand_pdb_id)

    @classmethod
    def get_ligand_pdb_path(cls, receptor_pdb_id, ligand_pdb_id):
        return cls.LIGAND_PDB_FILE_PATH.format(receptor_pdb_id, ligand_pdb_id)

    @classmethod
    def get_known_comlex_pdb_path(cls, receptor_pdb_id, ligand_pdb_id):
        return cls.KNOWN_COMPLEX_PDB_FILE_PATH.format(receptor_pdb_id, ligand_pdb_id)

    @classmethod
    def get_raptorx_matrix_file_path(cls, receptor_id, ligand_id):
        return cls.RAPTORX_MATRIX_FILE_PATH.format(receptor_id, ligand_id)