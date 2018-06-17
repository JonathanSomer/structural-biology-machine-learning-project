from Bio.PDB import PDBParser
from Constants import file_paths

class PDBParserWrapper():

    def __init__(self):
        self._parser = PDBParser()

class PatchDockResultComplexStructure(object):

    def __init__(self, receptor_struct, ligand_struct):
        self.receptor_struct = receptor_struct
        self.ligand_struct = ligand_struct

class PatchDockResultsParser(PDBParserWrapper):

    def get_structure(self, receptor_pdb_id, ligand_pdb_id, rank):

        pdb_path = file_paths.patchDock_result_path_format.format(receptor_pdb_id, ligand_pdb_id, rank)
        pdb_id = "{}-{}-PatchDockResult{}".format(receptor_pdb_id, ligand_pdb_id, rank)
        return self._parser.get_structure(pdb_id, pdb_path)

class KnownComplexParser(PDBParserWrapper):

    def get_structure(self, complex_pdb_id):
        pdb_path = file_paths.pdb_file_path_format(complex_pdb_id)
        return self._parser.get_structure(complex_pdb_id, pdb_path)
