from Bio.PDB import PDBParser
from Constants import FilePaths
from objects.patch_dock_results_complex_structure import *

class PDBParserWrapper():

    def __init__(self):
        self._parser = PDBParser()

    def get_structure(self, id, pdb_id):
        path = FilePaths.get_pdb_file_path(pdb_id)
        return StructureWrapper(self._parser.get_structure(id, path))

    def get_result_complex_structure(self, receptor_pdb_id, ligand_pdb_id, rank):
        receptor_path = FilePaths.get_patch_dock_result_receptor_path(receptor_pdb_id, ligand_pdb_id, rank)
        receptor_struct = self.get_structure(receptor_pdb_id, receptor_path)
        ligand_path = FilePaths.get_patch_dock_result_ligand_path(receptor_pdb_id, ligand_pdb_id, rank)
        ligand_struct = self.get_structure(ligand_pdb_id, ligand_path)
        return PatchDockResultComplexStructure(receptor_struct, ligand_struct)
