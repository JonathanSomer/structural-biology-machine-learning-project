from Bio.PDB import PDBParser
from Constants import FilePaths
from objects.complex import *

class PDBParserWrapper():

    def __init__(self):
        self._parser = PDBParser()

    def get_structure_by_file_path(self, id, pdb_path):
        return self._parser.get_structure(id, pdb_path)

    def get_receptor_structure(self, id, receptor_pdb_id, ligand_pdb_id):
        path = FilePaths.get_receptor_pdb_path(receptor_pdb_id, ligand_pdb_id)
        return self.get_structure_by_file_path(id, path)

    def get_ligand_structure(self, id, receptor_pdb_id, ligand_pdb_id):
        path = FilePaths.get_ligand_pdb_path(receptor_pdb_id, ligand_pdb_id)
        return self.get_structure_by_file_path(id, path)

    def get_known_complex_structure(self, receptor_pdb_id, ligand_pdb_id):
        path = FilePaths.get_known_comlex_pdb_path(receptor_pdb_id, ligand_pdb_id)
        return KnownComplex(self.get_structure_by_file_path("known", path))

    def get_patch_dock_result_complex_structure(self, receptor_pdb_id, ligand_pdb_id, rank):
        receptor_path = FilePaths.get_patch_dock_result_receptor_path(receptor_pdb_id, ligand_pdb_id, rank)
        receptor_struct = self.get_structure_by_file_path(receptor_pdb_id, receptor_path)
        ligand_path = FilePaths.get_patch_dock_result_ligand_path(receptor_pdb_id, ligand_pdb_id, rank)
        ligand_struct = self.get_structure_by_file_path(ligand_pdb_id, ligand_path)
        return PatchDockResultComplexStructure(receptor_struct, ligand_struct)
