from Bio.PDB import PDBParser
from Constants import file_paths

class PDBParserWrapper():

    def __init__(self):
        self._parser = PDBParser()

    def get_structure(self, id, file):
        return self._parser.get_structure(id, file)

    def get_result_complex_structure(self, receptor_pdb_id, ligand_pdb_id, rank):
        receptor_path = file_paths.patchDock_result_receptor_path_format.format(receptor_pdb_id, ligand_pdb_id, rank)
        receptor_struct = self.get_structure(receptor_pdb_id, receptor_path)
        ligand_path = file_paths.patchDock_result_ligand_path_format.format(receptor_pdb_id, ligand_pdb_id, rank)
        ligand_struct = self.get_structure(ligand_pdb_id, ligand_path)
        return PatchDockResultComplexStructure(receptor_struct, ligand_struct)

class PatchDockResultComplexStructure(object):

    def __init__(self, receptor_struct, ligand_struct):
        self.receptor_struct = receptor_struct
        self.ligand_struct = ligand_struct

    def get_atoms(self):
        for atom in self.receptor_struct:
            yield atom
        for atom in self.ligand_struct:
            yield atom
