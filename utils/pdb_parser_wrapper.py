# comment out old code
'''
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
'''