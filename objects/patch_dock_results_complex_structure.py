from Bio.PDB import NeighborSearch
from collections import defaultdict
from Bio.PDB import Structure

class StructureWrapper():

    def __init__(self, struct):
        self._struct = struct

class PatchDockResultComplexStructure(object):

    def __init__(self, receptor_struct, ligand_struct):
        self.receptor_struct = receptor_struct
        self.ligand_struct = ligand_struct

    def get_aa_neighbors(self, radius=8):
        all_neighbors = self._get_all_neighbors(radius)

    def _get_all_neighbors(self, radius):
        all_atoms = self.get_all_atoms()
        neighbors_search = NeighborSearch(all_atoms)
        return neighbors_search.search_all(radius, level='R')

    def get_all_atoms(self):
        return [atom for atom in self.receptor_struct.get_atoms()] + \
               [atom for atom in self.ligand_struct.get_atoms()]

    def _get_receptor_ligand_neighbors(self, all_neighbors):
        receptor_ligand_neighbors = set()
        for neighborA, neighborB in all_neighbors:
            None


