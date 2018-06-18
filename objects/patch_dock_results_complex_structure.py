from Bio.PDB import NeighborSearch
from Bio.PDB import PPBuilder

class StructureWrapper():

    def __init__(self, struct):
        self._struct = struct

    def get_sequence(self):
        ppb = PPBuilder()
        return ppb.build_peptides(self._struct)[0].get_sequence()

class PatchDockResultComplexStructure(object):

    def __init__(self, receptor_struct, ligand_struct):
        self.receptor_struct = receptor_struct
        self.ligand_struct = ligand_struct

    def get_aa_neighbors(self, radius=8):
        all_neighbors = self._get_all_neighbors(radius)
        return self._get_all_receptor_ligand_neighbors_col(all_neighbors)

    def _get_all_neighbors(self, radius):
        all_atoms = self.get_all_atoms()
        neighbors_search = NeighborSearch(all_atoms)
        return neighbors_search.search_all(radius, level='R')

    def get_all_atoms(self):
        return [atom for atom in self.receptor_struct.get_atoms()] + \
               [atom for atom in self.ligand_struct.get_atoms()]

    def _get_all_receptor_ligand_neighbors_col(self, all_neighbors):
        receptor_ligand_neighbors_col = set()
        for neighborA, neighborB in all_neighbors:
            if self._is_neighbors_are_from_different_protein(neighborA, neighborB):
                receptor_col, ligand_col = self._get_receptor_ligand_col(neighborA, neighborB)
                receptor_ligand_neighbors_col.add((receptor_col, ligand_col))
        return receptor_ligand_neighbors_col

    def _is_neighbors_are_from_different_protein(self, neighborA, neighborB):
        neighborA_prot_id, neighborB_prot_id = neighborA.full_id[0], neighborB.full_id[0] #full_id[0] == protein id
        is_neighborA_is_prot_res = neighborA.id[0] == ' ' #id[0] == ' ' iff atom is from protein, else from glucose, water..
        is_neighborB_is_prot_res = neighborB.id[0] == ' '
        if neighborA_prot_id != neighborB_prot_id and is_neighborA_is_prot_res and is_neighborB_is_prot_res:
            return True
        else:
            return False

    def _get_receptor_ligand_col(self, neighborA, neighborB):
        receptor_col, ligand_col = neighborA.id[1], neighborB.id[1]
        if neighborA.full_id[0] == self.ligand_struct.id:
            receptor_col, ligand_col = ligand_col, receptor_col
        return receptor_col, ligand_col

