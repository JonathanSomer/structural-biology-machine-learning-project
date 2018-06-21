from Bio.PDB import NeighborSearch
from Bio.PDB import PPBuilder
from utils import pdb_utils
from abc import ABCMeta

class AbstractComplex(object):
    __metaclass__ = ABCMeta

    def get_all_atoms(self):
        raise NotImplementedError("Abstract method")

    def get_residues(self):
        raise NotImplementedError("Abstract method")

    def get_all_residue_atoms(self):
        all_atoms = self.get_all_atoms()
        return [atom for atom in all_atoms if pdb_utils.is_protein_residue(atom.get_parent())]


class KnownComplex(AbstractComplex):

    def __init__(self, struct):
        self._struct = struct

    def get_all_atoms(self):
        return list(self._struct.get_atoms())

    def get_residues(self):
        for residue in self._struct.get_residues():
            if pdb_utils.is_protein_residue(residue):
                yield residue

class PatchDockResultComplexStructure(AbstractComplex):

    def __init__(self, receptor_struct, ligand_struct):
        self.receptor_struct = receptor_struct
        self.ligand_struct = ligand_struct

    def get_aa_neighbors(self, radius=8):
        # type: (self, int) -> set[tuple[int, int]]
        """
        Gets a set of neighboring receptor-ligand residues.
        :param radius: radius to search for neighbors (in Angstrom)
        :return: a set of tuples with the global ids of the receptor and ligand residues respectively
        """
        residue_neighbors = self.get_all_residue_neighbors(radius)
        aa_neighbors = set()
        for neighborA, neighborB in residue_neighbors:
            neighborA_prot_id, neighborB_prot_id = neighborA.full_id[0], neighborB.full_id[0]
            # skip neighbors if they come from the same protein or at least one of them is not a residue of a protein
            if (neighborA_prot_id == neighborB_prot_id
                    or not pdb_utils.is_protein_residue(neighborA)
                    or not pdb_utils.is_protein_residue(neighborB)):
                continue
            # swap neighbors if the first is from the ligand
            if neighborA_prot_id == self.ligand_struct.id:
                neighborA, neighborB = neighborB, neighborA
            aa_neighbors.add((pdb_utils.get_residue_global_id(neighborA), pdb_utils.get_residue_global_id(neighborB)))
        return aa_neighbors

    def get_all_residue_neighbors(self, radius):
        atoms = self.get_all_atoms()
        neighbors_search = NeighborSearch(atoms)
        return neighbors_search.search_all(radius, level='R')

    def get_all_atoms(self):
        return list(self.receptor_struct.get_atoms()) + list(self.ligand_struct.get_atoms())

    def get_residues(self):
        for receptor_residue in self.receptor_struct.get_residues():
            if pdb_utils.is_protein_residue(receptor_residue):
                yield receptor_residue
        for ligand_residue in self.ligand_struct.get_residues():
            if pdb_utils.is_protein_residue(ligand_residue):
                yield ligand_residue

     #todo
    def distance_from_known_complex(self):
        None
