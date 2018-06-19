from Bio.PDB import NeighborSearch
from Bio.PDB import PPBuilder
from utils import pdb_utils


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
        # type: (self, int) -> set[tuple[int, int]]
        """
        Gets a set of neighboring receptor-ligand residues.
        :param radius: radius to search for neighbors (in Angstrom)
        :return: a set of tuples with the global ids of the receptor and ligand residues respectively
        """

        residue_neighbors = self._get_all_residue_neighbors(radius)
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

    def _get_all_residue_neighbors(self, radius):
        atoms = self.get_all_atoms()
        neighbors_search = NeighborSearch(atoms)
        return neighbors_search.search_all(radius, level='R')

    def get_all_atoms(self):
        return list(self.receptor_struct.get_atoms()) + list(self.ligand_struct.get_atoms())
