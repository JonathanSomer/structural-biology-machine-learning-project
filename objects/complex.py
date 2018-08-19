from Bio.PDB import NeighborSearch
from abc import ABCMeta, abstractmethod
from enum import Enum
from Constants import *
from utils.pdb_utils import pdb_parser

NEIHGBOR_RADIUS = 5
LIGAND_STRUCT_ID = 'ligand'
RECEPTOR_STRUCT_ID = 'receptor'


class ComplexType(Enum):
    zdock_benchmark_bound = 1
    zdock_benchmark_unbound = 2
    patch_dock = 3


class Complex(object):
    __metaclass__ = ABCMeta

    def __init__(self):
        self._neighbours = None

    @property
    def complex_id(self):
        return self._complex_id

    @property
    def type(self):
        return self._type

    @property
    def receptor(self):
        return self._receptor

    @property
    def ligand(self):
        return self._ligand

    @abstractmethod
    def _init_complex(self):
        pass

    def _add_true_residue_indexes(self):
        '''
        adds for each residue (or ligand ot receptor) a 'true' index, i.e it's index when counting all residues
        starting from one in same order they appear in the pdb file.  the indices are global among chains of same protein
        '''
        for i, residue in enumerate(self.ligand.get_residues()):
            residue.true_index = i

        for i, residue in enumerate(self.receptor.get_residues()):
            residue.true_index = i

    def get_neighbours_residues(self):
        # type: () -> lst((int, int),...)
        '''
        :return: list of tuples (receptor_residue_index, ligand_residue_index) in which the euclidean distance
                 between them is at most $(NEIGHBOURS_RADIUS)
        '''
        if self._neighbours:
            return self._neighbours

        ligand_atoms = list(self.ligand.get_atoms())
        receptor_atoms = list(self.receptor.get_atoms())

        nb = NeighborSearch(ligand_atoms + receptor_atoms)
        all_neighbours = nb.search_all(NEIHGBOR_RADIUS, level='R')

        neighbor_indexes = []

        for residue_neighbor_A, residue_neighbor_B in all_neighbours:
            neighbor_A_struct_id, neighbor_B_struct_id = residue_neighbor_A.get_full_id()[0], \
                                                         residue_neighbor_B.get_full_id()[0]

            if neighbor_A_struct_id == neighbor_B_struct_id:  # means both from receptor or both from ligand
                continue

            if neighbor_A_struct_id == RECEPTOR_STRUCT_ID:
                receptor_residue, ligand_residue = residue_neighbor_A, residue_neighbor_B
            else:
                receptor_residue, ligand_residue = residue_neighbor_B, residue_neighbor_A

            neighbor_indexes.append((receptor_residue.true_index, ligand_residue.true_index))

        self._neighbours = neighbor_indexes
        return self._neighbours


class BenchmarkComplex(Complex):

    def __init__(self, complex_id, type=ComplexType.zdock_benchmark_bound):
        super(BenchmarkComplex, self).__init__()
        self._complex_id = complex_id
        self._type = type
        self._ligand, self._receptor = self._init_complex()
        self._add_true_residue_indexes()
        self._neighbours = None

    def _init_complex(self):
        bound = self.type == ComplexType.zdock_benchmark_bound
        ligand_pdb_file_path = get_zdock_benchmark_pdb_path(self.complex_id, ligand=True, bound=bound)
        receptor_pdb_file_path = get_zdock_benchmark_pdb_path(self.complex_id, ligand=False, bound=bound)
        ligand = pdb_parser.get_structure(LIGAND_STRUCT_ID, ligand_pdb_file_path)
        receptor = pdb_parser.get_structure(RECEPTOR_STRUCT_ID, receptor_pdb_file_path)
        return ligand, receptor


class PatchDockComplex(Complex):

    def __init__(self, complex_id, rank):
        super(PatchDockComplex, self).__init__()
        self._complex_id = complex_id
        self._type = ComplexType.patch_dock
        self.original_rank = rank
        self._ligand_chain_ids, self._receptor_chain_ids = self._infer_r_l_chain_ids_from_banchmark_complex()
        self._ligand, self._receptor = self._init_complex()
        self._add_true_residue_indexes()

    def _infer_r_l_chain_ids_from_banchmark_complex(self):
        benchmark_complex = BenchmarkComplex(self._complex_id, type=ComplexType.zdock_benchmark_unbound)
        receptor_chain_ids = [chain.get_id() for chain in list(benchmark_complex.receptor.get_chains())]
        ligand_chain_ids = [chain.get_id() for chain in list(benchmark_complex.ligand.get_chains())]

        return ligand_chain_ids, receptor_chain_ids

    def _init_complex(self):
        patch_dock_complex_path = get_patchdock_ranked_complex_pdb_path(self._complex_id, self.original_rank)
        ligand = pdb_parser.get_structure(LIGAND_STRUCT_ID, patch_dock_complex_path)
        receptor = pdb_parser.get_structure(RECEPTOR_STRUCT_ID, patch_dock_complex_path)

        for chain_id in self._receptor_chain_ids:
            self._remove_chain_from_struct(ligand, chain_id)

        for chain_id in self._ligand_chain_ids:
            self._remove_chain_from_struct(receptor, chain_id)

        return ligand, receptor

    def _remove_chain_from_struct(self, struct, chain_id):
        first_model = struct.get_list()[0]  # struct children are models
        first_model.detach_child(chain_id)  # model children are chains
