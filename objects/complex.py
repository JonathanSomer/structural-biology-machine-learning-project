from Bio.PDB import NeighborSearch
from Bio.PDB import PPBuilder
from utils import pdb_utils
from abc import ABCMeta, abstractmethod
from enum import Enum
from Constants import *
from utils.pdb_utils import pdb_parser


class ComplexType(Enum):
    zdock_benchmark_bound = 1
    zdock_benchmark_unbound = 2
    patch_dock = 3


class Complex(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def _init_complex(self):
        pass


class BenchmarkComplex(Complex):
    def __init__(self, complex_id, type=ComplexType.zdock_benchmark_bound):
        self.complex_id = complex_id
        self.type = type

        self.ligand, self.receptor = self._init_complex()

    def _init_complex(self):
        bound = self.type == ComplexType.zdock_benchmark_bound
        ligand_pdb_file_path = get_zdock_benchmark_pdb_path(self.complex_id, ligand=True, bound=bound)
        receptor_pdb_file_path = get_zdock_benchmark_pdb_path(self.complex_id, ligand=False, bound=bound)
        ligand = pdb_parser.get_structure_by_file_path(self.complex_id, ligand_pdb_file_path)
        receptor = pdb_parser.get_structure_by_file_path(self.complex_id, receptor_pdb_file_path)
        return ligand, receptor


class PatchDockComplex(Complex):
    def __init__(self, complex_id, rank, ligand_chains, receptor_chains):
        self.complex_id = complex_id
        self.ligand_chains, self.receptor_chains = ligand_chains, receptor_chains
        self.original_rank = rank
        self.type = ComplexType.patch_dock

        self.ligand, self.receptor = self._init_complex()
        self.capri_score = self.calculate_capri_score()
        self.raptor_score = self.calculate_raptor_score()

    def _init_complex(self):
        raise NotImplementedError("Should implement this method")

    def calculate_capri_score(self):
        raise NotImplementedError("Should implement this method")

    def calculate_raptor_score(self):
        raise NotImplementedError("Should implement this method")
