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

    @abstractmethod
    def get_neighbouring_residues(self):
        pass


class BenchmarkComplex(Complex):
    def __init__(self, complex_id, type=ComplexType.zdock_benchmark_bound):
        self._complex_id = complex_id
        self._type = type

        self._ligand, self._receptor = self._init_complex()

    def _init_complex(self):
        bound = self.type == ComplexType.zdock_benchmark_bound
        ligand_pdb_file_path = get_zdock_benchmark_pdb_path(self.complex_id, ligand=True, bound=bound)
        receptor_pdb_file_path = get_zdock_benchmark_pdb_path(self.complex_id, ligand=False, bound=bound)
        ligand = pdb_parser.get_structure(self.complex_id, ligand_pdb_file_path)
        receptor = pdb_parser.get_structure(self.complex_id, receptor_pdb_file_path)
        return ligand, receptor

    def get_neighbouring_residues(self):
        raise NotImplementedError("Should implement this method")


class PatchDockComplex(Complex):
    def __init__(self, complex_id, rank):
        self._complex_id = complex_id
        self._type = ComplexType.patch_dock

        self.ligand_chain_ids, self.receptor_chain_ids = self.get_chains()
        self.original_rank = rank
        self._ligand, self._receptor = self._init_complex()
        # self.capri_score = self.calculate_capri_score()
        # self.raptor_score = self.calculate_raptor_score()

    def get_chains(self):
        benchmark_complex = BenchmarkComplex(self._complex_id, type=ComplexType.zdock_benchmark_unbound)
        receptor_chain_ids = [chain.get_id() for chain in list(benchmark_complex.receptor.get_chains())]
        ligand_chain_ids = [chain.get_id() for chain in list(benchmark_complex.ligand.get_chains())]

        return ligand_chain_ids, receptor_chain_ids

    def _init_complex(self):
        patch_dock_complex_path = get_patchdock_ranked_complex_pdb_path(self._complex_id, self.original_rank)

        patch_dock_ligand = pdb_parser.get_structure(self._complex_id, patch_dock_complex_path)
        patch_dock_receptor = pdb_parser.get_structure(self._complex_id, patch_dock_complex_path)

        for chain_id in self.receptor_chain_ids:
            patch_dock_ligand.get_list()[0].detach_child(chain_id)

        for chain_id in self.ligand_chain_ids:
            patch_dock_receptor.get_list()[0].detach_child(chain_id)

        return patch_dock_ligand, patch_dock_receptor

    def get_neighbouring_residues(self):
        raise NotImplementedError("Should implement this method")

    def calculate_capri_score(self):
        raise NotImplementedError("Should implement this method")

    def calculate_raptor_score(self):
        raise NotImplementedError("Should implement this method")
