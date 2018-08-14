from Bio.PDB import NeighborSearch
from Bio.PDB import PPBuilder
from utils import pdb_utils
from abc import ABCMeta, abstractmethod
from enum import Enum
from Constants import *
from utils.pdb_utils import pdb_parser


NEIHGBOR_RADIUS = 5

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

    def get_neighbouring_residues(self):
        for i, residue in enumerate(self.ligand.get_residues()):
            residue.true_index = i

        for i, residue in enumerate(self.receptor.get_residues()):
            residue.true_index = i
        
        ligand_atoms = list(self.ligand.get_atoms())
        receptor_atoms = list(self.receptor.get_atoms())

        nb = NeighborSearch(ligand_atoms + receptor_atoms)
        search_results = nb.search_all(NEIHGBOR_RADIUS, level='R')
        
        neighbor_indexes = []

        for residue_1, residue_2 in search_results:
            type_1, type_2 = residue_1.get_full_id()[0], residue_2.get_full_id()[0]
            
            if type_1 == type_2:
                continue

            if type_1 == 'ligand':
                receptor_residue, ligand_residue = residue_2, residue_1
            else:
                receptor_residue, ligand_residue = residue_1, residue_2

            neighbor_indexes.append((receptor_residue.true_index, ligand_residue.true_index))

        return neighbor_indexes

class BenchmarkComplex(Complex):
    def __init__(self, complex_id, type=ComplexType.zdock_benchmark_bound):
        self._complex_id = complex_id
        self._type = type

        self._ligand, self._receptor = self._init_complex()

    def _init_complex(self):
        bound = self.type == ComplexType.zdock_benchmark_bound
        ligand_pdb_file_path = get_zdock_benchmark_pdb_path(self.complex_id, ligand=True, bound=bound)
        receptor_pdb_file_path = get_zdock_benchmark_pdb_path(self.complex_id, ligand=False, bound=bound)
        ligand = pdb_parser.get_structure('ligand', ligand_pdb_file_path)
        receptor = pdb_parser.get_structure('receptor', receptor_pdb_file_path)
        return ligand, receptor


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

        patch_dock_ligand = pdb_parser.get_structure('ligand', patch_dock_complex_path)
        patch_dock_receptor = pdb_parser.get_structure('receptor', patch_dock_complex_path)

        for chain_id in self.receptor_chain_ids:
            patch_dock_ligand.get_list()[0].detach_child(chain_id)

        for chain_id in self.ligand_chain_ids:
            patch_dock_receptor.get_list()[0].detach_child(chain_id)

        return patch_dock_ligand, patch_dock_receptor

    def calculate_capri_score(self):
        raise NotImplementedError("Should implement this method")

    def calculate_raptor_score(self):
        raise NotImplementedError("Should implement this method")
