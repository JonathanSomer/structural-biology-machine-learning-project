from Bio.PDB import NeighborSearch
from abc import ABCMeta, abstractmethod
from enum import Enum
from Constants import *
from utils.pdb_utils import pdb_parser
import re
import json

import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

warnings.simplefilter('ignore', PDBConstructionWarning)

NEIHGBOR_RADIUS = 5
LIGAND_STRUCT_ID = 'ligand'
RECEPTOR_STRUCT_ID = 'receptor'


class ComplexType(Enum):
    zdock_benchmark_bound = 1
    zdock_benchmark_unbound = 2
    patch_dock = 3


class Complex(object):
    __metaclass__ = ABCMeta

    def __init__(self, complex_id):
        self._complex_id = complex_id
        if not self._is_cached():
            self._cache_complex()
        cache = self._get_complex_cache()
        self._receptor_sequence, self._ligand_sequence = cache['r_seq'], cache['l_seq']
        self._neighbours = [tuple(nb) for nb in cache['nb5']]

    @property
    def complex_id(self):
        return self._complex_id

    @property
    def type(self):
        return self._type

    @property
    def receptor(self):
        return self._lazy_init('_receptor', self._init_complex)

    @property
    def ligand(self):
        return self._lazy_init('_ligand', self._init_complex)

    @property
    def receptor_sequence(self):
        return self._receptor_sequence

    @property
    def ligand_sequence(self):
        return self._ligand_sequence

    @staticmethod
    def get_structure_sequence(struct):
        # type: (Structure) -> str
        """
        Gets the structure sequence using PPBuilder
        :param struct: Structure object
        :return: struct sequence
        """
        ppb = PPBuilder()
        return ''.join([str(pp.get_sequence()) for pp in ppb.build_peptides(struct)])

    @abstractmethod
    def _init_complex(self):
        pass

    def _lazy_init(self, attr_str, init_fn):
        if not hasattr(self, attr_str) or getattr(self, attr_str, None) is None:
            init_fn(self)
        return getattr(self, attr_str)

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

    def _cache_complex(self):
        c_path = self._get_cache_path()
        cache = {
            "complex_id": self.complex_id,
            "l_seq": Complex.get_structure_sequence(self.ligand),
            "r_seq": Complex.get_structure_sequence(self.receptor),
            "nb5": self.get_neighbours_residues()
        }
        with open(c_path, 'w') as f:
            json.dump(cache, f)

    def _get_complex_cache(self):
        with open(self._get_cache_path(), 'r') as f:
            return json.load(f)

    def _is_cached(self):
        return os.path.isfile(self._get_cache_path())

    @abstractmethod
    def _get_cache_path(self):
        raise NotImplementedError("abs method")

class BenchmarkComplex(Complex):

    def __init__(self, complex_id, type=ComplexType.zdock_benchmark_bound):
        super(BenchmarkComplex, self).__init__(complex_id)
        self._type = type
        self._add_true_residue_indexes()

    def _init_complex(self):
        bound = self.type == ComplexType.zdock_benchmark_bound
        ligand_pdb_file_path = get_zdock_benchmark_pdb_path(self.complex_id, ligand=True, bound=bound)
        receptor_pdb_file_path = get_zdock_benchmark_pdb_path(self.complex_id, ligand=False, bound=bound)
        ligand = pdb_parser.get_structure(LIGAND_STRUCT_ID, ligand_pdb_file_path)
        receptor = pdb_parser.get_structure(RECEPTOR_STRUCT_ID, receptor_pdb_file_path)
        self._ligand, self._receptor = ligand, receptor

    def _get_cache_path(self):
        bound = self.type == ComplexType.zdock_benchmark_bound
        return get_zdock_benchmark_cache_path(self.complex_id, bound)


class PatchDockComplex(Complex):

    def __init__(self, complex_id, rank):
        super(PatchDockComplex, self).__init__(complex_id)
        self._type = ComplexType.patch_dock
        self.original_rank = rank
<<<<<<< HEAD
        self._ligand, self._receptor = self._init_complex()
        self.add_true_residue_indexes()
        self.init_patch_dock_score_components()
        # self.capri_score = self.calculate_capri_score()
        # self.raptor_score = self.calculate_raptor_score()

    def init_patch_dock_score_components(self):
        with open(get_patchdock_complex_score_file_path(self.complex_id), "r") as f:
            for line in f:
                pattern = re.compile("^\s+" + str(self.original_rank) + "\s\|.+")
                if pattern.match(line):
                    components = [x.strip() for x in line.split('|')]
                    # rank | score | pen.  | Area    | as1   | as2   | as12  | ACE     | hydroph | Energy  |cluster| dist. || Ligand Transformation
                    indexes = [1, 2, 3, 7]
                    self._score_components = [float(components[i]) for i in indexes]

    @property
    def score_components(self):
        return self._score_components

    def get_chains(self):
        benchmark_complex = BenchmarkComplex(self._complex_id, type=ComplexType.zdock_benchmark_unbound)
        receptor_chain_ids = [chain.get_id() for chain in list(benchmark_complex.receptor.get_chains())]
        ligand_chain_ids = [chain.get_id() for chain in list(benchmark_complex.ligand.get_chains())]

        return ligand_chain_ids, receptor_chain_ids
=======
        self._add_true_residue_indexes()
>>>>>>> 936cbacce422d94d78ac34973cbd20ba2ecaa4fc

    def _init_complex(self):
        # type: () -> None
        def _infer_r_l_chain_ids_from_benchmark_complex(self):
            benchmark_complex = BenchmarkComplex(self._complex_id, type=ComplexType.zdock_benchmark_unbound)
            receptor_chain_ids = [chain.get_id() for chain in list(benchmark_complex.receptor.get_chains())]
            ligand_chain_ids = [chain.get_id() for chain in list(benchmark_complex.ligand.get_chains())]
            return ligand_chain_ids, receptor_chain_ids

        def _remove_chain_from_struct(struct, chain_id):
            first_model = struct.get_list()[0]  # struct children are models
            first_model.detach_child(chain_id)  # model children are chains

        _ligand_chain_ids, _receptor_chain_ids = _infer_r_l_chain_ids_from_benchmark_complex(self)
        patch_dock_complex_path = get_patchdock_ranked_complex_pdb_path(self._complex_id, self.original_rank)
        ligand = pdb_parser.get_structure(LIGAND_STRUCT_ID, patch_dock_complex_path)
        receptor = pdb_parser.get_structure(RECEPTOR_STRUCT_ID, patch_dock_complex_path)

        for chain_id in self._receptor_chain_ids:
            _remove_chain_from_struct(ligand, chain_id)

        for chain_id in self._ligand_chain_ids:
            _remove_chain_from_struct(receptor, chain_id)

        self._ligand, self._receptor = ligand, receptor

    def _get_cache_path(self):
        return get_patchdock_ranked_complex_cache_path(self.complex_id, self.original_rank)
