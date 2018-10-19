import json
import re
from abc import ABCMeta, abstractmethod

from Bio.PDB import PDBParser, NeighborSearch, PPBuilder, Residue
from enum import Enum

from Constants import *

LIGAND_STRUCT_ID = 'ligand'
RECEPTOR_STRUCT_ID = 'receptor'
N_PATCH_DOCK_SCORE_COMPONENTS = 4

pdb_parser = PDBParser(QUIET=True)


class ComplexType(Enum):
    zdock_benchmark_bound = 1
    zdock_benchmark_unbound = 2
    patch_dock = 3


class Complex(object):
    __metaclass__ = ABCMeta

    def __init__(self, complex_id, reprocess=False):
        self._complex_id = complex_id
        self._neighbours = None
        if (not self._is_processed()) or reprocess:
            self._process_complex()
        data = self._load_processed_data()
        self._receptor_sequence, self._ligand_sequence = data['r_seq'], data['l_seq']
        for k in data.keys():
            if k.startswith('nb'):
                radius = int(k[2:])
                setattr(self, NEIGHBOUR_ATTR_TEMPLATE % radius, data[k])

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

    def get_neighbours_residues(self, neighbor_radius=NEIGHBOUR_DEFAULT_RADIUS):
        # type: () -> List[Tuple[int, int]]
        """
        :return: list of tuples (receptor_residue_index, ligand_residue_index) in which the euclidean distance
                 between them is at most $(NEIGHBOURS_RADIUS)
        """

        def is_protein_residue(residue):
            # type: (Residue) -> bool
            """
            Checks if the residue comes from a protein
            :param residue: the given residue
            :return: is from protein
            """
            return residue.id[0] == ' '

        def add_true_residue_indexes(struct):
            """
            adds for each residue in struct add a 'true' index, i.e it's index when counting all residues
            starting from one in same order they appear in the pdb file.
            the indices are global among chains of same protein
            """
            true_index = 0
            for residue in struct.get_residues():
                if not is_protein_residue(residue):
                    continue
                residue.true_index = true_index
                true_index += 1

        if neighbor_radius not in NEIGHBOUR_RADII:
            raise ValueError("radius value is invalid")

        neighbour_attr = NEIGHBOUR_ATTR_TEMPLATE % neighbor_radius
        if hasattr(self, neighbour_attr) and getattr(self, neighbour_attr) is not None:
            return getattr(self, neighbour_attr)

        add_true_residue_indexes(self.receptor)
        add_true_residue_indexes(self.ligand)
        # TODO: change this hardcoded shit
        if neighbor_radius == 8:
            ligand_atoms = [res['CA'] for res in self.ligand.get_residues() if res.has_id('CA')]
        else:
            ligand_atoms = list(self.ligand.get_atoms())
        if neighbor_radius == 8:
            receptor_atoms = [res['CA'] for res in self.receptor.get_residues() if res.has_id('CA')]
        else:
            receptor_atoms = list(self.receptor.get_atoms())

        nb = NeighborSearch(ligand_atoms + receptor_atoms)
        all_neighbours = nb.search_all(neighbor_radius, level='R')

        neighbor_indexes = []

        for residue_neighbor_A, residue_neighbor_B in all_neighbours:
            neighbor_A_struct_id, neighbor_B_struct_id = residue_neighbor_A.get_full_id()[0], \
                                                         residue_neighbor_B.get_full_id()[0]

            if neighbor_A_struct_id == neighbor_B_struct_id:  # means both from receptor or both from ligand
                continue
            # one of the residues is not animo-acid
            if not is_protein_residue(residue_neighbor_A) or not is_protein_residue(residue_neighbor_B):
                continue

            if neighbor_A_struct_id == RECEPTOR_STRUCT_ID:
                receptor_residue, ligand_residue = residue_neighbor_A, residue_neighbor_B
            else:
                receptor_residue, ligand_residue = residue_neighbor_B, residue_neighbor_A

            # import pdb; pdb.set_trace()
            neighbor_indexes.append((receptor_residue.true_index, ligand_residue.true_index))

        setattr(self, neighbour_attr, neighbor_indexes)
        return getattr(self, neighbour_attr)

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
            init_fn()
        return getattr(self, attr_str)

    def _process_complex(self):
        data_path = self._get_processed_path()
        data = {
            "complex_id": self.complex_id,
            "l_seq": Complex.get_structure_sequence(self.ligand),
            "r_seq": Complex.get_structure_sequence(self.receptor),
            "nb5": self.get_neighbours_residues(5),
            "nb8": self.get_neighbours_residues(8)  # should not be used!
        }
        data_dir_path = os.path.dirname(data_path)
        if not os.path.exists(data_dir_path):
            os.makedirs(data_dir_path)
        with open(data_path, 'w') as f:
            json.dump(data, f, sort_keys=True, indent=4)

    def _load_processed_data(self):
        with open(self._get_processed_path(), 'r') as f:
            return json.load(f)

    def _is_processed(self):
        return os.path.isfile(self._get_processed_path())

    @abstractmethod
    def _get_processed_path(self):
        raise NotImplementedError("abs method")


class BenchmarkComplex(Complex):

    def __init__(self, complex_id, type=ComplexType.zdock_benchmark_bound, reprocess=False):
        self._type = type
        super(BenchmarkComplex, self).__init__(complex_id, reprocess)

    def _init_complex(self):
        bound = self.type == ComplexType.zdock_benchmark_bound
        ligand_pdb_file_path = get_zdock_benchmark_pdb_path(self.complex_id, ligand=True, bound=bound)
        receptor_pdb_file_path = get_zdock_benchmark_pdb_path(self.complex_id, ligand=False, bound=bound)
        ligand = pdb_parser.get_structure(LIGAND_STRUCT_ID, ligand_pdb_file_path)
        receptor = pdb_parser.get_structure(RECEPTOR_STRUCT_ID, receptor_pdb_file_path)
        self._ligand, self._receptor = ligand, receptor

    def _get_processed_path(self):
        bound = self.type == ComplexType.zdock_benchmark_bound
        return get_zdock_benchmark_processed_data_path(self.complex_id, bound)

"""
LEGACY CODE
class PatchDockComplex(Complex):

    def __init__(self, complex_id, rank, reprocess=False):
        self.original_rank = rank
        self._type = ComplexType.patch_dock
        super(PatchDockComplex, self).__init__(complex_id, reprocess)
        self.init_patch_dock_score_components()

    def _init_complex(self):
        # type: () -> None
        patch_dock_ligand_path = get_patchdock_ranked_complex_pdb_path(self._complex_id, self.original_rank,
                                                                       ligand=True)
        patch_dock_receptor_path = get_patchdock_ranked_complex_pdb_path(self._complex_id, self.original_rank,
                                                                         ligand=False)
        ligand = pdb_parser.get_structure(LIGAND_STRUCT_ID, patch_dock_ligand_path)
        receptor = pdb_parser.get_structure(RECEPTOR_STRUCT_ID, patch_dock_receptor_path)

        self._ligand, self._receptor = ligand, receptor

    def _get_processed_path(self):
        return get_patchdock_ranked_complex_processed_data_path(self.complex_id, self.original_rank)

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
"""
