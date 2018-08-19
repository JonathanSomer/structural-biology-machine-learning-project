from Bio import pairwise2
from utils.pdb_utils import get_structure_sequence

def get_fnat_score(estimated_complex, banchmark_complex):
    return float(len(_intersect_neighbours_between_complexes(estimated_complex, banchmark_complex)) /
                 float(len(banchmark_complex.get_neighbours_residues())))

def _intersect_neighbours_between_complexes(c1, c2):
    neighbours_c1, neighbours_c2 = set(c1.get_neighbours_residues()), set(c2.get_neighbours_residues())
    receptor_pos_map = _get_position_map_between_structs(c1.receptor, c2.receptor)
    ligand_pos_map = _get_position_map_between_structs(c1.ligand, c2.ligand)
    joint_neighbors = [(neighbor_1_pos, neighbor_2_pos) for (neighbor_1_pos, neighbor_2_pos) in neighbours_c1
                       if (receptor_pos_map.get(neighbor_1_pos, -1), ligand_pos_map.get(neighbor_2_pos, -1))
                       in neighbours_c2]
    return joint_neighbors

def _get_position_map_between_structs(struct1, struct2):
    '''
    :return: a map from residues indices of struct1 to indices in struct2 corresponding to their pairwise alignment
    '''
    alignments = pairwise2.align.globalxx(
        get_structure_sequence(struct1),
        get_structure_sequence(struct2))
    return _get_position_map_from_alignment(*alignments[0])


def _get_position_map_from_alignment(align1, align2, score, begin, end):
    """format_alignment(align1, align2, score, begin, end) -> string
    Format the alignment prettily into a list of tuples of the identical residues.
    """
    index_map = {}
    i, align1_index, align2_index = 0, 0, 0
    while i < len(align1):
        if align1[i] != "-" and align2[i] != "-":
            index_map[align1_index] = align2_index
        if align1[i] != "-":
            align1_index += 1
        if align2[i] != "-":
            align2_index += 1
        i += 1
    return index_map

