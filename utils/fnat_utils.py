from Bio import pairwise2
from objects.complex import *

def get_fnat_score(estimated_complex, banchmark_complex):
    return float(len(_intersect_neighbours_between_complexes(estimated_complex, banchmark_complex)) /
                 float(len(banchmark_complex.get_neighbours_residues())))

def _intersect_neighbours_between_complexes(c1, c2):
    neighbours_c1, neighbours_c2 = set(c1.get_neighbours_residues()), set(c2.get_neighbours_residues())
    receptor_pos_map = get_position_map_between_sequences(c1.receptor_sequence, c2.receptor_sequence)
    ligand_pos_map = get_position_map_between_sequences(c1.ligand_sequence, c2.ligand_sequence)
    joint_neighbors = [(neighbor_1_pos, neighbor_2_pos) for (neighbor_1_pos, neighbor_2_pos) in neighbours_c1
                       if (receptor_pos_map.get(neighbor_1_pos, -1), ligand_pos_map.get(neighbor_2_pos, -1))
                       in neighbours_c2]
    return joint_neighbors

def get_position_map_between_sequences(seq1, seq2):
    '''
    :return: a map from residues indices of seq1 to indices in seq2 corresponding to their pairwise alignment
    '''
    alignments = pairwise2.align.globalxx(seq1, seq2)
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

def get_neighbours_for_bound_complex(bound_complex):
    unbound_complex = BenchmarkComplex(bound_complex.complex_id, type=ComplexType.zdock_benchmark_unbound)

    receptor_map = get_position_map_between_sequences(bound_complex.receptor_sequence,
                                                                 unbound_complex.receptor_sequence)
    ligand_map = get_position_map_between_sequences(bound_complex.ligand_sequence,
                                                               unbound_complex.ligand_sequence)
    # change from bound indices to unbound indices
    neighbours = [(receptor_map[receptor_id], ligand_map[ligand_id]) for receptor_id, ligand_id in bound_complex.get_neighbours_residues() if
                  receptor_id in receptor_map and ligand_id in ligand_map]
    return tuple(zip(*neighbours))
