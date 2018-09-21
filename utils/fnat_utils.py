from Bio import pairwise2
from objects.complex import *

seq_maps_cache = {}

def get_fnat_score(estimated_complex, benchmark_complex, n_radius=5):
    numerator = len(_intersect_neighbours_between_complexes(estimated_complex, benchmark_complex, n_radius))
    denomenator = float(len(benchmark_complex.get_neighbours_residues(n_radius)))
    return float(numerator / denomenator)


def _intersect_neighbours_between_complexes(c1, c2,  n_radius=5):
    residues_c1 = [tuple(pos) for pos in c1.get_neighbours_residues(n_radius)]
    residues_c2 = [tuple(pos) for pos in c2.get_neighbours_residues(n_radius)]
    neighbours_c1, neighbours_c2 = set(residues_c1), set(residues_c2)

    receptor_pos_map = get_position_map_between_sequences(c1.receptor_sequence, c2.receptor_sequence)
    ligand_pos_map = get_position_map_between_sequences(c1.ligand_sequence, c2.ligand_sequence)

    # map c1 to c2 alignment
    neighbours_c1_mapped = set()
    for pos1, pos2 in neighbours_c1:
        if pos1 in receptor_pos_map and pos2 in ligand_pos_map:
            new_pos = (receptor_pos_map[pos1], ligand_pos_map[pos2])
            neighbours_c1_mapped.add(new_pos)

    # intersection of mapped c1 and c2 neighbours
    return neighbours_c1_mapped.intersection(neighbours_c2)


def get_position_map_between_sequences(seq1, seq2):
    """
    :return: a map from residue indices of seq1 to indices in seq2 corresponding to their pairwise alignment
    """
    if (seq1, seq2) not in seq_maps_cache:
        alignments = pairwise2.align.globalxx(seq1, seq2)
        from Bio.pairwise2 import format_alignment
        seq_maps_cache[(seq1, seq2)] = _get_position_map_from_alignment(*alignments[0])
    return seq_maps_cache[(seq1, seq2)]


def _get_position_map_from_alignment(align1, align2, score, begin, end):
    """
    format_alignment(align1, align2, score, begin, end) -> string
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
    neighbours = [(receptor_map[receptor_id], ligand_map[ligand_id]) for receptor_id, ligand_id in
                  bound_complex.get_neighbours_residues() if
                  receptor_id in receptor_map and ligand_id in ligand_map]
    return tuple(zip(*neighbours))
