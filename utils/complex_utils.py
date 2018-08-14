from Bio import pairwise2

from utils.pdb_utils import get_structure_sequence
from objects import complex

def intersect_neighbours_between_complexes(c1, c2):
    neighbours_c1, neighbours_c2 = set(c1.get_neighbouring_residues()), set(c2.get_neighbouring_residues())
    print(neighbours_c1)
    print(neighbours_c2)
    receptor_map = get_joint_positions_from_structs(c1.receptor, c2.receptor)
    ligand_map = get_joint_positions_from_structs(c1.ligand, c2.ligand)
    joint_neighbors = [(i,j) for (i,j) in neighbours_c1
                       if (receptor_map.get(i, -1), ligand_map.get(j, -1)) in neighbours_c2]

    print(neighbours_c1.intersection(neighbours_c2))
    print(joint_neighbors)
    return joint_neighbors


def get_fnat_score(estimated_complex, banchmark_complex):
    return float(len(intersect_neighbours_between_complexes(estimated_complex, banchmark_complex))/
                 float(len(banchmark_complex.get_neighbouring_residues())))

def get_joint_positions_from_structs(struct1, struct2):
    """get paths to two pdb files
    return a list oa list of tuples of the identical residues.
    """
    alignments = pairwise2.align.globalxx(
        get_structure_sequence(struct1),
        get_structure_sequence(struct2))
    return get_joint_positions_from_alignment(*alignments[0])


def get_joint_positions_from_alignment(align1, align2, score, begin, end):
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


