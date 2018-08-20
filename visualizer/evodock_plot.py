import random
import matplotlib.pyplot as plt
from utils import fnat_utils

import numpy as np

from objects import complex as cmp
from utils import pdb_utils, raptorx_utils


def plot_rank_to_fnat(complexes):
    # type: (List[complex.Complex]) -> None
    plt.figure(0)
    new_ranks = range(0, len(complexes))
    original_ranks = [complex.original_rank for complex in complexes]
    capri_scores = [random.random() for complex in complexes]
    plt.scatter(new_ranks, capri_scores, c='b')
    plt.scatter(original_ranks, capri_scores, c='r')
    plt.show()


def plot_raptor_to_fnat(complexes, raptor_scores):
    # type: (List[cmp.Complex]) -> None
    plt.figure(1)
    capri_scores = [random.random() for complex in complexes]
    plt.scatter(raptor_scores, capri_scores, c='b')
    plt.show()


def plot_raptor_values_bar(bound_complexes, trim=0.01):
    complex_axis, general_avgs, interaction_avgs, outer_avgs = [], [], [], []
    interaction_big_prob_count, outer_big_prob_count = [], []
    for complex in bound_complexes:
        complex_id = complex.complex_id
        print complex_id
        unbound_complex = cmp.BenchmarkComplex(complex.complex_id, type=cmp.ComplexType.zdock_benchmark_unbound)
        receptor_len = len(unbound_complex.receptor_sequence)
        ligand_len = len(unbound_complex.ligand_sequence)
        rapt_mat = raptorx_utils.get_raptorx_matrix(complex_id, desired_shape=(receptor_len, ligand_len))

        receptor_map = fnat_utils.get_position_map_between_sequences(complex.receptor_sequence,
                                                                     unbound_complex.receptor_sequence)
        ligand_map = fnat_utils.get_position_map_between_sequences(complex.ligand_sequence,
                                                                   unbound_complex.ligand_sequence)
        neighbours = complex.get_neighbouring_residues()
        # change from bound indices to unbound indices
        neighbours = [(receptor_map[receptor_id], ligand_map[ligand_id]) for receptor_id, ligand_id in neighbours if
                      receptor_id in receptor_map and ligand_id in ligand_map]

        neighbours_indices = tuple(zip(*neighbours))
        mask = np.zeros(rapt_mat.shape, dtype=bool)  # np.ones_like(a,dtype=bool)
        mask[neighbours_indices] = True

        complex_axis.append(complex.complex_id)
        general_avgs.append(np.average(rapt_mat))
        interaction_avgs.append(np.average(rapt_mat[mask]))
        outer_avgs.append(np.average(rapt_mat[~mask]))
        interaction_big_prob_count.append(np.sum(rapt_mat[mask] >= trim) / float(rapt_mat[mask].size))
        outer_big_prob_count.append(np.sum(rapt_mat[~mask] >= trim) / float(rapt_mat[~mask].size))
    X = complex_axis
    _X = np.arange(len(X))
    plt.figure(0)
    plt.bar(_X - 0.2, interaction_avgs, 0.2)
    plt.bar(_X + 0.0, outer_avgs, 0.2)
    plt.bar(_X + 0.2, general_avgs, 0.2)
    plt.xticks(_X, X)  # set labels manually

    plt.figure(1)
    plt.bar(_X - 0.2, interaction_big_prob_count, 0.2)
    plt.bar(_X + 0.0, outer_big_prob_count, 0.2)
    plt.xticks(_X, X)  # set labels manually

    plt.show()
