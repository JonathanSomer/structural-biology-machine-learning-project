import random
import matplotlib.pyplot as plt
from utils import fnat_utils

import numpy as np
import scipy.stats as stats

from objects import complex as cmp
from objects.pipeline_handler import ResultsHelper
from utils import pdb_utils, raptorx_utils
from Constants import *

def plot_rank_to_fnat(complexes):
    # type: (List[complex.Complex]) -> None
    plt.figure(0)
    new_ranks = range(0, len(complexes))
    original_ranks = [complex.original_rank for complex in complexes]
    capri_scores = [random.random() for complex in complexes]
    plt.scatter(new_ranks, capri_scores, c='b')
    plt.scatter(original_ranks, capri_scores, c='r')
    plt.show()


def plot_raptor_to_fnat(result_helper):
    # type: (ResultsHelper) -> None
    """
    Scatter plot of (raptor score, fnat score) points for each patchdock results in result_helper
    Also adds a regression line
    :param result_helper:
    :return:
    """
    raptor_scores = np.array(result_helper.get_all_ranked_expectation_scores()).flatten()
    fnat_scores = np.array(result_helper.get_all_fnat_scores_after_reranking()).flatten()
    # regression line
    slope, intercept, r_value, p_value, std_err = stats.linregress(raptor_scores, fnat_scores)
    line = slope * raptor_scores + intercept

    plt.scatter(raptor_scores, fnat_scores, c='b')
    plt.plot(raptor_scores, line, 'r', label='y={:.3f}x+{:.3f} (R2={:.2f})'.format(slope, intercept, r_value))
    plt.legend(fontsize=9)
    plt.show()

def plot_average_raptor_score_in_interaction_site_vs_outside_site(bound_complexes, trim=0.01):
    complex_ids, average_for_all_residues, average_for_residues_in_interaction, average_for_residues_not_in_interaction = [], [], [], []
    interaction_big_prob_count, outer_big_prob_count = [], []

    complex_ids_raptor_score_is_larger_in_binding_site = []

    for complex in bound_complexes:
        complex_id = complex.complex_id
        print(complex_id)
        unbound_complex = cmp.BenchmarkComplex(complex.complex_id, type=cmp.ComplexType.zdock_benchmark_unbound)
        receptor_len = len(unbound_complex.receptor_sequence)
        ligand_len = len(unbound_complex.ligand_sequence)
        rapt_mat = raptorx_utils.get_raptorx_matrix(complex_id, desired_shape=(receptor_len, ligand_len))

        receptor_map = fnat_utils.get_position_map_between_sequences(complex.receptor_sequence,
                                                                     unbound_complex.receptor_sequence)
        ligand_map = fnat_utils.get_position_map_between_sequences(complex.ligand_sequence,
                                                                   unbound_complex.ligand_sequence)
        neighbours = complex.get_neighbours_residues()
        # change from bound indices to unbound indices
        neighbours = [(receptor_map[receptor_id], ligand_map[ligand_id]) for receptor_id, ligand_id in neighbours if
                      receptor_id in receptor_map and ligand_id in ligand_map]

        neighbours_indices = tuple(zip(*neighbours))
        mask_allow_closer_than_5_angstrom = np.zeros(rapt_mat.shape, dtype=bool)  # np.ones_like(a,dtype=bool)
        mask_allow_closer_than_5_angstrom[neighbours_indices] = True

        complex_ids.append(complex.complex_id)
        average_for_all_residues.append(np.average(rapt_mat))
        average_for_residues_in_interaction.append(np.average(rapt_mat[mask_allow_closer_than_5_angstrom]))
        average_for_residues_not_in_interaction.append(np.average(rapt_mat[~mask_allow_closer_than_5_angstrom]))

        if average_for_residues_in_interaction[-1] > average_for_residues_not_in_interaction[-1]:
            print("HA", complex_id)
            complex_ids_raptor_score_is_larger_in_binding_site.append(complex_id)

    X = complex_ids
    _X = np.arange(len(X))
    plt.figure(0)
    plt.bar(_X - 0.2, average_for_residues_in_interaction, 0.2)
    plt.bar(_X + 0.0, average_for_residues_not_in_interaction, 0.2)
    plt.bar(_X + 0.2, average_for_all_residues, 0.2)
    plt.xticks(_X, X)  # set labels manually

    plt.show()

    complex_ids_raptor_score_is_NOT_larger_in_binding_site = [c for c in complex_ids if c not in complex_ids_raptor_score_is_larger_in_binding_site]
    print("no improvement: ", complex_ids_raptor_score_is_NOT_larger_in_binding_site)
    print("improvement: ", complex_ids_raptor_score_is_larger_in_binding_site)
    print("all: ", complex_ids)

    plt.pie(x=[len(complex_ids_raptor_score_is_NOT_larger_in_binding_site), len(complex_ids_raptor_score_is_larger_in_binding_site)],
            labels=['no improvement', 'improvement'])

    plt.show()


# from visualizer.evodock_plot import *
# plot_raptor_values_bar_on_batch()
def plot_raptor_values_bar_on_batch():
    complex_ids = TRAIN_COMPLEX_IDS
    bound_complexes = []
    for complex_id in complex_ids:
        bound_complexes.append(cmp.BenchmarkComplex(complex_id=complex_id, type=cmp.ComplexType.zdock_benchmark_bound))

    plot_average_raptor_score_in_interaction_site_vs_outside_site(bound_complexes)

