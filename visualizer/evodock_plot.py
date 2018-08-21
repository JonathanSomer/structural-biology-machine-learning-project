import random
import matplotlib.pyplot as plt
from utils import fnat_utils

import numpy as np
import scipy.stats as stats

from objects import complex as cmp
from objects.pipeline_handler import ResultsHelper
from utils import pdb_utils, raptorx_utils
from Constants import *

def plot_rank_to_fnat(result_helper, accumulate=False):
    # type: (ResultsHelper) -> None
    """
    Scatter plot of (rank, fnat score) points for top10 patchdock results in result_helper before and after reranking
    :param result_helper: ResultsHelper with the data to plot
    :return: None
    """
    def process_scores(scores):
        if accumulate:
            scores = np.maximum.accumulate(scores, axis=1)
        return np.average(scores, axis=0)

    top = 10
    x = np.arange(10) + 1
    original_fnat_scores = process_scores(result_helper.get_all_fnat_scores(False, top))
    reranked_fnat_scores = process_scores(result_helper.get_all_fnat_scores(True, top))

    plt.ylim(ymin=0.0, ymax=max(np.max(original_fnat_scores), np.max(reranked_fnat_scores)) + 0.05)
    plt.scatter(x, original_fnat_scores, c='b')
    plt.scatter(x, reranked_fnat_scores, c='r', marker='X')
    plt.show()


def plot_raptor_to_fnat(result_helper):
    # type: (ResultsHelper) -> None
    """
    Scatter plot of (raptor score, fnat score) points for each patchdock results in result_helper
    Also adds a regression line
    :param result_helper: ResultsHelper with the data to plot
    :return: None
    """
    raptor_scores = np.array(result_helper.get_all_ranked_expectation_scores()).flatten()
    fnat_scores = np.array(result_helper.get_all_fnat_scores(True)).flatten()
    # regression line
    slope, intercept, r_value, p_value, std_err = stats.linregress(raptor_scores, fnat_scores)
    line = slope * raptor_scores + intercept

    plt.scatter(raptor_scores, fnat_scores, c='b')
    plt.plot(raptor_scores, line, 'r', label='y={:.3f}x+{:.3f} (R2={:.2f})'.format(slope, intercept, r_value))
    plt.legend(fontsize=9)
    plt.show()

# from visualizer.evodock_plot import *
# plot_average_raptor_score_in_binding_site_vs_not()
def plot_average_raptor_score_in_binding_site_vs_not(trim=0.01):
    bound_complexes = [cmp.BenchmarkComplex(complex_id=complex_id, type=cmp.ComplexType.zdock_benchmark_bound) for complex_id in TRAIN_COMPLEX_IDS]

    complex_ids = np.array([bound_complex.complex_id for bound_complex in bound_complexes])
    average_raptor_scores_for_neighbors, average_raptor_scores_for_non_neighbors = [], []

    for bound_complex in bound_complexes:
        complex_id = bound_complex.complex_id

        raptor_matrix = raptorx_utils.get_raptor_matrix_for_bound_complex(bound_complex)
        neighbor_indexes = fnat_utils.get_neighbours_for_bound_complex(bound_complex)

        neighbor_mask = np.zeros(raptor_matrix.shape, dtype=bool)
        neighbor_mask[neighbor_indexes] = True

        average_raptor_scores_for_neighbors.append(np.average(raptor_matrix[neighbor_mask]))
        average_raptor_scores_for_non_neighbors.append(np.average(raptor_matrix[~neighbor_mask]))

    # import pdb; pdb.set_trace()
    indexes = np.subtract(average_raptor_scores_for_non_neighbors, average_raptor_scores_for_neighbors).argsort()

    complex_ids = complex_ids[indexes]
    average_raptor_scores_for_neighbors = np.array(average_raptor_scores_for_neighbors)[indexes]
    average_raptor_scores_for_non_neighbors = np.array(average_raptor_scores_for_non_neighbors)[indexes]

    # import pdb; pdb.set_trace()
    _X = np.arange(len(complex_ids))
    plt.figure(0)
    plt.bar(_X - 0.2, average_raptor_scores_for_neighbors, 0.2)
    plt.bar(_X + 0.0, average_raptor_scores_for_non_neighbors, 0.2)
    plt.xticks(_X, complex_ids, rotation=45)  # set labels manually
    plt.legend(['average_raptor_score_for_neighbors', 'average_raptor_score_for_non_neighbors'],loc=1)
    plt.xlabel('Complex ID', fontsize=16)
    plt.ylabel('Average Raptor Score', fontsize=16)
    plt.show()

    improvement_complex_ids = complex_ids[np.greater(average_raptor_scores_for_neighbors, average_raptor_scores_for_non_neighbors)]
    no_improvement_complex_ids = [c for c in complex_ids if c not in improvement_complex_ids]

    explode = (0, 0.1)
    plt.pie(x=[len(no_improvement_complex_ids), len(improvement_complex_ids)], explode=explode, shadow=True, autopct='%1.1f%%')
    plt.legend(['average raptor score lower in binding site', 'average raptor score higher in binding site'], loc=2)
    plt.show()
