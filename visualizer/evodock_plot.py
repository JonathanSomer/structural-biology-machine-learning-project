gitimport random
import matplotlib.pyplot as plt
from utils import fnat_utils
from utils.fnat_utils import get_fnat_score
import numpy as np
import scipy.stats as stats

from objects.complex import BenchmarkComplex, ComplexType
from objects import complex as cmp
from objects.results_helper import ResultsHelper
from utils import pdb_utils, raptorx_utils
from Constants import *
from utils.capri_utils import FnatThresholds

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
    plt.xticks(_X, complex_ids)  # set labels manually

    plt.show()

    improvement_complex_ids = complex_ids[np.greater(average_raptor_scores_for_neighbors, average_raptor_scores_for_non_neighbors)]
    no_improvement_complex_ids = [c for c in complex_ids if c not in improvement_complex_ids]

    plt.pie(x=[len(no_improvement_complex_ids), len(improvement_complex_ids)],
            labels=['average raptor score lower in binding site', 'average raptor score higher in binding site'])

    plt.show()

def plot_fnat_above_treshold_per_complex(result_helper, treshold=FnatThresholds.Acceptable):
    fnat_before = np.array(result_helper.get_all_fnat_scores(after=False, top=TOP_RESULTS_COUNTS_FOR_CAPRI))
    fnat_after = np.array(result_helper.get_all_fnat_scores(after=True, top=TOP_RESULTS_COUNTS_FOR_CAPRI))

    above_before = np.array([len(scores[scores > treshold]) for scores in fnat_before])
    above_after = np.array([len(scores[scores > treshold]) for scores in fnat_after])

    sorted_indices = above_before.argsort()[::-1]
    above_before = above_before[sorted_indices]
    above_after = above_after[sorted_indices]

    plt.plot(above_after, 'bo')
    plt.plot(above_before, 'rx')
    plt.xlabel("complexes")
    plt.ylabel("acceptable count")
    plt.show()

def plot_max_patchdock_fnat_scores(result_helper, top=None):
    max_fnats_bound = np.array([max(fnats) for fnats
                                in result_helper.get_all_fnat_scores(after=False, top=top, bound=True)])
    max_fnats_unbound = np.array([max(fnats) for fnats
                                  in result_helper.get_all_fnat_scores(after=False, top=top, bound=False)])

    sorted_indices = max_fnats_bound.argsort()[::-1]
    max_fnats_bound = max_fnats_bound[sorted_indices]
    max_fnats_unbound = max_fnats_unbound[sorted_indices]

    plt.plot(max_fnats_bound, 'bo')
    plt.plot(max_fnats_unbound, 'rx')
    plt.plot([FnatThresholds.Acceptable for dot in max_fnats_bound]) #treshold line
    plt.show()

def plot_fnat_of_unbound(ids=TRAIN_COMPLEX_IDS):
    unbounds = [BenchmarkComplex(complex_id, type=ComplexType.zdock_benchmark_unbound)
                for complex_id in ids]
    bounds = [BenchmarkComplex(complex_id, type=ComplexType.zdock_benchmark_bound)
                for complex_id in ids]

    unbound_fnats = [get_fnat_score(unbound, bound) for unbound, bound in zip(unbounds, bounds)]
    unbound_fnats.sort(reverse=True)

    plt.plot([FnatThresholds.Acceptable for dot in unbound_fnats]) #treshold line
    plt.plot(unbound_fnats, 'bo')
    plt.ylabel("fnat of unbound")
    plt.xlabel("complexs")
    plt.show()
