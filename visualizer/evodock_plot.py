import matplotlib.pyplot as plt
from utils import fnat_utils, capri_utils, raptorx_utils
from utils.fnat_utils import get_fnat_score
import numpy as np
import scipy.stats as stats

from objects.complex import BenchmarkComplex, ComplexType
from objects import complex as cmp
from objects.results_helper import ResultsHelper
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
    original_fnat_scores = process_scores(result_helper.get_all_fnat_scores(after=False, top=top))
    reranked_fnat_scores = process_scores(result_helper.get_all_fnat_scores(after=True, top=top))

    plt.ylim(ymin=0.0, ymax=max(np.max(original_fnat_scores), np.max(reranked_fnat_scores)) + 0.05)
    plt.scatter(x, original_fnat_scores, c='b', label='Original Ranking')
    plt.scatter(x, reranked_fnat_scores, c='r', marker='X', label='Reranked')
    plt.legend(fontsize=9)
    plt.xlabel('Docking Rank')
    plt.ylabel('fnat')
    plt.show()


def plot_raptor_to_fnat(result_helper, standardize_fnat=True):
    # type: (ResultsHelper) -> None
    """
    Scatter plot of (raptor score, fnat score) points for each patchdock results in result_helper
    Also adds a regression line
    :param result_helper: ResultsHelper with the data to plot
    :param standardize_fnat: Should the fnat be standardize. Standardization is per complex.
    :return: None
    """

    def create_regression_line(x, y):
        # regression line
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
        line = slope * x + intercept
        plt.plot(x, line, 'r', label='y={:.3f}x+{:.3f} (R2={:.2f})'.format(slope, intercept, r_value))

    raptor_scores = np.log(result_helper.get_all_ranked_expectation_scores())
    fnat = result_helper.get_all_fnat_scores(after=True)
    if standardize_fnat:
        fnat = stats.zscore(fnat, axis=1)
        # remove rows with nan, this occurs if all the complex's results have the exact same fnat (usually 0.0)
        nan_mask = ~np.isnan(fnat).any(axis=1)
        fnat = fnat[nan_mask]
        raptor_scores = raptor_scores[nan_mask]
    raptor_scores_flat = raptor_scores.flatten()
    fnat_scores = fnat.flatten()

    plt.scatter(raptor_scores_flat, fnat_scores, c='b')
    create_regression_line(raptor_scores_flat, fnat_scores)
    # create_regression_line(raptor_scores[fnat_scores > 0], fnat_scores[fnat_scores > 0])
    plt.legend(fontsize=9)
    plt.xlabel('RaptorX Score')
    plt.ylabel('fnat' + 'normalized per complex' if standardize_fnat else '')
    plt.show()


def plot_average_raptor_score_in_binding_site_vs_not(trim=0.01):
    # todo: need to make this method use only surface resdius
    bound_complexes = [cmp.BenchmarkComplex(complex_id=complex_id, type=cmp.ComplexType.zdock_benchmark_bound) for
                       complex_id in get_all_training_complexes()]

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
    plt.legend(['average_raptor_score_for_neighbors', 'average_raptor_score_for_non_neighbors'], loc=1)
    plt.xlabel('Complex ID', fontsize=16)
    plt.ylabel('Average Raptor Score', fontsize=16)
    plt.show()

    improvement_complex_ids = complex_ids[
        np.greater(average_raptor_scores_for_neighbors, average_raptor_scores_for_non_neighbors)]
    no_improvement_complex_ids = [c for c in complex_ids if c not in improvement_complex_ids]

    explode = (0, 0.1)
    plt.pie(x=[len(no_improvement_complex_ids), len(improvement_complex_ids)], explode=explode, shadow=True,
            autopct='%1.1f%%')
    plt.legend(['average raptor score lower in binding site', 'average raptor score higher in binding site'], loc=2)
    plt.show()


def plot_fnat_above_threshold_per_complex(result_helper, threshold=capri_utils.FnatThresholds.Acceptable):
    fnat_before = np.array(result_helper.get_all_fnat_scores(after=False, top=TOP_RESULTS_COUNTS_FOR_CAPRI))
    fnat_after = np.array(result_helper.get_all_fnat_scores(after=True, top=TOP_RESULTS_COUNTS_FOR_CAPRI))

    above_before = np.array([len(scores[scores > threshold]) for scores in fnat_before])
    above_after = np.array([len(scores[scores > threshold]) for scores in fnat_after])

    sorted_indices = above_before.argsort()[::-1]
    above_before = above_before[sorted_indices]
    above_after = above_after[sorted_indices]

    plt.plot(above_after, 'bo')
    plt.plot(above_before, 'rx')
    plt.xlabel("complexes")
    plt.ylabel("acceptable count")
    plt.show()


def plot_max_patchdock_fnat_scores(result_helper, top=None):
    max_fnat = np.array([max(fnats) for fnats
                         in result_helper.get_all_fnat_scores(after=False, top=top)])

    sorted_indices = max_fnat.argsort()[::-1]
    max_fnat = max_fnat[sorted_indices]
    ids = np.array(result_helper.complex_ids)[sorted_indices]

    plt.plot(max_fnat, 'bo', label="bound")
    plt.plot([capri_utils.FnatThresholds.Acceptable for dot in max_fnat])  # treshold line
    plt.xticks(np.arange(len(ids)), ids, rotation=45)
    plt.xlabel("complex id")
    plt.ylabel("max fnat score for all patchdock results")
    plt.show()


def plot_fnat_of_unbound(ids):
    unbounds = [BenchmarkComplex(complex_id, type=ComplexType.zdock_benchmark_unbound)
                for complex_id in ids]
    bounds = [BenchmarkComplex(complex_id, type=ComplexType.zdock_benchmark_bound)
              for complex_id in ids]

    unbound_fnats = [get_fnat_score(unbound, bound) for unbound, bound in zip(unbounds, bounds)]
    unbound_fnats.sort(reverse=True)

    plt.plot([capri_utils.FnatThresholds.Acceptable for dot in unbound_fnats])  # treshold line
    plt.plot(unbound_fnats, 'bo')
    plt.ylabel("fnat of unbound")
    plt.xlabel("complexs")
    plt.show()


def raptor_to_fnat_plot_per_complex(result_helper):
    all_raptor_scores = result_helper.get_all_ranked_expectation_scores()
    all_fnat_scores = result_helper.get_all_fnat_scores(after=True)
    n_plots = len(result_helper.complex_ids)

    for i, raptor_scores, fnat_scores in zip(range(n_plots), all_raptor_scores, all_fnat_scores):
        raptor_scores_top10, fnat_scores_top_10 = raptor_scores[:10], fnat_scores[:10]

        # regression line
        slope, intercept, r_value, p_value, std_err = stats.linregress(raptor_scores, fnat_scores)
        line = slope * raptor_scores + intercept

        plt.subplot(2, n_plots / 2 + n_plots % 2, i + 1)

        # plt.plot([FnatThresholds.Acceptable for dot in range(max(np.array(raptor_scores_top10).astype(int)))]) #treshold line
        plt.scatter(raptor_scores, fnat_scores, c='b')
        plt.scatter(raptor_scores_top10, fnat_scores_top_10, c='r')
        plt.plot(raptor_scores, line, 'r',
                 label='y={:.3f}x+{:.3f} (R2={:.2f})'.format(slope, intercept, r_value))
        plt.legend(fontsize=9)
        plt.title(result_helper.complex_ids[i])

    plt.suptitle("fant over {} raptor scores".format(result_helper.reranker.scoring_method))
    plt.show()


def plot_improved_percentage(result_helper, top=[10]):
    # type: (ResultsHelper, List[int]) -> None
    """
    Plot a barplot with percentage of improvement in result reranking over all complexes
    :param result_helper: helper for computing results
    :param top: the amount of results to add to the calculation per complex
    :return: None
    """

    def autolabel(rects):
        """
        Attach a text label above each bar displaying its height
        """
        for rect in rects:
            height = rect.get_height()
            plt.text(rect.get_x() + rect.get_width() / 2., height + 0.1, "{:.2f}".format(height), ha='center',
                     va='bottom')

    if isinstance(top, int):
        top = [top]

    width = 1 / (len(top) + 1)
    x = ['better', 'equals', 'worse']
    ids = np.arange(len(x))
    for i, t in enumerate(top):
        original_score = result_helper.get_all_capri_scores(after=False, top=t)
        rerank_score = result_helper.get_all_capri_scores(after=True, top=t)
        total_comps = original_score.size

        # get x and y
        greater_count = np.sum(np.greater(rerank_score, original_score))
        equal_count = np.sum(np.equal(rerank_score, original_score))
        less_count = np.sum(np.less(rerank_score, original_score))
        y = 100 * np.array([greater_count, equal_count, less_count]) / total_comps

        # plot x and y in bar plot
        rects = plt.bar(ids + (i + 0.5 - len(top) / 2) * width, y, width=width, label="top {}".format(t))
        autolabel(rects)
    plt.title("Reranking impact per top results")
    plt.xlabel("reranking impact")
    plt.xticks(ids, x)
    plt.ylabel("percentage")
    plt.legend()
    plt.show()


def _plot_func_per_fnat_thresholds_before_vs_after(result_helper, func, top=10):
    # type: (ResultsHelper, Callable[[np.ndarray, capri_utils.FnatThresholds], float], List[int]) -> None
    """
    Plot a barplot with func values of before vs after per fnat threshold
    :param result_helper: helper for computing results
    :param func: function to perform
    :param top: the amount of results to add to the calculation per complex
    :return: None
    """

    def autolabel(rects):
        """
        Attach a text label above each bar displaying its height
        """
        for rect in rects:
            height = rect.get_height()
            plt.text(rect.get_x() + rect.get_width() / 2., height + 0.1, "{:.2f}".format(height), ha='center',
                     va='bottom')

    # if isinstance(top, int):
    #    top = [top]
    thresholds = list(capri_utils.FnatThresholds)
    width = 1 / 3
    x = [t.name for t in thresholds]
    ids = np.arange(len(x))

    original_fnats = result_helper.get_all_fnat_scores(after=False, top=top)
    rerank_fnats = result_helper.get_all_fnat_scores(after=True, top=top)
    y_before, y_after = [], []
    for i, threshold in enumerate(thresholds):
        # update y
        t_val = threshold.value
        y_before.append(func(original_fnats, t_val))
        y_after.append(func(rerank_fnats, t_val))

    # plot x and y in bar plot
    rects1 = plt.bar(ids - 0.5 * width, y_before, width=width, label="original ranking top {}".format(top))
    rects2 = plt.bar(ids + 0.5 * width, y_after, width=width, label="{} top {}".format(result_helper.reranker, top))
    autolabel(rects1)
    autolabel(rects2)
    plt.xticks(ids, x)
    plt.legend()


def plot_at_least_one_percentage_per_fnat_thresholds_before_vs_after(result_helper, top=10):
    # type: (ResultsHelper, int) -> None
    """
    Plot a barplot with percentage of complexes with at least one result above each threshold
    :param result_helper: helper for computing results
    :param top: the amount of results to add to the calculation per complex
    :return: None
    """

    def at_least_one_percentage(fnats, threshold):
        complex_num = fnats.shape[0]
        return 100 * np.sum(np.any(fnats >= threshold, axis=1)) / complex_num

    _plot_func_per_fnat_thresholds_before_vs_after(result_helper, at_least_one_percentage, top)
    plt.title("Reranking impact per top results")
    plt.xlabel("reranking impact")
    plt.ylabel("percentage")
    plt.show()
