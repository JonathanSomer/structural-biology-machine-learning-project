import numpy as np

from objects.complex import *
from utils.capri_utils import get_capri_score
from utils.fnat_utils import get_fnat_score


class ResultsHelper(object):

    def __init__(self, complex_ids, n_of_patchdock_results, reranker, ignore_failure=False, verbose=True):
        self.n_of_patchdock_results = n_of_patchdock_results
        self.reranker = reranker
        self.complex_helpers = {}
        for c_id in complex_ids:
            if verbose:
                print("ResultHelper Loading complex_id: %s" % c_id)
            if ignore_failure:

                try:
                    self.complex_helpers[c_id] = ComplexHelper(c_id, self.n_of_patchdock_results, self.reranker)
                except (OSError, IOError), e:
                    if verbose:
                        print("%s: %s" % (c_id, str(e)))
                    pass
            else:
                self.complex_helpers[c_id] = ComplexHelper(c_id, self.n_of_patchdock_results, self.reranker)
        self.complex_ids = list(self.complex_helpers.keys())

    def get_all_capri_scores_of_original_patchdock_ranking(self):
        return np.array([self.get_capri_score_of_original_patchdock_ranking(complex_id)
                         for complex_id in self.complex_ids])

    def get_capri_score_of_original_patchdock_ranking(self, complex_id):
        return self.complex_helpers[complex_id].get_capri_score_of_original_patchdock_ranking()

    def get_all_capri_scores_of_reranking(self):
        return np.array([self.get_capri_score_of_reranking(complex_id)
                         for complex_id in self.complex_ids])

    def get_capri_score_of_reranking(self, complex_id):
        return self.complex_helpers[complex_id].get_capri_score_of_reranking()

    def get_all_fnat_scores_reranking(self, after, top=None):
        return np.array([self.get_fnat_scores_reranking(complex_id, after, top)
                         for complex_id in self.complex_ids])

    def get_fnat_scores_reranking(self, complex_id, after, top=None):
        return self.complex_helpers[complex_id].get_fnat_scores_before_reranking(top)

    def get_all_ranked_expectation_scores(self):
        return np.array([self.get_ranked_expectation_scores(complex_id)
                         for complex_id in self.complex_ids])

    def get_ranked_expectation_scores(self, complex_id):
        return self.complex_helpers[complex_id].get_ranked_expectation_scores()


class ComplexHelper(object):

    def __init__(self, complex_id, n_of_patchdock_results, reranker):
        self.complex_id = complex_id
        self.n_of_patchdock_results = n_of_patchdock_results
        self.reranker = reranker
        self.original_ranked_complexes = self._get_patchdock_results_by_original_rank(complex_id,
                                                                                      self.n_of_patchdock_results)
        self.unbound_complex = BenchmarkComplex(complex_id, ComplexType.zdock_benchmark_unbound)
        self.bound_complex = BenchmarkComplex(complex_id, ComplexType.zdock_benchmark_bound)
        self.reranked_complexes = reranker.rerank(self.original_ranked_complexes, detailed_return=False)

    def _get_patchdock_results_by_original_rank(self, complex_id, n_of_patchdock_results):
        return [PatchDockComplex(complex_id, i + 1) for i in range(n_of_patchdock_results)]

    def get_capri_score_of_original_patchdock_ranking(self):
        return get_capri_score(self.original_ranked_complexes, self.bound_complex)

    def get_capri_score_of_reranking(self):
        return get_capri_score(self.reranked_complexes, self.bound_complex)

    def get_fnat_scores_reranking(self, after, top=None):
        scores = []
        complexes = self.reranked_complexes if after else self.original_ranked_complexes
        for i in range(top or self.n_of_patchdock_results):
            scores.append(get_fnat_score(complexes[i], self.bound_complex))
        return np.array(scores)

    def get_ranked_expectation_scores(self):
        detailed_rerank = self.reranker.rerank(self.original_ranked_complexes, detailed_return=True)
        return np.array([x[2] for x in detailed_rerank])
