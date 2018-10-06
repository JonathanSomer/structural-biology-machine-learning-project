import numpy as np

import utils.capri_utils as capri_utils
from objects.complex import *
from objects.processed_result import ComplexProcessedResult


class ResultsHelper(object):

    def __init__(self, complex_ids, n_of_patchdock_results, reranker,
                 ignore_failure=False, verbose=True, recache=False):
        self.n_of_patchdock_results = n_of_patchdock_results
        self.reranker = reranker
        self.complex_helpers = {}
        for c_id in complex_ids:
            if verbose:
                print("ResultHelper Loading complex_id: %s" % c_id)
            try:
                self.complex_helpers[c_id] = ComplexHelper(c_id, self.n_of_patchdock_results, self.reranker,
                                                           recache)
            except Exception as e:
                if verbose:
                    print("%s: %s" % (c_id, str(e)))
                if not ignore_failure:
                    raise
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

    def get_all_fnat_scores(self, after, top=None, bound=True):
        # type: (bool, Union[int, None]) -> np.ndarray
        """
        Gets the fnat scores before or after reranking (depending on after argument) for all complexes.
        :param after: Should order the scores by the reranking order or the original order
        :param top: return only the top results, by default None, which returns all scores
        :return: np.ndarray of fnat scores
        """
        return np.array([self.get_fnat_scores(complex_id, after, top, bound)
                         for complex_id in self.complex_ids])

    def get_fnat_scores(self, complex_id, after, top=None, bound=True):
        return self.complex_helpers[complex_id].get_fnat_scores(after, top, bound)

    def get_all_ranked_expectation_scores(self):
        return np.array([self.get_ranked_expectation_scores(complex_id)
                         for complex_id in self.complex_ids])

    def get_ranked_expectation_scores(self, complex_id):
        return self.complex_helpers[complex_id].get_ranked_expectation_scores()


class ComplexHelper(object):

    def __init__(self, complex_id, n_of_patchdock_results, reranker, recache=False):
        self.complex_id = complex_id
        self.n_of_patchdock_results = n_of_patchdock_results
        self.reranker = reranker
        self.recache = recache
        self.receptor_sequence = self.get_sequence_from_fasta(ligand=False)
        self.ligand_sequence = self.get_sequence_from_fasta(ligand=True)
        self.original_ranked_complexes = self._get_patchdock_results()
        self.detailed_reranked_complexes = reranker.rerank(self.original_ranked_complexes, detailed_return=True)
        self.reranked_complexes = [t[0] for t in self.detailed_reranked_complexes]

    def _get_patchdock_results(self):
        # type: (str, int) -> List[ComplexProcessedResult]
        with open(get_processed_data_json_path(self.complex_id)) as f:
            results_json = json.load(f)
        return [ComplexProcessedResult(self.complex_id, self.receptor_sequence, self.ligand_sequence, results_json[i])
                for i in range(self.n_of_patchdock_results)]

    def get_sequence_from_fasta(self, ligand):
        with open(get_complex_fasta_path(self.complex_id, ligand)) as f:
            next(f)
            return f.readline()

    def get_capri_score_of_original_patchdock_ranking(self):
        return self.get_capri_score(False)

    def get_capri_score_of_reranking(self):
        return self.get_capri_score(True)

    def get_capri_score(self, after):
        complexes = self.reranked_complexes if after else self.original_ranked_complexes
        top_10 = complexes[:10]
        x= sum([capri_utils.convert_fnat_to_capri(comp.get_fnat_score()) for comp in top_10])
        return x

    def get_fnat_scores(self, after, top=None, bound=True):
        complexes = self.reranked_complexes if after else self.original_ranked_complexes
        return np.array([complexes[i].get_fnat_score() for i in range(top or self.n_of_patchdock_results)])

    def get_ranked_expectation_scores(self):
        return np.array([t[2] for t in self.detailed_reranked_complexes])
