from objects import complex
from utils import pdb_utils, raptorx_utils
from abc import ABCMeta, abstractmethod

class Reranker(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def rerank(self, complexes, detailed_return=False):
        raise NotImplementedError("abstract class")

class RaptorxReranker(Reranker):

    def __init__(self, scoring_method, prob_trim=0.2):
        '''
        :param scoring_method: See utils.raptorx_utils.get_raptorx_score for details
        :param prob_trim: See utils.raptorx_utils.get_raptorx_score for details
        '''
        self.scoring_method = scoring_method
        self.prob_trim = prob_trim

    def rerank(self, complexes, detailed_return=False):
        # type: (self, List[complex.Complex]) -> Union[List[complex.Complex], List[complex.Complex, int, float]]
        """
        Reranks a list of Complexes by raptorx matrix and given scoring parameters.
        :param complexes: list of complexes to rerank
        :param detailed_return: returns a more detailed list with the score of each complex
        :return: reordered list by the score calculated using the conf params given in init
        """
        if not isinstance(complexes, list):
            raise TypeError("complexes argument must by of list type")
        if len(complexes) < 1:
            raise ValueError("List of complexes must contain at least one complex")
        comp = complexes[0]
        complex_id = comp.complex_id
        # Get raptorx matrix using the first complex data
        receptor_len = len(pdb_utils.get_structure_sequence(comp.receptor))
        ligand_len = len(pdb_utils.get_structure_sequence(comp.ligand))
        rapt_mat = raptorx_utils.get_raptorx_matrix(complex_id, desired_shape=(receptor_len, ligand_len))
        # Give each complex it's score
        ranks = []
        for i, res_complex in enumerate(complexes):
            neighbours = res_complex.get_neighbours_residues()
            score = float(raptorx_utils.get_raptorx_score(rapt_mat, neighbours, self.scoring_method, self.prob_trim))
            ranks.append((i, score))
        ranks = sorted(ranks, key=lambda rank: rank[1], reverse=True)
        if detailed_return:
            return [(complexes[i], i, score) for i, score in ranks]
        return [complexes[i] for i, score in ranks]
