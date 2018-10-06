import os

from enum import Enum
import numpy as np
from Constants import *
from abc import ABCMeta, abstractmethod
from objects.processed_result import ComplexProcessedResult


class RaptorXScoringMethod(Enum):
    """
    Enum for various scoring methods
    """
    log_likelihood = "sum of log likelihood"  # sum of log probablities
    likelihood = "sum of likelihood"  # multiply probabilities
    percentage = "cutoff by percentile"  # cutoff percentage
    sqrt_sum = "sum of sqrt"
    cbrt_sum = "sum of cube"
    norm = "norm"  # euclidean distance


class Reranker(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def rerank(self, complexes, detailed_return=False):
        raise NotImplementedError("abstract class")


class RaptorxReranker(Reranker):

    def __init__(self, scoring_method, prob_trim=0.2, method_arg=None):
        """
        :param scoring_method: See utils.raptorx_utils.get_raptorx_score for details
        :param prob_trim: See utils.raptorx_utils.get_raptorx_score for details
        """
        self.scoring_method = scoring_method
        self.prob_trim = prob_trim
        self._method_arg = method_arg

    def rerank(self, complexes, raptorx_matrix=None, detailed_return=False):
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
        if raptorx_matrix is None:
            raptorx_matrix = RaptorxReranker.get_raptorx_matrix(comp)
        # Give each complex it's score
        ranks = []
        for i, res_complex in enumerate(complexes):
            neighbours = res_complex.get_neighbours_residues()
            score = float(RaptorxReranker.get_raptorx_score(raptorx_matrix, neighbours, self.scoring_method,
                                                            self.prob_trim, self._method_arg))
            ranks.append((i, score))
        ranks = sorted(ranks, key=lambda rank: rank[1], reverse=True)
        if detailed_return:
            return [(complexes[i], i, score) for i, score in ranks]
        return [complexes[i] for i, score in ranks]

    @staticmethod
    def get_raptorx_matrix(res_complex, filepath=None):
        # type: (ComplexProcessedResult, str, bool) -> np.ndarray
        """
        Returns an numpy 2D ndarray of the score raptorx matrix with the given arguments.
        The first dimension is the receptor sequence and the second the ligand
        Can be given a complex id to look for or a file path of the matrix.
        :param res_complex: complex of the raptorx matrix, uses the conventional paths specified in Constants
        :param filepath: file path of the matrix. Ignored when complex_id is specified (default: None).
        :return: raptorx matrix
        """
        desired_shape = (len(res_complex.receptor_sequence), len(res_complex.ligand_sequence))
        # convert complex_id to file path
        if filepath is None:
            dirpath = get_raptorx_dir_path(res_complex.complex_id)
            matrix_file_extension = '.gcnn_inter'
            filepath_list = [file_name for file_name in os.listdir(dirpath) if
                             file_name.endswith(matrix_file_extension)]
            if len(filepath_list) == 0:
                raise IOError("Couldn't find any *%s file in %s" % (matrix_file_extension, dirpath))
            if len(filepath_list) > 1:
                raise IOError("More then one *%s file exists in %s" % (matrix_file_extension, dirpath))
            filepath = os.path.join(dirpath, filepath_list[0])
        # load raptorx matrix
        raptorx_mat = np.loadtxt(filepath, delimiter='\t')
        assert raptorx_mat.ndim == 2, "Invalid file format. Loaded matrix is not 2D."
        # check if matrix matches desired shape
        if desired_shape is not None and raptorx_mat.shape != desired_shape:
            if raptorx_mat.T.shape == desired_shape:
                raptorx_mat = raptorx_mat.T  # transpose matrix, we mixed receptor and ligand
            else:
                raise IndexError("File matrix shape %s and desired shape %s doesn't match"
                                 % (raptorx_mat.shape, desired_shape))
        return raptorx_mat

    @staticmethod
    def get_raptorx_score(raptorx_mat, neighbour_indices, method, trim, arg):
        # type: (np.ndarray, Iterable[Tuple[int, int]], RaptorXScoringMethod, float) -> float
        """
        Returns score based on raptorx matrix with the given score method.
        :param raptorx_mat: matrix with probabilities
        :param neighbour_indices: <receptor_index, ligand_index> indices for scoring
        :param method: score method as specified by the RaptorXScoringMethod enum
        :param trim: trim probabilities which are lower the trim value
        :return: interaction score
        """
        # assert method not in RaptorXScoringMethod.__members__, "method is not in RaptorXScoringMethod enum"
        # get new ndarray by list of incides
        neighbour_scores = raptorx_mat[tuple(zip(*neighbour_indices))]
        # trim results
        neighbour_scores = neighbour_scores[neighbour_scores >= trim]
        if method == RaptorXScoringMethod.log_likelihood:
            return np.sum(np.log(neighbour_scores))
        elif method == RaptorXScoringMethod.likelihood:
            return np.prod(neighbour_scores)
        elif method == RaptorXScoringMethod.percentage:
            cutoff_value = np.percentile(raptorx_mat.flatten(), arg)
            return len(neighbour_scores[neighbour_scores >= cutoff_value])
        elif method == RaptorXScoringMethod.sqrt_sum:
            return np.sum(np.sqrt(neighbour_scores))
        elif method == RaptorXScoringMethod.cbrt_sum:
            return np.sum(np.cbrt(neighbour_scores))
        elif method == RaptorXScoringMethod.norm:
            return np.linalg.norm(neighbour_scores)

