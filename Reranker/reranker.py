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
    sqrt_sum = "sum of sqrt"
    sum = "sum"
    cbrt_sum = "sum of cube"
    power_sum = "sum of {:.2f} power"
    norm = "norm"  # euclidean distance
    average = "average"
    count = "count"  # simply count the number of neighbours
    log_len_mul_sum = "log(len(nb lst)*sum(nb raptor scores)"


class Reranker(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def rerank(self, complexes, detailed_return=False):
        raise NotImplementedError("abstract class")


class RaptorxReranker(Reranker):

    def __init__(self, scoring_method, prob_trim=0.1, percentile_trim=0.8, method_arg=None, shuffle_raptor_matrix=False):
        """
        :param scoring_method: See utils.raptorx_utils.get_raptorx_score for details
        :param prob_trim: See utils.raptorx_utils.get_raptorx_score for details
        """
        self.prob_trim = prob_trim
        self.percentile_trim = percentile_trim
        self._method_arg = method_arg
        self.scoring_method = scoring_method
        self.shuffle_raptor_matrix= shuffle_raptor_matrix

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
            raptorx_matrix = RaptorxReranker.get_raptorx_matrix(comp, shuffle_raptor_matrix=self.shuffle_raptor_matrix)
        # Give each complex it's score
        ranks = []
        for i, res_complex in enumerate(complexes):
            neighbours = res_complex.get_neighbours_residues()
            score = float(self.get_raptorx_score(raptorx_matrix, neighbours))
            ranks.append((i, score))
        ranks = sorted(ranks, key=lambda rank: rank[1], reverse=True)
        if detailed_return:
            return [(complexes[i], i, score) for i, score in ranks]
        return [complexes[i] for i, score in ranks]

    @staticmethod
    def get_raptorx_matrix(res_complex, filepath=None, shuffle_raptor_matrix=False):
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

        if shuffle_raptor_matrix:
            shape = raptorx_mat.shape
            raptorx_mat = raptorx_mat.flatten()
            np.random.shuffle(raptorx_mat)
            raptorx_mat = np.reshape(raptorx_mat, shape)
        return raptorx_mat

    def get_raptorx_score(self, raptorx_mat, neighbour_indices):
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
        method = self.scoring_method
        neighbour_scores = raptorx_mat[tuple(zip(*neighbour_indices))]
        l = len(neighbour_scores)
        # trim results by probability and percentile
        neighbour_scores = neighbour_scores[neighbour_scores >= self.prob_trim]
        if neighbour_scores.size > 0:
            rapt_percentile = np.percentile(neighbour_scores, 100 * self.percentile_trim)
            neighbour_scores = neighbour_scores[neighbour_scores >= rapt_percentile]

        if method == RaptorXScoringMethod.log_likelihood:
            return np.sum(np.log(neighbour_scores))
        elif method == RaptorXScoringMethod.likelihood:
            return np.prod(neighbour_scores)
        elif method == RaptorXScoringMethod.sqrt_sum:
            return  np.sum(np.sqrt(neighbour_scores))
        elif method == RaptorXScoringMethod.cbrt_sum:
            return np.sum(np.cbrt(neighbour_scores))
        elif method == RaptorXScoringMethod.sum:
            return np.sum(neighbour_scores)
        elif method == RaptorXScoringMethod.power_sum:
            return np.sum(np.power(neighbour_scores, self._method_arg))
        elif method == RaptorXScoringMethod.norm:
            return np.linalg.norm(neighbour_scores)
        elif method == RaptorXScoringMethod.average:
            return 0.0 if neighbour_scores.size == 0 else np.average(neighbour_scores)
        elif method == RaptorXScoringMethod.count:  # ignore trims
            return raptorx_mat[tuple(zip(*neighbour_indices))].size
        elif method == RaptorXScoringMethod.log_len_mul_sum:
            return np.log(len(neighbour_scores)+1)*np.sum(neighbour_scores)

    def __str__(self):
        method = self.scoring_method
        method_str = method.value.format(self._method_arg) if method == RaptorXScoringMethod.power_sum else method.value
        return "reranker-method:{}".format(method_str)
