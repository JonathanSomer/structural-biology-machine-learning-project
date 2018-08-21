import os

from enum import Enum
import numpy as np
from objects import complex
from utils import pdb_utils
from Constants import get_raptorx_dir_path
from abc import ABCMeta, abstractmethod
from Reranker.regression_model import *
from objects.complex import *

class RaptorXScoringMethod(Enum):
    """
    Enum for various scoring methods
    """
    log_likelihood = 1  # sum of log probablities
    likelihood = 2  # multiply probabilities


class Reranker(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def rerank(self, complexes, detailed_return=False):
        raise NotImplementedError("abstract class")


class RaptorxReranker(Reranker):

    def __init__(self, scoring_method, prob_trim=0.2):
        """
        :param scoring_method: See utils.raptorx_utils.get_raptorx_score for details
        :param prob_trim: See utils.raptorx_utils.get_raptorx_score for details
        """
        self.scoring_method = scoring_method
        self.prob_trim = prob_trim

    @staticmethod
    def get_raptorx_matrix(complex_id, filepath=None, desired_shape=None):
        # type: (str, str, Tuple[int, int]) -> np.ndarray
        """
        Returns an numpy 2D ndarray of the score raptorx matrix with the given arguments.
        The first dimension is the receptor sequence and the second the ligand
        Can be given a complex id to look for or a file path of the matrix.
        :param complex_id: complex id of the raptorx matrix, uses the convensional paths specified in Constants
        :param filepath: file path of the matrix. Ignored when complex_id is specified (default: None).
        :param desired_shape: desired shape of the matrix. Ensures the matrix ligand and receptor didn't switch places.
        :return: raptorx matrix
        """
        # convert complex_id to file path
        if complex_id is not None:
            dirpath = get_raptorx_dir_path(complex_id)
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
    def get_raptorx_score(raptorx_mat, neighbour_indices, method, trim):
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
        receptor_len = len(comp.receptor_sequence)
        ligand_len = len(comp.ligand_sequence)
        rapt_mat = RaptorxReranker.get_raptorx_matrix(complex_id, desired_shape=(receptor_len, ligand_len))
        # Give each complex it's score
        ranks = []
        for i, res_complex in enumerate(complexes):
            neighbours = res_complex.get_neighbours_residues()
            score = float(RaptorxReranker.get_raptorx_score(rapt_mat, neighbours, self.scoring_method, self.prob_trim))
            ranks.append((i, score))
        ranks = sorted(ranks, key=lambda rank: rank[1], reverse=True)
        if detailed_return:
            return [(complexes[i], i, score) for i, score in ranks]
        return [complexes[i] for i, score in ranks]


# NOTE! complexes must be a subset of ACCEPTED_COMPLEXES
class SvmRegressionReranker(Reranker):

    def __init__(self, use_raptor=True):
        self._use_raptor = use_raptor
        self._classifier = NonZeroFnatClassifier()
        self._regressor = FnatRegressor()

    def rerank(self, complexes, detailed_return=False):
        complex_ids = np.unique([c.complex_id for c in complexes])
        assert len(complex_ids) == 1

        self._train_classifier(complex_ids, use_raptor=self._use_raptor)
        self._train_regressor(complex_ids, use_raptor=self._use_raptor)


        features = [[get_patch_dock_complex_features(c, include_raptor_score=self._use_raptor)] for c in complexes]
        print("feature len", len(features[0][0]))
        scores = [predict_fnat_from_features(f, self._classifier, self._regressor) for f in features]

        if detailed_return:
            return [(complexes[i], i, scores[i]) for i in np.argsort(scores)][::-1]
        return [complexes[i] for i in np.argsort(scores)][::-1]

    # PRIVATE METHODS:
    def _unison_shuffle(self, X, y):
        assert len(X) == len(y)
        p = np.random.permutation(len(X))
        return X[p], y[p]

    def _train_classifier(self, complex_ids, use_raptor=True):
        X_train, y_train, X_test, y_test = self._get_train_test_data(complex_ids, binary_labels=True, use_raptor=use_raptor)
        self._classifier.fit(X_train, y_train)

    def _train_regressor(self, complex_ids, use_raptor=True):
        X_train, y_train, X_test, y_test = self._get_train_test_data(complex_ids, binary_labels=False, use_raptor=use_raptor)
        self._regressor.fit(X_train, y_train)
        # print("regressor params:", self._regressor._regressor.coef_)

    def _get_train_test_data(self, complex_ids, binary_labels=True, use_raptor=True):        
        if binary_labels:
            X,y = load_features_and_binary_labels()
        else:
            X,y = load_features_and_continuous_labels(non_zero_data_only=False)

        if not use_raptor:
            X = self._remove_raptor_components(X)

        can_train_on = np.repeat(~np.isin(ACCEPTED_COMPLEXES, complex_ids), NUMBER_OF_TRANSFORMATIONS_PER_COMPLEX)
        X_train, y_train = self._unison_shuffle(X[can_train_on], y[can_train_on])
        X_test, y_test = self._unison_shuffle(X[~can_train_on], y[~can_train_on])

        if not binary_labels:
            has_positive_target = np.greater(y_train, 0)
            X_train, y_train = X_train[has_positive_target], y_train[has_positive_target]

        return X_train, y_train, X_test, y_test
    
    def _remove_raptor_components(self, X):
        return X[:, :N_PATCH_DOCK_SCORE_COMPONENTS]

