from enum import Enum
import os
import numpy as np
import Constants


class RaptorXScoringMethod(Enum):
    """
    Enum for various scoring methods
    """
    log_likelihood = 1  # sum of log probablities
    likelihood = 2  # multiply probabilities


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
        dirpath = Constants.get_raptorx_dir_path(complex_id)
        matrix_file_extension = '.gcnn_inter'
        filepath_list = [file_name for file_name in os.listdir(dirpath) if file_name.endswith(matrix_file_extension)]
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
    #assert method not in RaptorXScoringMethod.__members__, "method is not in RaptorXScoringMethod enum"
    # get new ndarray by list of incides
    neighbour_scores = raptorx_mat[tuple(zip(*neighbour_indices))]
    # trim results
    neighbour_scores = neighbour_scores[neighbour_scores >= trim]
    if method == RaptorXScoringMethod.log_likelihood:
        return np.sum(np.log(neighbour_scores))
    elif method == RaptorXScoringMethod.likelihood:
        return np.prod(neighbour_scores)
