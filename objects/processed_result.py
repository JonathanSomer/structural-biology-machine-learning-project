from Constants import NEIGHBOUR_RADII, NEIGHBOUR_ATTR_TEMPLATE, NEIGHBOUR_DEFAULT_RADIUS, FNAT_ATTR_TEMPLATE, get_raptorx_dir_path
from os import listdir, path
import numpy as np

class ComplexProcessedResult(object):

    def __init__(self, complex_id, receptor_sequence, ligand_sequence, attr_dict,
                 neighbor_radius=NEIGHBOUR_DEFAULT_RADIUS):
        self.__dict__ = {"_" + k: v for k, v in attr_dict.items()}
        self._complex_id = complex_id
        self._receptor_sequence = receptor_sequence
        self._ligand_sequence = ligand_sequence
        if neighbor_radius not in NEIGHBOUR_RADII:
            raise ValueError("radius value is invalid")
        self.neighbor_radius = neighbor_radius

    @property
    def complex_id(self):
        return self._complex_id

    @property
    def receptor_sequence(self):
        return self._receptor_sequence

    @property
    def ligand_sequence(self):
        return self._ligand_sequence

    def get_neighbours_residues(self):
        # type: () -> List[Tuple[int, int]]
        """
        :param: neighbor_radius - neighbour radius to use
        :return: list of tuples (receptor_residue_index, ligand_residue_index) in which the euclidean distance
                 between them is at most $(NEIGHBOURS_RADIUS)
        """
        return self._get_attr_by_template(NEIGHBOUR_ATTR_TEMPLATE)

    def get_fnat_score(self):
        # type: () -> int
        """
        Returns the fnat score for the given radius
        :param: neighbor_radius - neighbour radius to use
        :return: fnat score
        """
        return self._get_attr_by_template(FNAT_ATTR_TEMPLATE)

    def _get_attr_by_template(self, attr_template):
        rad_attr = attr_template % self.neighbor_radius
        return getattr(self, rad_attr)


    def _get_raptorx_matrix(self):
        # type: (ComplexProcessedResult, str, bool) -> np.ndarray
        """
        Returns an numpy 2D ndarray of the score raptorx matrix with the given arguments.
        The first dimension is the receptor sequence and the second the ligand
        Can be given a complex id to look for or a file path of the matrix.
        :return: raptorx matrix
        """
        desired_shape = (len(self.receptor_sequence), len(self.ligand_sequence))
        # convert complex_id to file path
        dirpath = get_raptorx_dir_path(self.complex_id)
        matrix_file_extension = '.gcnn_inter'
        filepath_list = [file_name for file_name in listdir(dirpath) if
                         file_name.endswith(matrix_file_extension)]
        if len(filepath_list) == 0:
            raise IOError("Couldn't find any *%s file in %s" % (matrix_file_extension, dirpath))
        if len(filepath_list) > 1:
            raise IOError("More then one *%s file exists in %s" % (matrix_file_extension, dirpath))
        filepath = path.join(dirpath, filepath_list[0])
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

    def _get_values_in_percentiles_range(self, arr, low_percentile, upper_percentile, include_upper=False):
        '''
        :param arr: array or list of number
        :param low_percentile: low bound percentile for range
        :param upper_percentile: upper bound percentile for range
        :param include_upper: if true upper bound included in range else not included
        :return: np array of numbers within that range. lower bound always included in range
        '''
        arr = np.array(arr)
        low_percentile_val = np.percentile(arr, low_percentile, 0)
        upper_percentile_val = np.percentile(arr, upper_percentile, 0)
        if include_upper:
            return arr[np.logical_and(arr >= low_percentile_val, arr <= upper_percentile_val)]
        return arr[np.logical_and(arr >= low_percentile_val, arr < upper_percentile_val)]


    def get_raptor_percentiles_vector(self, func=np.sum):
        '''
        :param func: func to apply over each percentile range
        :return: array of func result over each of the raptor nb values of percentile ranges
        '''
        # ranges = [(0, 36, False), (36, 68, False), (68, 84, False), (84, 92, False), (92, 96, False),
        #          (96, 98, False), (98, 100, True)]
        ranges = [(0, 38, False), (38, 70, False), (70, 86, False), (86, 94, False), (94, 98, False), (98, 100, True)]
        raptor_nb_vals = self._get_raptorx_matrix()[tuple(zip(*self.get_neighbours_residues()))]
        return [func(self._get_values_in_percentiles_range(raptor_nb_vals, low_percentile, upper_percentile, include_upper))
                for low_percentile, upper_percentile, include_upper in ranges]