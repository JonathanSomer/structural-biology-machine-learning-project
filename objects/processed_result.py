from os import listdir, path
import numpy as np
from Constants import *
import re

class ComplexProcessedResult(object):

    def __init__(self, complex_id, receptor_sequence, ligand_sequence, attr_dict, original_rank,
                 neighbor_radius=NEIGHBOUR_DEFAULT_RADIUS):
        self.__dict__ = {"_" + k: v for k, v in attr_dict.items()}
        self._complex_id = complex_id
        self._receptor_sequence = receptor_sequence
        self._ligand_sequence = ligand_sequence
        self._original_rank = original_rank
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

    @property
    def original_rank(self):
        return self._original_rank

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

    def get_patch_dock_score_components(self):
        with open(get_patchdock_complex_score_file_path(self.complex_id), "r") as f:
            for line in f:
                pattern = re.compile("^\s*" + str(self.original_rank) + "\s\|.+")
                if pattern.match(line):
                    components = [x.strip() for x in line.split('|')]
                    # rank | score | pen.  | Area    | as1   | as2   | as12  | ACE     | hydroph | Energy  |cluster| dist. || Ligand Transformation
                    indexes = [1, 2, 3, 7]
                    return [float(components[i]) for i in indexes]

        print("ERROR GETTING PATCH DOCK SCORE COMPONENTS --- NO REGEX MATCH WTF?")
        return []

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

    def get_nb_raptor_values(self):
        return self._get_raptorx_matrix()[tuple(zip(*self.get_neighbours_residues()))]

    def get_raptor_grouped_values_vector(self, func=np.sum):
        '''
        :param func: func to apply over each grouped raptor values
        :return: array of func result over each of the grouped raptor nb values
        '''
        group_slices = [slice(0, 2), slice(2, 6), slice(6, 14), slice(14, 30), slice(30, 62), slice(62, 1000)]
        sorted_nb_values = sorted(self.get_nb_raptor_values(), reverse=True)
        return [func(sorted_nb_values[group_slice])
                for group_slice in group_slices] + [len(sorted_nb_values)]


