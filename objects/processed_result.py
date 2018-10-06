from Constants import NEIGHBOUR_RADII, NEIGHBOUR_ATTR_TEMPLATE, NEIGHBOUR_DEFAULT_RADIUS, FNAT_ATTR_TEMPLATE


class ComplexProcessedResult(object):

    def __init__(self, complex_id, receptor_sequence, ligand_sequence, attr_dict):
        self.__dict__ = {"_" + k: v for k, v in attr_dict.items()}
        self._complex_id = complex_id
        self._receptor_sequence = receptor_sequence
        self._ligand_sequence = ligand_sequence

    @property
    def complex_id(self):
        return self._complex_id

    @property
    def receptor_sequence(self):
        return self._receptor_sequence

    @property
    def ligand_sequence(self):
        return self._ligand_sequence

    def get_neighbours_residues(self, neighbor_radius=NEIGHBOUR_DEFAULT_RADIUS):
        # type: (int) -> List[Tuple[int, int]]
        """
        :param: neighbor_radius - neighbour radius to use
        :return: list of tuples (receptor_residue_index, ligand_residue_index) in which the euclidean distance
                 between them is at most $(NEIGHBOURS_RADIUS)
        """
        return self._get_attr_by_radius(NEIGHBOUR_ATTR_TEMPLATE, neighbor_radius)

    def get_fnat_score(self, neighbor_radius=NEIGHBOUR_DEFAULT_RADIUS):
        # type: (int) -> int
        """
        Returns the fnat score for the given radius
        :param: neighbor_radius - neighbour radius to use
        :return: fnat score
        """
        return self._get_attr_by_radius(FNAT_ATTR_TEMPLATE, neighbor_radius)

    def _get_attr_by_radius(self, attr_template, neighbor_radius):
        if neighbor_radius not in NEIGHBOUR_RADII:
            raise ValueError("radius value is invalid")
        rad_attr = attr_template % neighbor_radius
        return getattr(self, rad_attr)
