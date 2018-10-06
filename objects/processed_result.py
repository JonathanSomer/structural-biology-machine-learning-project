from Constants import NEIGHBOUR_RADII, NEIGHBOUR_ATTR_TEMPLATE, NEIGHBOUR_DEFAULT_RADIUS, FNAT_ATTR_TEMPLATE


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
