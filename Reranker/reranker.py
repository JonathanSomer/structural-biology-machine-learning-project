from objects import complex
from utils import pdb_utils, raptorx_utils


def raptorx_reranker(complexes, scoring_method, prob_trim=0.2):
    # type: (List[complex.Complex], raptorx_utils.RaptorXScoringMethod, float) -> List[complex.Complex]
    """
    Reranks a list of Complexes by raptorx matrix and given scoring parameters.
    :param complexes: list of complexes to rerank
    :param scoring_method: See utils.raptorx_utils.get_raptorx_score for details
    :param prob_trim: See utils.raptorx_utils.get_raptorx_score for details
    :return: reordered list by the score calculated
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
        neighbours = res_complex.get_neighbouring_residues()
        score = raptorx_utils.get_raptorx_score(rapt_mat, neighbours, scoring_method, prob_trim)
        ranks.append((i, score))
    sorted(ranks, key=lambda rank: rank[1])
    return [complexes[i] for i, score in ranks]
