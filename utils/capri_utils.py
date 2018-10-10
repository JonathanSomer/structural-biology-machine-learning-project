from enum import Enum
from utils.fnat_utils import get_fnat_score


class CapriAcceptanceLevelScores(Enum):
    High = 1000000000
    Medium = 1000000
    Acceptable = 1000
    Incorrect = 1


class FnatThresholds(Enum):
    High = 0.5
    Medium = 0.3
    Acceptable = 0.1


# this will indicate how many of the first 10 results is of each acceptance criteria
# thousands for high, hundreds for medium, tens for acceptance, and ones for incorrect
def get_capri_score(ranked_complexes, benchmark_complex):
    top_10 = ranked_complexes[:10]
    return sum(_get_capri_score_for_estimated_complex(estimated_complex, benchmark_complex)
               for estimated_complex in top_10)


def convert_fnat_to_capri(fnat):
    if fnat > FnatThresholds.High.value:
        return CapriAcceptanceLevelScores.High.value
    elif fnat > FnatThresholds.Medium.value:
        return CapriAcceptanceLevelScores.Medium.value
    elif fnat > FnatThresholds.Acceptable.value:
        return CapriAcceptanceLevelScores.Acceptable.value
    else:
        return CapriAcceptanceLevelScores.Incorrect.value


def _get_capri_score_for_estimated_complex(estimated_complex, benchmark_complex):
    fnat = get_fnat_score(estimated_complex, benchmark_complex)
    return convert_fnat_to_capri(fnat)
