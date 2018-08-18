import objects.complex as complex
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

warnings.simplefilter('ignore', PDBConstructionWarning)

complex_id = '1ACB'

complex_bound = complex.BenchmarkComplex(complex_id)
complex_unbound = complex.BenchmarkComplex(complex_id, complex.ComplexType.zdock_benchmark_unbound)

# TODO: use complex_unbound to extract chains for ligand and receptor
ligand_chains, receptor_chains = None, None

complex_ranking = []
for i in range(1, 200):
	print(str(i))
	patch_dock_complex = complex.PatchDockComplex(complex_id, i)
	complex_ranking.append(patch_dock_complex)

# simple test for capri score
print(complex_ranking[0].capri_score)

# TODO: implement score_rank using DockQ~!
# TODO: rename score_rank for a more informative name
overall_preraptor_score = capri_score_ranking(complex_ranking)

# TODO: implement raptor_rerank and rename it!
complex_predictions_post_raptor = raptor_rerank(complex_id, complex_ranking)

overall_postraptor_score = capri_score_ranking(complex_predictions_post_raptor)

print(compare_score(overall_preraptor_score, overall_postraptor_score))
