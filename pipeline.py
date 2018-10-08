from Constants import *
from Reranker.reranker import RaptorxReranker, RaptorXScoringMethod
from objects.results_helper import ResultsHelper, ComplexHelper
from visualizer import evodock_plot as plot
# from Reranker.regression_model import *


complex_ids = get_all_training_complexes()
n_of_patchdock_results = 1000
reranker = RaptorxReranker(scoring_method=RaptorXScoringMethod.sqrt_sum, prob_trim=0.01, percentile_trim=0.8)

helper = ResultsHelper(complex_ids, n_of_patchdock_results, reranker, ignore_failure=True, verbose=False)

original_score = helper.get_all_capri_scores(after=False)
rerank_score = helper.get_all_capri_scores(after=True)

print("original ranking capri score: {}".format(original_score))
print("rerank capri score: {}".format(rerank_score))


fnat_before_ranking = helper.get_all_fnat_scores(after=False)
fnat_after_ranking = helper.get_all_fnat_scores(after=True)

# print("fnat scores before reranking: {}".format(fnat_before_ranking))
# print("fnat scores after reranking: {}".format(fnat_after_ranking))

raptor_scores = helper.get_all_ranked_expectation_scores()

# print("raptorx scores: {}".format(raptor_scores))
