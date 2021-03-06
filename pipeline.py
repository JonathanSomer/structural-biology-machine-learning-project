from Constants import *
from Reranker.reranker import RaptorxReranker, RaptorXScoringMethod
from objects.results_helper import ResultsHelper, ComplexHelper
from visualizer import evodock_plot as plot

# from Reranker.regression_model import *


complex_ids = BETWEEN_10_TO_40_ACCEPTED
n_of_patchdock_results = 1000
reranker = RaptorxReranker(scoring_method=RaptorXScoringMethod.log_len_mul_sum, prob_trim=0.01, percentile_trim=0.0)

helper = ResultsHelper(complex_ids, n_of_patchdock_results, reranker, ignore_failure=True, verbose=False)
plot.plot_avg_capri_per_fnat_thresholds_before_vs_after(helper, top=[10, 20])
plot.plot_improved_fnat_category_percentage(helper)
plot.plot_improved_capri_percentage(helper, top=[10, 20])
plot.plot_at_least_one_percentage_per_fnat_thresholds_before_vs_after(helper, top=[10, 20])
'''
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
'''
