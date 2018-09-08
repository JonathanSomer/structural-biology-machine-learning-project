from Constants import *
from Reranker.reranker import RaptorxReranker, RaptorXScoringMethod, SvmRegressionReranker
from objects.results_helper import ResultsHelper
from visualizer import evodock_plot as plot
from Reranker.regression_model import *


# complex_ids = ['1CGI'] #, '1ACB'
# complex_ids = ACCEPTED_COMPLEXES
complex_ids = TOP_RAPTOR_CORRELATION_IDS
n_of_patchdock_results = 200
# reranker = RaptorxReranker(scoring_method=RaptorXScoringMethod.log_likelihood, prob_trim=0.01)
reranker = SvmRegressionReranker(use_raptor=True, training_data_complex_ids=TOP_RAPTOR_CORRELATION_IDS, training_data_file_name=TOP_RAPTOR_CORRELATION_PICKLE)

helper = ResultsHelper(complex_ids, n_of_patchdock_results, reranker)

original_score = helper.get_all_capri_scores_of_original_patchdock_ranking()
rerank_score = helper.get_all_capri_scores_of_reranking()

print("original ranking capri score: {}".format(original_score))
print("rerank capri score: {}".format(rerank_score))


fnat_before_ranking = helper.get_all_fnat_scores(after=False)
fnat_after_ranking = helper.get_all_fnat_scores(after=True)

print("fnat scores before reranking: {}".format(fnat_before_ranking))
print("fnat scores after reranking: {}".format(fnat_after_ranking))

raptor_scores = helper.get_all_ranked_expectation_scores()

print("raptorx scores: {}".format(raptor_scores))
