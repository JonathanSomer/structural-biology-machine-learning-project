from Reranker.reranker import RaptorxReranker, RaptorXScoringMethod
from objects.pipeline_handler import ResultsHelper

train = ['1ACB', '1AVX', '1CGI', '1FQ1', '1FSK', '1I9R', '1IJK', '1IQD', '1KKL', '1KXQ', '1M10', '1MLC', '1NCA', '1QFW', '1R6Q', '1WEJ', '1ZM4', '2FD6', '2I25', '2JEL', '2VXT', '2W9E', '3EOA', '3HMX', '3L5W', '3MXW', '3RVW', '3V6Z', '4DN4', '4G6J', '4G6M']
with_problem = ['2Z0E', '2NZ8', '1JK9', '1F6M', '1K4C', '1JIW', '1NW9', '1BJ1', '1VFB']

complex_ids = ['1ACB', '1CGI']
n_of_patchdock_results = 200
reranker = RaptorxReranker(scoring_method=RaptorXScoringMethod.log_likelihood,
                           prob_trim=0.01)

helper = ResultsHelper(complex_ids, n_of_patchdock_results, reranker)

original_score = helper.get_all_capri_scores_of_original_patchdock_ranking()
rerank_score = helper.get_all_capri_scores_of_reranking()

print("original ranking capri score: {}".format(original_score))
print("rerank capri score: {}".format(rerank_score))


fnat_before_ranking = helper.get_all_fnat_scores_before_reranking()
fnat_after_ranking = helper.get_all_fnat_scores_after_reranking()

print("fnat scores before reranking: {}".format(fnat_before_ranking))
print("fnat scores after reranking: {}".format(fnat_after_ranking))

raptor_scores = helper.get_all_ranked_expectation_scores()

print("raptorx scores: {}".format(raptor_scores))