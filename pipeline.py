from Reranker.reranker import RaptorxReranker, RaptorXScoringMethod
from objects.pipeline_handler import ResultsHelper


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