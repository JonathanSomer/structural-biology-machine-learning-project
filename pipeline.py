from Reranker.reranker import RaptorxReranker
from utils.raptorx_utils import RaptorXScoringMethod
from objects.pipeline_handler import PipelineHandler

complex_id = '1ACB'
n_of_patchdock_results = 200
reranker = RaptorxReranker(scoring_method=RaptorXScoringMethod.log_likelihood,
                           prob_trim=0.01)

handler = PipelineHandler(complex_id, n_of_patchdock_results, reranker)

original_score = handler.get_capri_score_of_original_patchdock_ranking()
rerank_score = handler.get_capri_score_of_reranking()

print("original ranking capri score: {}".format(original_score))
print("rerank capri score: {}".format(rerank_score))


fnat_before_ranking = handler.get_fnat_scores_before_reranking()
fnat_after_ranking = handler.get_fnat_scores_after_reranking()

print("fnat scores before reranking: {}".format(fnat_before_ranking))
print("fnat scores after reranking: {}".format(fnat_after_ranking))

raptor_scores = handler.get_ranked_expectation_scores()

print("raptorx scores: {}".format(raptor_scores))