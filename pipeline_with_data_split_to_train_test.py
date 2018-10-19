from Constants import *
from Reranker.regression_reranker import SvmRegressionReranker
from objects.results_helper import ResultsHelper, ComplexHelper
from visualizer import evodock_plot as plot

TRAIN_SUBSET_RATIO = 0.8
N_TRANSFORMATIONS_PER_COMPLEX = 1000


# all_complex_ids = get_all_training_complexes()[:10]
all_complex_ids = get_all_training_complexes()

n_complexes_in_train_set = round(float(len(all_complex_ids))*TRAIN_SUBSET_RATIO)
train_complex_ids = all_complex_ids[:n_complexes_in_train_set]
test_complex_ids = all_complex_ids[n_complexes_in_train_set:]

reranker = SvmRegressionReranker(train_complex_ids)

helper = ResultsHelper(test_complex_ids, N_TRANSFORMATIONS_PER_COMPLEX, reranker, ignore_failure=True, verbose=False)
plot.plot_improved_fnat_category_percentage(helper)
plot.plot_improved_capri_percentage(helper, top=[10, 20, 50, 100])
plot.plot_at_least_one_percentage_per_fnat_thresholds_before_vs_after(helper, top=[10, 20, 50, 100])
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
