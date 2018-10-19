# from Constants import *
# from objects.complex import *
# import numpy as np
# from utils.raptorx_utils import *
# import warnings
# from Bio.PDB.PDBExceptions import PDBConstructionWarning
# from utils.fnat_utils import *
# from sklearn.linear_model import Ridge
# from sklearn import svm
# import pickle
# from reranker import Reranker

# warnings.simplefilter('ignore', PDBConstructionWarning)

# TRAIN_FEATURES_AND_LABELS_PICKLE_5 = 'features_and_labels_radius_5.pickle'
# TRAIN_FEATURES_AND_LABELS_PICKLE_8 = 'features_and_labels_radius_8.pickle'
# TOP_RAPTOR_CORRELATION_PICKLE = 'top_8_raptor_correlation.pickle'
# NON_ZERO_FNAT_CLASSIFIER_PICKLE = 'non_zero_fnat_classifier.pickle'
# FNAT_REGRESSOR_PICKLE = 'fnat_regressor.pickle'

# # from Reranker.regression_model import *
# # X,y =  get_features_and_labels(complex_ids=TOP_RAPTOR_CORRELATION_IDS, file_name=TOP_RAPTOR_CORRELATION_PICKLE)
# def get_features_and_labels(pickle_data=True, use_8_angstrom=True, re_cache=True, complex_ids=get_all_training_complexes(), file_name=TRAIN_FEATURES_AND_LABELS_PICKLE_8):

# 	X = []
# 	y = []

# 	for complex_id in complex_ids:

# 		benchmark_complex = BenchmarkComplex(complex_id, type=ComplexType.zdock_benchmark_bound, reprocess=re_cache)

# 		for rank in range(1, NUMBER_OF_TRANSFORMATIONS_PER_COMPLEX + 1):
# 			patch_dock_complex = PatchDockComplex(complex_id, rank, reprocess=re_cache)

# 			features = get_patch_dock_complex_features(patch_dock_complex)
# 			target = get_fnat_score(patch_dock_complex, benchmark_complex)

# 			X.append(features)
# 			y.append(target)
# 			print(complex_id, target)

# 	if pickle_data:
# 		data = { 'X' : X, 'y' : y }
# 		with open(get_file_from_ml_models_path(file_name), "wb") as f:
# 			pickle.dump(data, f)

# 	return X, y

# def get_patch_dock_complex_features(patch_dock_complex, include_raptor_score=True):
# 			features = patch_dock_complex.score_components
			
# 			raptor_matrix = get_raptorx_matrix(patch_dock_complex.complex_id)
# 			neighbor_indexes = patch_dock_complex.get_neighbours_residues()

# 			if include_raptor_score and len(features) <= N_PATCH_DOCK_SCORE_COMPONENTS:
# 				for method_idx in range(1, len(RaptorXScoringMethod) + 1):
# 					for trim in [0.01, 0.05, 0.1]:
# 						features.append(get_raptorx_score(raptor_matrix, neighbor_indexes, RaptorXScoringMethod(method_idx), trim))
# 			return features

# def load_features_and_continuous_labels(non_zero_data_only=True, file_name=TRAIN_FEATURES_AND_LABELS_PICKLE_8):
# 	with open(get_file_from_ml_models_path(file_name), "rb") as f:
# 		data = pickle.load(f)
# 		X = np.array(data['X'])
# 		y = np.array(data['y'])
# 		if non_zero_data_only:
# 			mask = ~np.equal(y, 0.0)
# 			X = X[mask]
# 			y = y[mask]

# 		return X, y 

# def load_features_and_binary_labels(file_name=TRAIN_FEATURES_AND_LABELS_PICKLE_8):
# 	with open(get_file_from_ml_models_path(file_name), "rb") as f:
# 		data = pickle.load(f)
# 		X = np.array(data['X'])
# 		y = ~np.equal(np.array(data['y']), 0.0)
# 		return X, y

# def predict_fnat_from_features(complex_features, non_zero_fnat_classifier, fnat_regressor):
# 	if non_zero_fnat_classifier.predict(complex_features):
# 		return fnat_regressor.predict(complex_features)
# 	else:
# 		return 0.0


# def _get_values_in_percentiles_range(arr, low_percentile, upper_percentile, include_upper=False):
# 	'''
# 	:param arr: array or list of number
# 	:param low_percentile: low bound percentile for range
# 	:param upper_percentile: upper bound percentile for range
# 	:param include_upper: if true upper bound included in range else not included
# 	:return: np array of numbers within that range. lower bound always included in range
# 	'''
# 	arr = np.array(arr)
# 	low_percentile_val = np.percentile(arr, low_percentile, 0)
# 	upper_percentile_val = np.percentile(arr, upper_percentile, 0)
# 	if include_upper:
# 		return arr[np.logical_and(arr >= low_percentile_val, arr <= upper_percentile_val)]
# 	return arr[np.logical_and(arr >= low_percentile_val, arr < upper_percentile_val)]

# def _get_raptor_percentiles_vector(raptor_vals, func=np.sum):
# 	'''
# 	:param raptor_vals: list of raptor scores of nb of specific result
# 	:param func: func to apply over each percentile range
# 	:return: array of func result over each of the percentile ranges
# 	'''
# 	# ranges = [(0, 36, False), (36, 68, False), (68, 84, False), (84, 92, False), (92, 96, False),
# 	#          (96, 98, False), (98, 100, True)]
# 	ranges = [(0, 38, False), (38, 70, False), (70, 86, False), (86, 94, False), (94, 98, False), (98, 100, True)]
# 	return [func(_get_values_in_percentiles_range(raptor_vals, low_percentile, upper_percentile, include_upper))
# 			for low_percentile, upper_percentile, include_upper in ranges]

# # NOTE! complexes must be a subset of ACCEPTED_COMPLEXES
# class SvmRegressionReranker(Reranker):

#     # alternative args:
#     # training_data_file_name=ACCEPTED_COMPLEXES, training_data_complex_ids=TRAIN_FEATURES_AND_LABELS_PICKLE_8
#     def __init__(self, use_raptor=True, training_data_file_name=TOP_RAPTOR_CORRELATION_PICKLE,
#                  training_data_complex_ids=TOP_RAPTOR_CORRELATION_IDS):
#         self._training_data_file_name = training_data_file_name
#         self._training_data_complex_ids = training_data_complex_ids
#         self._use_raptor = use_raptor
#         self._classifier = NonZeroFnatClassifier()
#         self._regressor = FnatRegressor()

#     def rerank(self, complexes, detailed_return=False):
#         complex_ids = np.unique([c.complex_id for c in complexes])
#         assert len(complex_ids) == 1

#         self._train_classifier(complex_ids, use_raptor=self._use_raptor)
#         self._train_regressor(complex_ids, use_raptor=self._use_raptor)

#         features = [[get_patch_dock_complex_features(c, include_raptor_score=self._use_raptor)] for c in complexes]
#         # print("feature len", len(features[0][0]))
#         scores = [predict_fnat_from_features(f, self._classifier, self._regressor) for f in features]

#         if detailed_return:
#             return [(complexes[i], i, scores[i]) for i in np.argsort(scores)][::-1]
#         return [complexes[i] for i in np.argsort(scores)][::-1]

#     # PRIVATE METHODS:
#     def _unison_shuffle(self, X, y):
#         assert len(X) == len(y)
#         p = np.random.permutation(len(X))
#         return X[p], y[p]

#     def _train_classifier(self, complex_ids, use_raptor=True):
#         X_train, y_train, X_test, y_test = self._get_train_test_data(complex_ids, binary_labels=True,
#                                                                      use_raptor=use_raptor)
#         self._classifier.fit(X_train, y_train)

#     def _train_regressor(self, complex_ids, use_raptor=True):
#         X_train, y_train, X_test, y_test = self._get_train_test_data(complex_ids, binary_labels=False,
#                                                                      use_raptor=use_raptor)
#         self._regressor.fit(X_train, y_train)
#         # print("regressor params:", self._regressor._regressor.coef_)

#     def _get_train_test_data(self, complex_ids, binary_labels=True, use_raptor=True):
#         if binary_labels:
#             X, y = load_features_and_binary_labels(file_name=self._training_data_file_name)
#         else:
#             X, y = load_features_and_continuous_labels(non_zero_data_only=False,
#                                                        file_name=self._training_data_file_name)

#         if not use_raptor:
#             X = self._remove_raptor_components(X)

#         can_train_on = np.repeat(~np.isin(self._training_data_complex_ids, complex_ids),
#                                  NUMBER_OF_TRANSFORMATIONS_PER_COMPLEX)
#         X_train, y_train = self._unison_shuffle(X[can_train_on], y[can_train_on])
#         X_test, y_test = self._unison_shuffle(X[~can_train_on], y[~can_train_on])

#         if not binary_labels:
#             has_positive_target = np.greater(y_train, 0)
#             X_train, y_train = X_train[has_positive_target], y_train[has_positive_target]

#         return X_train, y_train, X_test, y_test

#     def _remove_raptor_components(self, X):
#         return X[:, :N_PATCH_DOCK_SCORE_COMPONENTS]
