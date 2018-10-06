from Constants import *
from objects.complex import *
import numpy as np
from utils.raptorx_utils import *
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from utils.fnat_utils import *
from sklearn.linear_model import Ridge
from sklearn import svm
import pickle

warnings.simplefilter('ignore', PDBConstructionWarning)

TRAIN_FEATURES_AND_LABELS_PICKLE_5 = 'features_and_labels_radius_5.pickle'
TRAIN_FEATURES_AND_LABELS_PICKLE_8 = 'features_and_labels_radius_8.pickle'
TOP_RAPTOR_CORRELATION_PICKLE = 'top_8_raptor_correlation.pickle'
NON_ZERO_FNAT_CLASSIFIER_PICKLE = 'non_zero_fnat_classifier.pickle'
FNAT_REGRESSOR_PICKLE = 'fnat_regressor.pickle'

# from Reranker.regression_model import *
# X,y =  get_features_and_labels(complex_ids=TOP_RAPTOR_CORRELATION_IDS, file_name=TOP_RAPTOR_CORRELATION_PICKLE)
def get_features_and_labels(pickle_data=True, use_8_angstrom=True, re_cache=True, complex_ids=ACCEPTED_COMPLEXES, file_name=TRAIN_FEATURES_AND_LABELS_PICKLE_8):

	X = []
	y = []

	for complex_id in complex_ids:

		benchmark_complex = BenchmarkComplex(complex_id, type=ComplexType.zdock_benchmark_bound, reprocess=re_cache)

		for rank in range(1, NUMBER_OF_TRANSFORMATIONS_PER_COMPLEX + 1):
			patch_dock_complex = PatchDockComplex(complex_id, rank, reprocess=re_cache)

			features = get_patch_dock_complex_features(patch_dock_complex)
			target = get_fnat_score(patch_dock_complex, benchmark_complex)

			X.append(features)
			y.append(target)
			print(complex_id, target)

	if pickle_data:
		data = { 'X' : X, 'y' : y }
		with open(get_file_from_ml_models_path(file_name), "wb") as f:
			pickle.dump(data, f)

	return X, y

def get_patch_dock_complex_features(patch_dock_complex, include_raptor_score=True):
			features = patch_dock_complex.score_components
			
			raptor_matrix = get_raptorx_matrix(patch_dock_complex.complex_id)
			neighbor_indexes = patch_dock_complex.get_neighbours_residues()

			if include_raptor_score and len(features) <= N_PATCH_DOCK_SCORE_COMPONENTS:
				for method_idx in range(1, len(RaptorXScoringMethod) + 1):
					for trim in [0.01, 0.05, 0.1]:
						features.append(get_raptorx_score(raptor_matrix, neighbor_indexes, RaptorXScoringMethod(method_idx), trim))
			return features

def load_features_and_continuous_labels(non_zero_data_only=True, file_name=TRAIN_FEATURES_AND_LABELS_PICKLE_8):
	with open(get_file_from_ml_models_path(file_name), "rb") as f:
		data = pickle.load(f)
		X = np.array(data['X'])
		y = np.array(data['y'])
		if non_zero_data_only:
			mask = ~np.equal(y, 0.0)
			X = X[mask]
			y = y[mask]

		return X, y 

def load_features_and_binary_labels(file_name=TRAIN_FEATURES_AND_LABELS_PICKLE_8):
	with open(get_file_from_ml_models_path(file_name), "rb") as f:
		data = pickle.load(f)
		X = np.array(data['X'])
		y = ~np.equal(np.array(data['y']), 0.0)
		return X, y

class NonZeroFnatClassifier(object):
	def __init__(self):
		try:
			with open(get_file_from_ml_models_path(NON_ZERO_FNAT_CLASSIFIER_PICKLE), "rb") as f:
				self._clf = pickle.load(f)
		except:
			self._clf = svm.SVC(C=2.0, kernel='rbf')

	def fit(self, X, y):
		self._clf.fit(X, y)
		with open(get_file_from_ml_models_path(NON_ZERO_FNAT_CLASSIFIER_PICKLE), "wb") as f:
			pickle.dump(self._clf, f)


	def predict(self, X):
		return self._clf.predict(X)

	def score(self, X, y):
		return self._clf.score(X, y)

class FnatRegressor(object):
	def __init__(self):
		try:
			with open(get_file_from_ml_models_path(FNAT_REGRESSOR_PICKLE), "rb") as f:
				self._regressor = pickle.load(f)
		except:
			self._regressor = Ridge(alpha=4.0)

	def fit(self, X, y):
		self._regressor.fit(X, y)
		with open(get_file_from_ml_models_path(FNAT_REGRESSOR_PICKLE), "wb") as f:
			pickle.dump(self._regressor, f)


	def predict(self, X):
		return self._regressor.predict(X)

	def score(self, X, y):
		return self._regressor.score(X, y)

def predict_fnat_from_features(complex_features, non_zero_fnat_classifier, fnat_regressor):
	if non_zero_fnat_classifier.predict(complex_features):
		return fnat_regressor.predict(complex_features)
	else:
		return 0.0
