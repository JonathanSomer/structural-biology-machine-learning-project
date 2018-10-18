from Constants import *
from objects.complex import *
import numpy as np
from utils.raptorx_utils import *
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from utils.fnat_utils import *
import pickle
from Reranker.reranker import Reranker
from Reranker.non_zero_fnat_classifier import *
from Reranker.fnat_regressor import *
from Reranker.complex_processed_result_generator import *

NUMBER_OF_TRANSFORMATIONS_PER_COMPLEX = 1000

class SvmRegressionReranker(Reranker):
    def __init__(self, train_complex_ids):
        self._classifier = NonZeroFnatClassifier()
        self._regressor = FnatRegressor()

        self._train_classifier(train_complex_ids)
        self._train_regressor(train_complex_ids)
        print("Finished SvmRegressionReranker __init__")


    def rerank(self, complexes, detailed_return=False):
        features = [[self._generate_features_for_processed_complex(c)] for c in complexes]
        # print("feature len", len(features[0][0]))
        scores = [self._predict_fnat_from_features(f) for f in features]

        if detailed_return:
            return [(complexes[i], i, scores[i]) for i in np.argsort(scores)][::-1]
        return [complexes[i] for i in np.argsort(scores)][::-1]

    # PRIVATE METHODS:

    def _train_classifier(self, complex_ids):
        X, y = self._get_X_train_y_train(complex_ids)
        y = ~np.equal(np.array(y), 0.0)
        self._classifier.fit(X, y)

    def _train_regressor(self, complex_ids):
        X, y = self._get_X_train_y_train(complex_ids)
        mask = ~np.equal(y, 0.0)
        X = X[mask]
        y = y[mask]
        self._regressor.fit(X, y)
        # print("regressor params:", self._regressor._regressor.coef_)

    def _get_X_train_y_train(self, complex_ids):
        X = []; y = []
        processed_complexes_generator = ComplexProcessedResultGenerator()
        for complex_id in complex_ids:
            processed_complexes = processed_complexes_generator.generate(complex_id, NUMBER_OF_TRANSFORMATIONS_PER_COMPLEX)
            for processed_complex in processed_complexes:
                X.append(self._generate_features_for_processed_complex(processed_complex))
                y.append(processed_complex.get_fnat_score())

        return self._unison_shuffle(np.array(X), np.array(y))

    def _generate_features_for_processed_complex(self, processed_complex):
        features = []
        features.extend(processed_complex.get_patch_dock_score_components())
        # TODO : raptor feature
        return features

    def _predict_fnat_from_features(self, complex_features):
        if self._classifier.predict(complex_features): # non zero fnat classification
            return self._regressor.predict(complex_features)
        else:
            return 0.0

    def _unison_shuffle(self, X, y):
        assert len(X) == len(y)
        p = np.random.permutation(len(X))
        return X[p], y[p]
