from Constants import *
from objects.complex import *
import numpy as np
from utils.raptorx_utils import *
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from utils.fnat_utils import *

warnings.simplefilter('ignore', PDBConstructionWarning)


def get_features_and_labels(use_training_data=True):
	
	X = []
	y = []

	ids = TRAIN_COMPLEX_IDS if use_training_data  else TEST_COMPLEX_IDS

	for complex_id in ids:

		benchmark_complex = BenchmarkComplex(complex_id, type=ComplexType.zdock_benchmark_bound)

		for rank in range(1, NUMBER_OF_TRANSFORMATIONS_PER_COMPLEX + 1):
			pd_complex = PatchDockComplex(complex_id, rank)

			features = pd_complex.score_components
			
			raptor_matrix = get_raptorx_matrix(complex_id)
			neighbor_indexes = pd_complex.get_neighbours_residues()

			for method_idx in range(1, len(RaptorXScoringMethod) + 1):
				for trim in [0.01, 0.05, 0.1]:
					features.append(get_raptorx_score(raptor_matrix, neighbor_indexes, RaptorXScoringMethod(method_idx), trim))
			

			target = get_fnat_score(pd_complex, benchmark_complex)

			X.append(features)
			y.append(target)
			print(complex_id, rank, features, "---->", target)

	return X, y


