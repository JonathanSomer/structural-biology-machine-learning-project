from Constants import *
import json
import progressbar
from Reranker.complex_processed_result_generator import *

def main():
	generator = ComplexProcessedResultGenerator()

	map_complex_id_to_X = {}
	map_complex_id_to_y = {}

	for complex_id in progressbar.progressbar(get_all_training_complexes()):
		complexes = generator.generate(complex_id, 1000)
		
		map_complex_id_to_X[complex_id] = []
		map_complex_id_to_y[complex_id] = []

		for complex in complexes:
			try:
				map_complex_id_to_X[complex_id].append([complex.get_patch_dock_score()[0], complex.get_raptor_score()])
				map_complex_id_to_y[complex_id].append([complex.get_fnat_score()])
			except:
				continue

	with open(train_test_data_path(), "w") as file:
		json.dump({"map_complex_id_to_X" : map_complex_id_to_X, "map_complex_id_to_y" : map_complex_id_to_y}, file)

