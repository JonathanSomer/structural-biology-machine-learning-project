from Constants import *
from Constants import NEIGHBOUR_DEFAULT_RADIUS
import json
from objects.processed_result import ComplexProcessedResult


class ComplexProcessedResultGenerator(object):
	
	def generate(self, complex_id, n_transformations, neighbor_radius=NEIGHBOUR_DEFAULT_RADIUS):
		ligand_sequence = self._get_sequence_from_fasta(complex_id, is_ligand_sequence=True)
		receptor_sequence = self._get_sequence_from_fasta(complex_id, is_ligand_sequence=False)

		with open(get_processed_data_json_path(complex_id)) as f:
			results_json = json.load(f)

		complex_processed_results = []
		for original_rank in range(1, n_transformations + 1):
			c = ComplexProcessedResult(complex_id, receptor_sequence, ligand_sequence, results_json[original_rank-1], original_rank, neighbor_radius)
			complex_processed_results.append(c)

		return complex_processed_results

	def _get_sequence_from_fasta(self, complex_id, is_ligand_sequence):
		with open(get_complex_fasta_path(complex_id, is_ligand_sequence)) as f:
			return f.readlines()[1]
