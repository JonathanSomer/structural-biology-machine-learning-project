import csv
from Constants import FilePaths

def get_matrix_from_raptorx(receptor_id, ligand_id):
    filepath = FilePaths.get_raptorx_results_file_path(receptor_id, ligand_id)
    with open(filepath, 'rb') as f:
        matrix = []
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            row = [float(val) for val in row]
            matrix.append(row)
        return matrix
