import os.path

project_dir = os.path.dirname(os.path.abspath(__file__))

TOP_RESULTS_COUNTS_FOR_CAPRI = 10

BASE_DATA_PATH = os.path.join(project_dir, "data")

# ids for complexes with all data ready

ACCEPTED_COMPLEXES = ['1CGI', '1KXQ', '1ACB', '1AVX', '1M10', '2W9E']
TOP_RAPTOR_CORRELATION_IDS = ['1R6Q', '1ACB', '1CGI', '1IQD', '4DN4', '3RVW', '3V6Z', '1KXQ']

TRAIN_COMPLEX_IDS = ['1ACB', '1AVX', '1CGI', '1FSK', '1I9R', '1IJK', '1IQD', '1KXQ', '1M10', '1NCA',
                     '1QFW', '1R6Q', '2FD6', '2I25', '2JEL', '2VXT', '2W9E', '3EOA', '3HMX', '3L5W',
                     '3MXW', '3RVW', '3V6Z', '4DN4', '4G6M', '2Z0E', '1F6M', '1K4C', '1JIW', '1NW9',
                     '1BJ1', '1WEJ']

IDS_WITH_PROBLEM = ['2NZ8', '1JK9', '1MLC']
TEST_COMPLEX_IDS = ['4G6J', '1FQ1', '1ZM4']
NUMBER_OF_TRANSFORMATIONS_PER_COMPLEX = 200


def get_zdock_benchmark_pdb_path(complex_id, ligand=True, bound=True):
    base_path = os.path.join(BASE_DATA_PATH, complex_id, "benchmark")
    file_name = "{}_{}_{}.pdb".format(complex_id, 'l' if ligand else 'r', 'b' if bound else 'u')
    return os.path.join(base_path, file_name)


def get_patchdock_ranked_complex_pdb_path(complex_id, rank, ligand):
    base_path = os.path.join(BASE_DATA_PATH, complex_id, "patch_dock", "processed")
    file_name = "{}.patch_dock_output.{}.{}.pdb".format(complex_id, rank, 'l' if ligand else 'r')
    return os.path.join(base_path, file_name)


def get_patchdock_complex_score_file_path(complex_id):
    base_path = os.path.join(BASE_DATA_PATH, complex_id, "patch_dock")
    file_name = "{}.patch_dock_output".format(complex_id)
    return os.path.join(base_path, file_name)


def get_raptorx_dir_path(complex_id):
    return os.path.join(BASE_DATA_PATH, complex_id, "raptorx_results")


def get_patchdock_ranked_complex_cache_path(complex_id, rank):
    base_path = os.path.join(BASE_DATA_PATH, complex_id, "cache")
    file_name = "{}.cache.{}.json".format(complex_id, rank)
    return os.path.join(base_path, file_name)


def get_zdock_benchmark_cache_path(complex_id, bound=True):
    base_path = os.path.join(BASE_DATA_PATH, complex_id, "cache")
    file_name = "{}_{}_cache.json".format(complex_id, 'b' if bound else 'u')
    return os.path.join(base_path, file_name)


def get_file_from_ml_models_path(file_name):
    return os.path.join(BASE_DATA_PATH, "ml_models", file_name)
