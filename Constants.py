import os.path

project_dir = os.path.dirname(os.path.abspath(__file__))

BASE_DATA_PATH = r"C:\Users\catav\Documents\sim6\sadna\data\data_raptor_batch_with_raptor"


def get_zdock_benchmark_pdb_path(complex_id, ligand=True, bound=True):
    base_path = os.path.join(BASE_DATA_PATH, complex_id, "benchmark")
    file_name = "{}_{}_{}.pdb".format(complex_id, 'l' if ligand else 'r', 'b' if bound else 'u')
    return os.path.join(base_path, file_name)


def get_patchdock_ranked_complex_pdb_path(complex_id, rank):
    base_path = os.path.join(BASE_DATA_PATH, complex_id, "patch_dock")
    file_name = "{}.patch_dock_output.{}.pdb".format(complex_id, rank)
    return os.path.join(base_path, file_name)

def get_raptorx_dir_path(complex_id):
    return os.path.join(BASE_DATA_PATH, complex_id, "raptorx_results")

def get_patchdock_ranked_complex_cache_path(complex_id, rank):
    base_path = os.path.join(BASE_DATA_PATH, complex_id, "patch_dock")
    file_name = "{}.cache.{}.json".format(complex_id, rank)
    return os.path.join(base_path, file_name)

def get_zdock_benchmark_cache_path(complex_id, bound=True):
    base_path = os.path.join(BASE_DATA_PATH, complex_id, "benchmark")
    file_name = "{}_{}_cache.json".format(complex_id, 'b' if bound else 'u')
    return os.path.join(base_path, file_name)