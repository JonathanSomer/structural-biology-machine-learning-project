import os.path

project_dir = os.path.dirname(os.path.abspath(__file__))

BASE_DATA_PATH = os.path.join(project_dir, r"data")


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
