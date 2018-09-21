import os.path

project_dir = os.path.dirname(os.path.abspath(__file__))

TOP_RESULTS_COUNTS_FOR_CAPRI = 10

BASE_DATA_PATH = os.path.join(project_dir, "data")

LIGAND_SEPARATOR_HEADER = "LIGAND SEPARATOR\n"
LIGAND_HEADER = "LIGAND\n"
RECEPTOR_HEADER = "RECEPTOR\n"


TOP_RAPTOR_CORRELATION_IDS = ['1R6Q', '1ACB', '1CGI', '1IQD', '4DN4', '3RVW', '3V6Z', '1KXQ']

RIGID_COMLEXES = ['1AHW', '1BVK', '1DQJ', '1E6J', '1JPS', '1MLC', '1VFB', '1WEJ', '2FD6', '2I25', '2VIS', '2VXT', '2W9E', '3EOA', '3HMX', '3MXW', '3RVW', '4DN4', '4FQI', '4G6J', '4G6M', '4GXU', '1BJ1', '1FSK', '1I9R', '1IQD', '1K4C', '1KXQ', '1NCA', '1NSN', '2JEL', '1AVX', '1AY7', '1BUH', '1BVN', '1CLV', '1D6R', '1DFJ', '1E6E', '1EAW', '1EWY', '1EZU', '1F34', '1F51', '1FLE', '1GL1', '1GLA', '1GXD', '1HIA', '1JTD', '1JTG', '1JWH', '1MAH', '1OC0', '1OPH', '1OYV', '1PPE', '1R0R', '1TMQ', '1UDI', '1US7', '1WDW', '1YVB', '1Z5Y', '2A1A', '2A9K', '2ABZ', '2AYO', '2B42', '2GAF', '2J0T', '2MTA', '2O8V', '2OOB', '2OOR', '2OUL', '2PCC', '2SIC', '2SNI', '2UUY', '2YVJ', '3A4S', '3K75', '3LVK', '3PC8', '3SGQ', '3VLB', '4CPA', '4H03', '4HX3', '7CEI', '1A2K', '1AK4', '1AKJ', '1AZS', '1E96', '1EFN', '1EXB', '1FCC', '1FFW', '1FQJ', '1GCQ', '1GHQ', '1GPW', '1H9D', '1HCF', '1HE1', '1I4D', '1J2J', '1K74', '1KAC', '1KLU', '1KTZ', '1KXP', '1M27', '1ML0', '1OFU', '1PVH', '1QA9', '1RLB', '1RV6', '1S1Q', '1SBB', '1T6B', '1XD3', '1XU1', '1Z0K', '1ZHH', '1ZHI', '2A5T', '2AJF', '2B4J', '2BTF', '2FJU', '2G77', '2GTP', '2HLE', '2HQS', '2VDB', '2X9A', '3BIW', '3BP8', '3D5S', '3H2V', '3P57', '4M76']
MEDIUM_COMPLEXES = ['3EO1', '3G6D', '3HI6', '3L5W', '3V6Z', '1CGI', '1IJK', '1JIW', '1KKL', '1M10', '1NW9', '1R6Q', '1ZM4', '2NZ8', '2Z0E', '4FZA', '4IZ7', '4LW4', '1B6C', '1FC2', '1GP2', '1GRN', '1HE8', '1I2M', '1IB1', '1K5D', '1LFD', '1MQ8', '1N2C', '1SYX', '1WQ1', '1XQS', '2CFH', '2H7V', '2HRK', '2OZA', '3AAA', '3AAD', '3BX7', '3CPH', '3DAW', '3R9A', '3S9D', '3SZK', '4JCV']
DIFFICULT_COMPLEXES = ['1BGX', '2HMI', '1ACB', '1F6M', '1FQ1', '1JK9', '1JMO', '1JZD', '1PXV', '1ZLI', '2IDO', '2O3B', '2OT3', '3FN1', '3H11', '4GAM', '1ATN', '1BKD', '1DE4', '1E4K', '1EER', '1FAK', '1H1V', '1IBR', '1IRA', '1R8S', '1RKE', '1Y64', '2C0L', '2I9B', '2J7P', '3AAD', '3F1P', '3L89']

ALL_PROCESSED_COMPLEXES = ['1A2K', '1ACB', '1AHW', '1AK4', '1AKJ', '1ATN', '1AVX', '1AY7', '1AZS', '1B6C', '1BGX', '1BJ1', '1BKD', '1BUH', '1BVK', '1BVN', '1CGI', '1CLV', '1D6R', '1DE4', '1DFJ', '1DQJ', '1E4K', '1E6E', '1E6J', '1E96', '1EAW', '1EER', '1EFN', '1EWY', '1EZU', '1F34', '1F51', '1F6M', '1FAK', '1FC2', '1FCC', '1FFW', '1FLE', '1FQ1', '1FQJ', '1FSK', '1GCQ', '1GHQ', '1GL1', '1GLA', '1GP2', '1GPW', '1GRN', '1GXD', '1H1V', '1H9D', '1HCF', '1HE1', '1HE8', '1HIA', '1I2M', '1I4D', '1I9R', '1IB1', '1IBR', '1IJK', '1IQD', '1IRA', '1J2J', '1JIW', '1JK9', '1JMO', '1JPS', '1JTD', '1JTG', '1JWH', '1JZD', '1K4C', '1K5D', '1K74', '1KAC', '1KKL', '1KLU', '1KTZ', '1KXP', '1KXQ', '1LFD', '1M10', '1M27', '1MAH', '1ML0', '1MLC', '1MQ8', '1NCA', '1NSN', '1NW9', '1OC0', '1OFU', '1OPH', '1OYV', '1PPE', '1PVH', '1PXV', '1QA9', '1QFW', '1R0R', '1R6Q', '1R8S', '1RKE', '1RLB', '1RV6', '1S1Q', '1SBB', '1SYX', '1T6B', '1TMQ', '1UDI', '1US7', '1VFB', '1WEJ', '1WQ1', '1XD3', '1XQS', '1XU1', '1Y64', '1YVB', '1Z0K', '1Z5Y', '1ZHH', '1ZHI', '1ZLI', '1ZM4', '2A1A', '2A5T', '2A9K', '2ABZ', '2AJF', '2AYO', '2B42', '2B4J', '2BTF', '2C0L', '2CFH', '2FD6', '2FJU', '2G77', '2GAF', '2GTP', '2H7V', '2HLE', '2HMI', '2HQS', '2HRK', '2I25', '2I9B', '2IDO', '2J0T', '2J7P', '2JEL', '2MTA', '2NZ8', '2O3B', '2O8V', '2OOB', '2OOR', '2OT3', '2OUL', '2OZA', '2PCC', '2SIC', '2SNI', '2UUY', '2VDB', '2VIS', '2VXT', '2W9E', '2X9A', '2YVJ', '2Z0E', '3A4S', '3AAA', '3AAD', '3BIW', '3BP8', '3BX7', '3CPH', '3D5S', '3DAW', '3EO1', '3EOA', '3F1P', '3FN1', '3G6D', '3H11', '3H2V', '3HI6', '3HMX', '3K75', '3L5W', '3L89', '3LVK', '3MXW', '3P57', '3PC8', '3R9A', '3RVW', '3S9D', '3SGQ', '3SZK', '3V6Z', '3VLB', '4CPA', '4DN4', '4FQI', '4FZA', '4G6J', '4G6M', '4GAM', '4GXU', '4H03', '4HX3', '4IZ7', '4JCV', '4LW4', '4M76', '7CEI']

IDS_WITH_PROBLEM = ['2Z0E', '2NZ8', '1JK9', '1F6M', '1K4C', '1JIW', '1NW9', '1BJ1', '1MLC']

NUMBER_OF_TRANSFORMATIONS_PER_COMPLEX = 200


def get_zdock_benchmark_pdb_path(complex_id, ligand=True, bound=True):
    base_path = os.path.join(BASE_DATA_PATH, complex_id, "benchmark")
    file_name = "{}_{}_{}.pdb".format(complex_id, 'l' if ligand else 'r', 'b' if bound else 'u')
    return os.path.join(base_path, file_name)

def get_patchdock_raw_result_pdb_path(complex_id, rank):
    base_path = os.path.join(BASE_DATA_PATH, complex_id, "patch_dock")
    file_name = "{}.patch_dock_output.{}.pdb".format(complex_id, rank)
    return os.path.join(base_path, file_name)

def get_patchdock_ranked_complex_pdb_path(complex_id, rank, ligand):
    base_path = os.path.join(BASE_DATA_PATH, complex_id, "patch_dock")
    file_name = "{}.patch_dock_output.{}.{}.pdb".format(complex_id, rank, 'l' if ligand else 'r')
    return os.path.join(base_path, file_name)


def get_patchdock_complex_score_file_path(complex_id):
    base_path = os.path.join(BASE_DATA_PATH, complex_id, "patch_dock")
    file_name = "{}.patch_dock_output".format(complex_id)
    return os.path.join(base_path, file_name)

def get_raptorx_dir_path(complex_id):
    return os.path.join(BASE_DATA_PATH, complex_id, "raptorx_results")


def get_patchdock_ranked_complex_processed_data_path(complex_id, rank):
    base_path = os.path.join(BASE_DATA_PATH, complex_id)
    file_name = "{}.processed.{}.json".format(complex_id, rank)
    return os.path.join(base_path, file_name)


def get_zdock_benchmark_processed_data_path(complex_id, bound=True):
    base_path = os.path.join(BASE_DATA_PATH, complex_id)
    file_name = "{}_{}_processed.json".format(complex_id, 'b' if bound else 'u')
    return os.path.join(base_path, file_name)

def get_processed_data_json_path(complex_id):
    base_path = os.path.join(BASE_DATA_PATH, complex_id)
    file_name = "{}_processed.json".format(complex_id)
    return os.path.join(base_path, file_name)


def get_file_from_ml_models_path(file_name):
    return os.path.join(BASE_DATA_PATH, "ml_models", file_name)
