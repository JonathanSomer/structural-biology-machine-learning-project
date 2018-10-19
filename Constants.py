import os.path

project_dir = os.path.dirname(os.path.abspath(__file__))

TOP_RESULTS_COUNTS_FOR_CAPRI = 10

BASE_DATA_PATH = os.path.join(project_dir, "data")

LIGAND_SEPARATOR_HEADER = "LIGAND SEPARATOR\n"
LIGAND_HEADER = "LIGAND\n"
RECEPTOR_HEADER = "RECEPTOR\n"

TOP_RAPTOR_CORRELATION_IDS = ['1R6Q', '1ACB', '1CGI', '1IQD', '4DN4', '3RVW', '3V6Z', '1KXQ']

RIGID_COMLEXES = ['1AHW', '1BVK', '1DQJ', '1E6J', '1JPS', '1MLC', '1VFB', '1WEJ', '2FD6', '2I25', '2VIS', '2VXT',
                  '2W9E', '3EOA', '3HMX', '3MXW', '3RVW', '4DN4', '4FQI', '4G6J', '4G6M', '4GXU', '1BJ1', '1FSK',
                  '1I9R', '1IQD', '1K4C', '1KXQ', '1NCA', '1NSN', '2JEL', '1AVX', '1AY7', '1BUH', '1BVN', '1CLV',
                  '1D6R', '1DFJ', '1E6E', '1EAW', '1EWY', '1EZU', '1F34', '1F51', '1FLE', '1GL1', '1GLA', '1GXD',
                  '1HIA', '1JTD', '1JTG', '1JWH', '1MAH', '1OC0', '1OPH', '1OYV', '1PPE', '1R0R', '1TMQ', '1UDI',
                  '1US7', '1WDW', '1YVB', '1Z5Y', '2A1A', '2A9K', '2ABZ', '2AYO', '2B42', '2GAF', '2J0T', '2MTA',
                  '2O8V', '2OOB', '2OOR', '2OUL', '2PCC', '2SIC', '2SNI', '2UUY', '2YVJ', '3A4S', '3K75', '3LVK',
                  '3PC8', '3SGQ', '3VLB', '4CPA', '4H03', '4HX3', '7CEI', '1A2K', '1AK4', '1AKJ', '1AZS', '1E96',
                  '1EFN', '1EXB', '1FCC', '1FFW', '1FQJ', '1GCQ', '1GHQ', '1GPW', '1H9D', '1HCF', '1HE1', '1I4D',
                  '1J2J', '1K74', '1KAC', '1KLU', '1KTZ', '1KXP', '1M27', '1ML0', '1OFU', '1PVH', '1QA9', '1RLB',
                  '1RV6', '1S1Q', '1SBB', '1T6B', '1XD3', '1XU1', '1Z0K', '1ZHH', '1ZHI', '2A5T', '2AJF', '2B4J',
                  '2BTF', '2FJU', '2G77', '2GTP', '2HLE', '2HQS', '2VDB', '2X9A', '3BIW', '3BP8', '3D5S', '3H2V',
                  '3P57', '4M76']
MEDIUM_COMPLEXES = ['3EO1', '3G6D', '3HI6', '3L5W', '3V6Z', '1CGI', '1IJK', '1JIW', '1KKL', '1M10', '1NW9', '1R6Q',
                    '1ZM4', '2NZ8', '2Z0E', '4FZA', '4IZ7', '4LW4', '1B6C', '1FC2', '1GP2', '1GRN', '1HE8', '1I2M',
                    '1IB1', '1K5D', '1LFD', '1MQ8', '1N2C', '1SYX', '1WQ1', '1XQS', '2CFH', '2H7V', '2HRK', '2OZA',
                    '3AAA', '3AAD', '3BX7', '3CPH', '3DAW', '3R9A', '3S9D', '3SZK', '4JCV']
DIFFICULT_COMPLEXES = ['1BGX', '2HMI', '1ACB', '1F6M', '1FQ1', '1JK9', '1JMO', '1JZD', '1PXV', '1ZLI', '2IDO', '2O3B',
                       '2OT3', '3FN1', '3H11', '4GAM', '1ATN', '1BKD', '1DE4', '1E4K', '1EER', '1FAK', '1H1V', '1IBR',
                       '1IRA', '1R8S', '1RKE', '1Y64', '2C0L', '2I9B', '2J7P', '3AAD', '3F1P', '3L89']

IDS_WITH_PROBLEM = ['2Z0E', '2NZ8', '1JK9', '1F6M', '1K4C', '1JIW', '1NW9', '1BJ1', '1MLC']

NO_ACCEPTED = ['1OPH', '1JWH', '1OFU', '3MXW', '1BJ1', '3HMX', '1IQD', '1DQJ', '1GXD', '2BTF', '1ZHH', '1KTZ', '1QA9']
LESS_THEN_10_ACCEPTED = ['1KLU', '4G6J', '3RVW', '1FCC', '1SBB', '1E6J', '4M76', '1AZS', '1A2K', '1T6B', '3BIW', '3EOA', '2O8V', '2FD6', '1YVB', '4G6M', '1XD3', '1WEJ', '3K75', '2VXT', '1I9R', '4DN4', '1Z5Y', '2VDB', '1NCA', '1JPS', '1GLA', '1K4C', '1AKJ', '2JEL', '1K74', '1ZHI', '1PVH', '2OOR', '1RV6', '1HCF', '1F51', '1MAH', '1JTG', '1DFJ', '4HX3', '1BUH', '2GTP', '1TMQ']
BETWEEN_10_TO_40_ACCEPTED = ['2GAF', '2PCC', '1I4D', '1VFB', '1FSK', '1NSN', '1EZU', '3VLB', '1US7', '2I25', '2W9E', '1F34', '1KXP', '1UDI', '2A9K', '1EFN', '4H03', '1KXQ', '2B42', '1KAC', '1ML0', '3LVK', '1GHQ', '3H2V', '2B4J', '2AYO', '1RLB', '1Z0K', '2A5T', '1OC0', '1OYV', '1XU1', '1E6E', '1E96', '1GCQ', '2SNI', '1AVX', '2MTA', '1AK4', '2OUL', '2SIC', '2A1A', '1HIA', '2ABZ', '3P57', '1BVN', '1BVK', '1H9D', '1D6R', '3D5S', '2HLE', '2J0T', '2UUY', '1CLV', '1JTD', '1HE1', '1J2J']
ABOVE_40_ACCEPTED = ['4CPA', '1AY7', '1S1Q', '1GPW', '3SGQ', '7CEI', '1EAW', '2YVJ', '3A4S', '1M27', '1FLE', '1PPE', '3PC8', '1R0R', '1FFW', '2X9A', '1EWY', '1GL1', '2OOB']

NUMBER_OF_TRANSFORMATIONS_PER_COMPLEX = 200

NEIGHBOUR_DEFAULT_RADIUS = 5
NEIGHBOUR_RADII = [NEIGHBOUR_DEFAULT_RADIUS, 8]
NEIGHBOUR_ATTR_TEMPLATE = '_nb%d'
FNAT_ATTR_TEMPLATE = '_fnat%d'


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

def get_patchdock_complex_score_json_file_path(complex_id):
    base_path = os.path.join(BASE_DATA_PATH, complex_id, "patch_dock")
    file_name = "{}.patch_dock_output.json".format(complex_id)
    return os.path.join(base_path, file_name)


def get_raptorx_dir_path(complex_id):
    return os.path.join(BASE_DATA_PATH, complex_id, "raptorx_results")


def get_complex_fasta_path(complex_id, ligand):
    base_path = os.path.join(BASE_DATA_PATH, complex_id, "benchmark")
    file_name = "{}.fasta".format('ligand' if ligand else 'receptor')
    return os.path.join(base_path, file_name)


def get_processed_data_json_path(complex_id):
    base_path = os.path.join(BASE_DATA_PATH, complex_id)
    file_name = "{}_processed.json".format(complex_id)
    return os.path.join(base_path, file_name)


def get_zdock_benchmark_processed_data_path(complex_id, bound=True):
    base_path = os.path.join(BASE_DATA_PATH, complex_id)
    file_name = "{}_{}_processed.json".format(complex_id, 'b' if bound else 'u')
    return os.path.join(base_path, file_name)


def get_file_from_ml_models_path(file_name):
    return os.path.join(BASE_DATA_PATH, "ml_models", file_name)

def get_all_training_complexes(only_rigid=True):
    all_train = set([c_id for c_id in next(os.walk('data'))[1]
                    if os.path.isdir(get_raptorx_dir_path(c_id))
                    and os.listdir(get_raptorx_dir_path(c_id))
                    and os.path.isfile(get_processed_data_json_path(c_id))])
    if only_rigid:
        return list(all_train.intersection(RIGID_COMLEXES))
    else:
        return list(all_train)


"""
LEGACY paths
def get_patchdock_ranked_complex_processed_data_path(complex_id, rank):
    base_path = os.path.join(BASE_DATA_PATH, complex_id)
    file_name = "{}.processed.{}.json".format(complex_id, rank)
    return os.path.join(base_path, file_name)
"""
