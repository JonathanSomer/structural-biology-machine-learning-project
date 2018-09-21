import Constants
from objects.complex import PatchDockComplex, BenchmarkComplex
from utils.fnat_utils import get_fnat_score
import json
import os
from nova.transformation_getter import main as get_tranformations
from datetime import datetime

'''
!!! this script assumes benchmark has lignd/receptor headers
'''

def main(complex_ids, n_of_results, user_name, password, remove_results=True):

    for c_id in complex_ids:
        try:
            print("----------------Starts {} {}--------------".format(c_id, datetime.now()))
            create_patchdock_results(c_id, n_of_results, user_name, password)
            split_patch_dock_results(c_id, n_of_results)
            create_process_data_json(c_id, n_of_results)
            if remove_results: remove_patchdock_results(c_id, n_of_results)
            print("----------------Done {} {}--------------".format(c_id,  datetime.now()))
        except Exception as e:
            print(e.message)
            print("----------------Error with complex {}--------------".format(c_id))


def create_patchdock_results(complex_id, n_of_results, user_name, password):
    print("Creates {} patchdock results for {}".format(n_of_results, complex_id))
    get_tranformations(complex_code=complex_id, num_transformations=n_of_results,
                       username=user_name, password=password, test=False,
                       force=False, keep_transformations=True)
    print("Done creating results")

def split_patch_dock_results(complex_id, n_of_results):
    print("Splitting patchdock results")
    for rank in range(1, n_of_results + 1):
        raw_path = Constants.get_patchdock_raw_result_pdb_path(complex_id, rank)
        l_path = Constants.get_patchdock_ranked_complex_pdb_path(complex_id, rank, ligand=True)
        r_path = Constants.get_patchdock_ranked_complex_pdb_path(complex_id, rank, ligand=False)
        with open(raw_path, "r") as raw_f:
            r_s, l_s = raw_f.read().split(Constants.LIGAND_SEPARATOR_HEADER)
        with open(l_path, "w") as l_f:
            l_f.write(l_s)
        with open(r_path, "w") as r_f:
            r_f.write(r_s)

def create_process_data_json(complex_id, n_of_results):
    print("start creating process data json file...")
    bound = BenchmarkComplex(complex_id)
    res = []
    for rank in range(1, n_of_results + 1):
        if rank%100 == 0: print(rank)
        c = PatchDockComplex(complex_id, rank)
        j = {
            "rank": rank,
            "nb5": c.get_neighbours_residues(5),
            "fnat5": get_fnat_score(c, bound, 5),
            "nb8": c.get_neighbours_residues(8),
            "fnat8": get_fnat_score(c, bound, 8)
        }
        res.append(j)
    with open(Constants.get_processed_data_json_path(complex_id), 'w') as f:
        json.dump(res, f)
    print("Done")

def remove_patchdock_results(complex_id, n_of_results):
    for rank in range(1, n_of_results+1):
        if os.path.isfile(Constants.get_patchdock_raw_result_pdb_path(complex_id, rank)):
            os.remove(Constants.get_patchdock_raw_result_pdb_path(complex_id, rank))
        if os.path.isfile(Constants.get_patchdock_ranked_complex_pdb_path(complex_id, rank, ligand=True)):
            os.remove(Constants.get_patchdock_ranked_complex_pdb_path(complex_id, rank, ligand=True))
        if os.path.isfile(Constants.get_patchdock_ranked_complex_pdb_path(complex_id, rank, ligand=False)):
            os.remove(Constants.get_patchdock_ranked_complex_pdb_path(complex_id, rank, ligand=False))
        if os.path.isfile(Constants.get_patchdock_ranked_complex_processed_data_path(complex_id, rank)):
            os.remove(Constants.get_patchdock_ranked_complex_processed_data_path(complex_id, rank))
    if os.path.isfile(Constants.get_zdock_benchmark_processed_data_path(complex_id, bound=True)):
        os.remove(Constants.get_zdock_benchmark_processed_data_path(complex_id, bound=True))
    if os.path.isfile(Constants.get_zdock_benchmark_processed_data_path(complex_id, bound=False)):
        os.remove(Constants.get_zdock_benchmark_processed_data_path(complex_id, bound=False))


all_c_i = ['1A2K', '1ACB', '1AHW', '1AK4', '1AKJ', '1ATN', '1AVX', '1AY7', '1AZS', '1B6C', '1BGX', '1BJ1', '1BKD', '1BUH', '1BVK', '1BVN', '1CGI', '1CLV', '1D6R', '1DE4', '1DFJ', '1DQJ', '1E4K', '1E6E', '1E6J', '1E96', '1EAW', '1EER', '1EFN', '1EWY', '1EZU', '1F34', '1F51', '1F6M', '1FAK', '1FC2', '1FCC', '1FFW', '1FLE', '1FQ1', '1FQJ', '1FSK', '1GCQ', '1GHQ', '1GL1', '1GLA', '1GP2', '1GPW', '1GRN', '1GXD', '1H1V', '1H9D', '1HCF', '1HE1', '1HE8', '1HIA', '1I2M', '1I4D', '1I9R', '1IB1', '1IBR', '1IJK', '1IQD', '1IRA', '1J2J', '1JIW', '1JK9', '1JMO', '1JPS', '1JTD', '1JTG', '1JWH', '1JZD', '1K4C', '1K5D', '1K74', '1KAC', '1KKL', '1KLU', '1KTZ', '1KXP', '1KXQ', '1LFD', '1M10', '1M27', '1MAH', '1ML0', '1MLC', '1MQ8', '1NCA', '1NSN', '1NW9', '1OC0', '1OFU', '1OPH', '1OYV', '1PPE', '1PVH', '1PXV', '1QA9', '1QFW', '1R0R', '1R6Q', '1R8S', '1RKE', '1RLB', '1RV6', '1S1Q', '1SBB', '1SYX', '1T6B', '1TMQ', '1UDI', '1US7', '1VFB', '1WEJ', '1WQ1', '1XD3', '1XQS', '1XU1', '1Y64', '1YVB', '1Z0K', '1Z5Y', '1ZHH', '1ZHI', '1ZLI', '1ZM4', '2A1A', '2A5T', '2A9K', '2ABZ', '2AJF', '2AYO', '2B42', '2B4J', '2BTF', '2C0L', '2CFH', '2FD6', '2FJU', '2G77', '2GAF', '2GTP', '2H7V', '2HLE', '2HMI', '2HQS', '2HRK', '2I25', '2I9B', '2IDO', '2J0T', '2J7P', '2JEL', '2MTA', '2NZ8', '2O3B', '2O8V', '2OOB', '2OOR', '2OT3', '2OUL', '2OZA', '2PCC', '2SIC', '2SNI', '2UUY', '2VDB', '2VIS', '2VXT', '2W9E', '2X9A', '2YVJ', '2Z0E', '3A4S', '3AAA', '3AAD', '3BIW', '3BP8', '3BX7', '3CPH', '3D5S', '3DAW', '3EO1', '3EOA', '3F1P', '3FN1', '3G6D', '3H11', '3H2V', '3HI6', '3HMX', '3K75', '3L5W', '3L89', '3LVK', '3MXW', '3P57', '3PC8', '3R9A', '3RVW', '3S9D', '3SGQ', '3SZK', '3V6Z', '3VLB', '4CPA', '4DN4', '4FQI', '4FZA', '4G6J', '4G6M', '4GAM', '4GXU', '4H03', '4HX3', '4IZ7', '4JCV', '4LW4', '4M76', '7CEI']
not_processed = [c_i for c_i in all_c_i if not os.path.isfile(Constants.get_processed_data_json_path(c_i))]


