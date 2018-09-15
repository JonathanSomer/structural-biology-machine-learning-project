import Constants
import logging
import subprocess
from objects.complex import PatchDockComplex, BenchmarkComplex

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)-8s %(message)s',
                    filename='myapp.log',
                    filemode='w')

logger = logging.getLogger("logger")

def add_separator_headers_to_benchmark_files(complex_id):
    with open(Constants.get_zdock_benchmark_pdb_path(complex_id, ligand=True, bound=True), 'r+') as l_file:
        f_lines = l_file.readlines()
        if f_lines[0] != Constants.LIGAND_SEPARATOR_HEADER:
            header = Constants.LIGAND_SEPARATOR_HEADER + Constants.LIGAND_HEADER
            l_file.seek(0)
            logger.debug("adds header {} to ligand file of {}".format(repr(header), complex_id))
            l_file.write(header + "".join(f_lines))


    with open(Constants.get_zdock_benchmark_pdb_path(complex_id, ligand=False, bound=True), 'r+') as r_file:
        f_lines = r_file.readlines()
        if f_lines[0] != Constants.RECEPTOR_HEADER:
            r_file.seek(0)
            logger.debug("adds header {} to receptor file of {}".format(repr(Constants.RECEPTOR_HEADER), complex_id))
            r_file.write(Constants.RECEPTOR_HEADER + "".join(f_lines))


def create_patchdock_results(complex_id, n_of_results):
    logger.info("Creates {} patchdock results for {}".format(n_of_results, complex_id))
    subprocess.call("getPatchdockResults.sh", Constants.get_patchdock_complex_score_file_path(complex_id),
                    "1", str(n_of_results))
    logger.info("Done creating results")

def split_patch_dock_results(complex_id, n_of_results):
    logger.debug("Splitting patchdock results")
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
    res = []
    for rank in range(1, n_of_results + 1):
        l_path = Constants.get_patchdock_ranked_complex_pdb_path(complex_id, rank, ligand=True)
        r_path = Constants.get_patchdock_ranked_complex_pdb_path(complex_id, rank, ligand=False)
        c = PatchDockComplex(complex_id, complex_id)
        b = BenchmarkComplex(complex_id)
        j = {
            "rank": rank,
            "neighbors5": c.get_neighbours_residues(5),
            "neighbors8": c.get_neighbours_residues(8)

        }



split_patch_dock_results("1A2K", 3)
