import os
import shutil
import re

import numpy as np
import pandas as pd

from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.PDBExceptions import PDBConstructionException

import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

warnings.simplefilter('ignore', PDBConstructionWarning)

from objects import complex
from utils import pdb_utils, raptorx_utils
from Reranker import reranker
from visualizer import evodock_plot as plot

sanity_check = False

base_path = os.path.abspath(os.pardir)
if sanity_check:
    print base_path
    raise EnvironmentError('Please make sure this is the path to the root directory')
    exit()
benchmark_path = os.path.join(base_path, 'benchmark5')
data_path = os.path.join(base_path, 'data')
pdb_parser = PDBParser()


def get_benchmark_pdb_path(complex_id, ligand=True, bound=True):
    base_path = os.path.join(benchmark_path, 'structures')
    file_name = "{}_{}_{}.pdb".format(complex_id, 'l' if ligand else 'r', 'b' if bound else 'u')
    return os.path.join(base_path, file_name)


def copy_benchmark_pdb(complex_id, complex_benchmark_path, ligand=True, bound=True):
    pdb_path = get_benchmark_pdb_path(complex_id, ligand=ligand, bound=bound)
    if not os.path.exists(complex_benchmark_path):
        os.makedirs(complex_benchmark_path)
    complex_pdb_path = os.path.join(complex_benchmark_path, os.path.basename(pdb_path))
    if not os.path.exists(complex_pdb_path):
        shutil.copy(pdb_path, complex_benchmark_path)
    return complex_pdb_path


benchmark_metadata = pd.read_excel(os.path.join(benchmark_path, 'Table_BM5.xlsx'), header=2,
                                   skiprows=[0, 1, 3, 155, 201])

for complex in benchmark_metadata.itertuples(index=False, name='Complex'):
    # parse complex str to get complex_id.
    comp_re_search = re.search('([0-9a-zA-Z]*)_([a-zA-Z]*):([a-zA-Z]*)', complex[0])
    complex_id = comp_re_search.group(1)
    complex_path = os.path.join(data_path, complex_id)
    complex_benchmark_path = os.path.join(complex_path, 'benchmark')

    # copy pdbs
    ligand_u_path = copy_benchmark_pdb(complex_id, complex_benchmark_path, ligand=True, bound=False)
    receptor_u_path = copy_benchmark_pdb(complex_id, complex_benchmark_path, ligand=False, bound=False)
    ligand_b_path = copy_benchmark_pdb(complex_id, complex_benchmark_path, ligand=True, bound=True)
    receptor_b_path = copy_benchmark_pdb(complex_id, complex_benchmark_path, ligand=False, bound=True)

    ligand_u = pdb_parser.get_structure(complex_id + '_l_u', ligand_u_path)
    receptor_u = pdb_parser.get_structure(complex_id + '_r_u', receptor_u_path)

    ligand_u_seq = pdb_utils.get_structure_sequence(ligand_u)
    receptor_u_seq = pdb_utils.get_structure_sequence(receptor_u)

    with open(os.path.join(complex_benchmark_path, 'ligand.fasta'), 'w') as f:
        f.write('>ligand_%s\n' % complex_id)
        f.write(ligand_u_seq)

    with open(os.path.join(complex_benchmark_path, 'receptor.fasta'), 'w') as f:
        f.write('>receptor_%s\n' % complex_id)
        f.write(receptor_u_seq)
