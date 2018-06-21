from utils.pdb_parser_wrapper import PDBParserWrapper
from utils.pdb_utils import get_rmsd_of_super_imposed_complex

num_of_complexes = 100
receptor = "1bra"
ligand = "1aap"

parser = PDBParserWrapper()
known = parser.get_known_complex_structure(receptor, ligand)
distances_array = []

for i in range(1, num_of_complexes+1):
    structure = parser.get_patch_dock_result_complex_structure(receptor, ligand, i)
    rmsd_score = get_rmsd_of_super_imposed_complex(structure, known)
    distances_array.append((i,rmsd_score))

distances_array.sort(key=lambda x: x[1])
print(distances_array)
