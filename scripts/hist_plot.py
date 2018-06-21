from utils.pdb_parser_wrapper import PDBParserWrapper
from utils.raptorx_parser import get_matrix_from_raptorx
from matplotlib import pyplot as plt


parser = PDBParserWrapper()
struct = parser.get_patch_dock_result_complex_structure("1bra", "1aap", 66)

raptorx_mat = get_matrix_from_raptorx("1bra", "1aap")
neighbors = struct.get_aa_neighbors(radius=8)

neighbors_score = []
not_neighbors_score = []
for recrptor_index, row in enumerate(raptorx_mat):
    for ligand_index, score in enumerate(row):
        if (recrptor_index, ligand_index) in neighbors:
            neighbors_score.append(raptorx_mat[recrptor_index][ligand_index])
        else:
            not_neighbors_score.append(raptorx_mat[recrptor_index][ligand_index])

plt.hist(neighbors_score, cumulative=-1, color='red', histtype='step', normed=True)
plt.hist(not_neighbors_score, cumulative=-1, color='blue', histtype='step', normed=True)
plt.show()