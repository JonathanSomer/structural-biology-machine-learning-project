import os
from DockQ import DockQ
from argparse import ArgumentParser
import os.path
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO

# patch lig: A
# patch rec: X
# zdock lig: X
# zdock rec: A

# $ python capri_score_script.py 2B42 -patch_dock_ligand_chains A -patch_dock_receptor_chains X -zdock_benchmark_bound_ligand_chains X -zdock_benchmark_bound_receptor_chains A
# $ python capri_score_script.py 1AHW -patch_dock_ligand_chains A -patch_dock_receptor_chains L H -zdock_benchmark_bound_ligand_chains C -zdock_benchmark_bound_receptor_chains A B
# $ python capri_score_script.py 2B42 -patch_dock_ligand_chains A -patch_dock_receptor_chains X -zdock_benchmark_bound_ligand_chains B -zdock_benchmark_bound_receptor_chains A


NUMBER_OF_PATCH_DOCK_COMPLEXES = 1000
DATA_BASE_DIR = os.path.dirname(os.path.abspath(__file__)) + "/../data/"

def main():
	parser=ArgumentParser(description="Calculate CAPRI score for a PatchDock generated complex")
	parser.add_argument('complex_id',metavar='<complex_id>',type=str,nargs=1,help='4 letter complex code')
	parser.add_argument('-patch_dock_ligand_chains',metavar='patch_dock_ligand_chains', type=str,nargs='+', help='pdb chain order to group together partner 1')
	parser.add_argument('-patch_dock_receptor_chains',metavar='patch_dock_receptor_chains', type=str,nargs='+', help='pdb chain order to group together partner 2')
	parser.add_argument('-zdock_benchmark_bound_ligand_chains',metavar='zdock_benchmark_bound_ligand_chains', type=str,nargs='+', help='pdb chain order to group together from native partner 1')
	parser.add_argument('-zdock_benchmark_bound_receptor_chains',metavar='zdock_benchmark_bound_receptor_chains', type=str,nargs='+', help='pdb chain order to group together from native partner 2')
	args = parser.parse_args()

	complex_id = args.complex_id[0]
	BENCHMARK_DATA_BASE_DIR = DATA_BASE_DIR + complex_id + "/benchmark/"

	zdock_benchmark_bound_complex_path = BENCHMARK_DATA_BASE_DIR + complex_id + "_bound_complex.pdb"
	# if not os.path.isfile(zdock_benchmark_bound_complex_path):

	# 	io = PDBIO()

	# 	unbound_ligand = PDBParser().get_structure(complex_id, BENCHMARK_DATA_BASE_DIR + complex_id + "_l_u.pdb")
	# 	unbound_ligand_chains = unbound_ligand.get_chains()
	# 	unbound_ligand_chain_ids = [chain.get_id() for chain in unbound_ligand_chains]
	# 	# print(unbound_ligand_chain_ids)

	# 	bound_ligand = PDBParser().get_structure(complex_id, BENCHMARK_DATA_BASE_DIR + complex_id + "_l_b.pdb")
	# 	bound_ligand_chains = bound_ligand.get_chains()
	# 	bound_ligand_chain_ids = [chain.get_id() for chain in bound_ligand_chains]
	# 	# print(bound_ligand_chain_ids)

	# 	unbound_receptor = PDBParser().get_structure(complex_id, BENCHMARK_DATA_BASE_DIR + complex_id + "_r_u.pdb")
	# 	unbound_receptor_chains = unbound_receptor.get_chains()
	# 	unbound_receptor_chain_ids = [chain.get_id() for chain in unbound_receptor_chains]
	# 	# print(unbound_receptor_chain_ids)

	# 	bound_receptor = PDBParser().get_structure(complex_id, BENCHMARK_DATA_BASE_DIR + complex_id + "_r_b.pdb")
	# 	bound_receptor_chains = bound_receptor.get_chains()
	# 	bound_receptor_chain_ids = [chain.get_id() for chain in bound_receptor_chains]
	# 	# print(bound_receptor_chain_ids)

	# else:
	# 	print("as expected")


	os.chdir('./DockQ')
	for i in range(1,2):
		patch_dock_complex_path = DATA_BASE_DIR + complex_id + "/patch_dock_ranking/" + complex_id + "." + str(i) + ".pdb"
		
		print(patch_dock_complex_path)
		print(zdock_benchmark_bound_complex_path)


		cmd = ""
		cmd = cmd + "{} {} ".format(patch_dock_complex_path, zdock_benchmark_bound_complex_path)
		# cmd = cmd + "{} {} ".format(patch_dock_complex_path, patch_dock_complex_path)
		cmd = cmd + "-model_chain1 {} ".format(' '.join(args.patch_dock_ligand_chains))
		# print(' '.join(args.patch_dock_ligand_chains))
		cmd = cmd + "-model_chain2 {} ".format(' '.join(args.patch_dock_receptor_chains))
		cmd = cmd + "-native_chain1 {} ".format(' '.join(args.zdock_benchmark_bound_ligand_chains))
		# print(' '.join(args.zdock_benchmark_bound_ligand_chains))
		cmd = cmd + "-native_chain2 {} ".format(' '.join(args.zdock_benchmark_bound_receptor_chains))
		cmd = cmd + "-perm1 -perm2 -verbose"

		res = DockQ.main(cmd.split())
		print("CAPRI SCORE:" + res)

if __name__ == '__main__':
  main()
