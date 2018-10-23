import json
from Constants import *
import progressbar
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize


def main():

	threshold = 0.1
	N_PLOTS = 2
	with open(train_test_data_path(), "r") as file:
		data = json.load(file)

	map_complex_id_to_X = data["map_complex_id_to_X"]
	map_complex_id_to_y = data["map_complex_id_to_y"]

	# for key, arr in map_complex_id_to_X.items():
	# 	print(len(arr))

	X = []
	y = []

	last_indexes_patch = []
	last_indexes_raptor = []
	first_indexes_patch = []
	first_indexes_raptor = []
	mean_indexes_patch = []
	mean_indexes_raptor = []

	i = 1
	fig = plt.figure()
	fig.suptitle('Raptor & PatchDock Order Comparison')

	all_complexes = ['2GAF', '2PCC', '1I4D', '1VFB', '1FSK', '1NSN', '1EZU', '3VLB', '1US7', '2I25', '2W9E', '1F34', '1KXP', '1UDI', '2A9K', '1EFN', '4H03', '1KXQ', '2B42', '1KAC', '1ML0', '3LVK', '1GHQ', '3H2V', '2B4J', '2AYO', '1RLB', '1Z0K', '2A5T', '1OC0', '1OYV', '1XU1', '1E6E', '1E96', '1GCQ', '2SNI', '1AVX', '2MTA', '1AK4', '2OUL', '2SIC', '2A1A', '1HIA', '2ABZ', '3P57', '1BVN', '1BVK', '1H9D', '1D6R', '3D5S', '2HLE', '2J0T', '2UUY', '1CLV', '1JTD', '1HE1', '1J2J']
	all_complexes = get_all_training_complexes()

	for complex_id in progressbar.progressbar(all_complexes):
		complex_X = np.array(map_complex_id_to_X[complex_id])
		if len(complex_X) == 0:
			continue
		
		X = complex_X
		y = np.array(map_complex_id_to_y[complex_id]).squeeze()

		X_patch = X[:,0].squeeze()
		X_raptor = X[:,1].squeeze()

		y_patch = np.copy(y)
		y_raptor = np.copy(y)

		y_patch = y_patch[np.argsort(X_patch)][::-1]
		y_raptor = y_raptor[np.argsort(X_raptor)][::-1]

		y_patch = np.greater(y_patch, threshold)
		if sum(y_patch) ==0 :
			continue
		if len(y_patch) != 1000:
			continue
		patch_order = np.array([i for i in range(1,1001)])
		X_patch_neg = patch_order[~y_patch]
		X_patch_pos = patch_order[y_patch]


		y_raptor = np.greater(y_raptor, threshold)
		if sum(y_raptor) ==0 :
			continue
		raptor_order = np.array([i for i in range(1,1001)])
		X_raptor_neg = raptor_order[~y_raptor]
		X_raptor_pos = raptor_order[y_raptor]


		last_indexes_patch.append(np.max(np.argwhere(y_patch == True).squeeze()))
		last_indexes_raptor.append(np.max(np.argwhere(y_raptor == True).squeeze()))
		first_indexes_patch.append(np.min(np.argwhere(y_patch == True).squeeze()))
		first_indexes_raptor.append(np.min(np.argwhere(y_raptor == True).squeeze()))

		mean_indexes_patch.append(np.mean(np.argwhere(y_patch == True).squeeze()))
		mean_indexes_raptor.append(np.mean(np.argwhere(y_raptor == True).squeeze()))
		
		if i <= N_PLOTS*N_PLOTS:
			ax = fig.add_subplot(N_PLOTS,N_PLOTS,i)
			plt.plot(X_patch_neg, np.repeat(0, len(X_patch_neg)),'bs', label = "FNAT < {}".format(str(threshold))) 
			plt.plot(X_patch_pos, np.repeat(0, len(X_patch_pos)), 'r^', label = "FNAT >= {}".format(str(threshold)))
			plt.plot(X_raptor_neg, np.repeat(1, len(X_raptor_neg)),'bs') 
			plt.plot(X_raptor_pos, np.repeat(1, len(X_raptor_pos)), 'r^')
			plt.ylabel("0: PatchDock, 1: Raptor")
			plt.xlabel("Ranking 1-1000 (1 is highest)")

		i += 1

		
	ax.legend(bbox_to_anchor=(1.05, 0), loc='lower left', borderaxespad=0.)
	plt.show()

	print("LAST")
	print("-- mean:")
	print(np.mean(last_indexes_patch))
	print(np.mean(last_indexes_raptor))

	print("-- median:")
	print(np.median(last_indexes_patch))
	print(np.median(last_indexes_raptor))

	print("FIRST")
	print("-- mean:")
	print(np.mean(first_indexes_patch))
	print(np.mean(first_indexes_raptor))

	print("-- median:")
	print(np.median(first_indexes_patch))
	print(np.median(first_indexes_raptor))

	print("MEAN")
	print("-- mean:")
	print(np.mean(mean_indexes_patch))
	print(np.mean(mean_indexes_raptor))

	print("-- median:")
	print(np.median(mean_indexes_patch))
	print(np.median(mean_indexes_raptor))
	

