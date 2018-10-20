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

	i = 1
	fig = plt.figure()
	fig.suptitle('Raptor & PatchDock Scores + Accepted/Not-Accepted')

	for complex_id in progressbar.progressbar(get_all_training_complexes()):
		complex_X = np.array(map_complex_id_to_X[complex_id])
		if len(complex_X) == 0:
			continue

		complex_X[:,0] = (complex_X[:,0] - np.mean(complex_X[:,0])) / np.std(complex_X[:,0])
		complex_X[:,1] = (complex_X[:,1] - np.mean(complex_X[:,1])) / np.std(complex_X[:,1])
		
		X = complex_X
		y = np.array(map_complex_id_to_y[complex_id]).squeeze()
		
		y = np.greater(y, threshold)
		X = np.array(X)

		X_neg = X[~y]
		X_pos = X[y]

		ax = fig.add_subplot(N_PLOTS,N_PLOTS,i)
		# plt.plot(X_neg[:,0], np.repeat(1, len(X_neg)),'bs', X_pos[:,0], np.repeat(1, len(X_pos)), 'r^', X_neg[:,1], np.repeat(2, len(X_neg)),'bs', X_pos[:,1], np.repeat(2, len(X_pos)), 'r^')
		ax.set_title(complex_id)
		plt.plot(X_neg[:,0], X_neg[:,1], 'bs', label = "FNAT < {}".format(str(threshold)))
		plt.plot(X_pos[:,0], X_pos[:,1], 'r^', label = "FNAT >= {}".format(str(threshold)))
		ax.set_xlabel('PatchDock Score (Standardized)')
		ax.set_ylabel('Raptor Score (Standardized)')

		i += 1
		if i > N_PLOTS*N_PLOTS:
			break
		
	ax.legend(bbox_to_anchor=(1.05, 0), loc='lower left', borderaxespad=0.)
	plt.show()
	# y = ~np.equal(np.array(y), 0.0)



	

