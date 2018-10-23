import json
from Constants import *
import progressbar
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize
from sklearn.svm import SVC
from sklearn.metrics import roc_curve, auc

def main():

	threshold = 0.3
	N_PLOTS = 2
	with open(train_test_data_path(), "r") as file:
		data = json.load(file)

	map_complex_id_to_X = data["map_complex_id_to_X"]
	map_complex_id_to_y = data["map_complex_id_to_y"]

	# for key, arr in map_complex_id_to_X.items():
	# 	print(len(arr))

	X_train = []
	y_train = []

	X_test = []
	y_test = []

	all_complexes = get_all_training_complexes()
	# all_complexes = ['2GAF', '2PCC', '1I4D', '1VFB', '1FSK', '1NSN', '1EZU', '3VLB', '1US7', '2I25', '2W9E', '1F34', '1KXP', '1UDI', '2A9K', '1EFN', '4H03', '1KXQ', '2B42', '1KAC', '1ML0', '3LVK', '1GHQ', '3H2V', '2B4J', '2AYO', '1RLB', '1Z0K', '2A5T', '1OC0', '1OYV', '1XU1', '1E6E', '1E96', '1GCQ', '2SNI', '1AVX', '2MTA', '1AK4', '2OUL', '2SIC', '2A1A', '1HIA', '2ABZ', '3P57', '1BVN', '1BVK', '1H9D', '1D6R', '3D5S', '2HLE', '2J0T', '2UUY', '1CLV', '1JTD', '1HE1', '1J2J']
	train_complexes = all_complexes[:100]
	test_complexes = all_complexes[100:]

	for complex_id in progressbar.progressbar(train_complexes):
		complex_X = np.array(map_complex_id_to_X[complex_id])
		if len(complex_X) == 0:
			continue
		
		X_train.extend(complex_X)
		y_train.extend(np.array(map_complex_id_to_y[complex_id]).squeeze())

	for complex_id in progressbar.progressbar(test_complexes):
		complex_X = np.array(map_complex_id_to_X[complex_id])
		if len(complex_X) == 0:
			continue
		
		X_test.extend(complex_X)
		y_test.extend(np.array(map_complex_id_to_y[complex_id]).squeeze())


	y_test = np.greater(y_test, threshold)
	y_train = np.greater(y_train, threshold)

	clf = SVC(gamma='auto', probability=True )
	clf.fit(X_train, y_train) 


	import pdb; pdb.set_trace()
	fpr, tpr, _  = roc_curve(clf.predict_proba(X_test), y_test)
	roc_auc = auc(fpr, tpr)

	plt.figure()
	lw = 2
	plt.plot(fpr, tpr, color='darkorange', lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
	plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.title('Receiver operating characteristic example')
	plt.legend(loc="lower right")
	plt.show()

