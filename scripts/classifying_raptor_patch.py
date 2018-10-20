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

	fpr, tpr, _  = roc_curve(clf.predict_proba(X_test), y_test)
	roc_auc = auc(fpr, tpr)
	# print(clf.score(X_test, y_test))

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

