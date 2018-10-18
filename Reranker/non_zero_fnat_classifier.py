from sklearn import svm

class NonZeroFnatClassifier(object):
	def __init__(self):
		self._clf = svm.SVC(C=2.0, kernel='rbf')

	def fit(self, X, y):
		self._clf.fit(X, y)

	def predict(self, X):
		return self._clf.predict(X)

	def score(self, X, y):
		return self._clf.score(X, y)
