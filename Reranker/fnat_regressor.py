from sklearn.linear_model import Ridge

class FnatRegressor(object):
	def __init__(self):
		self._regressor = Ridge(alpha=4.0)

	def fit(self, X, y):
		self._regressor.fit(X, y)

	def predict(self, X):
		return self._regressor.predict(X)

	def score(self, X, y):
		return self._regressor.score(X, y)
