from util.approximate.base import Approximator

# This class provides a simple wrapper for the SkLearn decision tree.
class DecisionTree(Approximator):
    def __init__(self, *args, **kwargs):
        from sklearn.tree import DecisionTreeRegressor
        self.tree = DecisionTreeRegressor(*args, **kwargs)

    def _fit(self, x, y, *args, **kwargs):
        return self.tree.fit(x, y, *args, **kwargs)

    def _predict(self, *args, **kwargs):
        return self.tree.predict(*args, **kwargs)


if __name__ == "__main__":
    from util.approximate.testing import test_plot
    m = DecisionTree()
    p, x, y = test_plot(m, random=True, N=20)
    p.show()
