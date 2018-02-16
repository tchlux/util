from util.algorithms import VoronoiMesh

# Create an instance of an algorithm that can be ued to classify
# data. The "fit" method causes the classifier to construct the
# necessary regression models to predict each class.
class Classifier:
    def __init__(self, alg=VoronoiMesh):
        self.alg = alg
        self.models = []

    # Given x (matrix with instance rows and feature columns) and y
    # (class labels), optionally a stratified subsample size (the
    # number of random subsets to model for each class), and
    # optionally use feature weighting for fitting, generate a model
    # that can use the internal algorithm to make predictions at
    # previously unseen x points.
    def fit(self, x, y, stratification=None, feature_weighting=False):
        self.classes = sorted(set(y))
        # Build a model for each class
        for value in self.classes:
            response = np.where(y == value, 1, 0)
            self.models.append( self.alg() )
            self.models[-1].fit(x, response)

    # Given a new x (or set of x), calculate the value of each class
    # model for that x and return the (set of) predicted classes.
    def predict(self, x):
        classes = [self.classes[idx] for idx in
                   np.argmax([m(x) for m in self.models], axis=0)]
        return np.array(classes)
