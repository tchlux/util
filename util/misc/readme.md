<h1 align="center"><code>util.misc</code></h1>

The [island of misfit toys](https://www.youtube.com/watch?v=Gr6GbKciNCY). One day these functions will find a deserving home, until now they wait here.

## [`util.misc.multi_dim_analysis`](multi_dim_analysis.py)

##### [`def make_test_data`](multi_dim_analysis.py#L152)

Function for generating a multi-dimensional analysis test data set for analyzing the performance of continuous valued `R^d -> R` function approximation algorithms over varying values of `d` and varying ratios of training and testing data. This is `O(d!)` in complexity, so be weary of using more than 10 dimensional data. The testing is very thorough and should give useful insight into the relative performance of approximation techniques on given data.

##### [`class MDA_Iterator`](multi_dim_analysis.py#L25)

This class is designed to iterate over training and testing data constructed by `make_test_data`. At each iteration it yields:
  - `train` - Path to training data file.
  - `test`  - Path to testing data file.
  - `dim`   - The input dimension of the problem.
  - `names` - The names of the columns in the input file.
  - `num`   - The numeric identifier of the test file.


## [`util.misc.apriori`](apriori.py)

##### [`class AprioriTree`](apriori.py#L17)

Construct this class and use the [`mine_items`](apriori.py#L132) method to construct item sets. Item sets are increasing size sets of frequently occurring values in the given data. See the [example](apriori.py#L180) in code.


## [`util.misc.display`](display.py)

#### [`def matrix_to_img`](display.py#L9)

Given a matrix (as a `numpy.ndarray`) that has either shape 2 (`(H,W)` gray scale) or shape 3 (`(H,W,3)` color in RGB), convert that matrix to an image with [`PIL`](https://pillow.readthedocs.io/en/stable/). For shape 2 signed matrices, use a color mapping for negative, neutral, and positive colors.
