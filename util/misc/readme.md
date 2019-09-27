# util.misc

The [island of misfit toys](https://www.youtube.com/watch?v=Gr6GbKciNCY). One day these functions will find a deserving home, until now they wait here.

### [`util.misc.multi_dim_analysis`](multi_dim_analysis.py)

##### [`make_test_data`](multi_dim_analysis.py#L152)

Function for generating a multi-dimensional analysis test data set for analyzing the performance of continuous valued `R^d -> R` function approximation algorithms over varying values of `d` and varying ratios of training and testing data. This is `O(d!)` in complexity, so be weary of using more than 10 dimensional data. The testing is very thorough and should give useful insight into the relative performance of approximation techniques on given data.

##### [`MDA_Iterator`](multi_dim_analysis.py#L25)

This class is designed to iterate over training and testing data constructed by `make_test_data`. At each iteration it yields:
  - `train` - Path to training data file.
  - `test`  - Path to testing data file.
  - `dim`   - The input dimension of the problem.
  - `names` - The names of the columns in the input file.
  - `num`   - The numeric identifier of the test file.


### [`util.misc.apriori`](util/misc/apriori.py)
