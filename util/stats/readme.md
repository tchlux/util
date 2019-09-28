<h1 align="center"><code>util.stats</code></h1>

All sorts of useful statistics, mostly nonparametric. Compute [effect](difference.py#L78), fit [CDF's](distributions.py#L119), compute [metric principle components](metric_pca.py#L50), [rank probability](rank.py#L51), or quantify [sample statistics](samples.py#L46).

## [`util.stats.difference`](difference.py)

#### [`def effect`](difference.py#L78)

Compute the effect between two sequences. Use the following:
- `number` vs. `number` - Correlation coefficient between two sequences
- `category` vs. `number` - Compute `method` (mean) L1-norm difference between full distribution and conditional distributions for each possible category.
- `category` vs. `category` - Compute the `method` (mean) total difference between full distribution of other sequence given values of one sequence.

## [`util.stats.distributions`](distributions.py)

#### [`class Distribution`](distributions.py#L7)

A class for handling standard operators over distributions. It supports addition with another distribution, subtraction from another distribution (KS statistic), and multiplication by a constant.

#### [`def gauss_smooth`](distributions.py#L68)

Produce a new function that uses a normal distribution to weight nearby function values. Only works for 1-dimensional distributions.

#### [`def cdf_fit`](distributions.py#L119)

Fit a CDF (cumulative distribution function) to given data. Offers `fit` of `None` (piecewise constant), `linear` (interpolation), and `cubic` (interpolation) options.

#### [`def pdf_fit`](distributions.py#L280)

Fit a PDF (probability density function) to given data. Offers `fit` of `None` (piecewise constant), `linear` (interpolation), and `cubic` (interpolation) options. By default, it uses Gaussian smoothing over the CDF.

## [`util.stats.ks`](ks.py)

#### [`def ks_diff`](ks.py#L7)

Compute the [Kolmogorov-Smirnov](https://en.wikipedia.org/wiki/Kolmogorov–Smirnov_test) statistic between two provided functions, mathematically known as the max-norm difference. Has options for the method of computing the statistic, usually relies on optimization.

#### [`def ks_p_value`](ks.py#L54)

Given a KS statistic and the sample sizes of the two distributions, return the largest confidence with which the two distributions an be said to be the same. More technically, this is the largest confidence for which the null hypothesis (two distributions are same) cannot be refuted.

#### [`def normal_confidence`](ks.py#L87)

Return the "confidence" that the provided distribution is normal. Uses the KS statistic and the `ks_p_value` function to compute value.

## [`util.stats.metric_pca`](metric_pca.py)

#### [`def pca`](metric_pca.py#L79)

Use `sklearn` to compute the principle components of a matrix of row-vectors. Return the components and their magnitudes.

#### [`def mpca`](metric_pca.py#L50)

Compute the metric principle components of a set of points with associated values and a given metric. Scale all vectors between points by the magnitude difference between associated values. Use `sklearn` to compute the principle components of these metric between vectors. Return the vectors and their magnitudes.

## [`util.stats.plotting`](plotting.py)

#### [`def plot_percentiles`](plotting.py#L71)

Given a `Plot` object (from [`util.plot`](../plot#user-content-utilapproximate)), add in functions for a "percentile cloud" (shaded regions darker in middle) for stacked y values associated with each x value.

## [`util.stats.rank`](rank.py)

#### [`def insert_sorted`](rank.py#L6)

Insert a value into a sorted list of tuples (defaults to looking at the last element of the tuple). Maintains sorted order after insert.

#### [`def product`](rank.py#L22)

Compute the product of all values in an iterable.

#### [`def rank_probability`](rank.py#L51)

Given a target list of numbers (empirical distribution `A`) and a list of lists of numbers (empirical distributions `B1`, `B2`, ...), determine the probability that items from the target list will have rank `k` when selecting randomly (assuming order doesn't matter). AKA, return the probability that a random draw from (`A`, `B1`, `B2`, ... ) leaves the value from `A` in position `k`.

#### [`def performance_profile`](rank.py#L135)

Returns the proportion of items from a list are less than or equal to a value.

#### [`def data_profile`](rank.py#L143)

Computes the data profiles for a given set of performances. This is the likelihood that the performance is within some factor of the best possible performance.

## [`util.stats.samples`](samples.py)

#### [`def samples`](samples.py#L46)

Given optionally a sample size, max error, and confidence, compute whichever of these values was not provided and return it. If all three parameters are given, then return `True` or `False` to determine if that sample meets that criteria.

This estimate is provided distribution free, using only combinatorics. It will not be as strong of an estimate as could be made with a parametric form.
