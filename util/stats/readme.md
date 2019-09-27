<h1 align="center"><code>util.stats</code></h1>

Contains useful statistical functions for data analysis.

## [`util.stats.difference`](difference.py)

#### [`def effect`](difference.py#L78)

Compute the effect between two sequences. Use the following:
- `number` vs. `number` - Correlation coefficient between two sequences
- `category` vs. `number` - Compute `method` (mean) 1-norm difference between full distribution and conditional distributions for each possible category.
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

## [`util.stats.metric_pca`](metric_pca.py)

## [`util.stats.modes`](modes.py)

## [`util.stats.plotting`](plotting.py)

## [`util.stats.rank`](rank.py)

## [`util.stats.samples`](samples.py)

