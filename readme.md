<p align="center">
  <h1 align="center"><code>util</code></h1>
</p>

<p align="center">
A machine learning, optimization, and data science utilities
package.
</p>

This is a huge ongoing project that I use for all of my work. Pretty much every new algorithm I implement will start here. I do my best to keep it organized and readable. Once things are thoroughly polished (and presumed generally useful), I turn them into their own python
modules.

## INSTALLATION:

  Latest stable and tested release available by:

```bash
pip install https://github.com/tchlux/util/archive/1.0.17.zip
```

  Current cutting edge, (possibly unstable) code available by:

```bash
pip install git+https://github.com/tchlux/util.git
```

  Include in a [`requirements.txt`](https://pip.pypa.io/en/stable/user_guide/#requirements-files) file as:

```
util@https://github.com/tchlux/util/archive/1.0.17.zip
```

## DOCUMENTATION

The following list links to each module within the `util` package. This serves as the web-accessible documentation for the code. If this documentation seems wrong, the most up-to-date documentation is in the comment block preceding a function or class.

## Big Modules

#### [`util.approximate`](util/approximate#user-content-utilapproximate)

#### [`util.data`](util/data#user-content-utildata)

#### [`util.plot`](util/plot#user-content-utilplot)

## Medium Modules

#### [`util.optimize`](util/optimize#user-content-utiloptimize)

Numerous optimization utilities for both gradient-free and gradient-based minimization of arbitrary functions.

#### [`util.stats`](util/stats#user-content-utilstats)

All sorts of useful statistics, mostly nonparametric. Compute [effect](util/stats/difference.py#L78), fit [CDF's](util/stats/distributions.py#L119), compute [metric principle components](util/stats/metric_pca.py#L50), [rank probability](util/stats/rank.py#L51), or quantify [sample statistics](util/stats/samples.py#L46).

#### [`util.math`](util/math#user-content-utilmath)

Miscellaneous mathematics-related functions. The most notable codes are a [`Fraction`](util/math/fraction.py#L24) rational number implementation in pure python with unlimited precision, and codes that identify and construct [Fekete points](util/math/points.py#L106).

## Small Modules

#### [`util.random`](util/random#user-content-utilrandom)

#### [`util.parallel`](util/parallel#user-content-utilparallel)

#### [`util.system`](util/system#user-content-utilsystem)

#### [`util.decorators`](util/decorators#user-content-utildecorators)

#### [`util.misc`](util/misc#user-content-utilmisc)

Extensive [data splitting for validation](util/misc/multi_dim_analysis.py#L152), [Apriori tree](util/misc/apriori.py#L17), [image transformation](util/misc/image.py#L58), and [latex table generation](util/misc/paper.py#L39). The stuff that doesn't fit anywhere else ... yet.





