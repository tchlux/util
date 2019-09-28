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

#### [`util.math`](util/math#user-content-utilmath)

Miscellaneous mathematics-related functions. The most notable codes are a [`Fraction`](util/math/fraction.py#L24) rational number implementation in pure python with unlimited precision, and codes that identify and construct [Fekete points](util/math/points.py#L106).

#### [`util.optimize`](util/optimize#user-content-utiloptimize)

Numerous optimization utilities for both gradient-free ([`Adaptive Normal`](util/optimize/adaptive_normal.py#L3), [`DiRect`](util/optimize/direct.py#L77), [`random`](util/optimize/random.py#L3)) and gradient-based ([`SGD`](util/optimize/gradient_based.py#L18), [`ADAM`](util/optimize/gradient_based.py#L63), [`L_BFGS`](util/optimize/gradient_based.py#L7)) minimization of arbitrary functions.

#### [`util.stats`](util/stats#user-content-utilstats)

All sorts of useful statistics, mostly nonparametric. Compute [effect](util/stats/difference.py#L78), fit [CDF's](util/stats/distributions.py#L119), compute [metric principle components](util/stats/metric_pca.py#L50), [rank probability](util/stats/rank.py#L51), or quantify [sample statistics](util/stats/samples.py#L46).

## Small Modules

#### [`util.random`](util/random#user-content-utilrandom)

Functions for generating useful random sequences ([ranges](util/random/random.py#L16), [distributions](util/random/random.py#L88), [pairs](util/random/random.py#L107), [Latin designs](util/random/random.py#L138)).

#### [`util.parallel`](util/parallel#user-content-utilparallel)

Provides a parallelized implementation of the builtin Python [`map`](util/parallel/parallel.py#L84) function for easy drop-in parallelization.

#### [`util.system`](util/system#user-content-utilsystem)

Function [`run`](util/system/system.py#L146) is a (python+OS)-safe interface to command-line execution that cleanly handles errors. Class [`AtomicOpen`](util/system/system.py#L181) provides an atomic file operation class that uses system locking mechanisms to enforce atomic operations on files. Also provides a robust [`hash`](util/system/system.py#L1) function, easy-to-use [`save`](util/system/system.py#L24) / [`load`](util/system/system.py#L36) functions, and a [`Timer`](util/system/system.py#L235) object.

#### [`util.decorators`](util/decorators#user-content-utildecorators)

Miscellaneous useful function decorators (the `@<decorator>` on lines before `def` in code). Alias a function with [`same_as`](util/decorators/decorators.py#L24), [`cache`](util/decorators/decorators.py#L71) function calls, [`capture`](util/decorators/decorators.py#L423) outputs, and run functions in the [`background`](util/decorators/decorators.py#526).

#### [`util.misc`](util/misc#user-content-utilmisc)

Extensive [data splitting for validation](util/misc/multi_dim_analysis.py#L152), [Apriori tree](util/misc/apriori.py#L17), [image transformation](util/misc/image.py#L58), and [latex table generation](util/misc/paper.py#L39). The stuff that doesn't fit anywhere else ... yet.





