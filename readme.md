<p align="center">
  <h1 align="center">util</h1>
</p>

<p align="center">
A machine learning, optimization, and data science utilities
package. Almost all the useful code I write starts here.
</p>

<hr>

This is a huge ongoing project that I use for all of my work. Pretty
much every new algorithm I come up with will start here. I do my best
to keep it organized and readable. Once things are thoroughly polished
(and presumed useful), I turn them into their own python modules.

## INSTALLATION:

  Latest stable and tested release available by (sorry `util` is taken
  on PyPi):

```bash
pip install https://github.com/tchlux/util/archive/1.0.16.zip
```

  Current cutting edge, (possibly unstable) code available by:

```bash
pip install git+https://github.com/tchlux/util.git
```

## ABOUT

Contains various useful utilities.

### util.approximate

Classes for many different approximation algorithms. Contains wrappers for converting pure approximators (numeric outputs) to classifiers.

### util.data

Class "Data" behaves like a modified Pandas dataframe, but is written in pure python.

### util.decorators

Decorator function "same_as" makes the signature and documentation of a function copy another.
Decorator function "cache" generates a unique file name for (input,output) pairs from a function and stores the pair in a serialized file for faster re-execution.
Decorator function "stability_lock" uses a cache of (input,output) pairs to check if a function maintains the same behaviors after modifications are made.
Decorator function "timeout" uses a system call to cancel a python function (must respond to global interpreter lock) after a certain amount of time has elapsed.
Decorator function "type_check" performs (unpythonic) type checking of function inputs before executing a function.

### util.multi_dim_analysis

Function "make_test_data" splits a given data set into components that allow for detailed analysis of model performance with increasing dimension and amounts of training data.
Class "MDA_Iterator" simplifies the process of iterating through test cases generated by the function "make_test_data".

### util.optimize

Function "minimize" uses a meta-heuristic optimization algorithm to solve a global optimization problem given an arbitrary function.

### util.plot

Provides an extensive interface to HTML plotting through *plotly*. All documentation is within module, see documentation of submodule with "from util import plotly; help(plotly)" for more details.

### util.stats

Contains useful statistical functions for data analysis.

### util.system

Function "run" is a (python+OS)-safe interface to command-line execution that cleanly handles errors.
Class "AtomicOpen" provides an atomic file operation class that uses system locking mechanisms to enforce atomic operations on files.

