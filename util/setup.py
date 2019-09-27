#!python 

# Any python code that should be executed (as if __main__) during
# setup should go here to prepare the module on a new computer.

import warnings
class FailedCompilation(Warning): pass

# --------------------------------------------------------------------
try:
    import util.points
except:
    warnings.warn(FailedCompilation("The `util.points` Fekete module failed to compile. Check `fmodpy` configuration and attempt `import util.points`."))

# --------------------------------------------------------------------
try:
    import util.approximate.shepard
except:
    warnings.warn(FailedCompilation("The `util.approximate.shepard` module failed to compile. Check `fmodpy` configuration and attempt `import util.approximate.shepard`."))

# --------------------------------------------------------------------
try:
    import os
    from util.system import run
    this_directory = os.path.dirname(os.path.abspath(__file__))
    delaunay_path = os.path.join(this_directory, "approximate", "delaunay")
    code, out, err = run(["make"], cwd=delaunay_path)
    assert( code == 0 )
except:
    warnings.warn(FailedCompilation("The `util.approximate.delaunay` module failed to compile. Consider trying to manually `make` it after install."))

# --------------------------------------------------------------------
try:
    import util.approximate.voronoi
except:
    warnings.warn(FailedCompilation("The `util.approximate.voronoi` module failed to compile. Check `fmodpy` configuration and attempt `import util.approximate.voronoi`."))

# --------------------------------------------------------------------
try:
    import util.approximate.linear_shepard
except:
    warnings.warn(FailedCompilation("The `util.approximate.linear_shepard` module failed to compile. Check `fmodpy` configuration and attempt `import util.approximate.linear_shepard`."))

# --------------------------------------------------------------------
try:
    import util.approximate.modified_shepard
except:
    warnings.warn(FailedCompilation("The `util.approximate.modified_shepard` module failed to compile. Check `fmodpy` configuration and attempt `import util.approximate.modified_shepard`."))

# --------------------------------------------------------------------
try:
    import util.approximate.box_mesh
except:
    warnings.warn(FailedCompilation("The `util.approximate.box_mesh` module failed to compile. Check `fmodpy` configuration and attempt `import util.approximate.box_mesh`."))

