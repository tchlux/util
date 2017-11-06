#!/usr/local/bin/python3

class DependencyError(Exception): pass

# Try to import fmodpy (to prepare the fortran sources locally)
try:
    import fmodpy
except:
    raise(DependencyError("Missing python package 'fmodpy'.\n  pip install --user fmodpy"))

import os

# Prepare fortran source files
algorithms_path = os.path.join(os.path.dirname(os.getcwd()), "util", "algorithms")
algorithm_dirs = [a for a in os.listdir(algorithms_path)
                  if os.path.isdir(os.path.join(algorithms_path,a))]
print(algorithms_path)
for alg_name in algorithm_dirs:
    path = os.path.join(algorithms_path, alg_name)
    for f in os.listdir(path):
        if (alg_name in f) and (".f" in f[-4:]):
            source_file = os.path.join(path, f)
            if os.path.isfile(source_file): break
    else:
        print("Could not find source for '%s'."%alg_name)
        continue
    print(source_file)
    if (alg_name == "marspack"):
        fmodpy.wrap(source_file, output_directory=algorithms_path,
                    requested_funcs=["mars","fmod"])
    else:
        fmodpy.wrap(source_file, output_directory=algorithms_path)


