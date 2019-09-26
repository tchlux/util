# Try to import setuptools (if it fails, the user needs that package)
try: 
    from setuptools import setup, find_packages
except:
    # Custom error (in case user does not have setuptools)
    class DependencyError(Exception): pass
    raise(DependencyError("Missing python package 'setuptools'.\n  pip install --user setuptools"))

import os, sys

# Convenience function for reading information files
def read(f_name, empty_lines=False):
    text = []
    with open(os.path.join(package_about, f_name)) as f:
        for line in f:
            line = line.strip("\n")
            if (not empty_lines) and (len(line.strip()) == 0): continue
            if (len(line) > 0) and (line[0] == "%"): continue
            text.append(line)
    return text

# Go to the "about" directory in the package directory.
package_about = ""
package_name = read("package_name.txt")[0]
package_about = os.path.join(os.path.dirname(os.path.abspath(__file__)),package_name,"about")
package_requires = os.path.join(package_about, "requirements.txt")

if __name__ == "__main__":
    #      Read in the package description files     
    # ===============================================
    package = package_name
    version =read("version.txt")[0]
    description = read("description.txt")[0]
    keywords = read("keywords.txt")
    classifiers = read("classifiers.txt")
    name, email, git_username = read("author.txt")
    requirements = read("requirements.txt")
    # Handle the requirements.
    has_external_packages = any("://" in r for r in requirements)
    if has_external_packages:
        # Use "pip" to install requirements. This is more generally
        # accurate than trying to parse out URL's from requirements.txt.
        try:
            from pip._internal import main as pip_main
            pip_main(["install", "-r", package_requires])        
        except:
            external_packages = "\n  ".join(r for r in requirements if "://" in r)
            import warnings
            class ExternalRequirements(Warning): pass
            warnings.warn(ExternalRequirements(
                f"Could not install some external requirements:\n  {external_packages}"))
        # Remove the external requirments from the list.
        requirements = [r for r in requirements if "://" not in r]
    # Call "setup" to formally set up this module.
    setup(
        author = name,
        author_email = email,
        name=package,
        packages=find_packages(exclude=[]),
        include_package_data=True,
        install_requires=requirements,
        version=version,
        url = 'https://github.com/{git_username}/{package}'.format(
            git_username=git_username, package=package),
        download_url = 'https://github.com/{git_username}/{package}/archive/{version}.tar.gz'.format(
            git_username=git_username, package=package, version=version),
        description = description,
        scripts=["scripts/compile"],
        keywords = keywords,
        python_requires = '>=3.6',
        license='MIT',
        classifiers=classifiers
    )


# 
# Install specific versions of a package with something like:
# 
#   pip install https://github.com/tchlux/packager/archive/1.0.0.zip
# 
# 
