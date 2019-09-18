# Custom error (in case user does not have setuptools)
class DependencyError(Exception): pass

# Try to import setuptools (if it fails, the user needs that package)
try: 
    from setuptools import setup, find_packages
except:
    raise(DependencyError("Missing python package 'setuptools'.\n  pip install --user setuptools"))

import os
# Convenience function for reading information files
def read(f_name, dir_name="about"):
    text = []
    with open(os.path.join(dir_name, f_name)) as f:
        for line in f:
            line = line.strip()
            if (len(line) > 0) and (line[0] != "#"):
                text.append(line)
    return text

if __name__ == "__main__":
    #      Read in the package description files     
    # ===============================================
    package = os.path.basename(os.path.dirname(os.path.abspath(__file__)))
    version =read("version.txt")[0]
    description = read("description.txt")[0]
    requirements = read("requirements.txt")
    keywords = read("keywords.txt")
    classifiers = read("classifiers.txt")
    name, email, git_username = read("author.txt")

    # Translate the git requirements to proper requirements.
    dependency_links = [r for r in requirements if "github.com" in r]
    for r in dependency_links:
        try:    pkg_name = r.split("egg=")[1]
        except: raise(DependencyError("GitHub repositories must specify '#egg=<package-name>' at the end."))
        requirements[requirements.index(r)] = " ".join([pkg_name,"@","git+https://" + r.split("://")[1]])

    setup(
        author = name,
        author_email = email,
        name=package,
        packages=find_packages(exclude=[]),
        install_requires=requirements,
        version=version,
        description = description,
        keywords = keywords,
        classifiers=classifiers,
        include_package_data=True,
        url = 'https://github.com/{git_username}/{package}'.format(
            git_username=git_username, package=package),
        download_url = 'https://github.com/{git_username}/{package}/archive/{version}.tar.gz'.format(
            git_username=git_username, package=package, version=version),
        python_requires = '>=3.0',
        license='MIT',
    )


# When using Homebrew to install this library, the following
# additional packages may be required to make cv2 import correctly:
# 
#   brew install gtk+
#   brew install linuxbrew/xorg/libsm
# 
