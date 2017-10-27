from setuptools import setup, find_packages

with open("version.txt") as f:
    version = f.read().strip()

setup(
    author = 'Thomas C.H. Lux',
    author_email = 'tchlux@vt.edu',
    name='util',

    packages=find_packages(exclude=[]),
    install_requires=['numpy>=1.11'],
    version=version,
    url = 'https://github.com/tchlux/util',
    download_url = 'https://github.com/tchlux/util/archive/{version}.tar.gz'.format(version=version),
    description = 'A machine learning, optimization, and data science utilities package.',
    keywords = ['python', 'python3'],
    python_requires = '>=3.0',
    license='MIT',
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
]

)
