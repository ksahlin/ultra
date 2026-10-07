"""A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# using example setup file from https://github.com/pypa/sampleproject/blob/master/setup.py

from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(

    name='ultra_bioinformatics',  # Required
    version='0.3.0',  # Required
    description='Splice aligner of long transcriptomic reads to genome.',  # Required
    long_description=long_description,  # Optional
    long_description_content_type='text/markdown',
    url='https://github.com/ksahlin/uLTRA',  # Optional
    author='Kristoffer Sahlin',  # Optional
    author_email='ksahlin@math.su.se',  # Optional

    # The licence README.md has always declared. The text is in LICENSE.txt,
    # which had never actually been committed until it was added alongside
    # this.
    license='GPL-3.0-only',
    license_files=['LICENSE.txt'],

    # Classifiers help users find your project by categorizing it.
    #
    # For a list of valid classifiers, see
    # https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[  # Optional
        # Was '3 - Alpha', which it has not been for years. 14 releases on
        # PyPI between 2020-08 and 2023-05, 698 commits since 2019, peer
        # reviewed in Bioinformatics (Sahlin & Makinen 2021,
        # doi:10.1093/bioinformatics/btab540) with a citation request in the
        # README, and packaged by third parties into bioconda, biocontainers
        # and a Galaxy singularity image. People run it on real data and cite
        # it in papers. The 0.x version number is the usual argument for
        # Beta, but this classifier describes maturity of use, not semver.
        'Development Status :: 5 - Production/Stable',

        # setuptools >= 77 warns that this classifier is deprecated in favour
        # of the SPDX `license=` expression above, and will eventually make it
        # an error. It is kept because it is still what PyPI's licence filter
        # and most SBOM and licence-scanning tools read, while
        # `License-Expression` needs metadata 2.4 that older setuptools does
        # not emit. When that warning becomes an error, deleting this one line
        # is the whole fix.
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

        # Every version listed here was verified by running the full pipeline
        # over test/ and checking reads.sam byte for byte against the others:
        # 3.9.23, 3.10.21, 3.11.16, 3.12 and 3.13.15 all produce the identical
        # file. The floor is set by the dependencies, not the language -- the
        # sources contain no syntax newer than 3.4 and not one f-string, but
        # pysam and dill both declare Requires-Python >=3.9.
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
    ],

    keywords='Oxford Nanopore transcript long read error correction',  # Optional

    # You can just specify package directories manually here if your project is
    # simple. Or you can use find_packages().
    #
    # Alternatively, if you just want to distribute a single Python file, use
    # the `py_modules` argument instead as follows, which will expect a file
    # called `my_module.py` to exist:
    #
    #   py_modules=["my_module"],
    #
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),  # Required

    # If your package is for Python 2.7, and all versions of Python 3 starting with 3.4, write
    # Was '!=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, <4', which excluded Python
    # 3.0-3.3 and therefore ALLOWED Python 2.6 and 2.7 -- pip would install
    # this on a Python that cannot import it.
    python_requires='>=3.9',
    # This field lists other packages that your project depends on to run.
    # Any package you put here will be installed by pip when your project is
    # installed, so they must be valid existing projects.
    #
    # For an analysis of "install_requires" vs pip's requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['parasail',
                     'pysam',
                     'dill',
                     'intervaltree',
                     'gffutils',
                     'edlib'],  # Optional
    # dependency_links=[], # Optional
    # List additional groups of dependencies here (e.g. development
    # dependencies). Users will be able to install these using the "extras"
    # syntax, for example:
    #
    #   $ pip install sampleproject[dev]
    #
    # Similar to `install_requires` above, these must be valid existing
    # projects.
    # extras_require={  # Optional
    #     'dev': ['check-manifest'],
    #     'test': ['coverage'],
    # },

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # `pip` to create the appropriate form of executable for the target
    # platform.
    # entry_points={  # Optional
    #     'console_scripts': [
    #         'IsoCon=IsoCon.__main__()',
    #     ],
    # },
    scripts=['uLTRA'],
)