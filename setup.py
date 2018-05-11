# Always prefer setuptools over distutils
from setuptools import setup, find_packages
from codecs import open
from os import path
from segzoo import __version__

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='segzoo',
    version=__version__,
    description='System for turnkey analysis of semi-automated genome annotations',
    long_description=long_description,  # Optional
    long_description_content_type='text/markdown',
    url='https://bitbucket.org/hoffmanlab/segzoo',
    author='Marc Asenjo',
    author_email='asenjomarc@gmail.com',
    classifiers=[
        'Development Status :: 3 - Alpha',
        # 'Intended Audience :: Developers',
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6'
    ],
    packages=find_packages('.'),
    install_requires=['matplotlib', 'numpy', 'pandas', 'pysam', 'seaborn', 'snakemake', 'pybedtools'],  #'segtools'
    python_requires='>=3.5',
    package_data={  # Optional
        'segzoo': ['Snakefile'],
    },
    entry_points={  # Optional
        'console_scripts': [
            'segzoo=segzoo.main:main',
        ],
    }
)
