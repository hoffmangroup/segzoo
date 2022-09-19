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
    url='https://github.com/hoffmangroup/segzoo',
    author='Mickael Mendez',
    author_email='mickael.mendez@mail.utoronto.ca',
    classifiers=[
        'Development Status :: 3 - Alpha',
        # 'Intended Audience :: Developers',
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        'Programming Language :: Python :: 3.8'
    ],
    packages=find_packages('.'),
    install_requires=['seaborn', 'segtools', 'snakemake', 'pybedtools'],
    python_requires='==3.8.*',
    package_data={  # Optional
        'segzoo': ['Snakefile'],
    },
    entry_points={  # Optional
        'console_scripts': [
            'segzoo=segzoo.main:main',
        ],
    }

)
