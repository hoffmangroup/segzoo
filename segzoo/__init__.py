import argparse
import sys
import snakemake

from .version import __version__
from os import path


def main(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description='Download necessary files, run workflow and obtain results')
    parser.add_argument('--version', action='version', version=__version__)  # TODO un-hardcode version
    parser.add_argument('segmentation', help='.bed.gz file, the segmentation/annotation output from Segway')
    parser.add_argument('parameters', help='The params.params file used to obtain the gmtk-parameters')
    parser.add_argument('-o', '--outdir', default='outdir', help='Output directory to store all the results')
    parser.add_argument('--species', default='Homo_sapiens', help='Species of the genome used for the segmentation')
    parser.add_argument('--build', default='hg38', help='Build of the genome assembly used for the segmentation')
    parser.add_argument('--prefix', default=path.expanduser('~'),
                        help='Prefix where all the GGD data is going to be downloaded')

    parsed_args = parser.parse_args(args)

    print(parsed_args)

    # """Entry point for the application script"""
    here = path.abspath(path.dirname(__file__))
    snakemake.snakemake(path.join(here, "Snakefile"), config=vars(parsed_args))  # dryrun=True, printshellcmds=True)
