import argparse
import sys
import snakemake

from .version import __version__
from os import path


# TODO add build and species to the path of the output so that different runs in the same folder can be done
def main(args=sys.argv[1:]):
    here = path.abspath(path.dirname(__file__))
    default_prefix = path.dirname(path.dirname(path.dirname(path.dirname(here))))

    parser = argparse.ArgumentParser(description='Download necessary files, run workflow and obtain results',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version=__version__)
    parser.add_argument('segmentation', help='.bed.gz file, the segmentation/annotation output from Segway')
    parser.add_argument('--parameters', default=False,
                        help='The params.params file used to obtain the gmtk-parameters')
    parser.add_argument('-o', '--outdir', default='outdir', help='Output directory to store all the results')
    parser.add_argument('--species', default='Homo_sapiens', help='Species of the genome used for the segmentation')
    parser.add_argument('--build', default='hg38', help='Build of the genome assembly used for the segmentation')
    parser.add_argument('--prefix', default=default_prefix,
                        help='Prefix where all the GGD data is going to be downloaded')

    parsed_args = parser.parse_args(args)

    snakemake.snakemake(path.join(here, "Snakefile"), config=vars(parsed_args))  # dryrun=True, printshellcmds=True)
