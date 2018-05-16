import argparse
import sys
import snakemake

from .version import __version__
from os import path


def main(args=sys.argv[1:]):
    here = path.abspath(path.dirname(__file__))
    default_prefix = path.dirname(path.dirname(path.dirname(path.dirname(here))))

    description = '''
    Segzoo is a tool that allows to run various genomic analysis on a segmentation obtained by segway.
    The results of each analysis are made available as well as a summarizing visualization of the results.
    The tool will download all necessary data into a common directory and run all the analysis, storing
    the results in an output directory. All this information is then transformed into a set of tables
    that can be found in this same directory under the "data" folder, that are used to generate a final visualization.
    '''

    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version=__version__)
    parser.add_argument('segmentation', help='.bed.gz file, the segmentation/annotation output from Segway')
    parser.add_argument('--parameters', default=False,
                        help='The params.params file used to obtain the gmtk-parameters')
    parser.add_argument('-o', '--outdir', default='outdir', help='Output directory to store all the results')
    parser.add_argument('-j', default=1, metavar='CORES', type=int, help='Number of cores to use')
    parser.add_argument('--species', default='Homo_sapiens', help='Species of the genome used for the segmentation')
    parser.add_argument('--build', default='hg38', help='Build of the genome assembly used for the segmentation')
    parser.add_argument('--prefix', default=default_prefix,
                        help='Prefix where all the GGD data is going to be downloaded, followed by /share/ggd/SPECIES/BUILD')

    parsed_args = parser.parse_args(args)

    snakemake.snakemake(path.join(here, "Snakefile"), cores=parsed_args.j, config=vars(parsed_args))
    # dryrun=True, printshellcmds=True, printreason=True)
