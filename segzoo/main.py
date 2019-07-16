import argparse
import sys
import snakemake

from .version import __version__
from os import path


def main(args=sys.argv[1:]):
    here = path.abspath(path.dirname(__file__))
    default_prefix = path.dirname(path.dirname(path.dirname(path.dirname(here))))

    description = '''
    Run various genomic analysis on a segmentation obtained by segway.
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
                        help='Prefix where all the external data is going to be downloaded, followed by /share/ggd/SPECIES/BUILD')
    parser.add_argument('--download-only', action='store_true',
                        help='Execute only the rules that need internet connection, which store data in a shared directory')
    parser.add_argument('--mne', help='Allows specify an mne file to translate segment labels and track names on the shown on the figure')
    parser.add_argument('--normalize-gmtk', action='store_true', help='If set, normalize gmtk parameters column wise')
    parser.add_argument('--dendrogram', action='store_true', help='If set, perform hierarchical clustering of GMTK parameters table row-wise')
    parser.add_argument('--unlock', action='store_true', help='unlock directory (see snakemake doc)')

    parsed_args = parser.parse_args(args)

    if parsed_args.download_only:
        snakemake.snakemake(path.join(here, "Snakefile"), targets=["download_targets"], cores=parsed_args.j, config=vars(parsed_args))

    else:
        snakemake.snakemake(path.join(here, "Snakefile"), cores=parsed_args.j, config=vars(parsed_args), unlock=parsed_args.unlock)

    # printreason=True, dryrun=True, printshellcmds=True)
