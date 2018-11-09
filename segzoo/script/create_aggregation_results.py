import argparse
import sys

import pandas as pd

SPLICE_BASENAMES = ["5' flanking", "initial exon", "initial intron", "internal exon", "internal introns",
                    "terminal exon", "terminal intron", "3' flanking"]


# creates a data frame with segtools' output, computes the means of the values for each gene area and label,
# and writes it in a tab-delimited file


def parse_args(args):

    parser = argparse.ArgumentParser(description='Create a matrix from segtools-aggregation output')

    parser.add_argument('input', type=str, help='path to feature_aggregation.tab')
    parser.add_argument('output', type=str, help='path to store the table file, '
                                                 'including filename and extension (.tsv)')
    parsed_args = parser.parse_args(args)
    return parsed_args


def get_args():
    # allow to run the script from snakemake and cmd
    # if script is run from snakemake
    if 'snakemake' in globals():
        # convert snakemake.io.InputFiles object to str
        args = [str(snakemake.input), str(snakemake.output)]
        args = parse_args(args)
    else:
        args = parse_args(sys.argv[1:])

    return args


def first_col_with_matching_start(start, cols):
    """
    Return the first column from `cols` that starts with `start`.
    None otherwise.
    """
    for col in cols:
        if col.startswith(start):
            return col

if __name__ == '__main__':
    args = get_args()

    df = pd.read_table(args.input, skiprows=1, header=0)
    df.drop(columns=['offset', 'group'], inplace=True)

    component_df = df.groupby("component").mean().T

    # translate component names with `splice_basenames`
    col_translator = {}
    for base in SPLICE_BASENAMES:
        matching_col = first_col_with_matching_start(base, component_df.columns)

        if matching_col:
            col_translator[matching_col] = base

    component_df = component_df.rename(columns=col_translator)

    # save dataframe
    component_df[SPLICE_BASENAMES].to_csv(path_or_buf=args.output, sep='\t')
