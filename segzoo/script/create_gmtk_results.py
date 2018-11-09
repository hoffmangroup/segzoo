import argparse
import re
import sys

import pandas as pd


def unescape_trackname(name):
    """Returns a string where all the hexadicimal sub-strings are decoded.

    gmtk converts non-alphanumerical character to hexadecimal and prepend
    an underscore. For example the the following track name "abc+"
    is converted to "abc_2B".

    This function finds the occurences of "_XX" and replaces it with its decoded version.

    >>> some_name = 'endswithaplus_2B'
    >>> unescape_trackname(some_name)
    endswithaplus+
    """

    gmtk_hexa_pattern = re.compile('_[0-9a-fA-F]{2}')

    for found_pattern in gmtk_hexa_pattern.findall(name):
        hexa_string = found_pattern.replace('_', '')
        hexa_string_decoded = bytes.fromhex(hexa_string).decode('utf-8')

        name = name.replace(found_pattern, hexa_string_decoded)

    return name


def gmtk_means_to_df(filename):
    # Load all the lines in a list
    with open(filename, 'r') as file_handler:
        param_lines = [line.strip() for line in file_handler]

    # Here I assume that lines of interest start with the following pattern
    mean_seg_line_re = re.compile('\d*\s*mean_seg(\d+)_subseg\d+_([^ ]+) \d (.*)')
    mean_seg_lines = filter(mean_seg_line_re.match, param_lines)

    # Extract label, mean and track name information and store it in a list of tuple
    label_track_mean_tuples = [mean_seg_line_re.match(line).groups() for line in mean_seg_lines]

    gmtk_df = pd.DataFrame(label_track_mean_tuples, columns=['label', 'track', 'gmtk_mean'])

    # decode hexadecimal string
    gmtk_df.track = gmtk_df.track.map(unescape_trackname)

    # convert dataframe to a matrix form
    gmtk_matrix = gmtk_df.pivot(values='gmtk_mean', index='track', columns='label').astype(float)

    return gmtk_matrix


def parse_args(args):

    parser = argparse.ArgumentParser(description='Save GMTK mean parameters in tabular format from input.master '
                                                 'or params.params')

    parser.add_argument('input', type=str, help='path to a gmtk param file')
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


if __name__ == '__main__':
    args = get_args()

    df = gmtk_means_to_df(args.input)
    df.T.to_csv(args.output, sep='\t')
