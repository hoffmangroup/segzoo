#!/usr/bin/env python
# coding: utf-8

import sys
from os import path
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt

import pandas as pd
import seaborn as sns
import scipy

import scipy.cluster.hierarchy as sch

from collections import defaultdict

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.transforms import ScaledTranslation

from segzoo.gene_biotypes import BIOTYPES

mpl.use('Agg')

# VARIABLES AND CONSTANTS

# Table proportion modifiers
NUM_COMPONENTS = 8
GMTK_FACTOR = 1
MIX_FACTOR = 1.5
AGG_FACTOR = 1
OVERLAP_FACTOR = 1
OVERLAP_COLUMN_NUMBER = 2
ROW_CORRECTOR = 1

# Font scaling variables
FONT_SCALE = 1.5
sns.set(font_scale=FONT_SCALE)
# sns.set_context("poster")
LABEL_FONTSIZE = 20 * FONT_SCALE / 1.5
TITLE_FONTSIZE = 25 * FONT_SCALE / 1.5

# Table options and properties
TABLE_POS = "bottom"  # top / bottom / other to ommit
TABLE_HEIGHT = 1  # relative to the height of 2 rows from the mix matrix
TABLE_CONTENT = [['max', 'max', 'max', 'max', 'max', 65],
                 ['min', 'min', 'min', 'min', 'min', 35]]

# Color maps for the visualization
cmap_gmtk = sns.diverging_palette(220, 10, as_cmap=True)
cmap_mix = 'YlGn'
cmap_agg = 'Blues'
cmap_overlap = 'Reds'

def is_decimal_zero(num):
    """
    True if all num's decimals are 0.
    False otherwise

    >>> is_decimal_zero(12.0)
    True
    >>> is_decimal_zero(12.3)
    False
    """
    return int(num) == num


def gt_n_ints(num, n_digits=2):
    """
    True if the number num has more than n digits.
    False otherwise

    >>> gt_n_ints(12.123456)
    False
    >>> gt_n_ints(123.123456)
    True
    >>> gt_n_ints(1.23, 1)
    False
    """
    return len(str(int(num))) > n_digits


def human_format(num):
    """
    Shorten long numbers by replacing trailing intergers by a unit.

    >>> human_format(1_000)
    '1k'
    >>> human_format(1_234)
    '1.2k'
    >>> human_format(12_345)
    '12.3k'
    >>> human_format(123_456)
    '123k'
    >>> human_format(1_000_000)
    '1m'
    """
    magnit_chars = ['', 'k', 'm', 'g', 't', 'p']

    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0

    if gt_n_ints(num) or is_decimal_zero(num):
        return '{}{}'.format(int(num), magnit_chars[magnitude])
    else:
        return '{:.1f}{}'.format(num, magnit_chars[magnitude])

def pretty_number(n):
    """
    Add space every three digits from left to right

    >>> pretty_number(1000)
    '1 000'
    >>> pretty_number(100)
    '100'
    """
    return '{:,}'.format(int(n)).replace(',', ' ')

# Prepare the gmtk parameters in a DataFrame
def gmtk_parameters(args):
    def normalize_col(col):
        return (col-col.min())/(col.max()-col.min())

    df = pd.read_csv(args.gmtk, index_col=0, sep='\t')
    df.sort_index(inplace=True)
    if args.normalize_gmtk:
        df = df.apply(normalize_col, axis=0)
    return df, [df.max().max(), df.min().min()]


# Prepare nucleotide results in a Series format
def nucleotide(args):
    res_nuc_ann = pd.read_csv(args.nuc, index_col=0, sep='\t')['GC content'].round(2) * 100
    res_nuc_ann.sort_index(inplace=True)

    # Rename columns
    res_nuc_ann = res_nuc_ann.rename('GC content (%)')

    # Interpolation of the parameters to rescale them between vmin and vmax for the heatmap
    vmax = 65
    vmin = 35
    res_nuc_hm = res_nuc_ann.copy()
    res_nuc_hm = ((res_nuc_hm - vmin) / (vmax - vmin)).clip(0, 1)

    return res_nuc_hm, res_nuc_ann


# Prepare length_distribution results in a DataFrame
def length_distribution(args):
    # Preparing the annotation for the matrix, creating a new column called 'frac.segs'
    res_len_ann = pd.read_csv(args.len_dist, index_col=0, sep='\t')
    res_len_ann['frac.segs'] = (res_len_ann['num.segs'] / res_len_ann.loc['all']['num.segs']) * 100
    res_len_ann['frac.bp'] = res_len_ann['frac.bp'] * 100
    res_len_ann = res_len_ann.drop(['num.segs', 'num.bp'], axis=1).drop('all')

    # Rename columns
    res_len_ann.index = res_len_ann.index.map(int)  # labels need to be integers
    res_len_ann.sort_index(inplace=True)
    res_len_ann.columns = ['Mean length', 'Median length', 'Std length', 'Base pairs (%)', 'Segments (%)']

    res_len_hm = res_len_ann.copy()
    # Interpolation of the parameters to rescale them between 0 and 1
    for col in res_len_hm.columns:
        res_len_hm[col] = (res_len_hm[col] - res_len_hm[col].min()) / (res_len_hm[col].max() - res_len_hm[col].min())

    res_len_ann = res_len_ann.round(0)
    return res_len_hm, res_len_ann


# Prepare the mix matrix for the heatmap and its annotation, both in DataFrames
def mix_data_matrix(args):
    # Joining the matrices to create final heatmap and annotation
    res_nuc_hm, res_nuc_ann = nucleotide(args)
    res_len_hm, res_len_ann = length_distribution(args)

    res_ann = res_len_ann.join(res_nuc_ann)
    res_hm = res_len_hm.join(res_nuc_hm)

    return res_hm, res_ann

# Prepare the overlap results in Dataframe
def overlap(args):
    df = pd.read_csv(args.overlap, sep='\t', header=0, index_col=0)
    df.sort_index(inplace=True)
    df = df * 100
    df = df.apply(round).astype(int)
    return df

# Prepare the aggregation results in a dictionary of DataFrames by gene_biotype and return the maximum value
def aggregation(args):
    # Rename columns
    column_names = ["5' flanking", "initial exon", "initial intron", "internal exon", "internal introns", "terminal exon", "terminal intron", "3' flanking"]

    def to_percent(row):
        return (row / row.sum()).round(2) * 100

    df_dict = defaultdict()
    max_value = 0
    for biotype in BIOTYPES:
        filename = next(x for x in args.aggs if path.basename(path.dirname(x)) == biotype)
        biotype_df = pd.read_csv(filename, index_col=0, sep='\t').apply(to_percent, axis=1).fillna(0)
        biotype_df.sort_index(inplace=True)
        biotype_df.columns = column_names

        # Update max value
        max_value = max(biotype_df.values.max(), max_value)
        df_dict[biotype] = biotype_df

    return df_dict, max_value


def get_mne_ticklabels(filename, track_labels=[], label_labels=[]):
    """Parse mne file and return updated tracks and labels"""

    mne_df = pd.read_csv(filename, dtype=str, sep='\t')
    mne_df.sort_index(inplace=True)
    assert all(col in ['type', 'old', 'new'] for col in mne_df.columns)

    track_df = mne_df[mne_df.type == 'track']
    track_translator = dict(zip(track_df.old, track_df.new))
    track_labels = map(str, track_labels)
    new_tracks = [track_translator.get(old, old)
                  for old in track_labels]

    label_df = mne_df[mne_df.type == 'label']
    label_translator = dict(zip(label_df.old, label_df.new))
    label_labels = map(str, label_labels)
    new_labels = [label_translator.get(old, old)
                  for old in label_labels]

    return new_tracks, new_labels


# parse arguments
def parse_args(args):
    description = '''
    By running visualization.py directly, you are asked to specify the local results files as parameters.
    This enables fast generation of desired plots. 
    
    If you generated these results using Segzoo, they should be in the folder where your segzoo has been installed, 
    and under /outdir/results/
    
    In many cases, stats would be in your Segzoo environment, for example,
    segzoo_env/share/ggd/Homo_sapiens/hg38/rnaseq/gene_biotype/gene_biotype_stats
    
    If you run Segzoo, outfile would be your_segzoo_folder/outdir/plots/plot.png
    But you do not have to follow this convention.
    '''
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--gmtk', help='Gmtk parameter results produced by Segway')
    parser.add_argument('--normalize-gmtk', action='store_true', default=True,
                        help='If set, normalize gmtk parameters column wise')
    parser.add_argument('--dendrogram', action='store_true',
                        help='If set, perform hierarchical clustering of GMTK parameters table row-wise')
    parser.add_argument('--nuc', help='Nucleotide results file')
    parser.add_argument('--len_dist', help='Length distribution statistics')
    parser.add_argument('--overlap', help='The percentage of segments that overlap with a gene')
    parser.add_argument('--mne', help='Allows specify an mne file to translate segment '
                                      'labels and track names on the shown on the figure')
    parser.add_argument('--aggs', help='Aggregation results file')
    parser.add_argument('--stats', help='Gene biotype stats')
    parser.add_argument('--outfile', help='The path of the resulting visualization')
    return parser.parse_args(args)

if __name__ == '__main__':
    if 'snakemake' in dir():
        arg_list = []
        if snakemake.params.normalize_gmtk:
            arg_list.append('--normalize-gmtk')
        if snakemake.params.dendrogram:
            arg_list.append('--dendrogram')
        arg_list += ['--gmtk', snakemake.input.gmtk,
                     '--nuc', snakemake.input.nuc,
                     '--len_dist', snakemake.input.len_dist,
                     '--overlap', snakemake.input.olp,
                     '--mne', snakemake.input.mne,
                     '--aggs', snakemake.input.aggs,
                     '--stats', snakemake.input.stats,
                     '--outfile', snakemake.output.outfile
                     ]
        args = parse_args(arg_list)
    else:
        args = parse_args(sys.argv[1:])

    # Call the functions that obtain the results in DataFrames
    if args.gmtk:
        res_gmtk, gmtk_max_min = gmtk_parameters(args)
    else:
        res_gmtk = pd.DataFrame()
    res_mix_hm, res_mix_ann = mix_data_matrix(args)
    res_agg_dict, agg_vmax = aggregation(args)

    # Dimensioning variables
    DENDROGRAM_COL = 2
    GMTK_COL = res_gmtk.shape[1] * GMTK_FACTOR + 1
    MIX_COL = res_mix_hm.shape[1] * MIX_FACTOR + 1
    OVERLAP_COL = OVERLAP_COLUMN_NUMBER * OVERLAP_FACTOR + 1
    AGG_COL = len(BIOTYPES) * NUM_COMPONENTS * AGG_FACTOR + 1

    n_rows = res_mix_hm.shape[0] * ROW_CORRECTOR
    n_columns = DENDROGRAM_COL + GMTK_COL + MIX_COL + OVERLAP_COL + AGG_COL

    # Create grid with axes following the ratios desired for the dimensions
    figure, axes = plt.subplots(1, 5, figsize=(n_columns, n_rows),
                                gridspec_kw={"wspace": 8 / n_columns,
                                             "width_ratios": [DENDROGRAM_COL, GMTK_COL, MIX_COL, OVERLAP_COL, AGG_COL]})

    ax_dendrogram, ax_gmtk, ax_mix, ax_overlap, ax_agg = axes

    # Read labels from mne file
    if args.mne and not res_gmtk.empty:
        new_tracks, new_labels = get_mne_ticklabels(args['mne'], res_gmtk.columns, res_mix_hm.index)
    elif args.mne:
        new_tracks, new_labels = get_mne_ticklabels(args['mne'], [], res_mix_hm.index)
    else:
        new_tracks, new_labels = (res_gmtk.columns, res_mix_hm.index)

    # GMTK parameters
    if args.gmtk:
        if args.dendrogram:
            linkage_matrix = sch.linkage(res_gmtk, method='weighted')
            dendrogram = sch.dendrogram(linkage_matrix, ax=ax_dendrogram, orientation='left',
                                        color_threshold=0, above_threshold_color='k',
                                        leaf_font_size=LABEL_FONTSIZE)
            ax_dendrogram.spines['right'].set_visible(False)
            ax_dendrogram.spines['left'].set_visible(False)
            ax_dendrogram.spines['top'].set_visible(False)
            ax_dendrogram.spines['bottom'].set_visible(False)
            ax_dendrogram.set_facecolor((1, 1, 1))    # Set dendrogram background to white
            ax_dendrogram.set_xticklabels('')
            row_ordering = [int(item) for item in dendrogram['ivl']]
            row_ordering.reverse()
            res_gmtk = res_gmtk.loc[row_ordering]
        else:
            figure.delaxes(ax_dendrogram)

        divider_gmtk = make_axes_locatable(ax_gmtk)
        ax_gmtk_cbar = divider_gmtk.append_axes("right", size=0.35, pad=0.3)
        g_gmtk = sns.heatmap(res_gmtk, cmap=cmap_gmtk, ax=ax_gmtk, cbar_ax=ax_gmtk_cbar)
        cbar_gmtk = g_gmtk.collections[0].colorbar

        if args.normalize_gmtk:
            cbar_gmtk.set_ticks(gmtk_max_min)
            cbar_gmtk.ax.set_yticklabels(['col\nmax', 'col\nmin'], fontsize=LABEL_FONTSIZE)
        else:
            cbar_gmtk.ax.set_yticklabels(cbar_gmtk.ax.get_yticklabels(), fontsize=LABEL_FONTSIZE)

        # Setting titles and axis labels
        if not args.dendrogram:
            ax_gmtk.set_yticklabels(new_labels, rotation=0, fontsize=LABEL_FONTSIZE)  # put label names horizontally
        else:
            ax_gmtk.set_yticklabels('')
            ax_gmtk.set_ylabel('')
        ax_gmtk.set_xticklabels(new_tracks, rotation=90, fontsize=LABEL_FONTSIZE)
        ax_gmtk.set_title('GMTK parameters',
                          fontsize=TITLE_FONTSIZE,
                          position=(0, 1 + 0.6 / res_gmtk.shape[0] * FONT_SCALE / 1.5),
                          ha='left', va='bottom')
    else:
        figure.delaxes(ax_gmtk)
        figure.delaxes(ax_dendrogram)

    def sort_by_dendrogram(df):
        if args.gmtk and args.dendrogram:
            return df.loc[row_ordering]
        return df

    # Mix matrix
    divider_mix = make_axes_locatable(ax_mix)
    ax_mix_cbar = divider_mix.append_axes("right", size=0.35, pad=0.3)
    res_mix_ann = sort_by_dendrogram(res_mix_ann)
    res_mix_hm = sort_by_dendrogram(res_mix_hm)
    g_mix = sns.heatmap(res_mix_hm, annot=res_mix_ann.applymap(human_format), cbar=True, cmap=cmap_mix, vmin=0, vmax=1,
                        ax=ax_mix, cbar_ax=ax_mix_cbar, fmt='')
    cbar_mix = g_mix.collections[0].colorbar
    cbar_mix.set_ticks([0, 1])
    cbar_mix.ax.set_yticklabels(['low', 'high'], fontsize=LABEL_FONTSIZE)
    ax_mix.set_ylabel('')
    ax_mix.set_xticklabels(ax_mix.get_xticklabels(), rotation=90, fontsize=LABEL_FONTSIZE)

    # Setting axis labels for the mix matrix
    if args.gmtk:
        ax_mix.set_yticklabels([])
    else:
        ax_mix.set_yticklabels(new_labels, rotation=0, fontsize=LABEL_FONTSIZE)

    # Add min-max table
    mix_columns = res_mix_hm.shape[1]

    if TABLE_POS == "bottom":
        high_low_table = ax_mix.table(
            cellText=TABLE_CONTENT,
            cellColours=[[cbar_mix.cmap(0.99)] * mix_columns, [cbar_mix.cmap(0.01)] * mix_columns],
            bbox=[0, - (TABLE_HEIGHT + .25) / n_rows, 1, TABLE_HEIGHT / n_rows],  # [left,bottom,width,height]
            fontsize=LABEL_FONTSIZE,
            cellLoc='center')
        for j in range(mix_columns):
            high_low_table._cells[(0, j)]._text.set_color('white')   # TODO: do not access protected variables

        # Offset labels down to leave space for the table
        dx = 0
        dy = -(TABLE_HEIGHT + 0.25) * 55 / 72
        offset = ScaledTranslation(dx, dy, figure.dpi_scale_trans)

        for label in ax_mix.xaxis.get_majorticklabels():
            label.set_transform(label.get_transform() + offset)

    elif TABLE_POS == "top":
        high_low_table = ax_mix.table(
            cellText=TABLE_CONTENT,
            cellColours=[[cbar_mix.cmap(0.99)] * mix_columns, [cbar_mix.cmap(0.01)] * mix_columns],
            bbox=[0, 1.02, 1, TABLE_HEIGHT / n_rows],  # [left,bottom,width,height]
            fontsize=LABEL_FONTSIZE,
            cellLoc='center')
        # TODO: try to not iterate through every text. change colour by row
        for j in range(mix_columns):
            high_low_table._cells[(0, j)]._text.set_color('white')   # TODO: do not access protected variables
    
    # Overlap
    divider_overlap = make_axes_locatable(ax_overlap)
    ax_overlap_cbar = divider_overlap.append_axes("right", size=0.35, pad=0.3)
    overlap_hm = sort_by_dendrogram(overlap(args))
    g_overlap = sns.heatmap(overlap_hm, vmin=0, vmax=100, annot=True, cbar=True, fmt='.5g', yticklabels=False, cmap=cmap_overlap, ax=ax_overlap, cbar_ax=ax_overlap_cbar)
    
    cbar_overlap = g_overlap.collections[0].colorbar
    cbar_overlap.set_ticks([0, 100])
    cbar_overlap.ax.set_yticklabels(['0%', '100%'], fontsize=LABEL_FONTSIZE)
    
    ax_overlap.text(0, -0.6 * FONT_SCALE / 1.5, "Overlap", fontsize=TITLE_FONTSIZE, ha='left', va='bottom')
    ax_overlap.set_xticklabels(ax_overlap.get_xticklabels(), rotation=90, fontsize=LABEL_FONTSIZE)
    ax_overlap.set_title('Bases',
                          fontsize=LABEL_FONTSIZE,
                          position=(0, 1 + 0.6 / 10 * FONT_SCALE / 1.5),
                          ha='left', va='bottom')
    
    # Aggregation
    stats_df = pd.read_csv(args.stats, index_col=0, sep='\t')  # data stored when creating the gtf files
    divider_agg = make_axes_locatable(ax_agg)
    title_args = dict(fontsize=LABEL_FONTSIZE, position=(1.0, 1.0), ha='right', va='bottom')

    # Divide axes, plot heatmap and edit axis configuration for each biotype
    for biotype in BIOTYPES[1:]:
        ax_agg_aux = divider_agg.append_axes("right", size="100%", pad=0.3)
        sns.heatmap(sort_by_dendrogram(res_agg_dict[biotype]), annot=True, cbar=False, vmin=0, vmax=agg_vmax,
                    cmap=cmap_agg, ax=ax_agg_aux, fmt='.5g')
        ax_agg_aux.set_title('{} (n={})'.format(biotype, pretty_number(stats_df.loc[biotype, 'genes'])), **title_args)
        ax_agg_aux.set_yticklabels([])
        ax_agg_aux.set_xticklabels(ax_agg_aux.get_xticklabels(), rotation=90, fontsize=LABEL_FONTSIZE)

    if len(BIOTYPES) > 0:
        ax_agg_cbar = divider_agg.append_axes("right", size=0.35, pad=0.3)

        ax_agg.text(0, -0.6 * FONT_SCALE / 1.5, "Aggregation", fontsize=TITLE_FONTSIZE, ha='left', va='bottom')
        g_agg = sns.heatmap(sort_by_dendrogram(res_agg_dict[BIOTYPES[0]]), annot=True, cbar=True, vmin=0, vmax=agg_vmax,
                            cbar_ax=ax_agg_cbar,
                            cmap=cmap_agg, ax=ax_agg, fmt='.5g')
        ax_agg.set_title('{} (n={})'.format(BIOTYPES[0], pretty_number(stats_df.loc[BIOTYPES[0], 'genes'])), **title_args)
        ax_agg.set_yticklabels([])
        ax_agg.set_xticklabels(ax_agg.get_xticklabels(), rotation=90, fontsize=LABEL_FONTSIZE)

        # Edit the colorbar created by the first biotype
        cbar_agg = g_agg.collections[0].colorbar
        cbar_agg.set_ticks([0, agg_vmax])
        cbar_agg.ax.set_yticklabels(['0%', '{:.0f}%'.format(agg_vmax)],
                                    fontsize=LABEL_FONTSIZE)  # the format takes out decimals
    else:
        figure.delaxes(ax_agg)
    
    figure.savefig(args.outfile, bbox_inches='tight')
