#!/usr/bin/env python
# coding: utf-8

import sys
from os import path
import math
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt

import pandas as pd
import seaborn as sns

import scipy.cluster.hierarchy as sch

from collections import defaultdict

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.transforms import ScaledTranslation

from segzoo.gene_biotypes import BIOTYPES

mpl.use('Agg')

# VARIABLES AND CONSTANTS

# Heatmap proportion modifiers
NUM_COMPONENTS = 6
GMTK_FACTOR = 1
MIX_FACTOR = 1.5
AGG_FACTOR = 1
OVERLAP_FACTOR = 1
OVERLAP_COLUMN_NUMBER = 2
ROW_FACTOR = 1

# Font scaling variables
FONT_SCALE = 1.5
sns.set(font_scale=FONT_SCALE)

LABEL_FONTSIZE = 20 * FONT_SCALE / 1.5
MAIN_TITLE_FONTSIZE = 25 * FONT_SCALE / 1.5
SUBTITLE_FONTSIZE = 20 * FONT_SCALE / 1.5

# Set the vertical positions of the titles.
TOP_TITLE_Y = -0.6 * FONT_SCALE / 1.5
BOTTOM_TITLE_Y = -0.1 * FONT_SCALE / 1.5

plt.rc('xtick', labelsize=LABEL_FONTSIZE)
plt.rc('ytick', labelsize=LABEL_FONTSIZE)

YTICKLABELS = False
INDEX_AXES_X_COORD = -0.1
INDEX_AXES_Y_COORD = 1.12
INDEX_AXES_X_DATACOORD = -0.3

# Table options and properties
TABLE_POS = "bottom"  # top / bottom / other to omit
TABLE_HEIGHT = 1  # relative to the height of 2 rows from the mix matrix
MIX_TABLE_CONTENT = [['max', 'max', 'max', 'max', 'max', 65, 'max'],
                     ['min', 'min', 'min', 'min', 'min', 35, 'min']]

# Color maps for the visualization
cmap_gmtk = sns.diverging_palette(220, 10, as_cmap=True)
cmap_mix = 'YlGn'
cmap_agg = 'Blues'
cmap_overlap = 'Reds'


def is_decimal_zero(num):
    """
    True if all digits after num's decimal point are 0s.
    False otherwise

    >>> is_decimal_zero(12.0)
    True
    >>> is_decimal_zero(12.3)
    False
    """
    return int(num) == num


def more_than_n_digits(num, n_digits=2):
    """
    True if the number num has more than n digits.
    False otherwise

    >>> more_than_n_digits(12.123456)
    False
    >>> more_than_n_digits(123.123456)
    True
    >>> more_than_n_digits(1.23, 1)
    False
    """
    return len(str(int(num))) > n_digits


def human_format(num):
    """
    Shorten long numbers by replacing ending zeroes by a magnitude unit.

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

    if more_than_n_digits(num) or is_decimal_zero(num):
        return '{}{}'.format(int(num), magnit_chars[magnitude])
    else:
        return '{:.2f}{}'.format(num, magnit_chars[magnitude])


def prettify_number(n):
    """
    Add space every three digits from right to left

    >>> prettify_number(1000)
    '1 000'
    >>> prettify_number(100)
    '100'
    """
    return '{:,}'.format(int(n)).replace(',', ' ')


def gmtk_parameters(args):
    """
    Prepare the gmtk parameters in a DataFrame.
    """
    if not args.gmtk:
        return pd.DataFrame(), pd.DataFrame()

    def normalize_col(col):
        return (col - col.min()) / (col.max() - col.min())

    df = pd.read_csv(args.gmtk, index_col=0, sep='\t')
    df.sort_index(inplace=True)
    df.index = df.index.astype(str)

    table_content_df = pd.concat([
        df.max().apply(human_format),
        df.min().apply(human_format)],
        axis=1).T

    if args.normalize_gmtk:
        df = df.apply(normalize_col, axis=0)

    return df, table_content_df


def nucleotide(args):
    """
    Prepare nucleotide results in a Series format
    """
    res_nuc_ann = pd.read_csv(args.nuc, index_col=0, sep='\t')['GC content'].round(2) * 100
    res_nuc_ann.sort_index(inplace=True)
    res_nuc_ann = res_nuc_ann.rename('G+C content (%)')
    res_nuc_ann.index = res_nuc_ann.index.astype(str)

    # Interpolation of the parameters to rescale them between vmin and vmax for the heatmap
    vmax = 65
    vmin = 35
    res_nuc_hm = res_nuc_ann.copy()
    res_nuc_hm = ((res_nuc_hm - vmin) / (vmax - vmin)).clip(0, 1)

    return res_nuc_hm, res_nuc_ann


def length_distribution(args):
    """
    Prepare length_distribution results in a DataFrame
    """
    # Preparing the annotation for the matrix, creating a new column called 'frac.segs'
    res_len_ann = pd.read_csv(args.len_dist, index_col=0, sep='\t')
    res_len_ann['frac.segs'] = (res_len_ann['num.segs'] / res_len_ann.loc['all']['num.segs']) * 100
    res_len_ann['frac.bp'] = res_len_ann['frac.bp'] * 100
    res_len_ann = res_len_ann.drop(['num.segs', 'num.bp'], axis=1).drop('all')

    # Rename columns
    res_len_ann.index = res_len_ann.index.astype(str)
    res_len_ann.sort_index(inplace=True)
    res_len_ann.columns = ['Mean length', 'Median length', 'Std length', 'Bases (%)', 'Segments (%)']
    col_order = ['Segments (%)', 'Mean length', 'Median length', 'Std length', 'Bases (%)']
    res_len_ann = res_len_ann[col_order]
    res_len_hm = res_len_ann.copy()
    # Interpolation of the parameters to rescale them between 0 and 1
    for col in res_len_hm.columns:
        res_len_hm[col] = (res_len_hm[col] - res_len_hm[col].min()) / (res_len_hm[col].max() - res_len_hm[col].min())

    res_len_ann = res_len_ann.round(0)
    return res_len_hm, res_len_ann


def phastcons(args):
    """
    Prepare phastcons result in a DataFrame
    """
    res_ann = pd.read_csv(args.phastcons, index_col=0, sep='\t')
    res_ann.columns = ['Mean PhastCons', 'Max PhastCons']
    res_ann = res_ann[['Mean PhastCons']]
    res_ann.index = res_ann.index.astype(str)

    res_hm = (res_ann - res_ann.min()) / (res_ann.max() - res_ann.min())

    return res_hm, res_ann


def repeatmasker(args):
    """
    Prepare repeatmasker result in a DataFrame
    """
    df = pd.read_csv(args.repeatmasker, index_col=0, sep='\t')
    df.sort_index(inplace=True)
    df = df * 100
    df = df.apply(round).astype(int)
    df.columns = ['RepeatMasker']
    df.index = df.index.astype(str)
    return df


def mix_data_matrix(args):
    """
    Prepare the Mix DataFrame.
    """
    # Joining the matrices to create final heatmap and annotation
    res_nuc_hm, res_nuc_ann = nucleotide(args)
    res_len_hm, res_len_ann = length_distribution(args)
    res_phastcons_hm, res_phastcons_ann = phastcons(args)

    res_ann = res_len_ann.join(res_nuc_ann)
    res_ann = res_ann.join(res_phastcons_ann)
    res_ann.index = res_ann.index.astype(str)

    res_hm = res_len_hm.join(res_nuc_hm)
    res_hm = res_hm.join(res_phastcons_hm)
    res_hm.index = res_hm.index.astype(str)
    return res_hm, res_ann.applymap(human_format)


def overlap(args):
    """
    Prepare the Overlap DataFrame.
    """
    genic_df = genic_intergenic(args)
    repeatmasker_df = repeatmasker(args)

    res_df = genic_df.join(repeatmasker_df)
    res_df.index = res_df.index.astype(str)
    return res_df


def genic_intergenic(args):
    """
    Prepare the overlap results in Dataframe
    """
    df = pd.read_csv(args.overlap, sep='\t', header=0, index_col=0)
    df.sort_index(inplace=True)
    df = df * 100
    df = df.apply(round).astype(int)
    df.columns = ['Genic', 'Intergenic']
    df.index = df.index.astype(str)
    return df


def aggregation(args):
    """
    Prepare the aggregation results in a dictionary of DataFrames by gene_biotype and return the maximum value
    """
    # Rename columns
    column_names = ["5' flanking", "initial exon", "initial intron", "internal exon", "internal introns",
                    "terminal exon", "terminal intron", "3' flanking"]

    def to_percent(row):
        return (row / row.sum()).round(2) * 100

    df_dict = defaultdict()
    max_value = 0
    for biotype in BIOTYPES:
        filename = next(x for x in args.aggs if path.basename(path.dirname(x)) == biotype)
        biotype_df = pd.read_csv(filename, index_col=0, sep='\t').apply(to_percent, axis=1).fillna(0)
        biotype_df.sort_index(inplace=True)
        biotype_df.index = biotype_df.index.astype(str)
        biotype_df.columns = column_names

        # Update max value
        max_value = max(biotype_df.values.max(), max_value)
        df_dict[biotype] = biotype_df

    return df_dict, max_value


def get_mne_ticklabels(filename, track_labels=(), label_labels=()):
    """
    Parse mne file and return updated tracks and labels
    """
    if filename:
        mne_df = pd.read_csv(filename, dtype=str, sep='\t')
    else:
        mne_df = pd.DataFrame(columns=['type', 'old', 'new'])

    mne_df.sort_index(inplace=True)
    assert all(col in ['type', 'old', 'new'] for col in mne_df.columns)

    track_df = mne_df[mne_df.type == 'track']
    track_labels = list(map(str, track_labels))
    track_translator = dict(zip(track_df.old, track_df.new))
    # if tracks are missing in the mne file, add them in the translator with new value = old value
    track_translator = {old: track_translator.get(old, old) for old in track_labels}
    new_tracks = [track_translator[old] for old in track_labels]

    label_df = mne_df[mne_df.type == 'label']
    label_labels = list(map(str, label_labels))
    label_translator = dict(zip(label_df.old, label_df.new))
    # if labels are missing in the mne file, add them in the translator with new value = old value
    label_translator = {old: label_translator.get(old, old) for old in label_labels}
    new_labels = [label_translator[old] for old in label_labels]

    return new_tracks, new_labels, track_translator, label_translator


def calc_dendrogram_label_col(labels, zero_threshold=4, one_threshold=10, two_threshold=14, increment=4):
    """
    Return how many columns the invisible ax between dendrogram and gmtk parameters
    should have based on the length of the longest label in `labels`
    The defaults are set based on LABEL_FONTSIZE

    `increment` is the number of characters occupying one ax column

    >>> calc_dendrogram_label_col([12345])
    1
    >>> calc_dendrogram_label_col([1234567890123456])
    3
    """
    longest_label_len = max(len(str(item)) for item in labels)
    threshold_to_ncol = ((zero_threshold, 0), (one_threshold, 1), (two_threshold, 2))
    for threshold, ncol in threshold_to_ncol:
        if longest_label_len <= threshold:
            return ncol
    # Add one column for every increment increase in label length, rounding up
    # Add one to the numerator to make the visual look nicer
    return 2 + math.ceil((longest_label_len - two_threshold + 1) / increment)


def generate_table(ax, cbar, table_content, table_height, figure):
    """Generate table underneath a heatmap plot"""
    n_cols = len(table_content[0])
    high_low_table = ax.table(
        cellText=table_content,
        cellColours=[[cbar.cmap(0.99)] * n_cols, [cbar.cmap(0.01)] * n_cols],
        bbox=[0, - (TABLE_HEIGHT + .25) / table_height, 1, TABLE_HEIGHT / table_height],  # [left,bottom,width,height]
        fontsize=LABEL_FONTSIZE,
        cellLoc='center')
    for j in range(n_cols):
        high_low_table._cells[(0, j)]._text.set_color('white')  # TODO: do not access protected variables

    # Offset labels down to leave space for the table
    dx = 0
    dy = -(TABLE_HEIGHT + 0.25) * 55 / 72
    offset = ScaledTranslation(dx, dy, figure.dpi_scale_trans)

    for label in ax.xaxis.get_majorticklabels():
        label.set_transform(label.get_transform() + offset)


def main(args):

    # update environment variables
    global OVERLAP_COLUMN_NUMBER
    global NUM_COMPONENTS
    if args.repeatmasker:
        OVERLAP_COLUMN_NUMBER += 1
    if args.phastcons:
        NUM_COMPONENTS += 1

    # Call the functions that obtain the results in DataFrames
    res = {}

    res_gmtk, gmtk_table_content_df = gmtk_parameters(args)
    res['gmtk'] = res_gmtk

    res_mix_hm, res_mix_ann = mix_data_matrix(args)
    res['mix_hm'] = res_mix_hm
    res['mix_ann'] = res_mix_ann

    res['overlap'] = overlap(args)

    res_agg_dict, agg_vmax = aggregation(args)
    for biotype in BIOTYPES:
        res[f'agg_{biotype}'] = res_agg_dict[biotype]

    # Read labels from mne file
    new_tracks, new_labels, track_translator, label_translator = \
        get_mne_ticklabels(args.mne, res['gmtk'].columns, res['mix_hm'].index)

    if args.dendrogram:
        row_linkage_matrix = sch.linkage(res['gmtk'], method='weighted')
        row_dendrogram = sch.dendrogram(row_linkage_matrix, no_plot=True)
        row_ordering = [label_translator[str(l)] for l in row_dendrogram['leaves']][::-1]

        col_linkage_matrix = sch.linkage(res['gmtk'].transpose(), method='weighted')
        col_dendrogram = sch.dendrogram(col_linkage_matrix, no_plot=True)
        col_ordering = [res['gmtk'].columns[leaf_index] for leaf_index in col_dendrogram['leaves']]
        res['gmtk'].columns = res['gmtk'].columns.map(track_translator)
        res['gmtk'] = res['gmtk'][col_ordering]

    else:
        row_ordering = new_labels

    for key, df in res.items():
        df.index = df.index.map(label_translator)
        res[key] = df.loc[row_ordering]

    # Dimensioning variables
    dendrogram_col = 2
    dendrogram_labels_col = calc_dendrogram_label_col(new_labels) if args.dendrogram else 0
    # +1 for the color bar ?
    gmtk_col = res['gmtk'].shape[1] * GMTK_FACTOR + 1
    mix_col = res['mix_hm'].shape[1] * MIX_FACTOR + 1
    overlap_col = OVERLAP_COLUMN_NUMBER * OVERLAP_FACTOR + 1
    agg_col = len(BIOTYPES) * NUM_COMPONENTS * AGG_FACTOR + 1

    table_height = res['mix_hm'].shape[0] * ROW_FACTOR
    n_columns = dendrogram_col + dendrogram_labels_col + gmtk_col + mix_col + overlap_col + agg_col

    # If do not add space in between dendrogram and gmtk parameters, move the space ax to the left
    # to avoid adding extra space
    if dendrogram_labels_col:
        width_ratios = [dendrogram_col, dendrogram_labels_col - 1, gmtk_col, mix_col, overlap_col, agg_col]
    else:
        width_ratios = [dendrogram_col, gmtk_col, mix_col, overlap_col, agg_col]

    # Create grid with axes following the ratios desired for the dimensions
    figure, axes = plt.subplots(1, len(width_ratios),
                                figsize=(n_columns, table_height),
                                gridspec_kw={"wspace": 9 / n_columns,
                                             "width_ratios": width_ratios})

    # If the user does not choose to add space in between dendrogram and gmtk parameters, move the space ax to the left
    # to avoid adding extra space
    if dendrogram_labels_col:
        ax_dendrogram, ax_dendrogram_labels, ax_gmtk, ax_mix, ax_overlap, ax_agg = axes
    else:
        ax_dendrogram, ax_gmtk, ax_mix, ax_overlap, ax_agg = axes
        ax_dendrogram_labels = None

    # GMTK parameters
    if args.gmtk:
        if args.dendrogram:
            # Row-wise hierarchical clustering with dendrogram
            sch.dendrogram(row_linkage_matrix, ax=ax_dendrogram, orientation='left',
                           color_threshold=0, above_threshold_color='k',
                           leaf_font_size=int(LABEL_FONTSIZE),
                           labels=new_labels)

            ax_dendrogram.spines['right'].set_visible(False)
            ax_dendrogram.spines['left'].set_visible(False)
            ax_dendrogram.spines['top'].set_visible(False)
            ax_dendrogram.spines['bottom'].set_visible(False)
            ax_dendrogram.set_facecolor((1, 1, 1))  # Set dendrogram background to white
            ax_dendrogram.set_xticklabels('')
            ax_dendrogram.set_ylabel('Label')

            ax_dendrogram.text(INDEX_AXES_X_COORD, INDEX_AXES_Y_COORD, "a", fontsize=LABEL_FONTSIZE + 2, ha='right',
                               va='top', transform=ax_dendrogram.transAxes)

            if dendrogram_labels_col:
                ax_dendrogram_labels.set_facecolor('None')
                ax_dendrogram_labels.set_xticklabels('')
                ax_dendrogram_labels.set_yticklabels(())
        else:
            figure.delaxes(ax_dendrogram)

        divider_gmtk = make_axes_locatable(ax_gmtk)
        ax_gmtk_cbar = divider_gmtk.append_axes("right", size=0.35, pad=0.3)
        g_gmtk = sns.heatmap(res['gmtk'], cmap=cmap_gmtk, ax=ax_gmtk, cbar_ax=ax_gmtk_cbar)

        cbar_gmtk = g_gmtk.collections[0].colorbar
        cbar_gmtk.set_ticks([res['gmtk'].max().max(), res['gmtk'].min().min()])

        if args.normalize_gmtk:
            cbar_gmtk.ax.set_yticklabels(['col\nmax', 'col\nmin'])

            # Add min-max table for parameters matrix
            gmtk_table_content = gmtk_table_content_df.values.tolist()
            generate_table(ax_gmtk, cbar_gmtk, gmtk_table_content, table_height, figure)

        # Setting titles and axis labels
        if args.dendrogram:
            ax_gmtk.set_yticklabels(())
            ax_gmtk.set_ylabel('')
        else:
            ax_gmtk.set_ylabel('Label')
            ax_gmtk.tick_params(axis='y', rotation=0)
            ax_gmtk.text(INDEX_AXES_X_DATACOORD, INDEX_AXES_Y_COORD, "a", fontsize=LABEL_FONTSIZE + 2, ha='right',
                         va='top', transform=ax_gmtk.get_xaxis_transform())

        ax_gmtk.tick_params(axis='x', rotation=90)
        ax_gmtk.text(x=0, y=BOTTOM_TITLE_Y, s='Parameters', fontsize=MAIN_TITLE_FONTSIZE, ha='left', va='bottom')

    else:
        figure.delaxes(ax_gmtk)
        figure.delaxes(ax_dendrogram)

    # Mix matrix
    divider_mix = make_axes_locatable(ax_mix)
    ax_mix_cbar = divider_mix.append_axes("right", size=0.35, pad=0.3)
    g_mix = sns.heatmap(res['mix_hm'], annot=res['mix_ann'], cbar=True, cmap=cmap_mix, vmin=0, vmax=1,
                        ax=ax_mix, cbar_ax=ax_mix_cbar, fmt='')
    g_mix.text(INDEX_AXES_X_DATACOORD, INDEX_AXES_Y_COORD, "b", fontsize=LABEL_FONTSIZE + 2, ha='left',
               va='top', transform=g_mix.get_xaxis_transform())
    cbar_mix = g_mix.collections[0].colorbar
    cbar_mix.set_ticks([0, 1])
    cbar_mix.ax.set_yticklabels(['low', 'high'])

    ax_mix.tick_params(axis='x', rotation=90)
    if args.gmtk:
        ax_mix.set_ylabel('')
        ax_mix.set_yticklabels(())
    else:
        ax_mix.set_ylabel('Label')
        ax_mix.tick_params(axis='y', rotation=0)

    ax_mix.vlines(x=4, ymin=0, ymax=res['mix_hm'].shape[0], colors='black', lw=1)

    ax_mix.text(0, BOTTOM_TITLE_Y, "Segments", fontsize=SUBTITLE_FONTSIZE, ha='left', va='bottom')
    ax_mix.text(4, BOTTOM_TITLE_Y, "Bases", fontsize=SUBTITLE_FONTSIZE, ha='left', va='bottom')
    # Add min-max table for mix matrix
    generate_table(ax_mix, cbar_mix, MIX_TABLE_CONTENT, table_height, figure)

    # Overlap
    divider_overlap = make_axes_locatable(ax_overlap)
    ax_overlap_cbar = divider_overlap.append_axes("right", size=0.35, pad=0.3)
    g_overlap = sns.heatmap(res['overlap'], vmin=0, vmax=100, annot=True, cbar=True, fmt='.5g', yticklabels=YTICKLABELS,
                            cmap=cmap_overlap, ax=ax_overlap, cbar_ax=ax_overlap_cbar)

    cbar_overlap = g_overlap.collections[0].colorbar
    cbar_overlap.set_ticks([0, 100])
    cbar_overlap.ax.set_yticklabels(['0%', '100%'])

    ax_overlap.set_ylabel('')

    ax_overlap.text(INDEX_AXES_X_DATACOORD, INDEX_AXES_Y_COORD, "c", fontsize=LABEL_FONTSIZE + 2, ha='right',
                     va='top', transform=ax_overlap.get_xaxis_transform())
    ax_overlap.text(x=0, y=TOP_TITLE_Y, s="Overlap", fontsize=MAIN_TITLE_FONTSIZE, ha='left', va='bottom')
    ax_overlap.tick_params(axis='x', rotation=90)
    ax_overlap.text(x=0, y=BOTTOM_TITLE_Y, s='Bases', fontsize=SUBTITLE_FONTSIZE, ha='left', va='bottom')

    # Aggregation
    stats_df = pd.read_csv(args.stats, index_col=0, sep='\t')  # data stored when creating the gtf files
    divider_agg = make_axes_locatable(ax_agg)

    # We want one axis per biotype and the axis for the first biotype already exists: `ax_agg`.
    # Create the remaining axes.
    rem_axs = [divider_agg.append_axes("right", size="100%", pad=0.3)] * (len(BIOTYPES) - 1)
    ax_aggs = [ax_agg] + rem_axs

    # Divide axes, plot heatmap and edit axis configuration for each biotype
    for cur_ax, biotype in zip(ax_aggs, BIOTYPES):
        agg_df = res[f'agg_{biotype}']
        n_agg_columns = agg_df.shape[1]
        sns.heatmap(agg_df, annot=True, cbar=False, vmin=0, vmax=agg_vmax,
                    cmap=cmap_agg, ax=cur_ax, fmt='.5g', yticklabels=YTICKLABELS)
        n_gene = prettify_number(stats_df.loc[biotype, 'genes'])
        subtitle = f'{biotype} (n={n_gene})'
        cur_ax.text(x=n_agg_columns, y=BOTTOM_TITLE_Y, s=subtitle, ha='right', va='bottom')
        cur_ax.set_ylabel('')
        cur_ax.tick_params(axis='x', rotation=90)

    if len(BIOTYPES) > 0:
        # Add title
        ax_agg.text(x=0, y=TOP_TITLE_Y, s="Aggregation", fontsize=MAIN_TITLE_FONTSIZE, ha='left', va='bottom')

        # Add colorbar
        ax_agg_cbar = divider_agg.append_axes("right", size=0.35, pad=0.3)

        ax_agg.text(INDEX_AXES_X_DATACOORD, INDEX_AXES_Y_COORD, "d", fontsize=LABEL_FONTSIZE + 2, ha='right',
                    va='top', transform=ax_agg.get_xaxis_transform())
        ax_agg.set_ylabel('')

        # Edit the colorbar created by the first biotype
        cbar_agg = mpl.colorbar.ColorbarBase(ax_agg_cbar, cmap=cmap_agg, orientation='vertical')
        cbar_agg.set_ticks([0, 1])
        cbar_agg.ax.set_yticklabels(['0%', '{:.0f}%'.format(agg_vmax)])  # the format takes out decimals
    else:
        figure.delaxes(ax_agg)

    figure.savefig(args.outfile + '.png', bbox_inches='tight', dpi=350)
    figure.savefig(args.outfile + '.pdf', bbox_inches='tight')


def parse_args(args):
    """
    Parse command line arguments and return them as a Namespace object
    """

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
    parser.add_argument('--normalize-gmtk', action='store_true', help='If set, normalize gmtk parameters column wise')
    parser.add_argument('--dendrogram', action='store_true',
                        help='If set, perform hierarchical clustering of GMTK parameters table row-wise')
    parser.add_argument('--nuc', help='Nucleotide results file')
    parser.add_argument('--len_dist', help='Length distribution statistics')
    parser.add_argument('--genic', help='The percentage of segments that overlap with a gene')
    parser.add_argument('--mne', help='Allows specify an mne file to translate segment '
                                      'labels and track names on the shown on the figure')
    parser.add_argument('--aggs', help='Aggregation results file')
    parser.add_argument('--stats', help='Gene biotype stats')
    parser.add_argument('--phastcons', help='Phastcons result file')
    parser.add_argument('--repeatmasker', help='RepeatMasker result file')
    parser.add_argument('--outfile', help='The path of the resulting visualization, excluding file extension')
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
                     '--genic', snakemake.input.olp,
                     '--mne', snakemake.input.mne,
                     '--aggs', snakemake.input.aggs,
                     '--stats', snakemake.input.stats,
                     '--phastcons', snakemake.input.phastcons,
                     '--repeatmasker', snakemake.input.repeatmasker,
                     '--outfile', snakemake.params.outfile
                     ]
        args = parse_args(arg_list)
    else:
        args = parse_args(sys.argv[1:])
        args.aggs = args.aggs.split(',')

    main(args)


