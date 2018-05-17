import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from collections import defaultdict
from os import path
from mpl_toolkits.axes_grid1 import make_axes_locatable

from segzoo.gene_biotypes import BIOTYPES

# Prepare the gmtk parameters in a DataFrame
def gmtk_parameters():
    return pd.read_csv(snakemake.input.gmtk, index_col=0)


# Prepare nucleotide results in a Series format
def nucleotide():
    res_nuc_ann = pd.read_table(snakemake.input.nuc, index_col=0)['GC content'].round(2) * 100

    # Rename columns
    res_nuc_ann = res_nuc_ann.rename('GC content (%)')

    # Interpolation of the parameters to rescale them between vmin and vmax for the heatmap
    vmax = 0.65
    vmin = 0.35
    res_nuc_hm = res_nuc_ann.copy()
    res_nuc_hm = ((res_nuc_hm / 100 - vmin) / (vmax - vmin)).clip(0, 1)

    return res_nuc_hm, res_nuc_ann


# Prepare length_distribution results in a DataFrame
def length_distribution():
    # Preparing the annotation for the matrix, creating a new column called 'frac.segs'
    res_len_ann = pd.read_table(snakemake.input.len_dist, index_col=0)
    res_len_ann['frac.segs'] = (res_len_ann['num.segs'] / res_len_ann.loc['all']['num.segs']) * 100
    res_len_ann['frac.bp'] = res_len_ann['frac.bp'] * 100
    res_len_ann = res_len_ann.drop(['num.segs', 'num.bp'], axis=1).drop('all')

    # Rename columns
    res_len_ann.index = res_len_ann.index.map(int)  # labels need to be strings
    res_len_ann.columns = ['Mean length', 'Median length', 'st.dev length', 'Base pairs (%)', 'Segments (%)']

    res_len_hm = res_len_ann.copy()
    # Interpolation of the parameters to rescale them between 0 and 1
    for col in res_len_hm.columns:
        res_len_hm[col] = (res_len_hm[col] - res_len_hm[col].min()) / (res_len_hm[col].max() - res_len_hm[col].min())

    res_len_ann = res_len_ann.round(0)
    return res_len_hm, res_len_ann


# Prepare the mix matrix for the heatmap and its annotation, both in DataFrames
def mix_data_matrix():
    # Joining both matrix to create final heatmap and annotation
    res_nuc_hm, res_nuc_ann = nucleotide()
    res_len_hm, res_len_ann = length_distribution()

    res_ann = res_len_ann.join(res_nuc_ann)
    res_hm = res_len_hm.join(res_nuc_hm)

    return res_hm, res_ann


# Prepare the aggregation results in a dictionary of DataFrames by gene_biotype and return the maximum value
def aggregation():
    # Rename columns
    COLUMN_NAMES = ["5' flanking", "initial exon", "initial intron", "internal exon", "internal introns", \
                    "terminal exon", "terminal intron", "3' flanking"]

    def to_percent(row):
        return (row / row.sum()).round(2) * 100

    df_dict = defaultdict()
    max_value = 0
    for biotype in BIOTYPES:
        filename = next(x for x in snakemake.input.aggs if path.basename(path.dirname(x)) == biotype)
        biotype_df = pd.read_table(filename, index_col=0).apply(to_percent, axis=1).fillna(0)
        biotype_df.columns = COLUMN_NAMES

        # Update max value
        max_value = max(biotype_df.values.max(), max_value)
        df_dict[biotype] = biotype_df

    return df_dict, max_value


NUM_COMPONENTS = 8
GMTK_FACTOR = 1
MIX_FACTOR = 1.5
AGG_FACTOR = 1
ROW_CORRECTOR = 1

LABEL_FONTSIZE = 20
TITLE_FONTSIZE = 25

sns.set(font_scale=1.5)

# Color maps for the visualization
cmap_gmtk = sns.diverging_palette(220, 10, as_cmap=True)

# Call the functions that obtain the results in DataFrames
if snakemake.config['parameters']:
    res_gmtk = gmtk_parameters()
else:
    res_gmtk = pd.DataFrame()
res_mix_hm, res_mix_ann = mix_data_matrix()
res_agg_dict, agg_vmax = aggregation()

GMTK_COL = res_gmtk.shape[1] * GMTK_FACTOR + 1
MIX_COL = res_mix_hm.shape[1] * MIX_FACTOR + 1
AGG_COL = len(BIOTYPES) * NUM_COMPONENTS * AGG_FACTOR + 1

n_rows = res_mix_hm.shape[0] * ROW_CORRECTOR
n_columns = GMTK_COL + MIX_COL + AGG_COL

# Create grid with axes following the ratios desired for the dimensions
f, (ax_gmtk, ax_mix, ax_agg) = \
    plt.subplots(1, 3, figsize=(n_columns, n_rows), \
                 gridspec_kw={"wspace": 3.6 / n_columns, "width_ratios": [GMTK_COL, MIX_COL, AGG_COL]})

# GMTK parameters
if snakemake.config['parameters']:
    g_gmtk = sns.heatmap(res_gmtk, cmap=cmap_gmtk, cbar=True, vmin=0, vmax=1, ax=ax_gmtk)
    cbar_gmtk = g_gmtk.collections[0].colorbar
    cbar_gmtk.set_ticks([0, 1])
    cbar_gmtk.ax.set_yticklabels([0, 1], fontsize=LABEL_FONTSIZE)

    # Setting titles and axis labels
    ax_gmtk.set_yticklabels(ax_gmtk.get_yticklabels(), rotation=0, fontsize=LABEL_FONTSIZE)  # put label names horizontally
    ax_gmtk.set_xticklabels(ax_gmtk.get_xticklabels(), rotation=90, fontsize=LABEL_FONTSIZE)
    ax_gmtk.set_title('GMTK parameters', fontsize=TITLE_FONTSIZE, position=(0, 1 + 0.6 / res_gmtk.shape[0]), ha='left', va='bottom')
else:
    f.delaxes(ax_gmtk)

# Mix matrix
g_mix = sns.heatmap(res_mix_hm, annot=res_mix_ann, cbar=True, cmap='YlGn', vmin=0, vmax=1, ax=ax_mix, fmt='.5g')
cbar_mix = g_mix.collections[0].colorbar
cbar_mix.set_ticks([0, 1])
cbar_mix.ax.set_yticklabels(['low', 'high'], fontsize=LABEL_FONTSIZE)
ax_mix.set_ylabel('')
ax_mix.set_xticklabels(ax_mix.get_xticklabels(), rotation=90, fontsize=LABEL_FONTSIZE)

if snakemake.config['parameters']:
    ax_mix.set_yticklabels([])
else:
    ax_mix.set_yticklabels(ax_mix.get_yticklabels(), rotation=0, fontsize=LABEL_FONTSIZE)

# Aggregation
stats_df = pd.read_table(snakemake.input.stats, index_col=0)  # data stored when creating the gtf files
divider = make_axes_locatable(ax_agg)
title_args = dict(fontsize=LABEL_FONTSIZE, position=(1.0, 1.0), ha='right', va='bottom')

# Divide axes, plot heatmap and edit axis configuration for each biotype
for biotype in BIOTYPES[1:]:
    ax_agg_aux = divider.append_axes("right", size="100%", pad=0.3)
    sns.heatmap(res_agg_dict[biotype], annot=True, cbar=False, vmin=0, vmax=agg_vmax, cmap='Blues', ax=ax_agg_aux,
                fmt='.5g')
    ax_agg_aux.set_title('{} ({})'.format(biotype, stats_df.loc[biotype, 'genes']), **title_args)
    ax_agg_aux.set_yticklabels([])
    ax_agg_aux.set_xticklabels(ax_agg_aux.get_xticklabels(), rotation=90, fontsize=LABEL_FONTSIZE)

if len(BIOTYPES) > 0:
    ax_agg_cbar = divider.append_axes("right", size=0.35, pad=0.3)

    ax_agg.text(0, -0.6, "Aggregation", fontsize=TITLE_FONTSIZE, ha='left', va='bottom')
    g_agg = sns.heatmap(res_agg_dict[BIOTYPES[0]], annot=True, cbar=True, vmin=0, vmax=agg_vmax, cbar_ax=ax_agg_cbar,
                        cmap='Blues', ax=ax_agg, fmt='.5g')
    ax_agg.set_title('{} ({})'.format(BIOTYPES[0], stats_df.loc[BIOTYPES[0], 'genes']), **title_args)
    ax_agg.set_yticklabels([])
    ax_agg.set_xticklabels(ax_agg.get_xticklabels(), rotation=90, fontsize=LABEL_FONTSIZE)

    # Edit the colorbar created by the first biotype
    cbar_agg = g_agg.collections[0].colorbar
    cbar_agg.set_ticks([0, agg_vmax])
    cbar_agg.ax.set_yticklabels(['0%', '{:.0f}%'.format(agg_vmax)], fontsize=LABEL_FONTSIZE)  # the format takes out decimals
else:
    f.delaxes(ax_agg)

f.savefig(snakemake.output.outfile, bbox_inches='tight')
