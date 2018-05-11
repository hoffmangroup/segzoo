import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict
from os import path
from mpl_toolkits.axes_grid1 import make_axes_locatable


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
def aggregation(biotypes):
    # Rename columns
    COLUMN_NAMES = ["5' flanking", "initial exon", "initial intron", "internal exon", "internal introns", \
                    "terminal exon", "terminal intron", "3' flanking"]

    df_dict = defaultdict()
    max_value = 0
    for biotype in biotypes:
        filename = next(x for x in snakemake.input.aggs if path.basename(path.dirname(x)) == biotype)
        df_dict[biotype] = pd.read_table(filename, index_col=0)
        df_dict[biotype].columns = COLUMN_NAMES

        # Turn the results into a percentage of the labels's aggregation in a gene biotype
        for row in df_dict[biotype].index:
            df_dict[biotype].iloc[row] = (df_dict[biotype].iloc[row] / df_dict[biotype].iloc[row].sum()).round(2) * 100
        # Update max value
        max_value = max(df_dict[biotype].max().max(), max_value)

    return df_dict, max_value


# The biotypes wanted for the visualization, in order of appearance
BIOTYPES = [
    #     '3prime_overlapping_ncrna',
    #     'antisense',
    #     'IG_C_gene',
    #     'IG_C_pseudogene',
    #     'IG_D_gene',
    #     'IG_J_gene',
    #     'IG_J_pseudogene',
    #     'IG_V_gene',
    #     'IG_V_pseudogene',
    #     'known_ncrna',
    'lincRNA',
    #     'miRNA',
    #     'misc_RNA',
    #     'Mt_rRNA',
    #     'Mt_tRNA',
    #     'non_coding',
    #     'polymorphic_pseudogene',
    #     'processed_pseudogene',
    #     'processed_transcript',
    'protein_coding',
    #     'pseudogene',
    #     'rRNA',
    #     'sense_intronic',
    #     'sense_overlapping',
    #     'snoRNA',
    #     'snRNA',
    #     'TEC',
    #     'transcribed_processed_pseudogene',
    #     'transcribed_unitary_pseudogene',
    #     'transcribed_unprocessed_pseudogene',
    #     'translated_processed_pseudogene',
    #     'translated_unprocessed_pseudogene',
    #     'TR_C_gene',
    #     'TR_D_gene',
    #     'TR_J_gene',
    #     'TR_J_pseudogene',
    #     'TR_V_gene',
    #     'TR_V_pseudogene',
    #     'unitary_pseudogene',
    #     'unprocessed_pseudogene',
]

NUM_COMPONENTS = 8
GMTK_FACTOR = 1
MIX_FACTOR = 1.5
AGG_FACTOR = 1.8
ROW_CORRECTOR = 0.8

sns.set(font_scale=1.5)

# Color maps for the visualization
cmap_gmtk = sns.diverging_palette(220, 10, as_cmap=True)

# Call the functions that obtain the results in DataFrames
if snakemake.config['parameters']:
    res_gmtk = gmtk_parameters()
else:
    res_gmtk = pd.DataFrame()
res_mix_hm, res_mix_ann = mix_data_matrix()
res_agg_dict, agg_vmax = aggregation(BIOTYPES)

n_rows = res_mix_hm.shape[0] * ROW_CORRECTOR
n_columns = res_gmtk.shape[1] * GMTK_FACTOR + res_mix_hm.shape[1] * MIX_FACTOR + len(
    res_agg_dict) * NUM_COMPONENTS * AGG_FACTOR

# Create grid and plot the results in two or three separate plots
f, (ax_gmtk, ax_mix, ax_agg) = \
    plt.subplots(1, 3, figsize=(n_columns, n_rows), sharey=True, \
                 gridspec_kw={"wspace": 0.08,
                              "width_ratios": [res_gmtk.shape[1] * GMTK_FACTOR, res_mix_hm.shape[1] * MIX_FACTOR,
                                               len(res_agg_dict) * NUM_COMPONENTS * AGG_FACTOR]})

# GMTK parameters
if snakemake.config['parameters']:
    g_gmtk = sns.heatmap(res_gmtk, cmap=cmap_gmtk, cbar=True, vmin=0, vmax=1, ax=ax_gmtk)
    cbar_gmtk = g_gmtk.collections[0].colorbar
    cbar_gmtk.set_ticks([0, 1])

    # Setting titles and axis labels
    ax_gmtk.set_yticklabels(ax_gmtk.get_yticklabels(), rotation=0)  # put label names horizontally
    ax_gmtk.set_title('GMTK parameters', fontsize=25, position=(0, 1 + 0.5 / res_gmtk.shape[0]), ha='left', va='bottom')
else:
    f.delaxes(ax_gmtk)

# Mix matrix
g_mix = sns.heatmap(res_mix_hm, annot=res_mix_ann, cbar=True, cmap='YlGn', vmin=0, vmax=1, ax=ax_mix, fmt='.5g')
cbar_mix = g_mix.collections[0].colorbar
cbar_mix.set_ticks([0, 1])
cbar_mix.set_ticklabels(['low', 'high'])
ax_mix.set_ylabel('')

# Aggregation
ax_agg_list = [None] * len(BIOTYPES)  # Structure to store the divided axis
ax_agg_list[0] = ax_agg

divider = make_axes_locatable(ax_agg_list[0])
title_args = dict(fontsize=20, position=(1.0, 1.0), ha='right', va='bottom')

# Divide axes, plot heatmap and edit axis configuration for each biotype
for i in range(1, len(res_agg_dict)):
    ax_agg_list[i] = divider.append_axes("right", size="100%", pad=0.3)
    sns.heatmap(res_agg_dict[BIOTYPES[i]], annot=True, cbar=False, vmin=0, vmax=agg_vmax, cmap='Blues',
                ax=ax_agg_list[i], fmt='.5g')
    ax_agg_list[i].set_title(BIOTYPES[i], **title_args)
    ax_agg_list[i].set_yticklabels([])
    ax_agg_list[i].set_xticklabels(ax_agg_list[i].get_xticklabels(), rotation=90)
ax_agg_cbar = divider.append_axes("right", size="3%", pad=0.3)

g_agg = sns.heatmap(res_agg_dict[BIOTYPES[0]], annot=True, cbar=True, vmin=0, vmax=agg_vmax, cbar_ax=ax_agg_cbar,
                    cmap='Blues', ax=ax_agg_list[0], fmt='.5g')
ax_agg_list[0].text(0, -0.5, "Aggregation (%)", fontsize=25, ha='left', va='bottom')
ax_agg_list[0].set_title(BIOTYPES[0], **title_args)
ax_agg_list[i].set_xticklabels(ax_agg_list[i].get_xticklabels(), rotation=90)
# do not set yticklabels to [], because that will affect GMTK parameters, as it's set to yshare
# Edit the colorbar created by the first biotype
cbar_agg = g_agg.collections[0].colorbar
cbar_agg.set_ticks([0, agg_vmax])
cbar_agg.set_ticklabels(['0%', '{0}%'.format(agg_vmax)])

f.savefig(snakemake.output.outfile, bbox_inches='tight')
