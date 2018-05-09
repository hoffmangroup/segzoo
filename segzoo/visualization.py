import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
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


# Prepare the main matrix for the heatmap and its annotation, both in DataFrames
def mix_data_matrix():
    # Joining both matrix to create final heatmap and annotation
    res_nuc_hm, res_nuc_ann = nucleotide()
    res_len_hm, res_len_ann = length_distribution()

    res_ann = res_len_ann.join(res_nuc_ann)
    res_hm = res_len_hm.join(res_nuc_hm)

    return res_hm, res_ann


# Prepare the aggregation results for 'protein_coding' and 'lincRNA' and join them in a DataFrame
def aggregation():
    # TODO don't join
    pc_file = next(x for x in snakemake.input.aggs if path.basename(path.dirname(x)) == 'protein_coding')
    lrna_file = next(x for x in snakemake.input.aggs if path.basename(path.dirname(x)) == 'lincRNA')

    # Rename columns
    column_names = ["5' flanking", "initial exon", "initial intron", "internal exon", "internal introns", \
                    "terminal exon", "terminal intron", "3' flanking"]

    res_agg_pc = pd.read_table(pc_file, index_col=0)
    res_agg_pc.columns = column_names

    res_agg_linc = pd.read_table(lrna_file, index_col=0)
    res_agg_linc.columns = column_names

    res_agg = res_agg_pc.join(res_agg_linc, lsuffix='\n protein coding', rsuffix='\n lincRNA')

    # Turn the results into a percentage of the labels's aggregation in a gene biotype
    for row in res_agg.index:
        res_agg.iloc[row, :8] = res_agg.iloc[row, :8] / res_agg.iloc[row, :8].sum()
        res_agg.iloc[row, 8:] = res_agg.iloc[row, 8:] / res_agg.iloc[row, 8:].sum()

    res_agg_ann = res_agg.copy().round(2) * 100  # Turn into percentage
    # For the heatmap, map the maximum value of each gene_biotype to 1
    res_agg.iloc[:, :8] = res_agg.iloc[:, :8] / res_agg.iloc[:, :8].max().max()
    res_agg.iloc[:, 8:] = res_agg.iloc[:, 8:] / res_agg.iloc[:, 8:].max().max()

    return res_agg, res_agg_ann


sns.set(font_scale=1.5)

# Color maps for the visualization
cmap_gmtk = sns.diverging_palette(220, 10, as_cmap=True)
# cmap_mix = sns.diverging_palette(176, 360, as_cmap=True)
# cmap_agg = sns.diverging_palette(46, 360, as_cmap=True)
# cmap_mix = sns.light_palette('lightblue', as_cmap=True)
# cmap_agg = sns.light_palette('lightgreen', as_cmap=True)

if snakemake.config['parameters']:
    res_gmtk = gmtk_parameters()
else:
    res_gmtk = pd.DataFrame()
res_mix_hm, res_mix_ann = mix_data_matrix()
res_agg_hm, res_agg_ann = aggregation()

GMTK_FACTOR = 1
MIX_FACTOR = 1.5
AGG_FACTOR = 1.8
ROW_CORRECTOR = 0.8

n_rows = res_mix_hm.shape[0] * ROW_CORRECTOR
n_columns = res_gmtk.shape[1] * GMTK_FACTOR + res_mix_hm.shape[1] * MIX_FACTOR + res_agg_hm.shape[1] * AGG_FACTOR

# Create grid and plot the results in two or three separate plots
f, (ax_gmtk, ax_mix, ax_agg) = \
    plt.subplots(1, 3, figsize=(n_columns, n_rows), sharey=True, \
                 gridspec_kw={"wspace": 0.1,
                              "width_ratios": [res_gmtk.shape[1] * GMTK_FACTOR, res_mix_hm.shape[1] * MIX_FACTOR,
                                               res_agg_hm.shape[1] * AGG_FACTOR]})

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
divider = make_axes_locatable(ax_agg)
ax_agg2 = divider.append_axes("right", size="100%", pad=0.3)
ax_agg_cbar = divider.append_axes("right", size="3%", pad=0.3)

sns.heatmap(res_agg_hm.iloc[:, :8], annot=res_agg_ann.iloc[:, :8], cbar=False, cmap='Blues', ax=ax_agg, fmt='.5g')
g_agg = sns.heatmap(res_agg_hm.iloc[:, 8:], annot=res_agg_ann.iloc[:, 8:], cbar=True, vmin=0, vmax=1,
                    cbar_ax=ax_agg_cbar, cmap='Blues', ax=ax_agg2, fmt='.5g')

cbar_agg = g_agg.collections[0].colorbar
cbar_agg.set_ticks([0, 1])
cbar_agg.set_ticklabels(['low', 'high'])

ax_agg.text(0, -0.5, "Aggregation", fontsize=25, ha='left', va='bottom')

title_args = dict(fontsize=20, position=(1.0, 1.0), ha='right', va='bottom')
ax_agg.set_title('Protein coding', **title_args)
ax_agg2.set_title('lincRNA', **title_args)
ax_agg2.set_yticks([])

f.savefig(snakemake.output.outfile, bbox_inches='tight')
