import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from os import path

sns.set()

cmap = sns.diverging_palette(220, 10, as_cmap=True)


# Prepare nucleotide results
res_nuc = pd.read_table(snakemake.input.nuc, index_col=0)['GC_content']
res_nuc.index = res_nuc.index.map(str)

# TODO change directly from matrix:
# 2 nums after decimal point? Range for heatmap between 35% and 65%

# Prepare length_distribution results
# Preparing the annotation for the matrix, creating a new column called 'frac.segs'
res_len = pd.read_table(snakemake.input.len_dist, index_col=0)
res_len['frac.segs'] = (res_len['num.segs'] / res_len.loc['all']['num.segs'])
res_len = res_len.drop(['num.segs', 'num.bp'], axis=1)
res_len = res_len.drop('all')

# Joining both matrix to create final heatmap and annotation
res_ann = res_len.join(res_nuc)
res = res_ann.copy()

# Interpolation of the parameters to rescale them between 0 and 1
# ['mean.len', 'median.len','stdev.len', 'frac.segs', 'frac.bp', 'GC_content']
for col in res.columns:
    res[col] = (res[col] - res[col].min()) / (res[col].max()-res[col].min())

# Prepare the aggregation results for 'protein_coding' and 'lincRNA' and join them
pc_file = next(x for x in snakemake.input.aggs if path.basename(path.dirname(x)) == 'protein_coding')
lrna_file = next(x for x in snakemake.input.aggs if path.basename(path.dirname(x)) == 'lincRNA')

res_agg_pc = pd.read_table(pc_file, index_col=0)
res_agg_linc = pd.read_table(lrna_file, index_col=0)
res_agg_ann = res_agg_pc.join(res_agg_linc, lsuffix='_protein_coding', rsuffix='_lincRNA')
res_agg = res_agg_ann.copy()

for col in res_agg.columns:
    res_agg[col] = (res_agg[col] - res_agg[col].min()) / (res_agg[col].max()-res_agg[col].min())

# Create grid and plot the results in two separate plots
f, (ax1, ax2, ax3) =\
    plt.subplots(1, 3, figsize=(25, 5), sharey=True, gridspec_kw={"wspace": 0.05, "width_ratios": [0.5, 1, 2]})

# Prepare the gmtk parameters
if snakemake.config['parameters']:
    res_gmtk = pd.read_csv(snakemake.input.gmtk, index_col=0)
    sns.heatmap(res_gmtk, cbar=False, ax=ax1)
else:
    f.delaxes(ax1)

# TODO Change formats to strings and format the values in the matrix directly
sns.heatmap(res, annot=res_ann, cbar=False, vmax=1, vmin=0, cmap="YlGnBu", ax=ax2, fmt='.4g') 
sns.heatmap(res_agg, annot=res_agg_ann, cbar=False, cmap=cmap, ax=ax3, fmt='.5g')

ax1.set_title('GMTK_parameters')
ax2.set_title('Length Distribution + GC_Content')
ax2.set_ylabel('')
ax3.set_title('Aggregation (protein_coding + lincRNA)')

plt.gcf().subplots_adjust(bottom=0.5)
# plt.tight_layout()

f.savefig(snakemake.output.outfile)
