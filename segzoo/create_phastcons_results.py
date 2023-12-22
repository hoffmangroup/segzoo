import pandas as pd
import pybedtools

GENOME_FILE = snakemake.input.chromsizes
PHASTCONS_FILE = snakemake.input.phastcons
SEGMENTATION_FILE = snakemake.input.segmentation

bt = pybedtools.BedTool(SEGMENTATION_FILE)
res_bt = bt.map(PHASTCONS_FILE, c=4, o='collapse', g=GENOME_FILE)
res_df = res_bt.to_dataframe().rename(columns={
    'blockCount': 'collapsePhastcons'})
res_df = res_df[~(res_df.collapsePhastcons == '.')]

res = {
    'sum': [],
    'counts': [],
    'max': []}
for collapsed_scores in res_df.collapsePhastcons:
    collapsed_scores = list(map(float, collapsed_scores.split(',')))
    res['sum'].append(sum(collapsed_scores))
    res['max'].append(max(collapsed_scores))
    res['counts'].append(len(collapsed_scores))

res_df['sums'] = res['sum']
res_df['counts'] = res['counts']
res_df['maxs'] = res['max']

gb = res_df.groupby('name').apply(lambda x: x.sums.sum() / x.counts.sum())
df = pd.DataFrame(gb, columns=['phastcons'])
df['phastcons_max'] = res_df.groupby('name').apply(lambda x: x.maxs.sum() / x.counts.size)
df.to_csv(snakemake.output.outfile, sep='\t')
