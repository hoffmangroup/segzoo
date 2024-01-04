import pandas as pd
import pybedtools

REPEAT_FILE = snakemake.input.repeatmasker
SEGMENTATION_FILE = snakemake.input.segmentation
RES_NAMES = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'count_b',
             'base_b', 'len_a', 'frac']

bt = pybedtools.BedTool(SEGMENTATION_FILE)
res_bt = bt.coverage(REPEAT_FILE)
res_df = res_bt.to_dataframe(header=None, names=RES_NAMES)

gb = res_df.groupby('name').frac.mean()
df = gb.to_frame(name='repeatmasker')
df.to_csv(snakemake.output.outfile, sep='\t')
