import pandas as pd

# creates a data frame consisted of genic and intergenic overlap percentage with each segment label
# and writes it in a tab-delimited file

df = pd.read_csv(snakemake.input.infile, sep='\t', header=1, index_col=0)
repeat_percentage = df['repeat']/df['total']

percentage_df = repeat_percentage.to_frame(name='repeat')
percentage_df.to_csv(snakemake.output.outfile, sep='\t')
