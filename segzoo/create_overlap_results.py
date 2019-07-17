import pandas as pd

# creates a data frame consisted of genic and intergenic overlap percentage with each segment label
# and writes it in a tab-delimited file

df = pd.read_csv(snakemake.input.infile, sep='\t', header=1, index_col=0)
genic_percentage = df['gene']/df['total']
intergenic_percentage = df['none']/df['total']

percentage_df = pd.concat([genic_percentage, intergenic_percentage], axis=1, keys=['genic', 'intergenic'])
percentage_df.to_csv(snakemake.output.outfile, sep='\t')
