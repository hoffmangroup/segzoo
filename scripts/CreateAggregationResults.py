import pandas as pd
import numpy as np

#creates a data frame with segtools' output, computes the means of the values for each gene area and label,
#and writes it in a tab-delimited file

#[0:8] is the interval of groups to take if we only want to analyze the splicing
#[9:] would be for the translation

df = pd.read_table(snakemake.input[0], skiprows=1, header=0)
means = df.drop(columns=['offset', 'group']).groupby("component", sort=False).mean()[0:8].T
means.to_csv(path_or_buf=snakemake.output[0], sep='\t')