import pandas as pd

# creates a data frame with bedtools' output, adds up the parameters for all segments of the same label,
# and writes it in a tab-delimited file
# two rows are added to the table: GC_content and AT_content

usecols = ["4_usercol", "12_num_A", "13_num_C", "14_num_G", "15_num_T", "18_seq_len"]

df = pd.read_table(snakemake.input.infile, header=0, usecols=usecols)
df.columns = ["label", "num_A", "num_C", "num_G", "num_T", "length"]
sums = df.groupby("label").sum()
sums["GC content"] = sums["num_G"].add(sums["num_C"]).div(sums["length"])
sums["AT content"] = sums["num_A"].add(sums["num_T"]).div(sums["length"])
sums.to_csv(path_or_buf=snakemake.output.outfile, sep='\t')
