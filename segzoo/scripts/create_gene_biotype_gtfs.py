from pybedtools import BedTool
from collections import defaultdict
import os
from os.path import join, exists

GTF_GENE_FEATURE_INDEX = 2  # Position of the "feature" parameter in an interval (2 in GTF format)

gtf = BedTool(snakemake.input.gtf)

biotype_dict = defaultdict(list)
biotype_gene_dict = defaultdict(list)
results_log_dict = defaultdict(list)

# Run through all the intervals, adding to the pertinent dictionaries the intervals from each biotype
for interval in gtf:
    # This is necessary because the following line changes the order of attributes,
    # and this can lead to segtools crashing
    fields = interval.fields
    biotype = interval.attrs['gene_biotype']  # gene_type instead of gene_biotype in GH19
    biotype_dict[biotype].append(fields)
    if interval[GTF_GENE_FEATURE_INDEX] == 'gene':
        biotype_gene_dict[biotype].append(fields)

# Create all the files with all the intervals from each biotype, and add their sizes to the results log dictionary
for biotype, array in biotype_dict.items():
    if not exists(join(snakemake.params.outdir, "general", biotype)):
        os.makedirs(join(snakemake.params.outdir, "general", biotype))
    BedTool(array).moveto(join(snakemake.params.outdir, "general", biotype, snakemake.params.outfile))
    results_log_dict[biotype].append(len(array))

# Create the files with only the intervals that are genes from each biotype
# and add their sizes to the results log dictionary
for biotype, array in biotype_gene_dict.items():
    if not exists(join(snakemake.params.outdir, "gene", biotype)):
        os.makedirs(join(snakemake.params.outdir, "gene", biotype))
    BedTool(array).moveto(join(snakemake.params.outdir, "gene", biotype, snakemake.params.outfile))
    results_log_dict[biotype].append(len(array))

# Write down the resulting files' sizes in a new file in the "logs" folder
log_file = open(snakemake.log[0], 'w')
log_file.write("gene biotype\tannotations\tgenes\n")
for biotype, values in results_log_dict.items():
    log_file.write(biotype + "\t" + '{}'.format(values[0]) + "\t" + '{}'.format(values[1]) + "\n")
log_file.close()
