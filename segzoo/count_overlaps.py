#!/usr/bin/env python
# coding: utf-8

# In[17]:


# Compute the percentage of segments which overlap with a gene

from pybedtools import BedTool
import numpy as np

OUTFILE = snakemake.output.outfile
FEATURE_IDX = 2

SEGMENTATION_FILENAME = snakemake.input.seg
ANNOTATION_FILENAME = snakemake.input.ann


def filter_gene(interval):
    """
    Returns True if interval is a gene
    """
    return interval[FEATURE_IDX] == 'gene'


segmentation = BedTool(SEGMENTATION_FILENAME)
annotation = BedTool(ANNOTATION_FILENAME)

overlap = segmentation.intersect(annotation.filter(filter_gene), c=True)    # get the number of genes each interval hits

overlap_df = overlap.to_dataframe()    # blockCount is gene count

# 1 if this interval overlaps with a gene, 0 otherwise
overlap_df['overlap'] = np.where(overlap_df['blockCount'] > 0, 1, 0)

counts = overlap_df.groupby(['name'])['overlap'].sum()

percentage_overlap = counts / overlap_df.groupby('name').size() * 100
percentage_overlap.rename('Genic overlap (%)')

df = percentage_overlap.to_frame(name='Genic overlap (%)')

df.index.name = 'label'   # to be consistent with the naming convention
df.to_csv(OUTFILE, sep='\t')




