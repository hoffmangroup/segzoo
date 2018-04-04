from os.path import join, splitext
import sys

PREFIX = "/users/masenjoi/anaconda2/envs/snakemake-test/"
SPECIES = "Homo_sapiens"
BUILD = "hg38"

recipeGTF = "rnaseq"
recipeSEQ = "sequence"

fileGTF = BUILD + ".gtf"
fileSEQ = BUILD + ".fa"

BUILDPATH = os.path.join(PREFIX, "share", "ggd", SPECIES, BUILD)

RESULTS = "results"
RESULTFILE="results"
DATA = "data"
LOGS = "logs"
PLOTS = "plots"


SEGMENTATION = "/scratch/turnkey/data/2018-02-22/annotation/segway.bed.gz"
preprocessed_segmentation = os.path.join(DATA,"segmentation","segmentation.pkl.gz")

GMTK_PARAMS = "/scratch/turnkey/data/2018-02-22/train/params/params.params"


#OPTIONS
QUIET = ""
CLOBBER = "--clobber"
NOPLOT = ""


#Main rule that sets the targets to obtain
rule all:
	input:
		dynamic(join(DATA, "feature_distance", "gene_biotype", "{biotype}", "feature_distance.tab" )),
		#os.path.join(DATA, "feature_distance", "general", "feature_distance.tab"),
		gmtk=join(RESULTS, "gmtk_parameters", RESULTFILE),
		len_dist=join(RESULTS, "length_distribution", RESULTFILE),
		nuc=join(RESULTS, "nucleotide", RESULTFILE),
		agg=join(RESULTS, "aggregation", "general", RESULTFILE),
		aggs=dynamic(join(RESULTS, "aggregation", "gene_biotype", "{biotype}", RESULTFILE))
	output:
		join(PLOTS, "plot.png")
	script:
		"scripts/Visualization.py"

#DOWNLOADING RULES
#GGD conda installations into the anaconda environment folder
	
rule download_ggd_annotation:
	input:
	output:
		join(BUILDPATH, recipeGTF, fileGTF)
	shell:
		"conda install -p {PREFIX} -y {QUIET} {CLOBBER} -c ggd-alpha {BUILD}-gtf"

rule download_ggd_sequence:
	input:
	output:
		join(BUILDPATH, recipeSEQ, fileSEQ)
	shell:
		"conda install -p {PREFIX} -y {QUIET} {CLOBBER} -c ggd-alpha {BUILD}-sequence"

#PREPROCESSING RULES

#Preprocess segway annotation for future faster pasing
rule run_segtools_preprocess:
	input:
		SEGMENTATION
	output:
		preprocessed_segmentation
	params:
		outfile=splitext(splitext(preprocessed_segmentation)[0])[0]
	shell:
		"segtools-preprocess {QUIET} {CLOBBER} {SEGMENTATION} {params.outfile}"

#Parse the main gtf file to create one file per each gene biotype (uses pybedtools)
#TODO CHECK IF PARAMS CAN BE OMITTED FOR OUTPUTS
#TODO change 0s and 1s to references by name
rule create_gene_biotype_gtf:
	input:
		join(BUILDPATH, recipeGTF, fileGTF)
	output:
		dynamic(join(DATA, "gene_biotype_gtfs", "general", "{biotype}", fileGTF)),
		dynamic(join(DATA, "gene_biotype_gtfs", "gene", "{biotype}", fileGTF))
	log:
		join(LOGS, "create_gene_biotype_gtfs")
	params:
		outdir=join(DATA, "gene_biotype_gtfs"),
		outfile=fileGTF
	script:
		"scripts/CreateGeneBiotypeGtfs.py"


#ANALYZING RULES

#Segtools GMTK-Parameters execution
rule run_segtools_gmtk_parameters:
	input:
		GMTK_PARAMS
	output:
		csv=join(DATA, "gmtk_parameters", "gmtk_parameters.stats.csv"),
		result=join(RESULTS, "gmtk_parameters", RESULTFILE)
	params:
		outdir=join(DATA, "gmtk_parameters")
	shell:
		"segtools-gmtk-parameters {QUIET} {CLOBBER} --outdir {params.outdir} {GMTK_PARAMS};"
		"cp {output.csv} {output.result}" #TODO PLOT is necessary to obtain .csv and no variance is available


#Segtools Aggregation results obtention, general and by gene_biotype
rule run_segtools_aggregation_general:
	input:
		preprocessed_segmentation,
		gtf=join(BUILDPATH, recipeGTF, fileGTF)
	output:
		join(DATA, "aggregation", "general", "feature_aggregation.tab")
	params:
		outdir=join(DATA, "aggregation", "general")
	shell:
		"segtools-aggregation --mode=gene {QUIET} {CLOBBER} {NOPLOT} --outdir {params.outdir} {preprocessed_segmentation} {input.gtf}"

rule run_segtools_aggregation_gene_biotype:
	input:
		preprocessed_segmentation,
		gtf=join(DATA, "gene_biotype_gtfs", "general", "{biotype}", fileGTF)
	output:
		join(DATA, "aggregation", "gene_biotype", "{biotype}", "feature_aggregation.tab")
	params:
		outdir=join(DATA, "aggregation", "gene_biotype", "{biotype}")
	shell:
		"segtools-aggregation --mode=gene {QUIET} {CLOBBER} {NOPLOT} --outdir {params.outdir} {preprocessed_segmentation} {input.gtf}"

rule create_aggregation_results_general:
	input:
		join(DATA, "aggregation", "general", "feature_aggregation.tab")
	output:
		join(RESULTS, "aggregation", "general", "results")
	script:
		"scripts/CreateAggregationResults.py"

rule create_aggregation_results_gene_biotype:
	input:
		join(DATA, "aggregation", "gene_biotype", "{biotype}", "feature_aggregation.tab")
	output:
		join(RESULTS, "aggregation", "gene_biotype", "{biotype}", "results")
	script:
		"scripts/CreateAggregationResults.py"


#Segtools Length Distribution execution and results creation
#TODO Don't just copy, preprocess for easier obtention of results
rule run_segtools_length_distribution:
	input:
		preprocessed_segmentation,
	output:
		join(DATA, "length_distribution", "segment_sizes.tab"),
		join(RESULTS, "length_distribution", "results")
	params:
		outdir=join(DATA, "length_distribution")
	shell:
		"segtools-length-distribution {QUIET} {CLOBBER} {NOPLOT} --outdir {params.outdir} {preprocessed_segmentation};"
		"cp {output[0]} {output[1]}"


#Segtools Feature Distance execution for general file and gene_biotype
rule run_segtools_feature_distance_general:
	input:
		preprocessed_segmentation,
		gtf=join(BUILDPATH, recipeGTF, fileGTF)
	output:
		join(DATA, "feature_distance", "general", "feature_distance.tab"),
		outfile=join(DATA, "feature_distance", "general", "feature_distance_segments.tab")
	params:
		outdir=join(DATA, "feature_distance", "general")
	shell:
		"segtools-feature-distance {QUIET} {CLOBBER} {NOPLOT} --outdir {params.outdir} {preprocessed_segmentation} {input.gtf} > {output.outfile}"

rule run_segtools_feature_distance_gene_biotype:
	input:
		preprocessed_segmentation,
		gtf=join(DATA, "gene_biotype_gtfs", "gene", "{biotype}", fileGTF)
	output:
		join(DATA, "feature_distance", "gene_biotype", "{biotype}", "feature_distance.tab"),
		outfile=join(DATA, "feature_distance", "gene_biotype", "{biotype}", "feature_distance_segments.tab")
	params:
		outdir=join(DATA, "feature_distance", "gene_biotype", "{biotype}")
	shell:
		"segtools-feature-distance {QUIET} {CLOBBER} {NOPLOT} --outdir {params.outdir} {preprocessed_segmentation} {input.gtf} > {output.outfile}"


#Bedtools nucleotide results obtention
rule run_bedtools_nuc:
	input:
		SEGMENTATION,
		fasta=join(BUILDPATH, recipeSEQ, fileSEQ)
	output:
		join(DATA, "nucleotide", "bedtools_output")
	shell:
		"bedtools nuc -fi {input.fasta} -bed {SEGMENTATION} > {output}"

rule create_nucleotide_results:
	input:
		join(DATA, "nucleotide", "bedtools_output")
	output:
		join(RESULTS, "nucleotide", "results")
	script:
		"scripts/CreateNucleotideResults.py"



# Visualization

rule create_diagram:
	input:
	output:
		"dag.svg",
		"dag.png"
	shell:
		"snakemake --rulegraph all | dot -Tsvg > {PLOTS}/dag.svg"
