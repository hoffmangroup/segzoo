# Segzoo

## Introduction

Segzoo is a tool designed to automate various genomic analyses on segmentations obtained using Segway. It provides detailed results for each analysis and a comprehensive visualization summarizing the outcomes.

The tool has specific dependencies, including segtools, bedtools, and various Python packages. However, these dependencies will be automatically handled during the installation process.

## Installation and Quick Start

We recommend installing Segzoo using the [mamba](https://mamba.readthedocs.io/en/latest/mamba-installation.html#mamba-install) package manager.

1. To install Segzoo in a separate environment, open a terminal and execute `mamba create -c bioconda -n segzooenv segzoo -y`.
2. Once the installation is complete, activate the Segzoo environment by running `mamba activate segzooenv`.
3. To test Segzoo, download the [segmentation file](https://segway.hoffmanlab.org/2018/protocol/trackhub/hg38/segway.bed.gz) and the [GMTK parameters](https://segway.hoffmanlab.org/2018/protocol/params/params.params), and place them in a directory named, for example, `segzoo`.
4. After the files are in place, execute `segzoo segway.bed.gz --parameters params.params`.
5. After approximately 30 minutes, the resulting visualization will be stored in the `outdir/plots` folder within the current directory.

## Usage

To access the help and learn how to run Segzoo, execute `segzoo -h` or `segzoo --help`. Here are the available command-line arguments:

- `--version`: Check the currently installed version of Segzoo.
- `--parameters`: Specify a `params.params` file generated from Segway's training to obtain GMTK parameters in the final visualization. If not specified, GMTK parameters will not be displayed.
- `--prefix`: Specify the location where all necessary data, such as genome assembly, should be downloaded (default: the installation environment's directory).
- `-o` or `--outdir`: Specify the folder where all the results and the final visualization will be created (default: `outdir`).
- `-j`: Specify the number of cores to utilize (default: 1).
- `--species` and `--build`: Specify the species and build for which the segmentation was created (default: Homo_sapiens and hg38).
- `--download-only`: This option is designed to support cluster use. Running Segzoo with this argument will only execute the downloading rules of the pipeline and store the data using the specified prefix. Subsequently, runs on nodes without internet access can be performed by specifying the same prefix.
- `--mne`: Specify an `mne` file to translate segment labels and track names shown on the figure. Refer to the 'Using mne files' section for details.
- `--normalize-gmtk`: Allow row-wise normalization of GMTK parameters table.
- `--dendrogram`: Perform hierarchical clustering of GMTK parameters row-wise.

If you are interested in obtaining information on gene biotypes other than *protein coding and lincRNA*, which are the default, modify the `gene_biotypes.py` file in the installation folder of Segzoo accordingly. Similarly, the final visualization can be customized by modifying specific variables in `visualization.py`.

Once the `segzoo` command is executed, specifying the segmentation file and any desired optional arguments, the pipeline will commence. It will download all necessary data, run various analyses, and generate the final visualization. Please note that this execution may take some time.

## Results

Upon completion of the execution, a new directory will be created (default name: `outdir`). The following folders will be available:

- **data**: Contains the results for all the tools' analyses.
- **results**: Contains the processed result tables used in the visualization.
- **plots**: Contains the final visualization, which will resemble the example below:

![Plot](https://github.com/mmendez12/segzoo/blob/master/plots/plot.png)

In the visualization:
- The y-axis represents the labels of the segmentation for all the heatmaps.
- The x-axis displays the different results obtained for each of them.
- The left section showcases the learned parameters during the Segway training.
- Subsequently, a heatmap is displayed, with each column normalized to the color map's limits.
- The aggregation tables are presented in the specified order from `gene_biotypes.py`, potentially containing duplicates.
- The aggregation results for each label represent the percentage of counts in one component compared to all the idealized genes. Each row's values sum up to 100.
- The number of genes found for each biotype is provided after the biotype's name.

## Using MNE Files

The `mne` file can be utilized to translate segment labels and track names in the final figure. The file is tab-delimited and should contain three columns in any order:

- `old`: The original label or track name displayed when running `segzoo` with default parameters. Values in this column serve as keys in a Python dictionary or lookup table.
- `new`: Replace the `old` value with the corresponding `new` value from this column.
- `type`: Indicate whether the row should be used to translate a track or a label. This is especially useful when tracks and labels have the same `old` name.

The file header is mandatory and should include the three fields: old, new, and type.

Please note that only the tracks and labels defined in the `mne` file will be updated. Unused tracks and labels will remain unchanged. Here is an example of an `mne` file:

```plaintext
old    new    type
0      Quiescent    label
1      TSS    label
H3K4me3_robust_peaks    H3K4me3    track
```
