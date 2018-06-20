# Segzoo

## What is segzoo?

Segzoo is a tool that allows to automatically run various genomic analysis on a segmentation obtained with Segway.
The results of each analysis are made available as well as a summarizing visualization of the results.
The requirements for this tool include segtools, bedtools and python packages, but all of them are dependencies that will be taken care of during installation.

## Quick start

This quick start needs you to have anaconda already installed in your local computer (either with python 2 or 3).

1. Download the test [segmentation](https://segway.hoffmanlab.org/2018/protocol/trackhub/hg38/segway.bed.gz) and [GMTK parameters](https://segway.hoffmanlab.org/2018/protocol/params/params.params) and move them both in a directory called, for example, `segzoo`
2. Open a terminal in the mentioned directory and run `conda create -c bioconda -n python3 python=3.6 r-base r-latticeextra r-reshape2 r-cairo r-cluster bedtools -y`
3. After the last command has finished, run `source activate python3` followed by `pip install segzoo`
4. When finished, run `segzoo segway.bed.gz --parameters params.params`
5. After around 30 min, the resulting visualization will be stored in the current's directory `outdir/plots` folder


## How to install

Segzoo is a python 3 tool, so if you have python 2 installed it is highly recommended for you to install segzoo in a separate python 3 environment.
Although it can, this tool is not designed to run on a cluster node without internet access, so all the following steps should be done in a local computer.
To create such an environment run `conda create -n python3 python=3.6` where you can change the name of the environment, `python3`.
Accept all the installation steps.

Next, you need to activate this environment. Run `source activate python3` specifying the name of the environment you chose before.
Now that you already are in it, you can install segzoo. You can do that by running `pip install segzoo`,
which will require you to have bedtools already installed, as it's only in anaconda.
To install bedtools beforehand you can use `conda install -c bioconda bedtools`.

*Note*: currently it's being worked on uploading Segzoo to bioconda.
When this is finished it will be possible to install it just by using `conda install -c bioconda segzoo` which will take care of all the dependencies.

After accepting all installations, segzoo will be good to go!

## How to use

To access the help to know how to run segzoo you can run `segzoo -h` or `segzoo --help`. Here's a look at all possible arguments:

- `--version` to check the current version of segzoo installed
- `--parameters` to specify a params.params file resulting from segway's training to obtain GMTK parameters in the final visualization. If not specified, GMTK parameters won't show in the final visualization
- `--prefix` to specify where you want all needed data (like the genome assembly) to be downloaded (default: the installation environment's directory)
- `-o` or `--outdir` to specify the folder where all the results and the final visualization will be created (default: outdir)
- `-j` to specify the number of cores to use (default: 1)
- `--species` and `--build` specify the species and the build for which the segmentation was created (default: Homo_sapiens and hg38)
- `--download-only` is an option to support cluster use. Running Segzoo using this argument will only run the downloading rules of the pipeline, and store the data in using the specified prefix. After that, runs on nodes without internet access can be done by specifying that same prefix

If you are interested in obtaining information on different gene biotypes than *protein coding and lincRNA*, which are the default,
you can get to the installation folder of segzoo and modify the file `gene_biotypes.py` as you wish.
The same can be said for the final visualization, which can be altered by modifying some variables on top of `visualization.py`

After running the command `segzoo` by specifying the segmentation file and all the optional arguments that you want, the execution of the pipeline will begin.
All necessary data will be downloaded, tools will run the different analysis and the final visualization will be created. This execution may take some time.

## Results

After the execution has finished, the new directory will be created (**outdir** is the default name).
In the **data** folder you will be able to find the results for all the tools' analysis.
In **results** you will find the tables of processed results used in the visualization.
Finally, the visualization will be in the **plots** directory. It will look something like this:

![Plot](https://bitbucket.org/hoffmanlab/segzoo/raw/default/plots/plot.png)

The y-axis are the labels of the segmentation for all the heatmaps, while the x-axis are the different results obtained for each of them.

- In the left there are the learned parameters during the training of Segway.
- Next, a heatmap that has each different column normalized so that the maximum and minimum values are the limits of the color map used.
This applies to all but the GC content, which is normalized between 35% and 65% always. All this information is displayed in the table below
- The aggregation tables are shown in the same order as specified in `gene_biotypes.py`, and can contain duplicates
- The aggregation results displayed for each label are the percentage of counts in one component in comparison to all the idealized gene, so notice that each row adds up to 100
- The number of genes found for each biotype shown is specified after the biotype's name