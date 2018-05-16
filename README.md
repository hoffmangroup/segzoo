# Segzoo

## What is segzoo?

Segzoo is a tool that allows to run various genomic analysis on a segmentation obtained by segway. The results of each analysis are made available as well as a summarizing visualization of the results. The requirements for this tool include segtools, bedtools and python packages, but all of them are dependencies that will be treated during installation.

## How to install

Segzoo is a python 3 tool, so if you have python 2 installed it is highly recommended for you to install segzoo in a separate python 3 environment. To create such an environment run `conda create -n python3 python=3.6` where you can change the name of the environment, `python3`. Accept all the installation steps.

Next, you need to activate this environment. Run `source activate python3` specifying the name of the environment you chose before. Now that you already are in it, you can install segzoo. You can do that by running `pip install segzoo`, which will require you to have bedtools already installed, as it's only in anaconda. To install bedtools beforehand you can use `conda install -c bioconda bedtools`. Another option is to install it using `conda install -c bioconda segzoo` which will take care of all the dependencies (__WIP__).

After accepting all installations, segzoo will be good to go!

## How to use

To access the help to know how to run segzoo you can run `segzoo -h` or `segzoo --help`. Here's a look at the most important parameters:
- `--parameters` to specify a params.params file resulting from segway's training to obtain GMTK parameters in the final visualization.
- `--prefix` to specify where you want all needed data (like genome assembly) to be downloaded. The default is in your current environment *share* directory.
- `-o` or `--outdir` to specify the folder where all the results and the final visualization will be created

After running the command `segzoo` by specifying the segmentation file and all the optional arguments that you want, the execution of the pipeline will begin. All necessary data will be downloaded, tools will run the different analysis and the final visualization will be created. This execution may take some time.

## Results

After the execution has finished, the new directory will be created (**outdir** is the default name). In the **data** folder you will be able to find the results for all the tools' analysis. In **results** you will find the tables of processed results used in the visualization. In **logs** some relevant information about the run. Finally, the visualization will be in the **plots** directory. It will look something like this:

![Plot](https://bitbucket.org/hoffmanlab/segzoo/raw/default/plots/plot.png)

Y-axis are the labels of the segmentation for all the plots. As a note: the aggregation results displayed are the percentage of aggregations in one component in comparison to all the gene biotype, so notice that each row adds up to 100.
