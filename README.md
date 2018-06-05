# scKidneyTumors

This repository contains the various R scripts necessary to perform the analyses for [this paper]() as described in the supplementary materials and main text.  

The primary purpose of this code is to provide the interested and informed reader with additional documentation of the finest details of the methods used in the paper.  It is in principle possible to use these bits of code, together with the input data, to re-produce exactly the results found in the paper.  However, this would take a lot of work to match the file structure expected by the code and the exact meta-data formatting.

Nevertheless, there are several parts of the analysis that are fairly self contained and can probably be used on other data with very little effort and adaptation (the similarity analysis for instance).  While this is encouraged, we emphasise that the primary purpose of this code is to document the methods in the paper and as such should not be considered a supported piece of software.

## Description of files

The starting point for the analysis is single cell RNA-seq data generated with the Chromium 10X and mapped using [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger).  These data are pre-processed, quantified and clustered using [Seurat](http://satijalab.org/seurat/).  Genotyping and cluster similarity analysis is then performed on this processed data and the results saved to intermediate files.  Finally, the code in "makeFigures.R" is run to generate the results found in the paper.

The scripts found herein contribute to this process in the following way:

processSingleCellData.R - Contains all the functions necessary to pre-process, cluster and summarise the 10X data.

pipelineScript.R - Provides an example of how the functions in processSingleCellData.R are joined together to form a "pipeline" used to process the data.

genotyping.R - Provides the functions necessary to genotype individual cells from bulk DNA data.

similarity.R - Provides the functions used to train the logistic regression model used to infer cluster similarity.

postPipeline.R - Runs the genotyping and similarity code over the output of the pipeline script and stores the result in intermediate files.

makeFigures.R - Once all the above scripts have been run, uses the output to perform the analysis and generate tables and figures found in the paper.
