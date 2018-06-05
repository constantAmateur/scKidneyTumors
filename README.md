# scKidneyTumors

This repository contains the various R scripts necessary to perform the analyses for [this paper]() as described in the supplementary materials and main text.  The purpose of these files is to provide further specific documentation of the exact methods and software used.  Detailed documentation can be found in the comments in the files themselves and the include readme, but briefly:
 
processSingleCellData.R - Contains all the functions necessary to pre-process, cluster and summarise the 10X data.

pipelineScript.R - Provides an example of how the functions in processSingleCellData.R are joined together to form a "pipeline" used to process the data.

genotyping.R - Provides the functions necessary to genotype individual cells from bulk DNA data.

similarity.R - Provides the functions used to train the logistic regression model used to infer cluster similarity.

postPipeline.R - Runs the genotyping and similarity code over the output of the pipeline script and stores the result in intermediate files.

makeFigures.R - Once all the above scripts have been run, uses the output to perform the analysis and generate tables and figures found in the paper.
