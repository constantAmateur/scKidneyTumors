#This is an example script that will run the analysis pipeline for the cells identified in "barcodeLevelGroup_foetalEpitheliumAndVascular.tsv" in the top level of each 10X input directory (srcDirs)

##########
# Params #
##########

maxPC = 30
projectLabel = 'barcodeLevelGroup_foetalEpitheliumAndVascular'
plotDir = file.path('/lustre/scratch112/sanger/my4/Plots/KidneySC/analysis20171204',projectLabel)
srcDirs = c("~/scratch/HCA/Fetus16/cellranger131_count_23208_4834STDY7002875_GRCh38_2A_kidney_CD45pos","~/scratch/HCA/Fetus16/cellranger131_count_23208_4834STDY7002876_GRCh38_2B_kidney_CD45neg","~/scratch/HCA/Fetus17/cellranger131_count_4834STDY7002881_GRCh38_kidney","~/scratch/HCA/Fetus17/cellranger131_count_4834STDY7002885_GRCh38_kidney_CD45pos","~/scratch/HCA/Fetus17/cellranger131_count_4834STDY7002886_GRCh38_kidney_CD45neg")
manifestPaths = '~/Projects/KidneySC/manifest10X.csv'
markerTable = '~/Projects/Common/cellTypeMarkers.tsv'
dontLoad = FALSE

##################
# Pre-Processing #
##################

source('processSingleCellData.R')
channelDat = buildManifest(srcDirs,manifestPaths)
mGenes = unique(read.delim(markerTable,sep='\t',header=TRUE)$Gene)
#Load or create pre-processed count matricies
srat = buildSeurat(channelDat,projectLabel,plotDir,mtCut=0.2,geneCut=200,keepPois=TRUE,dontLoad=dontLoad,keepFile='barcodeLevelGroup_foetalEpitheliumAndVascular.tsv',saveRDS='sratV3.RDS')
#Make the batch regressed version
sratBC = batchRegression(srat,'Experiment',regressionMethod='combat',preserveZeros=TRUE,zeroCut=0,reNormalise=TRUE,dontLoad=dontLoad)
#Process both versions and save output
objects = list(srat)
if(sratBC@misc$batchCorrected)
  objects[[2]]=sratBC
for(object in objects){
  #Add extra meta-data
  object = calcMetadata(object,mitoCutHigh=0.1)
  #Feature selection
  object = selectFeatures(object,maxPC,storePC=maxPC,markerGenes=mGenes)
  #TSNE
  object = reduceDimensions(object)
  #Cluster
  object = clusterCells(object)
  #Find markers
  #object = callMarkerGenes(object)
  #Finally, plot
  object = makeStandardPlots(object,
                           markerTable = read.delim(markerTable,sep='\t',header=TRUE),
                           markerCellAnnotVars = c('res.1','Experiment','Organ','Location1','Location2','TissueDiseaseState','BiologicalRepNo','TechnicalRepNo','Sort'),
                           markerCellAnnotLabs = c('Cluster','Experiment','Organ','LeftRight','Location','TissueDiseaseState','BiologicalRepNo','TechnicalRepNo','FlowSort'),
                           markerGeneAnnotVars = c('Cell_type', 'Organ', 'Anatomical_region', 'Microscopic_region'),
                           markerGeneAnnotLabs = c('CellType','Tissue','Region1','Region2'))
}
