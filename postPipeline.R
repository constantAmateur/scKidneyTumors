#The extra things that are done after the pipeline has been run on all groups of cells.  Essentially this is just to do the similarity analysis and genotyping independently and store the results somewhere useful.
library(ggplot2)
library(glmnet)
library(reshape2)
library(ComplexHeatmap)


###############################
#*****************************#
#* Cell genotyping from bulk *#
#*****************************#
###############################

import(core,all=TRUE)
import(genome,all=TRUE)
import('genotyping.R',as='gt')

##########
# Params #
##########
mani = read.csv('~/Projects/KidneySC/manifest10X.csv')
mani$BAM = sprintf('~/HCA/remapV2/%s/outs/possorted_genome_bam.bam',mani$Sanger_study_ID)
plotDir='~/scratch/Plots/KidneySC/'
#Get global cluster IDs for all cells
baseDir = '~/scratch/Plots/KidneySC/analysis20171204/'
experiments = grep('barcodeLevelGroup',list.files(baseDir),value=TRUE)
experiments = grep('roximalTubular',experiments,value=TRUE,invert=TRUE)
experiments = setNames(file.path(baseDir,experiments,'batchCorrected'),gsub('barcodeLevelGroup_','',experiments))
#Exclude anything that doesn't have the meta.data object
experiments = experiments[file.exists(file.path(experiments,'metadata.RDS'))]
#For each, get the cluster IDs from the meta-data
clusters=c()
for(i in seq_along(experiments)){
  mDat = readRDS(file.path(experiments[i],'metadata.RDS'))$meta.data
  rownames(mDat) = paste0(mDat$Sanger_study_ID,'___',gsub('^[0-9]+_','',rownames(mDat)))
  clusters=c(clusters,setNames(paste0(names(experiments)[i],'___',mDat$res.1),rownames(mDat)))
}
snpFiles = c(RCC1='~/PipeLive/1711/PD37104a/PD37104a.caveman_c.snps.vcf.gz',
  RCC2='~/PipeLive/1711/PD37228c/PD37228c.caveman_c.snps.vcf.gz',
  PapRCC='~/PipeLive/1711/PD35918h/PD35918h.caveman_c.snps.vcf.gz',
  VHL_RCC='~/PipeLive/1711/PD36793a/PD36793a.caveman_c.snps.vcf.gz',
  Wilms1='~/PipeLive/1711/PD36165d/PD36165d.caveman_c.snps.vcf.gz',
  Wilms2='~/PipeLive/1711/PD37272a/PD37272a.caveman_c.snps.vcf.gz',
  Wilms3='~/PipeLive/1711/PD37276a/PD37276a.caveman_c.snps.vcf.gz',
  Wilms1NR = '~/PipeLive/1711/PD36165e/PD36165e.caveman_c.snps.vcf.gz'
  )
subFiles = gsub('\\.snps\\.','.annot.',snpFiles)
PDids = setNames(gsub('\\..*','',basename(snpFiles)),names(snpFiles))

#############
# Call Subs #
#############

#Get the targets to compare each genotype to
out=list()
for(i in seq_along(subFiles)){
  nom = names(subFiles)[i]
  subFile=subFiles[i]
  bams = setNames(mani$BAM[mani$Experiment==nom],mani$Sanger_study_ID[mani$Experiment==nom])
  if(nom=='Wilms1NR')
    bams = setNames(mani$BAM[mani$Experiment=='Wilms1'],mani$Sanger_study_ID[mani$Experiment=='Wilms1'])
  PDid = gsub('.*/(PD[0-9]+[a-z]).*','\\1',subFile)
  out[[i]] = gt$genoTypeBySubs(subFile,bams,sprintf('~/scratch/KidneySC/%s.subs',PDid),clusterIDs=clusters,nParallel=8,clonalAF=0.2)
  saveRDS(out[[i]],sprintf('~/scratch/KidneySC/%s.subs.RDS',PDid))
}

###########
# Call CN #
###########

#Define the CN segments for each sample
segments = list(RCC1="2,199481900,242985493,2,1,1,0
3,60197,84808933,2,1,1,0
3,84809350,197846280,2,1,2,0
9,209134,141068960,2,1,1,0
13,19455253,65486222,2,1,1,0
14,20433516,107283886,2,1,1,0",
RCC2="2,211334242,242985493,2,1,1,0
3,65982,90362964,2,1,1,0",
PapRCC="3,93505756,197811124,2,1,3,1
12,188285,133812333,2,1,3,1
17,9034,79998834,2,1,3,1",
VHL_RCC="2,55984,242985493,2,1,3,1
3,60596,83348972,2,1,1,0
3,83352258,197606877,2,1,3,1
5,850203,107560911,2,1,3,1
5,107562191,180687907,2,1,5,2
7,115401,159122682,2,1,3,1
9,203937,141068960,2,1,3,1
10,266373,135235890,2,1,3,1
11,196944,134944142,2,1,3,1
12,188285,133839356,2,1,3,1
13,19455957,114999838,2,1,3,1
15,20021973,102431166,2,1,3,1
16,84170,53627229,2,1,3,1
16,53628412,89997381,2,1,2,0
17,27074565,79998834,2,1,3,1
18,125371,78017073,2,1,2,0
19,226776,59097308,2,1,3,1
20,20000786,62954871,2,1,3,1
21,15345102,48101335,2,1,3,1",
Wilms1="11,78383801,134937738,2,1,1,0
12,203339,133839356,2,1,3,1
16,35272667,89973832,2,1,1,0",
Wilms2="2,96057011,242852391,2,1,2,0",
Wilms3="16,86084,89998157,2,1,4,1",
Wilms1NR="16,32156537,33767707,2,1,3,1"
)
gsegs = lapply(segments,function(segments) as.data.frame(t(sapply(lapply(strsplit(segments,'\n')[[1]],strsplit,','),`[[`,1))))
gsegs = lapply(gsegs,function(segments) GRanges(segments[,1],IRanges(as.numeric(segments[,2]),as.numeric(segments[,3])),total=as.integer(segments[,6]),minor=as.integer(segments[,7])))
#Define the q-value cut-off to use for each sample
qCut = 0.1
qCuts = c(RCC1=0.5,RCC2=0.1)
#Do the processing
out=list()
for(i in seq_along(segments)){
  nom = names(segments)[i]
  snpFile = snpFiles[nom]
  bams = setNames(mani$BAM[mani$Experiment==nom],mani$Sanger_study_ID[mani$Experiment==nom])
  PDid = gsub('.*/(PD[0-9]+[a-z]).*','\\1',snpFile)
  #Get the informative SNPs
  qq = qCut
  if(nom %in% names(qCuts))
    qq = qCuts[nom]
  hetSNPs = gt$getCN_SNPs(snpFile,segments[[nom]],qCut=qq)
  out[[i]] = gt$genoTypeByCN(hetSNPs,bams,sprintf('~/scratch/KidneySC/%s.segments',PDid),clusterIDs=clusters,nParallel=8,skipIfExists=FALSE)
  saveRDS(out[[i]],sprintf('~/scratch/KidneySC/%s.segments.RDS',PDid))
  null = gt$plotGenotypedCN(out[[i]]$cCnts,out[[i]]$clCnts,out[[i]]$sgCnts,sprintf('~/scratch/KidneySC/%s.segments',PDid))
}

##########################
#************************#
#* Cell type similarity *#
#************************#
##########################

############
# Preamble #
############

#############
# Functions
import(genome,all=TRUE)
import(core,all=TRUE)
import('~/Projects/KidneySC/Code/similarity.R',as='sim')
##########
# Params
plotDir='~/scratch/Plots/KidneySC/analysis20171204/'
tfs = read.table('~/Projects/Common/Homo_sapiens_transcription_factors_gene_list.txt',sep='\t',header=TRUE)
rtoc = readRDS('~/scratch/KidneySC/globalRawTableOfCounts.RDS')
#Genes to always exclude
hkGeneREGEX='^(EIF[0-9]|RPL[0-9]|RPS[0-9]|RPN1|POLR[0-9]|SNX[0-9]|HSP[AB][0-9]|H1FX|H2AF[VXYZ]|PRKA|NDUF[ABCSV]|ATP[0-9]|PSM[ABCDEFG][0-9]|UBA[0-9]|UBE[0-9]|USP[0-9]|TXN)'
coreExcludeGenes = unique(c(grep('\\.[0-9]+_',rownames(rtoc),value=TRUE), #Poorly characterised
                        grep('MALAT1',rownames(rtoc),value=TRUE), #Contamination
                        grep('^HB[BGMQDAZE][12]?_',rownames(rtoc),value=TRUE), #Contamination
                        grep('^MT-',rownames(rtoc),value=TRUE), #Mitochondria
                        grep(hkGeneREGEX,rownames(rtoc),value=TRUE) #Housekeeping genes
                        ))


##########################
# Train different models #
##########################

###################
# Foetal Clusters
trainDat = sim$loadTrainingData(file.path(plotDir,'barcodeLevelGroup_foetalEpitheliumAndVascularV3/batchCorrected/'))
trainDat$mDat$Trainer = trainDat$mDat$res.1
#Add in the MAST cells
extras = sim$loadTrainingData(file.path(plotDir,'barcodeLevelGroup_tumourImmune/batchCorrected/'))
#Get the cells and their class
cells = c(rownames(trainDat$mDat),rownames(extras$mDat[extras$mDat$res.1=='16',]))
classes = c(trainDat$mDat$res.1,rep('MAST',sum(extras$mDat$res.1=='16')))
#Any genes to exclude go here
excludeGenes=coreExcludeGenes
#Any genes that we wish to give extra weight should go here
includeGenes=c()
#Get and normalise the data
dat = Seurat::LogNormalize(rtoc[,cells])
dat = dat[(rownames(dat) %in% includeGenes) | (Matrix::rowSums(dat>0)>3 & !(rownames(dat)%in%excludeGenes)),]
dat = t(dat)
fitFoetalClusters = sim$multinomialFitCV(dat,classes,nParallel=10)
#Get the genes that are going to be informative
saveRDS(fitFoetalClusters,file.path(plotDir,'lrFoetalClustersV4.RDS'))
###################
# Normal Clusters
#trainDat = sim$loadTrainingData(file.path(plotDir,'barcodeLevelGroup_normalEpitheliumAndVascularV2/batchCorrected/'))
trainDat = sim$loadTrainingData(file.path(plotDir,'barcodeLevelGroup_normalEpitheliumAndVascularWithoutProximalTubularV2/batchCorrected/'))
#Label and merge clusters
#annotMap = c(EN17='G',EN0='PT1',EN7='PT1',EN6='PT2',EN5='PT2',EN4='PT2',EN9='PT2',EN3='PT2',EN18='PT3',EN19='D',EN21='C1',EN22='C2',EN20='P',EN15='U1/2',EN25='U3/4',EN23='M',EN24='F',EN16='GE',EN12='DV',EN10='AV')
#Drop the ambiguous clusters and merge the ones that should really be the same
#trainDat$mDat = trainDat$mDat[paste0('EN',trainDat$mDat$res.1) %in% names(annotMap),]
#trainDat$mDat$Trainer = annotMap[paste0('EN',trainDat$mDat$res.1)]
trainDat$mDat$Trainer = trainDat$mDat$res.1
#Add in the MAST cells
extras = sim$loadTrainingData(file.path(plotDir,'barcodeLevelGroup_tumourImmune/batchCorrected/'))
#Get the cells and their class
cells = c(rownames(trainDat$mDat),rownames(extras$mDat[extras$mDat$res.1=='16',]))
classes = c(trainDat$mDat$res.1,rep('MAST',sum(extras$mDat$res.1=='16')))
#Any genes to exclude go here
excludeGenes=coreExcludeGenes
#Any genes that we wish to give extra weight should go here
includeGenes=c()
#Get and normalise the data
dat = Seurat::LogNormalize(rtoc[,cells])
dat = dat[(rownames(dat) %in% includeGenes) | (Matrix::rowSums(dat>0)>3 & !(rownames(dat)%in%excludeGenes)),]
dat = t(dat)
fitNormalClusters = sim$multinomialFitCV(dat,classes,nParallel=10)
#Get the genes that are going to be informative
saveRDS(fitNormalClusters,file.path(plotDir,'lrNormalClustersWithControl.RDS'))
