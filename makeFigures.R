#' This will create the figures and tables for the paper.  As well as some other dross not worth publishing.

#############
# Libraries #
#############
library(TCGAbiolinks)
library(ComplexHeatmap)
library(GenomicRanges)
library(ggplot2)
library(cowplot)
library(reshape2)
library(gridExtra)
library(edgeR)
library(biomaRt)
library(monocle)
library(Matrix)
import('~/Projects/KidneySC/Code/similarity.R',as='sim')
import(core,all=TRUE)

#################
# Shared params #
#################
#' Shared objects used by all figures
#Manifest file path
fp_mani = '~/Projects/KidneySC/manifest10X.csv'
#Directory where output of 10X processing is stored
baseDir = '~/scratch/Plots/KidneySC/analysis20171204'
#Directory to store the output in
plotDir='~/scratch/Plots/KidneySC/Paper/Science_R1/'
#Path to global table of counts with all cells
fp_rtoc = '~/scratch/KidneySC/globalRawTableOfCounts.RDS'
#Key genes that define different cell types
fp_markers = '~/Projects/KidneySC/cellTypeMarkers.tsv'
#Pattern used to get BAM files from Sanger ID
fp_bams = '~/HCA/remapV2/%s/outs/possorted_genome_bam.bam'
#Path to list of TFs
fp_tfs = '~/Projects/Common/Homo_sapiens_transcription_factors_gene_list.txt'
#Table to convert from symbolic names to paper labels
fp_num = '~/Projects/KidneySC/FigureAndTableOverview.txt'
#Path where pre-computed genotyping objects are saved
fp_genotyping = '~/scratch/KidneySC/'
#Core compartments.  These are the ones that we really care about and we will use to define a unique mapping.  Within this group each cell must map to only one compartment.
coreComps=c("barcodeLevelGroup_foetalEpitheliumAndVascularV3", 
"barcodeLevelGroup_foetalImmuneV3",
"barcodeLevelGroup_normalEpitheliumAndVascularWithoutProximalTubularV2",
"barcodeLevelGroup_normalImmune", 
"barcodeLevelGroup_proximalTubularV2", 
"barcodeLevelGroup_rubbish",
"barcodeLevelGroup_tumourEpitheliumAndVascular", 
"barcodeLevelGroup_tumourImmune",
"barcodeLevelGroup_notKidney",
"barcodeLevelGroup_unassigned")
#Name mapping.  Give different names associated with each compartment 
compNames = data.frame(dirtyName = c(coreComps,
                                  "barcodeLevelGroup_normalEpitheliumAndVascular"
                                  ),
                       prettyName = c('Foetal Epithelium and Vascular',
                                      'Foetal Immune',
                                      'Normal Epithelium and Vascular without PT',
                                      'Normal Immune',
                                      'Normal Proximal Tubules',
                                      'Indistinct',
                                      'Tumour Epithelium and Vascular',
                                      'Tumour Immune',
                                      'Not Kidney',
                                      'Unassigned',
                                      'Normal Epithelium and Vascular'
                                      ),
                       clusterName = c('F',
                                       'IF',
                                       'N',
                                       'IN',
                                       'PT',
                                       'R',
                                       'T',
                                       'IT',
                                       'NK',
                                       'U',
                                       'EN'
                                       )
                       )
compNames$prettyComputerNames = gsub('\\s','_',compNames$prettyName)
compNames$shortDirtyNames = gsub('.*_','',compNames$dirtyName)
#Load them and prepare them
mani = read.csv(fp_mani)
#Global table of raw counts
rtoc = readRDS(fp_rtoc)
#Sub-figure labelling
subLabels = c('a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z')
#Colour and shape specification
diseaseCols = c(Tumour='#d8b365',Normal='#5ab4ac',NephrogenicRest='#f5f5f5')
develCols = c(Fetal='#af8dc3',Child='#f7f7f7',Adult='#7fbf7b')
#These just get recycled without being mapped to anything in particular
tsneCols = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
tsneCols = c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#b15928")
#tsneCols = c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666')
#Set the min/max for row normalised
rowNormLims = c(-3,3)
diseaseShapes = c(Tumour=24,Normal=22,NephrogenicRest=25)
gg_colour_hue = function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
#Reconstruct the mapping of DropletIDs to compartments
clustTab = list()
for(srcDir in unique(mani$path)){
  srcDir = gsub('/lustre/scratch112/sanger/my4','/lustre/scratch119/realdata/mdt1/team274/my4',srcDir)
  print(srcDir)
  #tgtFile = file.path(srcDir,'barcodesToFilter.tsv')
  #if(file.exists(tgtFile) & file.info(tgtFile)$size>0)
  #  tmp = read.table(tgtFile,sep='\t')
  #Get anything else
  filts=list()
  tgts = grep('barcodeLevelGroup_.*\\.tsv$',list.files(srcDir),value=TRUE)
  for(tgt in tgts){
    tgtFile = file.path(srcDir,tgt)
    if(file.info(tgtFile)$size>0)
      tmp = read.table(tgtFile,sep='\t',header=TRUE)
    filts[[tgt]]=tmp
  }
  clustTab[[srcDir]]=do.call(rbind,filts)
}
clustTab = do.call(rbind,clustTab)
#Filter to include only the core, unique, compartments
clustTab = clustTab[clustTab$Group %in% gsub('.*_','',coreComps),]
clustTab$Compartment = clustTab$ClusterID = clustTab$Cluster =NA
#Add in cluster IDs from the core compartments
for(i in seq_along(coreComps)){
  tgt = file.path(baseDir,coreComps[i],'batchCorrected','metadata.RDS')
  print(tgt)
  if(!file.exists(tgt))
    next
  mDat = readRDS(file.path(baseDir,coreComps[i],'batchCorrected','metadata.RDS'))$meta.data
  rownames(mDat) = paste0(mDat$Sanger_study_ID,'___',gsub('^[0-9]+_','',rownames(mDat)))
  m = match(rownames(mDat),clustTab$DropletID)
  clustTab$Compartment[m] = coreComps[i]
  clustTab$ClusterID[m] = mDat$res.1
  clustTab$Cluster[m] = paste0(gsub('.*_','',coreComps[i]),'___',clustTab$ClusterID[m])
}
#Make sure we have no graphics open
graphics.off()
#Make consistent output
set.seed(1)
#Aliases for immune cell clusters
immuneAliases =c(IF1='fMNP1',
  IF3='fNKT',
  IF5='fMNP2',
  IF6='fMk',
  IN0='NK1',
  IN1='NKT1',
  IN4='Th',
  IN6='MNP1',
  IN7='MNP2',
  IN8='MNP3',
  IN9='NKT2',
  IN10='8T',
  IN11='NP',
  IN13='B',
  IN14='NK2',
  IN15='MST',
  IN16='PDC',
  IT0='tTh1',
  IT2='t8T1',
  IT3='t8T2',
  IT4='tMNP1',
  IT6='tTr',
  IT7='tMNP2',
  IT8='tMNP3',
  IT9='tNKT',
  IT10='tNK1',
  IT11='tMNP4',
  IT12='tNK2',
  IT13='tNK3',
  IT14='tMNP5',
  IT15='tT',
  IT16='tMST1',
  IT17='tPDC',
  IT18='t8T3',
  IT19='tTh2',
  IT20='tP',
  IT21='tMST2',
  IT22='tE',
  IT23='tB',
  IT24='tNP')
#Make a gene frequency object from all maps in the data
geneFrac = list()
for(i in seq_len(nrow(compNames))){
  #Skip the uninformative ones 
  if(compNames$prettyName[i] %in% c('Indistinct','Not Kidney','Unassigned'))
    next
  print(compNames$prettyName[i])
  mDat = readRDS(file.path(baseDir,compNames$dirtyName[i],'batchCorrected','metadata.RDS'))
  rownames(mDat$meta.data) = paste0(mDat$meta.data$Sanger_study_ID,'___',gsub('^[0-9]+_','',rownames(mDat$meta.data)))
  geneFrac[[compNames$clusterName[i]]] = do.call(cbind,lapply(split(rownames(mDat$meta.data),paste0(compNames$clusterName[i],mDat$meta.data$res.1)),function(e) Matrix::rowSums(rtoc[,e,drop=FALSE]>0)/length(e)))
}
geneFrac = do.call(cbind,geneFrac)



#########################
# Params for similarity #
#########################
#The grouping (and ordering) of clusters
#Cluster annotation and grouping
tumourGroups = list(pRCC=15,
                  RCC=c(17,4,9,6,7,10,12,0,1),
                  WT=c(18,16,5,13),
                  MyoFibroblasts=c(14,8),
                  Vascular=c(3,11,2))
normGroups = list(Podocytes=11,
                  ProximalTubules=c(12,22,13),
                  LoopOfHenle=c(0,8),
                  DistalTubules=c(14),
                  CollectingDuct=c(16,17),
                  PelvicEpithelium=c(20),
                  UretericEpithelium=c(19,23,15,10),
                  Fibroblasts=c(21),
                  MyoFibroblasts=c(18),
                  Endothelium=c(9,4,2,7),
                  Other=c(1,3,5,6)
                  )
foetGroups = list(Endothelium=3,
                  MyoFibroblasts=9,
                  Fibroblasts=6,
                  UretericBud=c(11,4),
                  CapMesenchyme=c(7,10,1),
                  PrimitiveVesicle=c(8),
                  Ganglia=13,
                  Other=c(12,2,5)
                  )

############################
# Global Genotyping params #
############################
alpha=0.025
mani$BAM = sprintf(fp_bams,mani$Sanger_study_ID)
#Files defining the locations of all SNPs in each sample
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


#####################
# Params for tables #
#####################
#Path to the various annotation tables to merge
nonImmune = '~/Projects/KidneySC/Supplementary_Table1_17.csv'
foetalImmune = '~/Projects/KidneySC/FoetalImmuneAnnotation.csv'
normalImmune = '~/Projects/KidneySC/NormalImmuneAnnotation.csv'
tumourImmune = '~/Projects/KidneySC/TumourImmuneAnnotation.csv'

###############
# clusterInfo #
###############
#Create annotation table.  Used subsequently, that's why it's first.
#Load them and merge them
all = read.table(nonImmune,sep=',',header=TRUE)
#Add extra columns in
all$Genotype='Normal'
all$Genotype[all$Cluster_ID %in% paste0('T',c(0,1,4,6,7,9,10,12,15,17,18))]='Tumour'
all$Genotype[all$Cluster_ID=='T5']='Tumour+NR'
all$Genotype[all$Cluster_ID=='T13']='Mixed'
all$Genotype[all$Cluster_ID=='T16']='NR'
#Foetal immune
tmp = read.table(foetalImmune,sep=',',header=TRUE)
tmp = data.frame(Cluster_ID = paste0('IF',tmp$Cluster),
                 Alias = '-',
                 Data_set = tmp$Map,
                 Referring_to_Figure = '-',
                 Referring_to_Supplementary_Figure = '-',
                 Number_of_cells = 'TBC',
                 Category = 'Foetal_kidney_immune',
                 Cell_type1 = tmp$Cell_type,
                 Cell_type2 = '-',
                 Cell_type3 = '-',
                 Positive_marker_mRNA = gsub(',',';',gsub('Junk','-',tmp$Key_markers......)),
                 Genotype='-')
all = rbind(all,tmp)
#Normal immune
tmp = read.table(normalImmune,sep=',',header=TRUE)
tmp = data.frame(Cluster_ID = paste0('IN',tmp$Cluster),
                 Alias = '-',
                 Data_set = tmp$Map,
                 Referring_to_Figure = '-',
                 Referring_to_Supplementary_Figure = '-',
                 Number_of_cells = 'TBC',
                 Category = 'Normal_mature_kidney_immune',
                 Cell_type1 = tmp$Cell_type,
                 Cell_type2 = '-',
                 Cell_type3 = '-',
                 Positive_marker_mRNA = gsub(',',';',gsub('Junk','-',tmp$Key_markers......)),
                 Genotype='Normal')
all = rbind(all,tmp)
#Tumour immune
tmp = read.table(tumourImmune,sep=',',header=TRUE)
tmp = data.frame(Cluster_ID = paste0('IT',tmp$Cluster),
                 Alias = '-',
                 Data_set = tmp$Map,
                 Referring_to_Figure = '-',
                 Referring_to_Supplementary_Figure = '-',
                 Number_of_cells = 'TBC',
                 Category = 'Kidney_tumour_immune',
                 Cell_type1 = tmp$Cell_type,
                 Cell_type2 = '-',
                 Cell_type3 = '-',
                 Positive_marker_mRNA = gsub(',',';',gsub('Junk','-',tmp$Key_markers......)),
                 Genotype='Normal')
all = rbind(all,tmp)
#Fix up some errors in the source file
all$Data_set = gsub('foetalEpitheliumAndVascular.*','foetalEpitheliumAndVascularV3',all$Data_set)
all$Data_set = gsub('foetalImmune.*','foetalImmuneV3',all$Data_set)
all$Cell_type1 = gsub('\n','',all$Cell_type1)
all$Positive_marker_mRNA = gsub('\n','',all$Positive_marker_mRNA)
#Add in the aliases for the immune clusters
m = match(names(immuneAliases),all$Cluster_ID)
all$Alias[m] = immuneAliases
#Now for each, load the meta-data and get counts
for(nom in unique(all$Data_set)){
  print(nom)
  tgt = file.path(baseDir,nom,'batchCorrected','metadata.RDS')
  mDat = readRDS(tgt)$meta.data
  cnts = table(mDat$res.1)
  clusts = as.character(seq(min(as.numeric(names(cnts))),max(as.numeric(names(cnts)))))
  cnts = cnts[clusts]
  w = all$Data_set==nom
  all$Number_of_cells[w] = cnts[gsub('[A-Z]*','',all$Cluster_ID[w])]
  #Get counts by experiment
  cnts = table(mDat$res.1,mDat$Experiment)[clusts,]
  #Process one column at a time
  for(cnom in colnames(cnts)){
    tmp = paste0('Number_',cnom)
    if(!(tmp %in% colnames(all)))
      all[,tmp]=0
    all[w,tmp] = cnts[gsub('[A-Z]*','',all$Cluster_ID[w]),cnom]
  }
}
all$Referring_to_Supplementary_Figure='1'
#Fix mis-labelling
all$Cell_type2 =  gsub('Capillary','Cap',all$Cell_type2)
write.table(all,file.path(plotDir,'clusterInfo.tsv'),sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
#Extract markers of different compartments from this
markers = lapply(strsplit(all$Positive_marker_mRNA,';'),function(e) gsub('[+-]|\\s','',e))
markers = split(markers,grepl('Immune',all$Data_set))
markers = lapply(markers,function(e) unique(unlist(e,use.names=FALSE)))
markers = lapply(markers,function(e) sort(e[!(e %in% c('','Private'))]))
markers = data.frame(Gene = unlist(markers,use.names=FALSE),
                     MarkerOf = rep(c('EpitheliumOrVascular','Immune'),lengths(markers))
                     )
#Add extra information about some
markersClusterUsage = list()
markersClusterUsage$G = c('WT1','PODXL','PTPRO')
markersClusterUsage$PT1 = c('SLC22A8','SLC17A3','SLC22A7','SLC16A9','SLC7A13','SLC34A1','SLC13A3','VCAM1')
markersClusterUsage$PT2 = c('SLC22A8','SLC17A3','SLC22A7','SLC16A9','SLC7A13','SLC34A1','SLC13A3','VCAM1')
markersClusterUsage$PT3 = c('SLC22A8','SLC17A3','SLC22A7','SLC16A9','SLC7A13','SLC34A1','SLC13A3','VCAM1')
markersClusterUsage$H = c('SLC12A1','CLDN16')
markersClusterUsage$D = c('KCNJ1','SLC8A1','AVPR2','CLDN8','AQP2')
markersClusterUsage$C1 = c('CLCNKB','ATP6V0D2','SLC26A4','SLC4A1')
markersClusterUsage$C2 = c('CLCNKB','ATP6V0D2','SLC26A4','SLC4A1')
markersClusterUsage$P = c('SAA2','KRT23')
markersClusterUsage$U1 = c('S100P','DHRS2','UPK1A','UPK1B','KRT5','PVRL4','TP63')
markersClusterUsage$U2 = c('S100P','DHRS2','UPK1A','UPK1B','KRT5','PVRL4','TP63')
markersClusterUsage$U3 = c('S100P','DHRS2','UPK1A','UPK1B','KRT5','PVRL4','TP63')
markersClusterUsage$U4 = c('S100P','DHRS2','UPK1A','UPK1B','KRT5','PVRL4','TP63')
markersClusterUsage$F = c('MMP2','EMILIN1','SFRP2')
markersClusterUsage$M = c('ACTA2','PDGFRB')
markersClusterUsage$GE = c('PECAM1','PTPRB','KDR','VCAM1','AQP1','SEMA3G','CLDN5','SLC14A1','PLVAP')
markersClusterUsage$DV = c('PECAM1','PTPRB','KDR','VCAM1','AQP1','SEMA3G','CLDN5','SLC14A1','PLVAP')
markersClusterUsage$AV1 = c('PECAM1','PTPRB','KDR','VCAM1','AQP1','SEMA3G','CLDN5','SLC14A1','PLVAP')
markersClusterUsage$AV2 = c('PECAM1','PTPRB','KDR','VCAM1','AQP1','SEMA3G','CLDN5','SLC14A1','PLVAP')
markersClusterUsage$PV = c('WT1','PODXL','PTPRO')
markersClusterUsage$UB = c('CDH16','POU3F3','RET','GATA3','TFCP2L1','ELF3','HNF1B')
markersClusterUsage$CM = c('SIX1','SIX2','PAX2','CITED1')
markersClusterUsage$fF = c('MMP2','EMILIN1','SFRP2')
markersClusterUsage$fM = c('ACTA2','PDGFRB')
markersClusterUsage$fE = c('PECAM1','PTPRB','KDR','VCAM1','AQP1','SEMA3G','CLDN5','SLC14A1','PLVAP')
markersClusterUsage$fGa = c('SCG2','CHGB','STMN2')
markersClusterUsage = melt(markersClusterUsage)
colnames(markersClusterUsage) = c('Gene','Cluster')
write.table(markers,fp_markers,sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)



#######################
# normMap and foetMap #
#######################
#The fancy looking tSNE plots
tgts = c('F','N')
m = match(tgts,compNames$clusterName)
srcDirs = compNames$dirtyName[m]
srcDirs = file.path(baseDir,srcDirs)
cPrefix = setNames(tgts,srcDirs)
figNoms = setNames(c('foetMap.pdf','normMap.pdf'),srcDirs)
for(srcDir in srcDirs){
  message(srcDir)
  ############################################
  # Load detail for this cell and prep data
  toc = readRDS(file.path(srcDir,'batchCorrected','tableOfCounts.RDS'))
  mDat = readRDS(file.path(srcDir,'batchCorrected','metadata.RDS'))
  #Add a marker of developmental state
  mDat$meta.data$developmentalState = cut(x=mDat$meta.data$Age_donor.MonthsPostConception.,breaks=c(0,9,18*12,Inf),labels=c('Fetal','Child','Adult'))
  #Fix up CellIDs
  rownames(mDat$meta.data) = paste0(mDat$meta.data$Sanger_study_ID,'___',gsub('^[0-9]+_','',rownames(mDat$meta.data)))
  colnames(toc) = rownames(mDat$meta.data)
  #############
  # tSNE plot
  message("Making tSNE")
  df = mDat$meta.data
  tmp = split(df,df$res.1)
  #Work out where to draw the labels
  cLabs = lapply(tmp,function(e) c(x=median(e$tSNE1),y=median(e$tSNE2)))
  cLabs = cbind(res.1=names(cLabs),as.data.frame(do.call(rbind,cLabs)))
  #Set cluster ordering and ensure consistent levels
  lvls = as.character(seq(min(as.numeric(df$res.1)),max(as.numeric(df$res.1))))
  df$res.1 = factor(df$res.1,levels=lvls)
  cLabs$res.1 = factor(cLabs$res.1,levels=lvls)
  #Set the bandwidth to a fixed fraction of width and height
  bwidth = .05*c(sum(abs(range(df$tSNE1))),sum(abs(range(df$tSNE2))))
  #Manual colours
  manCols = ((seq_along(lvls)-1)%%length(tsneCols))+1
  #Adjust colours manually. Based on a six colour pallette.  Also define mask clusters
  toMask=c()
  if(grepl('tumourEpitheliumAndVascular',srcDir)){
    manCols[lvls %in% as.character(c(4,6))] = 3
    manCols[lvls %in% as.character(c(0,8))] = 6
    manCols[lvls %in% as.character(c(14,11,15))] = 5
    manCols = tsneCols[manCols]
  }else if(grepl('foetalEpitheliumAndVascular',srcDir)){
    manCols[lvls %in% as.character(c(10))] = 4
    manCols[lvls %in% as.character(c(9))] = 3
    manCols[lvls %in% as.character(c(6))] = 4
    #manCols[lvls %in% as.character(c(7))] = 4
    #manCols[lvls %in% as.character(c(12))] = 5
    toMask=c(0,12,2,5,4,1,10)
    manCols = tsneCols[manCols]
    manCols[11+1]='#009FE3'
    manCols[7+1]='#E6007E'
    manCols[8+1]='#CA9E67'
  }else if(grepl('normalEpitheliumAndVascular',srcDir)){
    #manCols = rep('#808080',length(manCols))
    #names(manCols) = seq(0,length(manCols)-1)
    #proxToDist = c(11,13,22,12,0,8,14,16,17,20,10,15,23,19)
    #manCols[as.character(proxToDist)]= gg_colour_hue(length(proxToDist))
    manCols = tsneCols[manCols]
    toMask = c(1,3,5,6,8)
  }
  #Make the masked ones grey
  manCols[toMask+1] = '#808080'
  #This looks weird (and is weird), but we need to trick ggplot into letting us specify manual alpha for both geoms somehow
  df$isMasked = as.numeric(df$res.1 %in% toMask)
  df$isMasked2 = factor(as.numeric(df$isMasked+2))
  df$isMasked = factor(df$isMasked)
  #Make cluster names with prefix
  cLabs$ClusterID = paste0(cPrefix[srcDir],cLabs$res.1)
  gg = ggplot(df,aes(tSNE1,tSNE2,colour=res.1)) +
    #geom_tile(data=rastMap,aes(x,y,fill=res.1),alpha=0.5,inherit.aes=FALSE)+
    #geom_contour(data=minContDat,aes(x,y,z=z,colour=res.1),alpha=0.3,bins=1,inherit.aes=FALSE) +
    stat_density2d(aes(x=tSNE1,y=tSNE2,colour=res.1,alpha=isMasked),h=bwidth) +
    geom_point(size=1.0,aes(shape=TissueDiseaseState,alpha=isMasked2),stroke=0.1) +
    scale_shape_manual(values=diseaseShapes,breaks=names(diseaseShapes))+
    scale_colour_manual(values=manCols) +
    scale_alpha_manual(values=c(`0`=0.3,`1`=0.1,`2`=1.0,`3`=0.25)) + 
    #Randomise default colours
    #scale_colour_manual(values=sample(ggColours(length(unique(df$res.1))))) +
    geom_text(data=cLabs,aes(x,y,label=ClusterID),inherit.aes=FALSE,size=10) +
    ggtitle(gsub('.*_','',srcDir)) +
    guides(colour=FALSE,fill=FALSE,shape=FALSE,alpha=FALSE)
  pdf(file.path(plotDir,figNoms[srcDir]),width=14,height=14)
  plot(gg)
  dev.off()
}

########
# vegf #
########
#Make barplots for vegf figure
gns = c('CD68','ACKR1','FLT4','FLT1','KDR','VEGFC','VEGFA')
clusts = c('cR2','tE1','tE2')
#Get the simple ones
df = geneFrac[match(gns,gsub('_.*','',rownames(geneFrac))),all$Cluster_ID[match(clusts,all$Alias)]]
colnames(df) = clusts
#Add the tumour macrophages
tmp = all$Cluster_ID[grep('tMNP',all$Alias)]
mDat = readRDS(file.path(baseDir,compNames$dirtyName[compNames$clusterName=='IT'],'batchCorrected','metadata.RDS'))
#Fix up CellIDs
rownames(mDat$meta.data) = paste0(mDat$meta.data$Sanger_study_ID,'___',gsub('^[0-9]+_','',rownames(mDat$meta.data)))
tmp = rtoc[match(gns,gsub('_.*','',rownames(rtoc))),rownames(mDat$meta.data)[mDat$meta.data$res.1 %in% gsub('IT','',tmp)]]
df = cbind(df,'tMNP'=Matrix::rowSums(tmp>0)/ncol(tmp))
#Melt it and plot
df = melt(df)
df[,1] =as.character(df[,1])
df[,2] =as.character(df[,2])
colnames(df) = c('Gene','Cluster','Frac')
df$Symbol = gsub('_.*','',df$Gene)
#Set order
df$Symbol = factor(df$Symbol,level=rev(gns))
#Make the barplots
gg = ggplot(df,aes(Symbol,Frac)) +
  geom_bar(aes(fill=Symbol),stat='identity') +
  ylim(0,1) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  facet_grid(~Cluster)
pdf(file.path(plotDir,'vegf.pdf'))
plot(gg)
dev.off()

 

####################################
# Special heatmaps for main figure #
####################################
tgts = c('N','F')
m = match(tgts,compNames$clusterName)
srcDirs = compNames$dirtyName[m]
srcDirs = file.path(baseDir,srcDirs)
mDats = lapply(srcDirs,function(e) {
               #Define the different clusters and markers in order
               mDat = readRDS(file.path(e,'batchCorrected','metadata.RDS'))
               #Fix up CellIDs
               rownames(mDat$meta.data) = paste0(mDat$meta.data$Sanger_study_ID,'___',gsub('^[0-9]+_','',rownames(mDat$meta.data)))
               mDat$meta.data})
markerList = list(normNeph = list(
                                  genes = c("TP63", "PVRL4", "KRT5", "UPK1B", "UPK1A", "DHRS2", "S100P", "KRT23", "SAA2", "SLC4A1", "SLC26A4", "ATP6V0D2", "CLCNKB", "AQP2", "CLDN8", "AVPR2", "SLC8A1", "KCNJ1", "CLDN16", "SLC12A1", "SLC13A3", "SLC34A1", "SLC7A13", "SLC16A9", "SLC22A7", "SLC17A3", "SLC22A8", "PTPRO", "PODXL", "WT1"),
                                  clusters = c('G','PT1','PT2','PT3','H','D','C1','C2','P','U1','U2','U3','U4'),
                                  mDat = mDats[[1]],
                                  prefix = 'N'
                                  ),
                  normOther = list(
                                   genes = c("PLVAP", "SLC14A1", "CLDN5", "SEMA3G", "AQP1", "VCAM1", "KDR", "PTPRB", "PECAM1", "PDGFRB", "ACTA2", "SFRP2", "EMILIN1", "MMP2"),
                                   clusters = c('F','M','GE','DV','AV1','AV2'),
                                   mDat = mDats[[1]],
                                   prefix = 'N'
                                   ),
                  devNeph = list(
                                 genes = c("CITED1", "PAX2", "SIX2", "SIX1", "HNF1B", "ELF3", "TFCP2L1", "GATA3", "RET", "POU3F3", "CDH16", "PTPRO", "PODXL", "WT1"),
                                 clusters = c('UB','CM','PV'),
                                 mDat = mDats[[2]],
                                 prefix = 'F'
                                 ),
                  devOther = list(
                                  genes = c("STMN2","CHGB","SCG2", "PLVAP", "SLC14A1", "CLDN5", "SEMA3G", "AQP1", "VCAM1", "KDR", "PTPRB", "PECAM1", "PDGFRB", "ACTA2", "SFRP2", "EMILIN1", "MMP2"),
                                  clusters = c('fE','fM','fF','fGa'),
                                  mDat = mDats[[2]],
                                  prefix='F'
                                  )
                  )
for(nom in names(markerList)){
  gns = markerList[[nom]]$genes
  clusts = markerList[[nom]]$clusters
  mDat = markerList[[nom]]$mDat
  prefix = markerList[[nom]]$prefix
  #Get the fixed up gene names
  gns = rownames(rtoc)[match(gns,gsub('_.*','',rownames(rtoc)))]
  #Get the data for these
  tmp = rtoc[gns,rownames(mDat)]
  tmp = do.call(cbind,lapply(split(rownames(mDat),mDat$res.1),function(e)  Matrix::rowSums(tmp[,e,drop=FALSE]>0)/length(e)))
  #Fix up names
  colnames(tmp) = all$Alias[match(paste0(prefix,colnames(tmp)),all$Cluster_ID)]
  tmp = tmp[,clusts]
  rownames(tmp) = gsub('_.*','',rownames(tmp))
  #Define the ordering of the columns
  #Row-norm and use the default heatmap colour
  tmp = t(scale(t(tmp)))
  tmp[tmp<rowNormLims[1]]=rowNormLims[1]
  tmp[tmp>rowNormLims[2]]=rowNormLims[2]
  #Drop pure NA rows
  w = apply(tmp,1,function(e) any(!is.na(e)))
  tmp = tmp[w,]
  gg_hm = Heatmap(tmp,
          #col = circlize::colorRamp2(c(-3,0,3),c('black','green','purple')),
          col = circlize::colorRamp2(c(rowNormLims[1],0,rowNormLims[2]), c("blue", "#EEEEEE", "red")),
          name='z-score of cells in cluster\nexpressing gene',
          row_title='Genes',
          row_names_gp = gpar(fontsize=8),
          column_title='Clusters',
          cluster_columns=FALSE,
          cluster_rows=FALSE,
          column_order = clusts,
          row_order = gsub('_.*','',rev(gns))
          #column_order = colnames(tmp)[order(as.numeric(gsub('[A-Za-z]+','',colnames(tmp))))]
          )
  pdf(file.path(plotDir,paste0(nom,'Marks.pdf')))
  draw(gg_hm)
  dev.off()
}


#############
# fullMap*  #
#############
#tSNE maps and marker heatmaps
#Which compartments to make?
tgts = c('T','F','EN','N','IT','IF','IN')
#tgts = c('F','IF','N','IN','T','IT','PT','EN','R')
m = match(tgts,compNames$clusterName)
srcDirs = file.path(baseDir,compNames$dirtyName[m])
figNoms = setNames(compNames$prettyName[m],srcDirs)
cPrefix = setNames(tgts,srcDirs)
figPath = setNames(paste0('fullMap',tgts,'.pdf'),srcDirs)
for(srcDir in srcDirs){
  pdf(file.path(plotDir,figPath[srcDir]),width=32,height=18)
  message(srcDir)
  ############################################
  # Load detail for this cell and prep data
  toc = readRDS(file.path(srcDir,'batchCorrected','tableOfCounts.RDS'))
  mDat = readRDS(file.path(srcDir,'batchCorrected','metadata.RDS'))
  #Add a marker of developmental state
  mDat$meta.data$developmentalState = cut(x=mDat$meta.data$Age_donor.MonthsPostConception.,breaks=c(0,9,18*12,Inf),labels=c('Fetal','Child','Adult'))
  #Fix up CellIDs
  rownames(mDat$meta.data) = paste0(mDat$meta.data$Sanger_study_ID,'___',gsub('^[0-9]+_','',rownames(mDat$meta.data)))
  colnames(toc) = rownames(mDat$meta.data)
  #############
  # tSNE plot
  message("Making tSNE")
  df = mDat$meta.data
  tmp = split(df,df$res.1)
  #Work out where to draw the labels
  cLabs = lapply(tmp,function(e) c(x=median(e$tSNE1),y=median(e$tSNE2)))
  cLabs = cbind(res.1=names(cLabs),as.data.frame(do.call(rbind,cLabs)))
  #Set cluster ordering and ensure consistent levels
  lvls = as.character(seq(min(as.numeric(df$res.1)),max(as.numeric(df$res.1))))
  df$res.1 = factor(df$res.1,levels=lvls)
  cLabs$res.1 = factor(cLabs$res.1,levels=lvls)
  #Set the bandwidth to a fixed fraction of width and height
  bwidth = .05*c(sum(abs(range(df$tSNE1))),sum(abs(range(df$tSNE2))))
  #Manual colours
  manCols = ((seq_along(lvls)-1)%%length(tsneCols))+1
  #Adjust colours manually. Based on a six colour pallette.  Also define mask clusters
  if(grepl('tumourEpitheliumAndVascular',srcDir)){
    manCols[lvls %in% as.character(c(4,6))] = 3
    manCols[lvls %in% as.character(c(0,8))] = 6
    manCols[lvls %in% as.character(c(14,11,15))] = 5
    manCols = tsneCols[manCols]
  }else if(grepl('foetalEpitheliumAndVascular',srcDir)){
    manCols[lvls %in% as.character(c(10))] = 4
    manCols[lvls %in% as.character(c(9))] = 3
    manCols[lvls %in% as.character(c(6))] = 4
    #manCols[lvls %in% as.character(c(7))] = 4
    #manCols[lvls %in% as.character(c(12))] = 5
    manCols = tsneCols[manCols]
    manCols[11+1]='#009FE3'
    manCols[7+1]='#E6007E'
    manCols[8+1]='#CA9E67'
  }else if(grepl('normalEpitheliumAndVascular',srcDir)){
    #manCols = rep('#808080',length(manCols))
    #names(manCols) = seq(0,length(manCols)-1)
    #proxToDist = c(11,13,22,12,0,8,14,16,17,20,10,15,23,19)
    #manCols[as.character(proxToDist)]= gg_colour_hue(length(proxToDist))
    manCols = tsneCols[manCols]
  }else{
    manCols = tsneCols[manCols]
  }
  #Make cluster names with prefix
  cLabs$ClusterID = paste0(cPrefix[srcDir],cLabs$res.1)
  m = match(cLabs$ClusterID,all$Cluster_ID)
  w = !is.na(m) & all$Alias[m]!='-'
  cLabs$ClusterID[w] = sprintf('%s (%s)',cLabs$ClusterID[w],all$Alias[m[w]])
  gg_tsne = ggplot(df,aes(tSNE1,tSNE2,colour=res.1)) +
    stat_density2d(aes(x=tSNE1,y=tSNE2,colour=res.1),alpha=0.3,h=bwidth) +
    geom_point(size=1.0,aes(shape=TissueDiseaseState),stroke=0.1) +
    scale_shape_manual(values=diseaseShapes,breaks=names(diseaseShapes))+
    scale_colour_manual(values=manCols) +
    geom_text(data=cLabs,aes(x,y,label=ClusterID),inherit.aes=FALSE,size=10) +
    ggtitle(figNoms[srcDir]) +
    guides(colour=FALSE,fill=FALSE,shape=FALSE)
  #pdf(file.path(plotDir,sprintf('tSNE_%s.pdf',gsub('.*_','',srcDir))),width=14,height=14)
  #plot(gg)
  #dev.off()
  #Make a heatmap of the key markers
  gMarks = markers[markers$Gene %in% gsub('_.*','',rownames(rtoc)),]
  #Exclude the immune from epithelium and visa-versa
  if(grepl('I',cPrefix[srcDir])){
    gMarks = gMarks[gMarks$MarkerOf=='Immune',]
  }else if(cPrefix[srcDir]!='R'){
    #Include both for rubbish
    gMarks = gMarks[gMarks$MarkerOf!='Immune',]
  }
  #Get the extra algorithmic markers
  extraMarkers = lapply(mDat$clusterSummaries,function(e) head(rownames(e$markerGeneDat)[!e$markerGeneDat$isKeyGene],n=5))
  #Add specific markers to the specific cluster.  Keep a record of these
  specMarks = c()
  #for(i in seq_along(extraMarkers)){
  #  nom = names(extraMarkers)[i]
  #  if(grepl('\\(.+\\)',cLabs[nom,'ClusterID'])){
  #    tClust = gsub('.*\\((.+)\\)$','\\1',cLabs[nom,'ClusterID'])
  #    tmp = markersClusterUsage$Gene[markersClusterUsage$Cluster==tClust]
  #    extraMarkers[[i]] = c(extraMarkers[[i]],tmp)
  #    specMarks = c(specMarks,tmp)
  #  }
  #}
  #Any that don't have a special assigment, add as an extra category
  extraMarkers$literature = gMarks[!(gMarks[,1] %in% specMarks),1]
  #Unlist and make split
  ss = rep(names(extraMarkers),lengths(extraMarkers))
  extraMarkers = unlist(extraMarkers,use.names=FALSE)
  #Add the ENSEMBL ID if we're missing it
  w = grepl('_',extraMarkers)
  extraMarkers[!w] = rownames(rtoc)[match(extraMarkers[!w],gsub('_.*','',rownames(rtoc)))]
  #Get the data for these
  tmp = rtoc[extraMarkers,rownames(mDat$meta.data)]
  tmp = do.call(cbind,lapply(split(rownames(mDat$meta.data),mDat$meta.data$res.1),function(e)  Matrix::rowSums(tmp[,e,drop=FALSE]>0)/length(e)))
  #Fix up names
  colnames(tmp) = cLabs[colnames(tmp),'ClusterID']
  rownames(tmp) = gsub('_.*','',rownames(tmp))
  #Add an indicater if it's a lit markers
  rownames(tmp)[rownames(tmp) %in% markersClusterUsage$Gene] = paste0(rownames(tmp)[rownames(tmp) %in% markersClusterUsage$Gene],'*')
  ss[ss%in%cLabs$res.1] = cLabs[ss[ss%in%cLabs$res.1],'ClusterID']
  #Define the ordering of the columns
  #First, order numerically
  o = order(as.numeric(gsub('[A-Z]+([0-9]+).*','\\1',colnames(tmp))))
  tmp = tmp[,o]
  #Now do custom ordering for special compartments
  if(cPrefix[srcDir]=='N'){
    o = c('G','PT1','PT2','PT3','H','D','C1','C2','P','U1','U2','U3','U4','F','M','GE','DV','AV1','AV2')
    o = match(paste0(o,')'),gsub('.*\\(','',colnames(tmp)))
    o = c(o,seq_len(ncol(tmp))[-o])
  }else if(cPrefix[srcDir]=='F'){
    o = c('UB','CM','PV','fE','fM','fF','fGa')
    o = match(paste0(o,')'),gsub('.*\\(','',colnames(tmp)))
    o = c(o,seq_len(ncol(tmp))[-o])
  }else{
    #Otherwise just keep numeric ordering
    o = seq_len(ncol(tmp))
  }
  tmp = tmp[,o]
  #Order genes in the same way, then stick lit markers at end
  o = order(ifelse(ss%in%colnames(tmp),match(ss,colnames(tmp)),ncol(tmp)+1))
  ss = ss[o]
  tmp = tmp[o,]
  #Row-norm and use the default heatmap colour
  tmp = t(scale(t(tmp)))
  tmp[tmp<rowNormLims[1]]=rowNormLims[1]
  tmp[tmp>rowNormLims[2]]=rowNormLims[2]
  #Drop pure NA rows
  w = apply(tmp,1,function(e) any(!is.na(e)))
  tmp = tmp[w,]
  ss = ss[w]
  #Make split names fit
  ss = gsub(' ','\n',ss)
  #Make sure the order is maintained
  ss  = factor(ss,levels=unique(ss))
  gg_hm = Heatmap(tmp,
          #col = circlize::colorRamp2(c(-3,0,3),c('black','green','purple')),
          col = circlize::colorRamp2(c(rowNormLims[1],0,rowNormLims[2]), c("blue", "#EEEEEE", "red")),
          name='z-score of cells in cluster\nexpressing gene',
          row_title='Genes',
          row_names_gp = gpar(fontsize=8),
          column_title='Clusters',
          cluster_columns=FALSE,
          cluster_rows=TRUE,
          split=ss
          #column_order = colnames(tmp)[order(as.numeric(gsub('[A-Za-z]+','',colnames(tmp))))]
          )
  gA = ggplotGrob(gg_tsne)
  gB = grid.grabExpr(draw(gg_hm))
  maxHeight = grid::unit.pmax(gA$heights, gB$heights)
  gA$heights =  as.list(maxHeight)
  gB$heights =  as.list(maxHeight)
  gridExtra::grid.arrange(gA,gB,nrow=1,newpage=FALSE)
  dev.off()
}

#################################
# foetPseudotime and foetTabTFs #
#################################
#First need to do the two different DE analysis
##################################
# DE analysis of UB versus CM/PV
srcDir = file.path(baseDir,'/barcodeLevelGroup_foetalEpitheliumAndVascularV3/batchCorrected/')
mDat = readRDS(file.path(srcDir,'metadata.RDS'))
toc = readRDS(file.path(srcDir,'tableOfCounts.RDS'))
rownames(mDat$meta.data) = paste0(mDat$meta.data$Sanger_study_ID,'___',gsub('^[0-9]+_','',rownames(mDat$meta.data)))
colnames(toc) = rownames(mDat$meta.data)
#Define which clusters we're interested in
#tgtClusters = as.character(c(11,8,7))
tgtClusters = as.character(c(11,8,7))
w = which(mDat$meta.data$res.1 %in% tgtClusters)
#Make a monocle object with target cells
pd = new('AnnotatedDataFrame',data=mDat$meta.data[w,])
#Create a gene object thing
tmp = data.frame(EnsGeneID=gsub('.*_','',rownames(rtoc)),
                 Symbol = gsub('_.*','',rownames(rtoc)),
                 row.names=rownames(rtoc))
tmp$gene_short_name=tmp$Symbol
fd = new("AnnotatedDataFrame",data=tmp)
inDat = rtoc[,rownames(mDat$meta.data)[w]]
#devGenes = c("PAX2_ENSG00000075891","HNF1B_ENSG00000275410","PROM2_ENSG00000155066","GATA3_ENSG00000107485","SIX1_ENSG00000126778","SIX2_ENSG00000170577","WT1-AS_ENSG00000183242","WT1_ENSG00000184937")
#devGenes = c('WT1_ENSG00000184937','HNF1B_ENSG00000275410','WT1-AS_ENSG00000183242','GATA3_ENSG00000107485','SIX1_ENSG00000126778','PODXL_ENSG00000128567')
#devGenes = unique(c(c('WT1_ENSG00000184937','HNF1B_ENSG00000275410','WT1-AS_ENSG00000183242','GATA3_ENSG00000107485','SIX1_ENSG00000126778','PODXL_ENSG00000128567'),unlist(lapply(mDat$clusterSummaries[tgtClusters],function(e) rownames(e$markerGeneDat[!e$markerGeneDat$isKeyGene,])[seq(min(50,sum(!e$markerGeneDat$isKeyGene)))]),use.names=FALSE)))
#Do the processing
mcle = newCellDataSet(inDat,phenoData=pd,featureData=fd,expressionFamily=negbinomial.size())
#Get rid of uninformative genes
keepGenes = rownames(exprs(mcle))[Matrix::rowSums(exprs(mcle)>0)>3]
mcle = mcle[keepGenes,]
mcle = estimateSizeFactors(mcle)
mcle = estimateDispersions(mcle)
#Do the differential expression
pData(mcle)$DevType = ifelse(pData(mcle)$res.1=='11','UB','CM/PV')
deGenes = differentialGeneTest(mcle,fullModelFormulaStr = "~DevType")
tfs = read.table(fp_tfs,sep='\t',header=TRUE)
tfDE = deGenes[deGenes$EnsGeneID %in% tfs$Ensembl.ID,]
#Save the list of DE genes and the expression
UB_DE = tfDE
UB_mcle = mcle
##Get the expression matrix to add expression in two groups
#exprMat = Matrix::t(Matrix::t(exprs(mcle[rownames(tfDE),]))/sizeFactors(mcle))
#tmp = lapply(split(seq(ncol(exprMat)),pData(mcle)$DevType),function(e) Matrix::rowMeans(exprMat[,e,drop=FALSE]))
#tmp = do.call(cbind,tmp)
#tfDE = cbind(tfDE,tmp[rownames(tfDE),])
#tfDE = tfDE[order(tfDE$qval),]
##Write out the table
#write.table(tfDE[,c(3,4,5,6,8,9)],file.path(plotDir,'Table6Left.tsv'),sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)
###########################
# Pseudotime of CM and PV
srcDir = file.path(baseDir,'/barcodeLevelGroup_foetalEpitheliumAndVascularV3/batchCorrected/')
mDat = readRDS(file.path(srcDir,'metadata.RDS'))
toc = readRDS(file.path(srcDir,'tableOfCounts.RDS'))
rownames(mDat$meta.data) = paste0(mDat$meta.data$Sanger_study_ID,'___',gsub('^[0-9]+_','',rownames(mDat$meta.data)))
colnames(toc) = rownames(mDat$meta.data)
#Define which clusters we're interested in
tgtClusters = tgtClusters[tgtClusters!=11]
w = which(mDat$meta.data$res.1 %in% tgtClusters)
#Make a monocle object with target cells
pd = new('AnnotatedDataFrame',data=mDat$meta.data[w,])
#Create a gene object thing
tmp = data.frame(EnsGeneID=gsub('.*_','',rownames(rtoc)),
                 Symbol = gsub('_.*','',rownames(rtoc)),
                 row.names=rownames(rtoc))
tmp$gene_short_name=tmp$Symbol
fd = new("AnnotatedDataFrame",data=tmp)
inDat = rtoc[,rownames(mDat$meta.data)[w]]
#Do the processing
mcle = newCellDataSet(inDat,phenoData=pd,featureData=fd,expressionFamily=negbinomial.size())
#Get rid of uninformative genes
keepGenes = rownames(exprs(mcle))[Matrix::rowSums(exprs(mcle)>0)>3]
#But make sure we keep the ones that are DE
keepGenes = unique(c(keepGenes,rownames(UB_DE[UB_DE$qval < 0.01,])))
mcle = mcle[keepGenes,]
mcle = estimateSizeFactors(mcle)
mcle = estimateDispersions(mcle)
#Get highly expressed, high dispersion genes to use as base set of variable genes
dispTab = dispersionTable(mcle)
varGenes= dispTab$gene_id[dispTab$mean_expression>0.1 & dispTab$dispersion_empirical >= dispTab$dispersion_fit]
mcle = setOrderingFilter(mcle,varGenes)
#Now we're up to doing the pseudotime inference
mcle = reduceDimension(mcle,verbose=TRUE,max_components=2,norm_method='log',pseudo_expr=1)
mcle = orderCells(mcle)
##Plot the result
#pdf(file.path(plotDir,'Figure2d.pdf'))
#plot(plot_cell_trajectory(mcle,color_by='res.1'))
#dev.off()
#Save the co-ordinates in a text file
#pData(mcle)$DDR1 = mcle@reducedDimS[1,rownames(pData(mcle))]
#pData(mcle)$DDR2 = mcle@reducedDimS[2,rownames(pData(mcle))]
#write.table(pData(mcle),file.path(plotDir,'kidDevPseudoTimeCoords.tsv'),sep='\t',row.names=TRUE,col.names=TRUE,quote=FALSE)
##Alternative plot on the tSNE
#mDat$meta.data$Pseudotime=NA
#mDat$meta.data$Pseudotime[match(rownames(pData(mcle)),rownames(mDat$meta.data))]=pData(mcle)$Pseudotime
##
##Drop the blobby part of plot
#gg = ggplot(mDat$meta.data[mDat$meta.data$res.1 != '0' & !is.na(mDat$meta.data$Pseudotime),],aes(tSNE1,tSNE2,fill=Pseudotime)) +
#  geom_point(shape=21,colour='black',size=10) +
#  scale_fill_gradient(low='white',high='black') +
#  ggtitle("Pseudotime ordering of developing Nephron")
#pdf(file.path(plotDir,'kidDevPseudoTime_tSNE.pdf'))
#plot(gg)
#dev.off()
#write.table(mDat$meta.data,file.path(plotDir,'developingNephronPseudotime_tSNE.tsv'),sep='\t',row.names=TRUE,col.names=TRUE,quote=FALSE)
#Find DE genes along this trajectory
modelStr = '~sm.ns(Pseudotime,df=3)'
ptimeDE = differentialGeneTest(mcle,fullModelFormulaStr = modelStr,cores=8)
#Use harshest multiple hypothesis correction
ptimeDE$qval = p.adjust(ptimeDE$pval,method='bonferroni')
ptimeDE = ptimeDE[order(ptimeDE$qval),]
#Get the reduced TF matrix
tfs = read.table(fp_tfs,sep='\t',header=TRUE)
tfDE = ptimeDE[ptimeDE$EnsGeneID %in% tfs$Ensembl.ID,]
#######################
# Now combine the two
#The heatmap basically just shows smoothed and scaled expression data, where the expression data is calculated as exprs(dat)/Size_Factor then passed through vstExprs to normalise
#We want to use the true DE mcle object to pull out all genes (those DE either in the pseudotime or true DE), and scale and cluster them as one block.  Then we will split and plot differently but with the same colour scheme.
oldDE = read.table('~/Table6.tsv',sep='\t',header=TRUE)
rownames(oldDE) = paste0(oldDE$Symbol,'_',oldDE$EnsGeneID)
oldDE$qval = p.adjust(oldDE$pval,method='bonferroni')
oldGenes = rownames(oldDE)[oldDE$qval<0.01]
#Work out which genes are DE
shared = table(c(rownames(tfDE),rownames(UB_DE)))
shared = names(shared)[shared==2]
qval = pmin(p.adjust(tfDE[shared,'pval'],method='bonferroni'),p.adjust(UB_DE[shared,'pval'],method='bonferroni'))
names(qval) = shared
newGenes = shared[qval<0.01]
genesToPlot = newGenes
#First get the smoothed data 
cds = mcle[genesToPlot,]
ptData = data.frame(Pseudotime = seq(min(pData(cds)$Pseudotime),max(pData(cds)$Pseudotime),length.out=100))
mat = genSmoothCurves(cds,cores = 1,trend_formula = modelStr,relative_expr=TRUE,new_data = ptData)
#Normalise data
mat = vstExprs(cds,expr_matrix=mat)
#Make the DE into a homogeneous block
#Now get the data for cluster 11 only
cds = UB_mcle[genesToPlot,pData(UB_mcle)$res.1=='11']
tmp = exprs(cds)
tmp = t(t(tmp)/sizeFactors(UB_mcle)[colnames(tmp)])
tmp = vstExprs(cds,expr_matrix=tmp)
#Construct the final matrix
mat = cbind(outer(rowMeans(tmp),rep(1,10)),mat)
#Set NAs to 0
mat[is.na(mat)]=0
#mat = cbind(tmp,mat)
#Scale
mat = t(scale(t(mat)))
mat[mat< rowNormLims[1]] = rowNormLims[1]
mat[mat> rowNormLims[2]] = rowNormLims[2]
#Order the UB cells
#hc = hclust(dist(t(mat[,seq(ncol(tmp))])))
#mat = mat[,c(hc$order,seq(ncol(tmp)+1,ncol(mat)))]
#Finally, construct the heatmap with the cool colours
bks = seq(-3.1, 3.1, by = 0.1)
hmcols = monocle:::blue2green2red(length(bks) - 1)
rownames(mat) = gsub('_.*','',rownames(mat))
pdf(file.path(plotDir,'foetPseudotime.pdf'),width=14,height=14)
ph_res = pheatmap::pheatmap(mat, useRaster = TRUE, cluster_cols = FALSE,
    cluster_rows = TRUE, show_rownames = TRUE,
    show_colnames = FALSE, cutree_rows = 3,
    treeheight_row = 20,
    breaks = bks, fontsize = 6, color = hmcols, silent = TRUE,
    filename = NA)
grid::grid.rect(gp = grid::gpar("fill", col = NA))
grid::grid.draw(ph_res$gtable)
dev.off()
######################################
# Create a scored table to write out
##Get the expression matrix to add expression in two groups
exprMat = Matrix::t(Matrix::t(exprs(UB_mcle[rownames(UB_DE),]))/sizeFactors(UB_mcle))
pData(UB_mcle)$DevType2 = ifelse(pData(UB_mcle)$res.1=='11',
                                 'UB',
                                 ifelse(pData(UB_mcle)$res.1=='7',
                                        'CM',
                                        'PV'
                                        )
                                 )
tmp = lapply(split(seq(ncol(exprMat)),pData(UB_mcle)$DevType2),function(e) Matrix::rowMeans(exprMat[,e,drop=FALSE]))
tmp = do.call(cbind,tmp)
#Now get the statistics
out = data.frame(EnsGeneID = gsub('.*_','',shared),
                 Symbol = gsub('_.*','',shared),
                 PseudoTimePval = tfDE[shared,'pval'],
                 DE_Pval = UB_DE[shared,'pval'],
                 qval = qval[shared]
                 )
#Add in expression
out = cbind(out,tmp[shared,c('UB','CM','PV')])
out = out[order(out$qval),]
#Write out the table
write.table(out,file.path(plotDir,'foetTabTFs.tsv'),sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)

################################################
# wilmsPseudotime and wilmsTFs and wilmsTabTFs #
################################################
# Pseudotime of Wilms and NR 
#Wilms + NR cells
srcDir = file.path(baseDir,'barcodeLevelGroup_tumourEpitheliumAndVascular/batchCorrected/')
mDat = readRDS(file.path(srcDir,'metadata.RDS'))
toc = readRDS(file.path(srcDir,'tableOfCounts.RDS'))
rownames(mDat$meta.data) = paste0(mDat$meta.data$Sanger_study_ID,'___',gsub('^[0-9]+_','',rownames(mDat$meta.data)))
colnames(toc) = rownames(mDat$meta.data)
#All cells that could be tumour or NR
tgtClusters = as.character(c(5,18,16))
w = which(mDat$meta.data$res.1 %in% tgtClusters & grepl('Wilms',mDat$meta.data$Experiment))
labs = paste0(mDat$meta.data$Experiment[w],'_',mDat$meta.data$TissueDiseaseState[w])
pDat = cbind(mDat$meta.data[w,],cellType=labs)
#Finalise and group everything
pd = new('AnnotatedDataFrame',data=pDat)
#Create a gene object thing
tmp = data.frame(EnsGeneID=gsub('.*_','',rownames(rtoc)),
                 Symbol = gsub('_.*','',rownames(rtoc)),
                 row.names=rownames(rtoc))
tmp$gene_short_name=tmp$Symbol
fd = new("AnnotatedDataFrame",data=tmp)
inDat = rtoc[,rownames(pDat)]
#Do the processing
mcle = newCellDataSet(inDat,phenoData=pd,featureData=fd,expressionFamily=negbinomial.size())
#Get rid of uninformative genes
keepGenes = rownames(exprs(mcle))[Matrix::rowSums(exprs(mcle)>0)>3]
mcle = mcle[keepGenes,]
mcle = estimateSizeFactors(mcle)
mcle = estimateDispersions(mcle)
#Get highly expressed, high dispersion genes to use as base set of variable genes
dispTab = dispersionTable(mcle)
varGenes= dispTab$gene_id[dispTab$mean_expression>0.1 & dispTab$dispersion_empirical >= dispTab$dispersion_fit]
mcle = setOrderingFilter(mcle,varGenes)
#Now we're up to doing the pseudotime inference
mcle = reduceDimension(mcle,verbose=TRUE,max_components=2,norm_method='log',pseudo_expr=1)
mcle = orderCells(mcle)
#Set the root state properly
mcle = orderCells(mcle,root_state=3)
#Make the result fancy
df = pData(mcle)
df$DDR1 = mcle@reducedDimS[1,rownames(pData(mcle))]
df$DDR2 = mcle@reducedDimS[2,rownames(pData(mcle))]
#Add each cells similarity to UB and PV tumour pops
testDat = sim$loadTrainingData(file.path(baseDir,'barcodeLevelGroup_tumourEpitheliumAndVascular/batchCorrected/'))
dat = Seurat::LogNormalize(rtoc[,colnames(testDat$toc)])
dat = t(dat)
#Load the trained models
fits = readRDS(file.path(baseDir,'lrFoetalClustersV4.RDS'))
#Predict on tumour data.
preds = list()
for(mark in c('8','11')){
  message(sprintf("Predicting probabilities for cluster %s",mark))
  preds[[mark]] = predict(fits[[mark]],newx=dat[,rownames(fits[[mark]]$glmnet.fit$beta)],s='lambda.1se',newoffset=rep(0,nrow(dat)))
}
#Now make a matrix
pp = do.call(cbind,preds)
colnames(pp) = c('PV','UB')
pp = (1+exp(-pp))**-1
df = cbind(df,pp[rownames(df),])
#Create a combined colour for each point
#Take the two probabilities (pr(UB),pr(PV)) and project them onto a line joining (1,0) and (0,1)
df$jointProb = ((1-df$UB)+df$PV)/2
#Estimate the point desity
dens = MASS::kde2d(df$DDR1,df$DDR2)
gr = data.frame(with(dens, expand.grid(x,y)), as.vector(dens$z))
names(gr) = c("xgr", "ygr", "zgr")
# Fit a model
mod = loess(zgr~xgr*ygr, data=gr)
df$dens = predict(mod, newdata=data.frame(xgr=df$DDR1, ygr=df$DDR2))
#Transparency setter.  Get the three centroids (these values come from k-mediods clustering with k=3)
centroids = matrix(c(6.7,-0.6,-2.2,4.7,-4.5,-2.3),nrow=3,byrow=TRUE)
#Get distance from centroid
df$distFromClust = apply(sqrt(apply(centroids,1,function(e) rowSums((t(t(df[,c('DDR1','DDR2')])-e))**2))),1,min)
#Manually add jitter based on distance from cluster
width=height=0.8
dd = exp(-df$distFromClust)
dd = dd/max(dd)
dx = rnorm(nrow(df),0,dd*width)
dx[abs(dx)>width] = sign(dx)[abs(dx)>width] * width
dy = rnorm(nrow(df),0,dd*height)
dy[abs(dy)>height] = sign(dy)[abs(dy)>height] * height
#Create jittered points
df$x = df$DDR1 + dx
df$y = df$DDR2 + dy
pdf(file.path(plotDir,'wilmsPseudotime.pdf'),width=14,height=14)
gg = ggplot(df,aes(-y,x,colour=jointProb,shape=TissueDiseaseState)) +
  #Add un-jittered points in black at bottom
  geom_point(aes(x=-DDR2,y=DDR1),inherit.aes=FALSE,alpha=1/10,size=3.0) +
  #Do two plots to control order of plotting
  geom_point(data=df[df$TissueDiseaseState=='NephrogenicRest',],size=6.0) +
  geom_point(data=df[df$TissueDiseaseState=='Tumour',],size=6.0) +
  scale_colour_gradientn(colours=c('#009FE3','#CA9E67')) +
  scale_shape_manual(values=c('NephrogenicRest'=0,'Tumour'=3)) +
  labs(colour='Probability \nUB to PV')
plot(gg)
dev.off()
#Do TF thing
branchDE = BEAM(mcle, branch_point = 1, cores = 8)
branchDE = branchDE[order(branchDE$qval),]
tfs = read.table(fp_tfs,sep='\t',header=TRUE)
tfDE = branchDE[branchDE$EnsGeneID %in% tfs$Ensembl.ID,]
pdf(file.path(plotDir,'wilmsTFs.pdf'))
plot_genes_branched_heatmap(mcle[row.names(tfDE[tfDE$qval<.01,]),],num_clusters=3,branch_states=c(1,2),branch_labels=c('Tumour','NR'),use_gene_short_name=TRUE,show_rownames=TRUE)
dev.off()
#Add in the gene expression to the list
exprMat = Matrix::t(Matrix::t(exprs(mcle[rownames(tfDE),]))/sizeFactors(mcle))
#Split by pseudotime bins
ptBins = split(rownames(pData(mcle)),pData(mcle)$State)
tmp = lapply(ptBins,function(e) Matrix::rowMeans(exprMat[,e,drop=FALSE]))
tmp = do.call(cbind,tmp)
colnames(tmp) = c('PV_NephRest','PV_Tumor','UB')
tfDE = cbind(tfDE,tmp)
#Record the list of the significant TFs
write.table(tfDE[,c(3,4,5,6,9,10,11)],file.path(plotDir,'wilmsTabTFs.tsv'),sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)
#phyper(sum(tfDE[match(gsub('.*_','',genesToPlot),tfDE$EnsGeneID),'qval']<.01,na.rm=TRUE),sum(!is.na(tfDE[match(gsub('.*_','',genesToPlot),tfDE$EnsGeneID),'qval'])),nrow(tfDE)-sum(!is.na(tfDE[match(gsub('.*_','',genesToPlot),tfDE$EnsGeneID),'qval'])),sum(tfDE$qval<.01))



###################
# wilmsSimilarity #
###################
# Similarity of tumours to foetal clusters
#Assumes predictors have already been trained elsewhere
testDat = sim$loadTrainingData(file.path(baseDir,'barcodeLevelGroup_tumourEpitheliumAndVascular/batchCorrected/'))
dat = Seurat::LogNormalize(rtoc[,colnames(testDat$toc)])
dat = t(dat)
#Load the trained models
fits = readRDS(file.path(baseDir,'lrFoetalClustersV4.RDS'))
#Predict on tumour data.
#pdf(file.path(plotDir,'SupFigure4a.pdf'))
preds = list()
for(mark in names(fits)){
  message(sprintf("Predicting probabilities for cluster %s",mark))
  if(mark=='13')
    next
  preds[[mark]] = predict(fits[[mark]],newx=dat[,rownames(fits[[mark]]$glmnet.fit$beta)],s='lambda.1se',newoffset=rep(0,nrow(dat)))
  #Make a plot of the logits in the new space
  testDat$mDat$logits[match(rownames(preds[[mark]]),rownames(testDat$mDat))] = preds[[mark]][,1]
  #Truncate logits at +/- 5
  m = abs(testDat$mDat$logits)>5
  testDat$mDat$logits[m]=5*sign(testDat$mDat$logits)[m]
  gg = ggplot(testDat$mDat,aes(tSNE1,tSNE2,fill=logits)) +
    geom_point(shape=21,colour='black',size=0.8,stroke=1/10) + 
    scale_fill_gradientn(colours=c('blue','white','red'),
                           values = scales::rescale(c(min(testDat$mDat$logits),0,max(testDat$mDat$logits))),
                           limits = range(testDat$mDat$logits)) +
    labs(colour=sprintf('LogOddsOf_F%s',mark)) +
    ggtitle(sprintf("Predicted log odds of Foetal cluster F%s",mark))
  #plot(gg)
}
#dev.off()
#Now make a matrix
pp = do.call(cbind,preds)
colnames(pp) = names(preds)
#clustToGroups = setNames(rep(names(tumourGroups),lengths(tumourGroups)),unlist(tumourGroups))
#clusterPreds = apply(pp,2,function(e) sapply(split(e,clustToGroups[testDat$mDat$res.1]),mean))
#Alternative
clusterPreds = apply(pp,2,function(e) sapply(split(e,testDat$mDat$res.1),mean))
rownames(clusterPreds) = paste0('T',rownames(clusterPreds))
colnames(clusterPreds) = paste0('F',colnames(clusterPreds))
#Convert to percentages
clusterPreds = (1+exp(-clusterPreds))**-1
#Simple grouped thing where we just re-order manually
m = match(rownames(clusterPreds),paste0('T',unlist(tumourGroups,use.name=FALSE)))
row_splitter = factor(rep(names(tumourGroups),lengths(tumourGroups))[m],levels=names(tumourGroups))
gg = Heatmap(clusterPreds,
             name='Predicted\nSimilarity',
             column_title='Foetal Clusters',
             row_title='Tumour Clusters',
             show_row_names=TRUE,
             split=row_splitter,
             row_order=match(paste0('T',unlist(tumourGroups,use.names=FALSE)),rownames(clusterPreds)),
             column_order=match(paste0('F',unlist(foetGroups,use.names=FALSE)),colnames(clusterPreds)),
             cluster_rows=FALSE,
             cluster_columns=FALSE
             )
pdf(file.path(plotDir,'wilmsSimilarity.pdf'),width=14,height=14)
draw(gg)
dev.off()

#################
# rccSimilarity #
#################
# Similarity of tumours to normal clusters 
#Assumes predictors have already been trained elsewhere
#Load test data (tumour)
testDat = sim$loadTrainingData(file.path(baseDir,'barcodeLevelGroup_tumourEpitheliumAndVascular/batchCorrected/'))
dat = Seurat::LogNormalize(rtoc[,colnames(testDat$toc)])
dat = t(dat)
#Load the trained models
#fits = readRDS(file.path(baseDir,'lrNormalClustersV2.RDS'))
fits = readRDS(file.path(baseDir,'lrNormalClustersWithControl.RDS'))
#pdf(file.path(plotDir,'SupFigure4b.pdf'))
preds = list()
for(mark in names(fits)){
  message(sprintf("Predicting probabilities for cluster %s",mark))
  preds[[mark]] = predict(fits[[mark]],newx=dat[,rownames(fits[[mark]]$glmnet.fit$beta)],s='lambda.1se',newoffset=rep(0,nrow(dat)))
  if(all(preds[[mark]]==0))
    next
  #Make a plot of the logits in the new space
  testDat$mDat$logits[match(rownames(preds[[mark]]),rownames(testDat$mDat))] = preds[[mark]][,1]
  #Truncate logits at +/- 5
  m = abs(testDat$mDat$logits)>5
  testDat$mDat$logits[m]=5*sign(testDat$mDat$logits)[m]
  gg = ggplot(testDat$mDat,aes(tSNE1,tSNE2,fill=logits)) +
    geom_point(shape=21,colour='black',size=0.8,stroke=1/10) + 
    scale_fill_gradientn(colours=c('blue','white','red'),
                           values = scales::rescale(c(min(testDat$mDat$logits),0,max(testDat$mDat$logits))),
                           limits = range(testDat$mDat$logits)) +
    labs(colour=sprintf('LogOddsOf_N%s',mark)) +
    ggtitle(sprintf("Predicted log odds of Normal cluster N%s",mark))
  #plot(gg)
}
#dev.off()
#Now make a matrix
pp = do.call(cbind,preds)
colnames(pp) = names(preds)
#clustToGroups = setNames(rep(names(tumourGroups),lengths(tumourGroups)),unlist(tumourGroups))
#clusterPreds = apply(pp,2,function(e) sapply(split(e,clustToGroups[testDat$mDat$res.1]),mean))
#Alternative
clusterPreds = apply(pp,2,function(e) sapply(split(e,testDat$mDat$res.1),mean))
rownames(clusterPreds) = paste0('T',rownames(clusterPreds))
colnames(clusterPreds) = paste0('N',colnames(clusterPreds))
#Convert to percentages
clusterPreds = (1+exp(-clusterPreds))**-1
#Simple grouped thing where we just re-order manually
m = match(rownames(clusterPreds),paste0('T',unlist(tumourGroups,use.name=FALSE)))
row_splitter = factor(rep(names(tumourGroups),lengths(tumourGroups))[m],levels=names(tumourGroups))
gg = Heatmap(clusterPreds,
             name='Predicted\nSimilarity',
             column_title='Normal Clusters',
             row_title='Tumour Clusters',
             show_row_names=TRUE,
             split=row_splitter,
             row_order=match(paste0('T',unlist(tumourGroups,use.names=FALSE)),rownames(clusterPreds)),
             column_order=match(paste0('N',unlist(normGroups,use.names=FALSE)),colnames(clusterPreds)),
             cluster_rows=FALSE,
             cluster_columns=FALSE
             )
pdf(file.path(plotDir,'rccSimilarity.pdf'),width=14,height=14)
draw(gg)
dev.off()

####################################################
# TCGA boxplots for many markers and globalMarkers #
####################################################
###############################################
# Find markers of specific foetal populations
#Cut-offs for selection of markers
inOverOut = 0.2
outMax = 0.1
#Get the marker populations
foetalMarks = list()
allDat = list()
for(clustNum in grep('^F[0-9]+$',colnames(geneFrac),value=TRUE)){
  clustNum = gsub('F','',clustNum)
  print(clustNum)
  #Define which columns to not look at when filtering out-of-group
  oog = c(grep('T',colnames(geneFrac),value=TRUE),paste0('F',clustNum))
  nephGroups = c(0,1,4,7,8,10,11)
  if(clustNum %in% c('7','8','11'))
    oog = unique(c(oog,paste0('F',nephGroups)))
  inGroup = geneFrac[,paste0('F',clustNum)]
  outGroup = apply(geneFrac[,-match(oog,colnames(geneFrac)),drop=FALSE],1,max)
  nephGroup = apply(geneFrac[,paste0('F',nephGroups[nephGroups!=as.integer(clustNum)])],1,max)
  tumGroup = geneFrac[,c(paste0('T',c(5,16)),paste0('F',c(11,7,8)))]
  dat = data.frame(inGroup=inGroup,
                   outGroup=outGroup,
                   nephGroup=nephGroup)
  dat = cbind(dat,tumGroup)
  allDat[[paste0('F',clustNum)]] = dat[order(dat$inGroup-dat$outGroup,decreasing=TRUE),]
}
#Pick the markers based on the sort
foetalMarks = lapply(allDat,function(e) e[e$inGroup > e$outGroup+inOverOut & e$outGroup<outMax,])
#Further require that the ones where we've excluded the rest of the nephron at least be highest in the part of the nephron we're saying they're specific for.
for(m in c('F7','F8','F11')){
  foetalMarks[[m]] = foetalMarks[[m]][foetalMarks[[m]]$inGroup > foetalMarks[[m]]$nephGroup,]
}
fMarks = setNames(unlist(lapply(foetalMarks,rownames),use.names=FALSE),rep(names(foetalMarks),sapply(foetalMarks,nrow)))
fMarks = fMarks[!duplicated(fMarks)]
finalMarkers = data.frame(gene = gsub('_.*','',fMarks),
                          markType = '+',
                          cellType = ifelse(names(fMarks)=='F7',
                                        'CM',
                                        ifelse(names(fMarks)=='F8',
                                               'PV',
                                               ifelse(names(fMarks)=='F11',
                                                      'UB',
                                                      names(fMarks)
                                                      )
                                               )
                                        ))
###################################################
# Save globally specific statistics for UB and PV
gStats = rbind(allDat$F8,allDat$F11)
gStats$Symbol = gsub('_.*','',rownames(gStats))
gStats$EnsemblID = gsub('.*_','',c(rownames(allDat$F8),rownames(allDat$F11)))
rownames(gStats)=NULL
gStats$MarkerOf = rep(c('PV','UB'),c(nrow(allDat$F8),nrow(allDat$F11)))
gStats$MeetsCriteria = FALSE
gStats$MeetsCriteria[paste0(gStats$Symbol,'_',gStats$EnsemblID) %in% rownames(foetalMarks$F8) & gStats$MarkerOf=='PV']=TRUE
gStats$MeetsCriteria[paste0(gStats$Symbol,'_',gStats$EnsemblID) %in% rownames(foetalMarks$F11) & gStats$MarkerOf=='UB']=TRUE
#Sort to put markers at top
gStats = gStats[order(gStats$MeetsCriteria,gStats$MarkerOf,decreasing=TRUE),]
#Keep key columns and re-order
gStats = gStats[,c(9,10,11,1,2,3,12)]
colnames(gStats)[4:6] = c('MarkerClusterFrequency','MaxOutOfClusterFrequency','MaxFetalNephFrequencyExcludingMarker')
#########################
# And of PT populations 
allDat = list()
for(clustNum in grep('^N[0-9]+$',colnames(geneFrac),value=TRUE)){
  clustNum = gsub('N','',clustNum)
  print(clustNum)
  #Define which columns to not look at when filtering out-of-group
  oog = c(grep('T',colnames(geneFrac),value=TRUE),grep('^EN[0-9]+',colnames(geneFrac),value=TRUE),grep('^PT[0-9]+',colnames(geneFrac),value=TRUE),paste0('N',clustNum))
  inGroup = apply(geneFrac[,paste0('N',clustNum),drop=FALSE],1,max)
  outGroup = apply(geneFrac[,-match(oog,colnames(geneFrac)),drop=FALSE],1,max)
  tumGroup = geneFrac[,grep('^T',colnames(geneFrac))]
  dat = data.frame(inGroup=inGroup,
                   outGroup=outGroup
                   )
  dat = cbind(dat,tumGroup)
  allDat[[paste0('N',clustNum)]] = dat[order(dat$inGroup-dat$outGroup,decreasing=TRUE),]
}
#Pick the markers based on the sort
normMarks = lapply(allDat,function(e) e[(e$inGroup > e$outGroup+inOverOut & e$outGroup<outMax),])
nMarks = setNames(unlist(lapply(normMarks,rownames),use.names=FALSE),rep(names(normMarks),sapply(normMarks,nrow)))
nMarks = nMarks[!duplicated(nMarks)]
finalMarkers = rbind(finalMarkers,
                     data.frame(gene = gsub('_.*','',nMarks),
                                markType = '+',
                                cellType = names(nMarks))
                     )
################################
# Add literature based markers 
finalMarkers = rbind(finalMarkers,
                     data.frame(gene = c('SLC17A3','SLC7A13','VCAM1','SLC12A1','KCNJ1','SLC4A1','SLC26A4','PDGFRB','ACTA2','SEMA3G','SLC14A1','PLVAP','EMILIN1','SFRP2','MMP2','WT1','CA9','MET','SMARCB1','TPSAB1','TPSB2'),
                                markType = '+',
                                cellType = c('PT','PT3','PT1','LoH','DT','CD-B','CD-A','M','M','GE','DVR','AVR','F','F','F','WT','RCC','RCC','RT','MAST','MAST'))
                     )
finalMarkers = rbind(finalMarkers,
                     data.frame(gene = c('EGFR',
'NTRK1',
'NTRK3',
'TFE3',
'TFEB',
'BCOR',
'FGFR1'
),
                                markType='+',
                                cellType='MISC'))
#########################
# Finalise marker table 
#Clean up duplicates
finalMarkers = finalMarkers[!(finalMarkers[,1] %in% finalMarkers[duplicated(finalMarkers[,1]),1]) | !(finalMarkers$cellType %in% c('EN18','N12')),]
#Exclude crappy genes
finalMarkers = finalMarkers[grep('^RP[0-9]',finalMarkers[,1],invert=TRUE),]
#Add ENSEMBL ID
finalMarkers$ENS_ID = gsub('.*_','',rownames(rtoc))[match(finalMarkers$gene,gsub('_.*','',rownames(rtoc)))]
#### Validation uses 10X bulk
####Want to show that each marker is DE in foetal versus normal and DE (higher in tumour) in normal versus Wilms
####Similarly, show that the markers picked for RCC are DE (or not DE) in normal versus tumour in the way we would expect
####Finally, show the same pattern of DE using public, bulk RNA
####Find and load all the bulk RNA samples we have
###all10X = read.table('~/Projects/Common/all10X.txt',sep='\t',header=TRUE)
###all10X = all10X[all10X$CGP_RNA_ID!='',]
####Drop the ones that we don't actually have data for
###all10X = all10X[file.exists(file.path('~/PipeLive/1700/',all10X$CGP_RNA_ID,paste0(all10X$CGP_RNA_ID,'.count.gz'))),]
####Now construct a table of counts from each of them
###bulkRNA = unique(all10X$CGP_RNA_ID)
###fnoms = file.path('~/PipeLive/1700/',bulkRNA,paste0(bulkRNA,'.count.gz'))
###toc = lapply(fnoms,function(e) read.table(e,sep='\t'))
###noms = toc[[1]][,1]
###toc = do.call(cbind,lapply(toc,function(e) e[,2]))
###rownames(toc)=noms
###colnames(toc)=bulkRNA
###toc = toc[grep('ENSG',rownames(toc)),]
####Convert to sensible gene names and remove useless things
###mart =useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
###gns = rownames(toc)
###geneNames = getBM( attributes= c("ensembl_gene_id","hgnc_symbol"),filters="ensembl_gene_id",values=gns,mart= mart)
####Discard anything without a unique symbol
###symbs = geneNames[match(rownames(toc),geneNames[,1]),2]
###symbs[table(symbs)[symbs]!=1]=NA
###o = !is.na(symbs) & symbs!=''
###toc = toc[o,]
###rownames(toc) = paste0(symbs[o],'_',rownames(toc))
###geneDat = data.frame(row.names=rownames(toc),
###                     ENS_ID = gsub('.*_','',rownames(toc)),
###                     SYMB = gsub('_.*','',rownames(toc)))
####Create the sample manifest
###samps = all10X[match(colnames(toc),all10X$CGP_RNA_ID),]
####Now initiate the edgeR object
###samps$Group = paste0(gsub('[0-9]$','',samps$Experiment),'_',samps$TissueDiseaseState)
###dgl = DGEList(toc,samples=samps,genes=geneDat,group=samps$Group)
###dgl = calcNormFactors(dgl)
####Estimate dispersion
###design = model.matrix(~Group,data=samps)
###dgl = estimateDisp(dgl,design)
####Create a TPM object
###geneLens = read.table('/lustre/scratch117/casm/team274/references/ENSEMBL_genome_annotation_files/Human/Human_ENSEMBL_V75_gene_average_length.txt',sep='\t',header=TRUE)
###tpm = toc/geneLens[match(gsub('.*_','',rownames(toc)),geneLens[,1]),2]
###tpm = t(t(tpm)/colSums(tpm)*1e6)
###mat = tpm[match(finalMarkers$ENS_ID,gsub('.*_','',rownames(tpm))),]
###mat = t(t(mat)/dgl$samples[colnames(mat),'norm.factors'])
###ss = finalMarkers$cellType
###rownames(mat) = gsub('_.*','',rownames(mat))
###w = !is.na(mat[,1])
###mat = mat[w,]
###ss = ss[w]
###mat = log2(mat)
###mat[!is.finite(mat)]=NA
####Heatmap showing expression of these markers
###pdf(file.path(plotDir,'NativeBulk.pdf'),width=14,height=14)
###cols = setNames(seq(length(unique(dgl$samples[colnames(mat),'Group']))),unique(dgl$samples[colnames(mat),'Group']))
###topAnn = HeatmapAnnotation(dgl$samples[colnames(mat),'Group',drop=FALSE],col=list(Group=cols))
###Heatmap(mat,top_annotation=topAnn,split=ss,cluster_column=TRUE,cluster_rows=TRUE)
###dev.off()
####For each of the chosen markers, check if it's DE between tumour and normal.
##########################################
# Processing of the TARGET mRNA-seq data
srcDir = '~/scratch/KidneySC/TARGET_Bulk'
subTypes = c('RT','WT')
cnts = list()
for(subType in subTypes){
  print(subType)
  files = list.files(file.path(srcDir,subType))
  tmani = read.delim(file.path(srcDir,subType,grep('sdrf.txt$',files,value=TRUE)),sep='\t',header=TRUE)
  if('Derived.Array.Data.File.1' %in% colnames(tmani)){
    tmani = tmani[match(files,tmani$Derived.Array.Data.File.1),]
  }else{
    tmani = tmani[match(files,tmani$Derived.Array.Data.File),]
  }
  w = which(!is.na(tmani[,1]))
  samps = files[w]
  tmani = tmani[w,]
  #Load the actual data
  #CCSK uses cufflinks FPKM values which are a little fucked up.  Do our best to make it consistent
  dat = lapply(samps,function(e) {
               tmp = read.delim(file.path(srcDir,subType,e),sep='\t',header=TRUE)
               tmp[,c('gene','RPKM')]
                    })
  #Keep only the common genes
  cGenes = table(unlist(lapply(dat,function(e) e[,1])))
  cGenes = names(cGenes)[cGenes==length(dat)]
  dat = do.call(cbind,lapply(dat,function(e) e[match(cGenes,e[,1]),2]))
  rownames(dat) = cGenes
  cnts[[subType]] = list(dat=dat,tmani=tmani)
}
#Standardise the manifest
manifest = do.call(rbind,lapply(cnts,function(e) {
                                data.frame(submitter_id = e$tmani$Extract.Name,
                                           cases = e$tmani$Sample.Name,
                                           tissue.definition = e$tmani$Description,
                                           age_at_diagnosis = NA,
                                           disease = e$tmani$Characteristics.DiseaseState.)}))
#Keep detailed manifest
bigManifest = lapply(cnts,function(e) e$tmani)
#Collect a core set of genes
for(i in seq_along(cnts)){
  if(i==1)
    geneDat = rownames(cnts[[i]]$dat)
  geneDat = geneDat[geneDat%in%rownames(cnts[[i]]$dat)]
}
#Keep only these genes
cnts = do.call(cbind,lapply(cnts,function(e) e$dat[geneDat,]))
colnames(cnts) = manifest$submitter_id
rownames(manifest) = manifest$submitter_id
####################################################################
# Load and process what data we can from the TCGA biolinks service
datFile = 'TCGA_data.RDS'
#Fetch data from TCGA directly
projects = c('TARGET-RT','TARGET-WT','TCGA-KICH','TCGA-KIRP','TCGA-KIRC')
if(file.exists(file.path(baseDir,datFile))){
  dat = readRDS(file.path(baseDir,datFile))
}else{
  #Prepare to get all the data
  queries = lapply(projects,GDCquery,
           data.category = 'Transcriptome Profiling',
           data.type = 'Gene Expression Quantification',
           workflow.type = 'HTSeq - FPKM',
           legacy = FALSE)
  #Actually download it
  lapply(queries,GDCdownload)
  #Now load it and the associated meta-data
  dat = lapply(queries,GDCprepare,summarizedExperiment=TRUE)
  saveRDS(dat,file.path(baseDir,datFile))
}
names(dat) = projects
#Flatten out format
#The gene information is the same across all projects
geneDat = dat[[1]]@rowRanges
#Record the maximal manifest information
bigManifest = c(bigManifest,lapply(dat,function(e) e@colData))
#Make a simple manifest with just the core information
tmp = lapply(dat,function(e){
                 tmp = as.data.frame(e@colData)
                 #The TARGET samples
                 if(ncol(tmp)==42){
                   data.frame(row.names=rownames(tmp),
                              submitter_id = tmp$barcode,
                              cases = tmp$sample,
                              tissue.definition=tmp$definition,
                              age_at_diagnosis=tmp$age_at_diagnosis,
                              disease = tmp$tumor.definition)
                 }else{
                   data.frame(row.names=rownames(tmp),
                              submitter_id = tmp$barcode,
                              cases = tmp$sample,
                              tissue.definition = tmp$definition,
                              age_at_diagnosis = tmp$age_at_diagnosis,
                              disease = unlist(tmp$disease_type,use.names=FALSE))
                 }
           })
tmp = do.call(rbind,tmp)
rownames(tmp) = gsub('.*\\.','',rownames(tmp))
#Construct the big table of expression
tCnts = do.call(cbind,lapply(dat,function(e) e@assays[[1]]))
#Label everything
rownames(tCnts) = names(geneDat)
colnames(tCnts) = rownames(tmp)
#Keep only the genes for which we have data everywhere
geneDat = geneDat[names(geneDat) %in% rownames(cnts)]
tCnts = tCnts[match(names(geneDat),rownames(tCnts)),]
#Merge with existing cnts and manifest
cnts = cbind(cnts[match(names(geneDat),rownames(cnts)),],tCnts)
manifest = rbind(manifest,tmp)
################
# Cleanup data
#Anything that occurs multiple times, keep the first one
manifest = manifest[!duplicated(manifest$submitter_id),]
#Drop anything that isn't primary
manifest$tissue.definition[manifest$tissue.definition=='Primary solid Tumor']='Primary Solid Tumor'
manifest = manifest[manifest$tissue.definition %in% c('Primary Solid Tumor','Solid Tissue Normal'),]
#Simplify the disease classification
diseaseMap = c('Kidney Chromophobe'='KICH',
               'Kidney Renal Clear Cell Carcinoma'='KIRC',
               'Kidney Renal Papillary Cell Carcinoma'='KIRP',
               'Rhabdoid tumor (kidney) (RT)'='RT',
               'Rhabdoid Tumor of the Kidney'='RT',
               'Childhood Kidney Wilms Tumor'='WT',
               'Clear Cell Sarcoma of the Kidney'='CCSK',
               'Wilms tumor (WT)'='WT')
manifest$disease = diseaseMap[manifest$disease]
#Make extra columns
manifest$DiseaseState = ifelse(grepl('Normal',manifest$tissue.definition),'Normal','Tumour')
manifest$Group = ifelse(manifest$DiseaseState=='Tumour',manifest$disease,'Normal')
#Now we need to discard things so the manifest is concordant with the count table
cnts = cnts[,match(manifest$submitter_id,colnames(cnts))]
#Convert to TPM
cnts = t(t(cnts)/colSums(cnts))*1e6
##########################################
# Refine the grouping beyond the default
#Split normals into Paed and Adult
o=with(manifest,which((is.na(age_at_diagnosis) | age_at_diagnosis < 10*365) & Group=='Normal'))
manifest$Group[o]='PaedNormal'
##Add in the translocation RCCs
#transLocIDs = c("TCGA.AK.3456.01A.02R.1325.07",
#"TCGA.BP.4756.01A.01R.1289.07",
#"TCGA.BP.4758.01A.01R.1289.07",
#"TCGA.B8.5546.01A.01R.1541.07",
#"TCGA.B0.5705.01A.11R.1541.07",
#"TCGA.CJ.5681.01A.11R.1541.07",
#"TCGA.A3.3313.01A.02R.1325.07")
#transLocIDs = gsub('\\.','-',transLocIDs)
#w = which(manifest$submitter_id %in% transLocIDs)
#manifest$Group[w] = 'Transloc'
#manifest$disease[w] = 'Transloc'
#####################
# Add foetal counts
fCnts = read.delim('~/scratch/KidneySC/BulkFoetalExpression/foetalCounts.tsv',sep='\t',header=TRUE)
#Keep only those genes that we have TCGA data for
fCnts = fCnts[gsub('\\..*','',fCnts$gene_id) %in% names(geneDat),]
#Normalise them by gene length
geneLens = read.table('/lustre/scratch117/casm/team274/references/ENSEMBL_genome_annotation_files/Human/Human_ENSEMBL_V75_gene_average_length.txt',sep='\t',header=TRUE)
m = match(gsub('\\..*','',fCnts$gene_id),geneLens[,1])
#Drop the non-numeric columns
fGeneDat = fCnts[,1:10]
fCnts = fCnts[,-(1:15)]
#Convert to TPM
fCnts = fCnts/geneLens[m,2]
fCnts = t(t(fCnts)/colSums(fCnts))*1e6
#Add them to the counts
cnts = cbind(cnts,fCnts[match(names(geneDat),gsub('\\..*','',fGeneDat$gene_id)),c('Kidney_1','Kidney_2')])
#Add them to manifest
fMani = data.frame(row.names=colnames(fCnts),
                   submitter_id = colnames(fCnts),
                   cases = colnames(fCnts),
                   tissue.definition='Solid Tissue Normal',
                   age_at_diagnosis=NA,
                   disease=NA,
                   DiseaseState = 'Normal',
                   Group = 'Foetal')
manifest = rbind(manifest,fMani[c('Kidney_1','Kidney_2'),])
##############
# Add in CMN
cmn = read.delim('~/Projects/Common/CMN_filtered_STAR2_ENSv75_TPM_ENS_IDs.tsv',sep='\t',header=TRUE)
colnames(cmn) = paste0('CMN.',colnames(cmn))
#Keep only these samples, as per Sam's email detailing which are really CMN
cmn = cmn[,gsub('CMN.','',colnames(cmn)) %in% c("PR37199a","PR37200a","PR37202a","PR37204a","PR37208a","PR37210a","PR37212a","PR37213a","PR37214a","PR37216a","PR37218a","PR37219a","PR37221a","PR37222a","PR37223a","PR37224a","PR37225a","PR37205a")]
#Add in our data
cnts = cbind(cnts,cmn[match(names(geneDat),rownames(cmn)),])
cnts = t(t(cnts)/colSums(cnts,na.rm=TRUE))*1e6
#Add CMN to the manifest
cmnMani = data.frame(row.names=colnames(cmn),
                     submitter_id = colnames(cmn),
                     cases = colnames(cmn),
                     tissue.definition='Primary solid Tumor',
                     age_at_diagnosis=NA,
                     disease='CMN',
                     DiseaseState='Tumour',
                     Group = 'CMN'
                     )
manifest = rbind(manifest,cmnMani)
####################
# Add pRCC subtype
pap = read.table('~/Projects/KidneySC/papillary_subtypes.txt',sep='\t',header=TRUE)
pap$ID = gsub('\\..*','',pap$id_2)
manifest$subtype = pap$type[match(gsub('-01[AB]','',manifest$cases),pap$ID)]
#Modify group
extra = ifelse(is.na(manifest$subtype),
                        '',
                        ifelse(grepl('Type 1',manifest$subtype),
                               '_1',
                               ifelse(grepl('Type 2',manifest$subtype),
                                      '_2',
                                      '_?'
                                      )
                               )
                        )
manifest$Group = paste0(manifest$Group,extra)
#Save this data, so we can use it again if needed
saveRDS(list(manifest=manifest,cnts=cnts),file.path(plotDir,'bulkData.RDS'))
############################
#### CMN and RT inference
####Pick the samples to use and order them
###m = match(colnames(cnts),rownames(manifest))
###p = split(seq_along(m),manifest$Group[m])
####Pick 10 at random or the whole list
###set.seed(1)
###p = lapply(p,function(e) sample(e,min(65,length(e))))
####Order categories and drop un-necessary ones
###p = p[c('CMN','RT','WT','KICH','KIRC','KIRP','Normal')]
####Decide which genes to use
###g = finalMarkers[finalMarkers$cellType %in% c('UB','PV','F13','F','RCC','WT','RT','F6','F1','F10','F12','F3','F6','F9','CD-A','CD-B',grep('^[FN][0-9]+',finalMarkers$cellType,value=TRUE)),]
####Further refine and re-name
###g$cellType[g$cellType %in% c('F3')] = 'FVasc'
###g$cellType[g$cellType%in% c('F6','N21')] = 'F'
###g$cellType[g$cellType=='F13'] = 'Gang'
###g$cellType[g$cellType %in% c('N10','N15','N23')] = 'Ureter'
###g$cellType[g$cellType %in% c('N11')] = 'Glom'
###g$cellType[g$cellType %in% c('N12')] = 'PT3'
###g$cellType[g$cellType %in% c('N14')] = 'Dist'
###g$cellType[g$cellType %in% c('N18')] = 'Mess'
###g$cellType[g$cellType %in% c('N20')] = 'Pelvis'
###g$cellType[g$cellType %in% c('N8')] = 'LoH'
###g$cellType[g$cellType %in% c('N4','N9')] = 'Vasc'
###g$cellType[g$cellType %in% c('CD-A','CD-B','N16','N17')] = 'CD'
####g$cellType[g$cellType %in% c('RCC','WT','RT')] = 'Canon'
###m = match(g$ENS_ID,names(geneDat))
###ss = g$cellType
###w = !is.na(m)
###ss = ss[w]
###mat = log2(cnts[m[w],unlist(p)])
###rownames(mat) = g$gene[w]
####Set the -Inf to NA
###mat[!is.finite(mat)]=NA
####Now do the plotting
###cols = c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494')
###cols = setNames(cols,names(p))
###topAnn = HeatmapAnnotation(data.frame(Group=rep(names(p),lengths(p))),col=list(Group=cols))
###mat = t(scale(t(mat)))
####Take the average of the scaled values
###mat = do.call(rbind,lapply(split(seq(nrow(mat)),g$cellType[w]),function(e) apply(mat[e,,drop=FALSE],2,mean,na.rm=TRUE)))
###ss = rownames(mat)
###ss = rep(names(p),lengths(p))
####Heirarchical cluster within groups to order
####mat = mat[,unlist(lapply(split(seq(ncol(mat)),rep(names(p),lengths(p)))[names(p)],function(e) e[hclust(dist(t(mat[,e])))$order]))]
####What limimts on the heatmap
###hLims = max(abs(quantile(mat,c(0.01,.99),na.rm=TRUE)))
###cols = circlize::colorRamp2(c(-hLims,0,hLims),c('purple','black','yellow'))
###pdf(file.path(plotDir,'CMN_RT.pdf'),width=10,height=7)
###Heatmap(t(mat),
###        col = cols,
###        name='ScaledExpression',
###        #top_annotation=topAnn,
###        split=ss,
###        cluster_column=FALSE,
###        cluster_rows=FALSE,
###        show_column_names=TRUE,
###        show_row_names=FALSE
###        )
###dev.off()
#############
#Create data for plot
df = log2(cnts)
#Keep the ones we care about
m = match(finalMarkers$ENS_ID,names(geneDat))
df = as.data.frame(df[m,])
df$geneLab = rownames(df)
df$cellType = factor(finalMarkers$cellType,levels=unique(finalMarkers$cellType))
#Order as in finalMarkers
df$GeneSymbol = factor(geneDat$external_gene_name[m],levels=unique(geneDat$external_gene_name[m]))
df = df[!is.na(df$GeneSymbol),]
#Create table to write out
tab = df
m = match(colnames(tab),rownames(manifest))
colnames(tab) = gsub('CMN_','',gsub('NA_','',paste0(manifest$Group[m],'_',ifelse(manifest$disease[m]==manifest$Group[m],'NA',manifest$disease[m]),'_',colnames(tab))))
write.table(tab,file.path(plotDir,'ValidationMarkerTable.tsv'),sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE)
df = melt(df,id.vars=c('GeneSymbol','geneLab','cellType'),variable.name='SampleID',value.name='FPKM',stringsAsFactors=FALSE)
df = dplyr::inner_join(df,manifest,by=c(SampleID='submitter_id'))
#Do some renaming to make things consistent
df = df[,-match('disease',colnames(df))]
colnames(df) = gsub('Group','disease',colnames(df))
#Set ordering
df$disease = factor(df$disease,levels=c('Foetal','PaedNormal','Normal','CMN','RT','WT','KICH','KIRP','KIRC'))
#Exclude NAs and extreme outliers
df = df[!is.na(df$FPKM),]
df = df[abs(df$FPKM)<10,]
#Make the plot
gg = ggplot(df,aes(x=GeneSymbol,y=FPKM,colour=disease)) + 
  geom_boxplot(outlier.colour=NA,fill='white') +
  #geom_point(position=position_jitterdodge(),alpha=1/20,size=.1) +
  #Attempt to manually make the guides by sensible.
  #scale_colour_manual(breaks=c('TRUE','FALSE',unique(df$disease)),
  #                    values=seq_along(unique(df$disease)),
  #                    limits=c(unique(df$disease))) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(face='italic'))+
  ylab('log2(TPM)') +
  facet_wrap(~cellType+GeneSymbol,scales='free_x')
pdf(file.path(plotDir,'BulkValidationBoxplots.pdf'),width=28,height=28)
plot(gg)
dev.off()
##Heatmap by these genes
#mat = log2(cnts)
#w = which(names(geneDat) %in% finalMarkers$ENS_ID)
#mat = mat[w,]
#rownames(mat) = geneDat$external_gene_name[w]
#mat = as.matrix(mat)
#mat[!is.finite(mat)]=NA
##Group genes
#ss = sapply(split(finalMarkers$cellType,finalMarkers$ENS_ID)[names(geneDat)[w]],paste,collapse=',')
#cols = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6')
#cols = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd')
#pdf(file.path(plotDir,'BulkValidationHeatmap.pdf'),width=28,height=28)
#cols = setNames(cols,unique(manifest$Group))
#topAnn = HeatmapAnnotation(manifest[colnames(mat),'Group',drop=FALSE],col=list(Group=cols))
#Heatmap(t(scale(t(mat))),top_annotation=topAnn,split=ss,cluster_column=TRUE)
#dev.off()
#####################################################################
# Finally, do statistical tests for significance of key comparisons
# x = log2(cnts) of gene for group 1
# y = log2(cnts) of gene for group 2
#DE if wilcox.test(x,y,paired=FALSE)<0.01 and abs(median(x)-median(y))>1
#Enrichment of SLC17A3 in ccRCC and pRCC versus chrRCC
nom = grep('SLC17A3_',rownames(rtoc),value=TRUE)
nom = gsub('.*_','',nom)
x = log2(cnts[nom,manifest$Group=='KICH'])
tsts=list()
for(tgt in c('KIRP','KIRC')){
  y = log2(cnts[nom,manifest$Group==tgt])
  tsts[[tgt]] = wilcox.test(x,y,paired=FALSE)
  tsts[[tgt]]$shift = median(x)-median(y)
}
#Same for VCAM1
nom = grep('VCAM1_',rownames(rtoc),value=TRUE)
nom = gsub('.*_','',nom)
x = log2(cnts[nom,manifest$Group=='KICH'])
tsts=list()
for(tgt in c('KIRP','KIRC')){
  y = log2(cnts[nom,manifest$Group==tgt])
  tsts[[tgt]] = wilcox.test(x,y,paired=FALSE)
  tsts[[tgt]]$shift = median(x)-median(y)
}
#Depletion of SLC7A13 in RCC versus normal
nom = grep('SLC7A13_',rownames(rtoc),value=TRUE)
nom = gsub('.*_','',nom)
x = log2(cnts[nom,manifest$Group=='Normal'])
tsts=list()
for(tgt in c('KICH','KIRP','KIRC')){
  y = log2(cnts[nom,manifest$Group==tgt])
  tsts[[tgt]] = wilcox.test(x,y,paired=FALSE)
  tsts[[tgt]]$shift = median(x)-median(y)
}
#Test for depletion of a bunch of other nephron markers in ccRCC and pRCC
tsts=list()
noms = c('SLC12A1_','KCNJ1_','SLC4A1_','SLC26A4_')
for(nom in noms){
  nom = grep(nom,rownames(rtoc),value=TRUE)
  nom = gsub('.*_','',nom)
  x = log2(cnts[nom,manifest$Group=='Normal'])
  for(tgt in c('KIRP','KIRC')){
    y = log2(cnts[nom,manifest$Group==tgt])
    tsts[[paste0(nom,'_',tgt)]] = wilcox.test(x,y,paired=FALSE)
    tsts[[paste0(nom,'_',tgt)]]$shift = median(x)-median(y)
  }
  ##Make boxplots
  #w = manifest$Group %in% c('KIRC','KIRP','KICH','Normal')
  #boxplot(split(log2(cnts[nom,w]),manifest$Group[w]))
}

#######################
# wilmsBulkValidation #
#######################
bulk = readRDS(file.path(plotDir,'bulkData.RDS'))
mat = log2(bulk$cnts)
gns = c('WT1','CA9','SMARCB1','NMU','RSPO1')
gns = rownames(rtoc)[match(gns,gsub('_.*','',rownames(rtoc)))]
df = mat[gsub('.*_','',gns),]
df = melt(df)
df[,1] = as.character(df[,1])
df[,2] = as.character(df[,2])
colnames(df) = c('Gene','Sample','TPM')
m = match(df$Sample,rownames(bulk$manifest))
df$Group = manifest$Group[m]
#Fix names and set ordering
nomMap = c('RT'='MRT','CMN'='CMN','WT'='Wilms','Foetal'='Fetal','PaedNormal'='Paed','Normal'='Adult','KIRC'='ccRCC','KIRP'='pRCC','KICH'='chrRCC')
df$Group=factor(nomMap[df$Group],levels=as.character(nomMap))
df$Symbol = gsub('_.*','',gns[match(df$Gene,gsub('.*_','',gns))])
df$Symbol = factor(df$Symbol,levels=gsub('_.*','',gns))
#Make the plots
gg = ggplot(df,aes(x=Group,y=TPM)) + 
  stat_boxplot(geom='errorbar') +
  geom_boxplot(outlier.colour=NA,fill='white',width=0.5) +
  #geom_point(position=position_jitterdodge(),alpha=1/20,size=.1) +
  #Attempt to manually make the guides by sensible.
  #scale_colour_manual(breaks=c('TRUE','FALSE',unique(df$disease)),
  #                    values=seq_along(unique(df$disease)),
  #                    limits=c(unique(df$disease))) +
  #theme(axis.title.x=element_blank(),
  #      axis.text.x=element_blank(),
  #      axis.ticks.x=element_blank(),
  #      strip.text = element_text(face='italic'))+
  ylab('log2(TPM)') +
  facet_grid(~Symbol,scales='free_x')
pdf(file.path(plotDir,'wilmsBulkValidation.pdf'),width=21,height=7)
plot(gg)
dev.off()

###########################################
# rccBulkValidation and rccSingleCellFrac #
###########################################
bulk = readRDS(file.path(plotDir,'bulkData.RDS'))
mat = log2(bulk$cnts)
gns = c('SLC17A3','SLC7A13','VCAM1','SLC12A1','KCNJ1','SLC4A1','SLC26A4')
gns = rownames(rtoc)[match(gns,gsub('_.*','',rownames(rtoc)))]
df = mat[gsub('.*_','',gns),]
df = melt(df)
df[,1] = as.character(df[,1])
df[,2] = as.character(df[,2])
colnames(df) = c('Gene','Sample','TPM')
m = match(df$Sample,rownames(bulk$manifest))
df$Group = manifest$Group[m]
#Exclude everything but the normal and adult tumours
#df = df[df$Group %in% c('KIRC','KIRP','KICH','Normal'),]
df = df[df$Group %in% c('KIRC','KIRP_1','KIRP_2','KIRP_?','KICH','Normal'),]
#Fix names and set ordering
#nomMap = c('Normal'='Adult','KIRC'='ccRCC','KIRP'='pRCC','KICH'='chrRCC')
nomMap = c('Normal'='Adult','KIRC'='ccRCC','KIRP_?'='pRCC','KIRP_1'='pRCC1','KIRP_2'='pRCC2','KICH'='chrRCC')
df$Group=factor(nomMap[df$Group],levels=as.character(nomMap))
df$Symbol = gsub('_.*','',gns[match(df$Gene,gsub('.*_','',gns))])
df$Symbol = factor(df$Symbol,levels=gsub('_.*','',gns))
#Make the plots
gg = ggplot(df,aes(x=Group,y=TPM)) + 
  stat_boxplot(geom='errorbar') +
  geom_boxplot(outlier.colour=NA,fill='white',width=0.5) +
  #geom_point(position=position_jitterdodge(),alpha=1/20,size=.1) +
  #Attempt to manually make the guides by sensible.
  #scale_colour_manual(breaks=c('TRUE','FALSE',unique(df$disease)),
  #                    values=seq_along(unique(df$disease)),
  #                    limits=c(unique(df$disease))) +
  #theme(axis.title.x=element_blank(),
  #      axis.text.x=element_blank(),
  #      axis.ticks.x=element_blank(),
  #      strip.text = element_text(face='italic'))+
  ylab('log2(TPM)') +
  facet_grid(~Symbol,scales='free_x')
pdf(file.path(plotDir,'rccBulkValidation.pdf'),width=21,height=7)
plot(gg)
dev.off()
####################################################
# Now add the single cell expression fraction part
df = geneFrac[gns,]
df = melt(df)
df[,1] = as.character(df[,1])
df[,2] = as.character(df[,2])
colnames(df) = c('Gene','Cluster','Frac')
#Get the relevant clusters and map them to groups
df = df[grepl('^[NT]',df$Cluster),]
#Thing to map the necessary clusters.  Others are dropped
clMap = c(T15='pR',
          T9='cR3',
          T4='cR2',
          T10='cR6',
          T12='cR7',
          T7='cR5',
          T6='cR4',
          T17='cR1',
          N13='PT1',
          N12='PT3',
          N22='PT2',
          N17='C2',
          N16='C1',
          N0='H',
          N14='D')
df = df[df$Cluster %in% names(clMap),]
df$Cluster = clMap[df$Cluster]
#Split into normal, cRCC and pRCC
df$Compartment = ifelse(grepl('^pR',df$Cluster),
                        'pRCC',
                        ifelse(grepl('^cR',df$Cluster),
                               'ccRCC',
                               'Normal'
                               )
                        )
#Fix some names
df$Symbol = gsub('_.*','',df$Gene)
#Set ordering
df$Symbol = factor(df$Symbol,levels=gsub('_.*','',gns))
df$Compartment = factor(df$Compartment,levels=c('Normal','pRCC','ccRCC'))
#First do the proximal tubul
gg = ggplot(df,aes(Compartment,Frac)) +
  geom_label(aes(label=Cluster)) +
  ylim(0,1) + 
  facet_grid(~Symbol,scale='free_x')
pdf(file.path(plotDir,'rccSingleCellFrac.pdf'),width=21,height=7)
plot(gg)
dev.off()

######################
# vascBulkValidation #
######################
bulk = readRDS(file.path(plotDir,'bulkData.RDS'))
mat = log2(bulk$cnts)
gns = c('PDGFRB','ACTA2','SEMA3G','SLC14A1','PLVAP','EMILIN1','SFRP2','MMP2')
gns = rownames(rtoc)[match(gns,gsub('_.*','',rownames(rtoc)))]
gnGroups = c('Myofibroblasts','Myofibroblasts','Endothelium','Endothelium','Endothelium','Fibroblast','Fibroblast','Fibroblast')
df = mat[gsub('.*_','',gns),]
df = melt(df)
df[,1] = as.character(df[,1])
df[,2] = as.character(df[,2])
colnames(df) = c('Gene','Sample','TPM')
m = match(df$Sample,rownames(bulk$manifest))
df$Group = manifest$Group[m]
#Exclude everything but the normal and sampled tumours tumours
df = df[df$Group %in% c('KIRC','KIRP','WT','Normal'),]
#Fix names and set ordering
nomMap = c('Normal'='Normal','WT'='Wilms','KIRC'='ccRCC','KIRP'='pRCC')
df$Group=factor(nomMap[df$Group],levels=as.character(nomMap))
df$Symbol = gsub('_.*','',gns[match(df$Gene,gsub('.*_','',gns))])
#df$Symbol = factor(df$Symbol,levels=gsub('_.*','',gns))
df$CellType = factor(gnGroups[match(df$Symbol,gsub('_.*','',gns))],levels=unique(gnGroups))
#Split into three plots, one for each cell type
ggs = list()
for(cellType in levels(df$CellType)){
  sdf = df[df$CellType==cellType,]
  sdf$Symbol = factor(sdf$Symbol,levels=gsub('_.*','',gns[gsub('_.*','',gns) %in% sdf$Symbol]))
  #Make the plots
  gg = ggplot(sdf,aes(x=Group,y=TPM)) + 
    stat_boxplot(geom='errorbar') +
    geom_boxplot(outlier.colour=NA,fill='white',width=0.5) +
    #geom_point(position=position_jitterdodge(),alpha=1/20,size=.1) +
    #Attempt to manually make the guides by sensible.
    #scale_colour_manual(breaks=c('TRUE','FALSE',unique(df$disease)),
    #                    values=seq_along(unique(df$disease)),
    #                    limits=c(unique(df$disease))) +
    #theme(axis.title.x=element_blank(),
    #      axis.text.x=element_blank(),
    #      axis.ticks.x=element_blank(),
    #      strip.text = element_text(face='italic'))+
    ylab('log2(TPM)') +
    xlab('')+
    facet_grid(CellType~Symbol,scales='free_x')
  ggs[[cellType]] = gg
}
#Make into a sensible layout and print
pdf(file.path(plotDir,'vascBulkValidation.pdf'),width=10,height=14)
grid.arrange(grobs = list(ggs[[1]],ggs[[2]],ggs[[3]]),widths=c(2,1),layout_matrix=rbind(c(1,NA),c(2,2),c(3,3)))
dev.off()





####################
# ptBulkValidation #
####################
bulk = readRDS(file.path(plotDir,'bulkData.RDS'))
mat = log2(bulk$cnts)
gns = c("SLC17A3_ENSG00000124564", "VCAM1_ENSG00000162692", "SLC7A13_ENSG00000164893",
"NMU_ENSG00000109255", "RSPO1_ENSG00000169218")
ss = c('PT1','PT1','PT1','UB','PV')
mat = mat[match(gsub('.*_','',gns),rownames(mat)),]
w = bulk$manifest$Group %in% c('KIRC','KIRP','WT','Normal')
mat = mat[,w]
mat[!is.finite(mat)]=NA
mat = t(mat)
#Truncate scale (or don't) and fill in
#lims = quantile(mat,c(0.01,0.99),na.rm=TRUE)
#mat[is.na(mat) | mat< lims[1]] = lims[1]
#mat[mat> lims[2]] =lims[2]
mat[!is.finite(mat)] = min(mat,na.rm=TRUE)
#mat = scale(mat)
top = bulk$manifest$Group[w]
top[top=='KIRC']='ccRCC'
top[top=='KIRP']='pRCC'
top[top=='WT']='Wilms'
colnames(mat) = gsub('_.*','',gns)
pdf(file.path(plotDir,'ptBulkValidation.pdf'),width=14,height=14)
draw(Heatmap(mat,
        col = circlize::colorRamp2(c(-5,0,10),c('purple','black','yellow')),
        name = 'log2(TPM)',
        split=top,
        show_row_names=FALSE,
        cluster_column=FALSE))
dev.off()



###############
# ageAnalysis #
###############
#Split data into tissues
neph = c('F1','F4','F7','F8','F10','F11',paste0('N',c(0,8,11,12,13,14,22)),paste0('T',c(4,5,6,7,9,10,12,15,16,17,18)))
subsets = list(neph=neph)
#Split samples into age bins
tmp = unique(mani$label)
m = match(tmp,mani$label)
labelAgeMap = paste0(ifelse(mani$TissueDiseaseState[m]!='Normal','Tumour','Normal'),'_',
                     ifelse(mani$Age_donor.MonthsPostConception.[m] < 9,
                            'Foetal',
                            ifelse(mani$Age_donor.MonthsPostConception.[m] < 100,
                                   'Child',
                                   ifelse(mani$Age_donor.MonthsPostConception.[m] < 225,
                                          'Adolescent',
                                          'Adult'
                                          )
                                   )
                            )
                     )
names(labelAgeMap) = tmp
#Want to basically create a geneFrac matrix where we exclude children
geneFracSinKids = list()
clustCountsSinKids = c()
for(i in seq_len(nrow(compNames))){
  #Skip the uninformative ones 
  if(compNames$prettyName[i] %in% c('Indistinct','Not Kidney','Unassigned'))
    next
  print(compNames$prettyName[i])
  mDat = readRDS(file.path(baseDir,compNames$dirtyName[i],'batchCorrected','metadata.RDS'))
  rownames(mDat$meta.data) = paste0(mDat$meta.data$Sanger_study_ID,'___',gsub('^[0-9]+_','',rownames(mDat$meta.data)))
  #Work out which ones to exclude.  Always exclude Adolescent.
  blackList = rownames(mDat$meta.data)[grepl('Adolescent',labelAgeMap[mDat$meta.data$label]) | grepl('Child',labelAgeMap[mDat$meta.data$label])]
  tmp = split(rownames(mDat$meta.data),paste0(compNames$clusterName[i],mDat$meta.data$res.1))
  geneFracSinKids[[compNames$clusterName[i]]] = do.call(cbind,lapply(tmp,function(e) Matrix::rowSums(rtoc[,e[!(e%in%blackList)],drop=FALSE]>0)/length(e[!(e%in%blackList)])))
  clustCountsSinKids = c(clustCountsSinKids,sapply(tmp,function(e) sum(!(e%in%blackList))))
}
geneFracSinKids = do.call(cbind,geneFracSinKids)
#And another without adults
geneFracSinAdults = list()
clustCountsSinAdults = c()
for(i in seq_len(nrow(compNames))){
  #Skip the uninformative ones 
  if(compNames$prettyName[i] %in% c('Indistinct','Not Kidney','Unassigned'))
    next
  print(compNames$prettyName[i])
  mDat = readRDS(file.path(baseDir,compNames$dirtyName[i],'batchCorrected','metadata.RDS'))
  rownames(mDat$meta.data) = paste0(mDat$meta.data$Sanger_study_ID,'___',gsub('^[0-9]+_','',rownames(mDat$meta.data)))
  #Work out which ones to exclude.  Always exclude Adolescent.
  blackList = rownames(mDat$meta.data)[grepl('Adolescent',labelAgeMap[mDat$meta.data$label]) | grepl('Adult',labelAgeMap[mDat$meta.data$label])]
  tmp = split(rownames(mDat$meta.data),paste0(compNames$clusterName[i],mDat$meta.data$res.1))
  geneFracSinAdults[[compNames$clusterName[i]]] = do.call(cbind,lapply(tmp,function(e) Matrix::rowSums(rtoc[,e[!(e%in%blackList)],drop=FALSE]>0)/length(e[!(e%in%blackList)])))
  clustCountsSinAdults = c(clustCountsSinAdults,sapply(tmp,function(e) sum(!(e%in%blackList))))
}
geneFracSinAdults = do.call(cbind,geneFracSinAdults)
#Make a merged count
tmp = unique(names(clustCountsSinKids),names(clustCountsSinAdults))
clustCounts = rep(0,length(tmp))
names(clustCounts) = tmp
clustCounts[names(clustCountsSinKids)] = clustCountsSinKids
#Now for the non-foetal ones, add in the childrens counts too
tmp = clustCountsSinAdults
tmp = tmp[grep('F',names(tmp),invert=TRUE)]
clustCounts[names(tmp)] = clustCounts[names(tmp)]+tmp
#Do the comparisons of the nephron versus other stages
inOverOut = 0.2
noiseRate = 0.1
pCut = 0.05
#Bake the bulk validation into the selection process
bulkMin = 1
bulkRat = 1
#Prepare the bulk data
bulk = readRDS(file.path(plotDir,'bulkData.RDS'))
fExp = rowMeans(bulk$cnts[,bulk$manifest$Group=='Foetal'])
cExp = rowMeans(bulk$cnts[,bulk$manifest$Group=='PaedNormal'])
aExp = rowMeans(bulk$cnts[,bulk$manifest$Group=='Normal'])
nExp = rowMeans(bulk$cnts[,grep('Normal',bulk$manifest$Group)])
cfExp = rowMeans(bulk$cnts[,bulk$manifest$Group%in%c('PaedNormal','Foetal')])
afExp = rowMeans(bulk$cnts[,bulk$manifest$Group%in%c('Normal','Foetal')])
markers = list()
for(i in seq_along(subsets)){
  nom = names(subsets)[i]
  print(nom)
  #This is easier if we just do it one at a time
  ##################
  # High in Foetus #
  ##################
  #Get the high fraction
  wHigh  = grepl('F',subsets[[nom]])
  fClust = apply(geneFrac[,subsets[[nom]][wHigh],drop=FALSE],1,which.max)
  fClust = subsets[[nom]][wHigh][fClust]
  fMax = apply(geneFrac[,subsets[[nom]][wHigh],drop=FALSE],1,max)
  ##################
  # Foetus v Child
  #Get the low fraction
  #Exclude adults from consideration
  lowDat = geneFracSinAdults
  cCnts = clustCountsSinAdults
  #Include the other foetal clusters
  oog = grep('F',colnames(lowDat),value=TRUE)
  oog = oog[!(oog %in% subsets[[nom]][wHigh])]
  #And everything normal
  oog = c(oog,grep('^I?[FT]',colnames(lowDat),value=TRUE,invert=TRUE))
  #Get the data
  pVals = pbinom(t(t(lowDat[,oog,drop=FALSE])*cCnts[oog])-1,outer(rep(1,nrow(lowDat)),cCnts[oog]),noiseRate,lower.tail=FALSE)
  w = apply(pVals,1,which.min)
  pVals = pVals[cbind(seq(nrow(pVals)),w)]
  oMax = (lowDat[,oog,drop=FALSE])[cbind(seq(nrow(lowDat)),w)]
  #Work out which ones to keep
  df = data.frame(fMax=fMax,
                  oMax=oMax,
                  oPval=pVals,
                  fCnt = cCnts[fClust],
                  oCnt=cCnts[oog[w]],
                  fClust = fClust,
                  oClust = oog[w],
                  bulkHigh = fExp[gsub('.*_','',names(fMax))],
                  bulkLow = cExp[gsub('.*_','',names(fMax))])
  df = df[order(df$fMax-df$oMax,decreasing=TRUE),]
  markers[[paste0(nom,'_FvC')]] = df
  ##################
  # Foetus v Adult
  #Exclude adults from consideration
  lowDat = geneFracSinKids
  cCnts = clustCountsSinKids
  #Include the other foetal clusters
  oog = grep('F',colnames(lowDat),value=TRUE)
  oog = oog[!(oog %in% subsets[[nom]][wHigh])]
  #And everything normal
  oog = c(oog,grep('^I?[FT]',colnames(lowDat),value=TRUE,invert=TRUE))
  #Get the data
  pVals = pbinom(t(t(lowDat[,oog,drop=FALSE])*cCnts[oog])-1,outer(rep(1,nrow(lowDat)),cCnts[oog]),noiseRate,lower.tail=FALSE)
  w = apply(pVals,1,which.min)
  pVals = pVals[cbind(seq(nrow(pVals)),w)]
  oMax = (lowDat[,oog,drop=FALSE])[cbind(seq(nrow(lowDat)),w)]
  #Work out which ones to keep
  df = data.frame(fMax=fMax,
                  oMax=oMax,
                  oPval=pVals,
                  fCnt = cCnts[fClust],
                  oCnt=cCnts[oog[w]],
                  fClust = fClust,
                  oClust = oog[w],
                  bulkHigh = fExp[gsub('.*_','',names(fMax))],
                  bulkLow = aExp[gsub('.*_','',names(fMax))]
                  )
  df = df[order(df$fMax-df$oMax,decreasing=TRUE),]
  markers[[paste0(nom,'_FvA')]] = df
  ######################
  # Foetus v Non-Foetus
  #And the low fraction
  t1 = geneFracSinKids
  colnames(t1) = paste0('Adult_',colnames(t1))
  t2 = geneFracSinAdults
  colnames(t2) = paste0('Child_',colnames(t2))
  cCnts = c(setNames(clustCountsSinKids,colnames(t1)),setNames(clustCountsSinAdults,colnames(t2)))
  lowDat = cbind(t1,t2)
  #Include the other foetal clusters
  oog = grep('F',colnames(lowDat),value=TRUE)
  oog = oog[!(gsub('(Child|Adult)_','',oog) %in% subsets[[nom]][wHigh])]
  #And everything normal
  oog = c(oog,grep('_I?[FT]',colnames(lowDat),value=TRUE,invert=TRUE))
  #Get the data
  pVals = pbinom(t(t(lowDat[,oog,drop=FALSE])*cCnts[oog])-1,outer(rep(1,nrow(lowDat)),cCnts[oog]),noiseRate,lower.tail=FALSE)
  w = apply(pVals,1,which.min)
  pVals = pVals[cbind(seq(nrow(pVals)),w)]
  oMax = (lowDat[,oog,drop=FALSE])[cbind(seq(nrow(lowDat)),w)]
  #Work out which ones to keep
  df = data.frame(fMax=fMax,
                  oMax=oMax,
                  oPval=pVals,
                  fCnt = cCnts[fClust],
                  oCnt=cCnts[oog[w]],
                  fClust = fClust,
                  oClust = oog[w],
                  bulkHigh = fExp[gsub('.*_','',names(fMax))],
                  bulkLow = nExp[gsub('.*_','',names(fMax))]
                  )
  df = df[order(df$fMax-df$oMax,decreasing=TRUE),]
  markers[[paste0(nom,'_FvAC')]] = df
  #################
  # High in Child #
  #################
  #Get the high fraction
  wHigh  = grepl('N',subsets[[nom]])
  fClust = apply(geneFracSinAdults[,subsets[[nom]][wHigh],drop=FALSE],1,which.max)
  fClust = subsets[[nom]][wHigh][fClust]
  fMax = apply(geneFracSinAdults[,subsets[[nom]][wHigh],drop=FALSE],1,max)
  ##################
  # Child v Foetus
  #And the low fraction
  #Need to include the adult from the target clusters as well, so get that first
  lowDat = geneFracSinAdults
  cCnts = clustCountsSinAdults
  #Include the other normal clusters
  oog = grep('^N',colnames(lowDat),value=TRUE)
  oog = oog[!(oog %in% subsets[[nom]][wHigh])]
  #Include immune normal and the two foetal groups
  oog = c(oog,grep('IN',colnames(lowDat),value=TRUE),grep('F',colnames(lowDat),value=TRUE))
  #Get the data
  pVals = pbinom(t(t(lowDat[,oog,drop=FALSE])*cCnts[oog])-1,outer(rep(1,nrow(lowDat)),cCnts[oog]),noiseRate,lower.tail=FALSE)
  w = apply(pVals,1,which.min)
  pVals = pVals[cbind(seq(nrow(pVals)),w)]
  oMax = (lowDat[,oog,drop=FALSE])[cbind(seq(nrow(lowDat)),w)]
  #Work out which ones to keep
  df = data.frame(fMax=fMax,
                  oMax=oMax,
                  oPval=pVals,
                  fCnt = cCnts[fClust],
                  oCnt=cCnts[oog[w]],
                  fClust = fClust,
                  oClust = oog[w],
                  bulkHigh = cExp[gsub('.*_','',names(fMax))],
                  bulkLow = fExp[gsub('.*_','',names(fMax))]
                  )
  df = df[order(df$fMax-df$oMax,decreasing=TRUE),]
  markers[[paste0(nom,'_CvF')]] = df
  ##################
  # Child v Adult
  #Exclude adults from consideration
  t1 = geneFracSinKids
  colnames(t1) = paste0('Adult_',colnames(t1))
  t2 = geneFracSinAdults
  colnames(t2) = paste0('Child_',colnames(t2))
  cCnts = c(setNames(clustCountsSinKids,colnames(t1)),setNames(clustCountsSinAdults,colnames(t2)))
  lowDat = cbind(t1,t2)
  #Include all the other clusters from child
  oog = grep('Child_N',colnames(lowDat),value=TRUE)
  oog = oog[!(oog %in% paste0('Child_',subsets[[nom]][wHigh]))]
  #And all the normal clusters from adult
  oog = c(oog,grep('Adult_N',colnames(lowDat),value=TRUE))
  #And all the normal immune
  oog = c(oog,grep('_IN',colnames(lowDat),value=TRUE))
  #Get the data
  pVals = pbinom(t(t(lowDat[,oog,drop=FALSE])*cCnts[oog])-1,outer(rep(1,nrow(lowDat)),cCnts[oog]),noiseRate,lower.tail=FALSE)
  w = apply(pVals,1,which.min)
  pVals = pVals[cbind(seq(nrow(pVals)),w)]
  oMax = (lowDat[,oog,drop=FALSE])[cbind(seq(nrow(lowDat)),w)]
  #Work out which ones to keep
  df = data.frame(fMax=fMax,
                  oMax=oMax,
                  oPval=pVals,
                  fCnt = cCnts[fClust],
                  oCnt=cCnts[oog[w]],
                  fClust = fClust,
                  oClust = oog[w],
                  bulkHigh = cExp[gsub('.*_','',names(fMax))],
                  bulkLow = aExp[gsub('.*_','',names(fMax))]
                  )
  df = df[order(df$fMax-df$oMax,decreasing=TRUE),]
  markers[[paste0(nom,'_CvA')]] = df
  ######################
  # Child v Adult+Foetus
  #Exclude adults from consideration
  t1 = geneFracSinKids
  colnames(t1) = paste0('Adult_',colnames(t1))
  t2 = geneFracSinAdults
  colnames(t2) = paste0('Child_',colnames(t2))
  cCnts = c(setNames(clustCountsSinKids,colnames(t1)),setNames(clustCountsSinAdults,colnames(t2)))
  lowDat = cbind(t1,t2)
  #Include all the other clusters from child
  oog = grep('Child_N',colnames(lowDat),value=TRUE)
  oog = oog[!(oog %in% paste0('Child_',subsets[[nom]][wHigh]))]
  #And all the normal clusters from adult
  oog = c(oog,grep('Adult_N',colnames(lowDat),value=TRUE))
  #And all the normal immune
  oog = c(oog,grep('_IN',colnames(lowDat),value=TRUE))
  #And all the foetal samples
  oog = c(oog,grep('F',colnames(lowDat),value=TRUE))
  #Get the data
  pVals = pbinom(t(t(lowDat[,oog,drop=FALSE])*cCnts[oog])-1,outer(rep(1,nrow(lowDat)),cCnts[oog]),noiseRate,lower.tail=FALSE)
  w = apply(pVals,1,which.min)
  pVals = pVals[cbind(seq(nrow(pVals)),w)]
  oMax = (lowDat[,oog,drop=FALSE])[cbind(seq(nrow(lowDat)),w)]
  #Work out which ones to keep
  df = data.frame(fMax=fMax,
                  oMax=oMax,
                  oPval=pVals,
                  fCnt = cCnts[fClust],
                  oCnt=cCnts[oog[w]],
                  fClust = fClust,
                  oClust = oog[w],
                  bulkHigh = cExp[gsub('.*_','',names(fMax))],
                  bulkLow = afExp[gsub('.*_','',names(fMax))]
                  )
  df = df[order(df$fMax-df$oMax,decreasing=TRUE),]
  markers[[paste0(nom,'_CvFA')]] = df
  #################
  # High in Adult #
  #################
  #Get the high fraction
  wHigh  = grepl('N',subsets[[nom]])
  fClust = apply(geneFracSinKids[,subsets[[nom]][wHigh],drop=FALSE],1,which.max)
  fClust = subsets[[nom]][wHigh][fClust]
  fMax = apply(geneFracSinKids[,subsets[[nom]][wHigh],drop=FALSE],1,max)
  ##################
  # Adult v Foetus
  #And the low fraction
  #Need to include the adult from the target clusters as well, so get that first
  lowDat = geneFracSinKids
  cCnts = clustCountsSinKids
  #Include the other normal clusters
  oog = grep('^N',colnames(lowDat),value=TRUE)
  oog = oog[!(oog %in% subsets[[nom]][wHigh])]
  #Include immune normal and the two foetal groups
  oog = c(oog,grep('IN',colnames(lowDat),value=TRUE),grep('F',colnames(lowDat),value=TRUE))
  #Get the data
  pVals = pbinom(t(t(lowDat[,oog,drop=FALSE])*cCnts[oog])-1,outer(rep(1,nrow(lowDat)),cCnts[oog]),noiseRate,lower.tail=FALSE)
  w = apply(pVals,1,which.min)
  pVals = pVals[cbind(seq(nrow(pVals)),w)]
  oMax = (lowDat[,oog,drop=FALSE])[cbind(seq(nrow(lowDat)),w)]
  #Work out which ones to keep
  df = data.frame(fMax=fMax,
                  oMax=oMax,
                  oPval=pVals,
                  fCnt = cCnts[fClust],
                  oCnt=cCnts[oog[w]],
                  fClust = fClust,
                  oClust = oog[w],
                  bulkHigh = aExp[gsub('.*_','',names(fMax))],
                  bulkLow = fExp[gsub('.*_','',names(fMax))]
                  )
  df = df[order(df$fMax-df$oMax,decreasing=TRUE),]
  markers[[paste0(nom,'_AvF')]] = df
  ##################
  # Adult v Child
  t1 = geneFracSinKids
  colnames(t1) = paste0('Adult_',colnames(t1))
  t2 = geneFracSinAdults
  colnames(t2) = paste0('Child_',colnames(t2))
  cCnts = c(setNames(clustCountsSinKids,colnames(t1)),setNames(clustCountsSinAdults,colnames(t2)))
  lowDat = cbind(t1,t2)
  #Include all the other clusters from Adult
  oog = grep('Adult_N',colnames(lowDat),value=TRUE)
  oog = oog[!(oog %in% paste0('Adult_',subsets[[nom]][wHigh]))]
  #And all the normal clusters from adult
  oog = c(oog,grep('Child_N',colnames(lowDat),value=TRUE))
  #And all the normal immune
  oog = c(oog,grep('_IN',colnames(lowDat),value=TRUE))
  #Get the data
  pVals = pbinom(t(t(lowDat[,oog,drop=FALSE])*cCnts[oog])-1,outer(rep(1,nrow(lowDat)),cCnts[oog]),noiseRate,lower.tail=FALSE)
  w = apply(pVals,1,which.min)
  pVals = pVals[cbind(seq(nrow(pVals)),w)]
  oMax = (lowDat[,oog,drop=FALSE])[cbind(seq(nrow(lowDat)),w)]
  #Work out which ones to keep
  df = data.frame(fMax=fMax,
                  oMax=oMax,
                  oPval=pVals,
                  fCnt = cCnts[fClust],
                  oCnt=cCnts[oog[w]],
                  fClust = fClust,
                  oClust = oog[w],
                  bulkHigh = aExp[gsub('.*_','',names(fMax))],
                  bulkLow = cExp[gsub('.*_','',names(fMax))]
                  )
  df = df[order(df$fMax-df$oMax,decreasing=TRUE),]
  markers[[paste0(nom,'_AvC')]] = df
  ######################
  # Adult v Child+Foetus
  t1 = geneFracSinKids
  colnames(t1) = paste0('Adult_',colnames(t1))
  t2 = geneFracSinAdults
  colnames(t2) = paste0('Child_',colnames(t2))
  cCnts = c(setNames(clustCountsSinKids,colnames(t1)),setNames(clustCountsSinAdults,colnames(t2)))
  lowDat = cbind(t1,t2)
  #Include all the other clusters from Adult
  oog = grep('Adult_N',colnames(lowDat),value=TRUE)
  oog = oog[!(oog %in% paste0('Adult_',subsets[[nom]][wHigh]))]
  #And all the normal clusters from adult
  oog = c(oog,grep('Child_N',colnames(lowDat),value=TRUE))
  #And all the normal immune
  oog = c(oog,grep('_IN',colnames(lowDat),value=TRUE))
  #And all the foetal samples
  oog = c(oog,grep('F',colnames(lowDat),value=TRUE))
  #Get the data
  pVals = pbinom(t(t(lowDat[,oog,drop=FALSE])*cCnts[oog])-1,outer(rep(1,nrow(lowDat)),cCnts[oog]),noiseRate,lower.tail=FALSE)
  w = apply(pVals,1,which.min)
  pVals = pVals[cbind(seq(nrow(pVals)),w)]
  oMax = (lowDat[,oog,drop=FALSE])[cbind(seq(nrow(lowDat)),w)]
  #Work out which ones to keep
  df = data.frame(fMax=fMax,
                  oMax=oMax,
                  oPval=pVals,
                  fCnt = cCnts[fClust],
                  oCnt=cCnts[oog[w]],
                  fClust = fClust,
                  oClust = oog[w],
                  bulkHigh = aExp[gsub('.*_','',names(fMax))],
                  bulkLow = cfExp[gsub('.*_','',names(fMax))]
                  )
  df = df[order(df$fMax-df$oMax,decreasing=TRUE),]
  markers[[paste0(nom,'_AvFC')]] = df
}
#Now extract the significant genes
ageGenes = lapply(markers,function(e) e[e$fMax > e$oMax + inOverOut & e$oPval>pCut & !is.na(e$bulkHigh) & log2(e$bulkHigh)>bulkMin & !is.na(e$bulkLow) & log2(e$bulkHigh/e$bulkLow) > bulkRat,])
#Keep just Foetus v Normal and Adult v Foetus
ageGenes = ageGenes[c('neph_FvAC','neph_AvF')]
ss = rep(names(ageGenes),sapply(ageGenes,nrow))
genesToPlot = unlist(lapply(ageGenes,rownames))
w = grep('RP[0-9]',genesToPlot,invert=TRUE)
genesToPlot = genesToPlot[w]
ss = ss[w]
#Now get the expression matrix in each group
expAges = list()
clCnts = list()
for(lab in unique(gsub('[0-9]+','',unlist(subsets)))){
  print(lab)
  #Load the meta-data
  i = match(lab,compNames$clusterName)
  mDat = readRDS(file.path(baseDir,compNames$dirtyName[i],'batchCorrected','metadata.RDS'))
  rownames(mDat$meta.data) = paste0(mDat$meta.data$Sanger_study_ID,'___',gsub('^[0-9]+_','',rownames(mDat$meta.data)))
  mDat = mDat$meta.data
  #Create a new label/age grouping
  mDat$ageGroup = paste0(lab,mDat$res.1,'_',labelAgeMap[mDat$label])
  #Process one subset at a time
  for(nom in names(subsets)){
    print(nom)
    clusts = subsets[[nom]][gsub('[0-9]+','',subsets[[nom]])==lab]
    #Get the cells in this cluster and split them by age
    w = which(paste0(lab,mDat$res.1) %in% clusts)
    cells = split(rownames(mDat)[w],mDat$ageGroup[w])
    if(nom %in% names(expAges)){
      expAges[[nom]] = cbind(expAges[[nom]],do.call(cbind,lapply(cells,function(e) Matrix::rowSums(rtoc[,e,drop=FALSE]>0)/length(e))))
      clCnts[[nom]] = c(clCnts[[nom]],table(mDat$ageGroup[w]))
    }else{
      expAges[[nom]] = do.call(cbind,lapply(cells,function(e) Matrix::rowSums(rtoc[,e,drop=FALSE]>0)/length(e)))
      clCnts[[nom]] = table(mDat$ageGroup[w])
    }
  }
}
#Build the expression plots for the target genes
pdf(file.path(plotDir,'ageAnalysis.pdf'),width=14,height=28)
for(nom in names(subsets)){
  tmp = expAges[[nom]]
  tmp = tmp[genesToPlot[grep(nom,ss)],,drop=FALSE]
  tmp = tmp[,gsub('_.*','',colnames(tmp)) %in% subsets[[nom]]]
  #Keep the highly expressed in foetus one and order sensibly
  ord = strsplit(colnames(tmp),'_')
  #Define ordering of the three levels
  ord = lapply(ord,function(e) c(factor(e[2],levels=c('Normal','Tumour')),
                                 e[1],
                                 factor(e[3],levels=c('Foetal','Child','Adolescent','Adult'))))
  ord = order(sapply(ord,`[`,1),sapply(ord,`[`,2),sapply(ord,`[`,3))
  tmp = tmp[,ord]
  #Fix up names
  rownames(tmp) = gsub('_.*','',rownames(tmp))
  #Exclude adolescent group
  tmp = tmp[,grep('Adolescent',colnames(tmp),invert=TRUE)]
  #And Tumours, don't really need them
  tmp = tmp[,grep('^T',colnames(tmp),invert=TRUE)]
  #tmp = tmp[,grep('^F',colnames(tmp),invert=TRUE)]
  #Split columns and re-order
  cSplit = gsub('.*[0-9]_','',colnames(tmp))
  ord = lapply(split(colnames(tmp),factor(cSplit,levels=unique(cSplit))),function(e) if(length(e)==1){e}else{e[hclust(dist(t(tmp[,e,drop=FALSE])))$order]})
  tmp = tmp[,unlist(ord,use.names=FALSE),drop=FALSE]
  #Scale
  tmp = t(scale(t(tmp)))
  tmp[!is.finite(tmp)]=NA
  w = rowSums(is.na(tmp))
  tmp = tmp[w==0,]
  #Make foetal,child and adult bulk heatmaps
  tt = log2(bulk$cnts[gsub('.*_','',genesToPlot[match(rownames(tmp),gsub('_.*','',genesToPlot))]),])
  tt[tt< -5]=-5
  tt[tt> 10]=10
  fAnn = anno_boxplot(tt[,bulk$manifest$Group=='Foetal'],which='row')
  cAnn = anno_boxplot(tt[,bulk$manifest$Group=='PaedNormal'],which='row')
  aAnn = anno_boxplot(tt[,bulk$manifest$Group=='Normal'],which='row')
  rightAnn = HeatmapAnnotation(foetal=fAnn,
                               child=cAnn,
                               adult=aAnn,
                               which='row',
                               show_annotation_name=TRUE,
                               annotation_name_side='top',
                               width=unit(8,'cm'))
  #Without row norm
  draw(Heatmap(tmp,
               cluster_columns=FALSE,
               split = ss[grep(nom,ss)][w==0],
               name=nom,
               column_title='scExpression')+rightAnn)
  for(i in seq(length(unique(ss)))){
    decorate_heatmap_body(nom,{
                          tmp = cumsum(lengths(ord))/sum(lengths(ord))
                          tmp = tmp[-length(tmp)]
                          for(i in seq_along(tmp))
                            grid.lines(c(tmp[i],tmp[i]),c(0,1),gp=gpar(col='black',lty=2))},slice=i)
  }
}
dev.off()

####################
# genotypeOverview #
####################
#Cluster level genotyping
bigdf=list()
for(PDlab in names(PDids)){
  print(PDlab)
  PDid=PDids[PDlab]
  #Load the Subs.  They always exist
  subs = readRDS(file.path(fp_genotyping,sprintf('%s.subs.RDS',PDid)))
  #And only include the clonal ones
  df = subs$cnts[subs$cnts$MUT_AF>0.2,]
  #Use the most up to date cluster definitions
  df$Cluster = clustTab$Cluster[match(df$DropletID,clustTab$DropletID)]
  df$ClusterID = gsub('.*___','',df$Cluster)
  df$Experiment = gsub('___.*','',df$Cluster)
  df = split(df,df$Cluster)
  df = lapply(df,function(e)
              data.frame(Cluster=unique(e$Cluster),
                         Good_depth=sum(e$Good_depth),
                         matCount = sum(e$AltCount),
                         patCount = sum(e$RefCount),
                         numSNPs = length(unique(paste0(e$Chr,':',e$Pos))),
                         numCells = length(unique(e$DropletID))
                         ))
  df = do.call(rbind,df)
  df$MAF = df$matCount/df$Good_depth
  df$MAF_Low = ifelse(df$matCount==0,0,qbeta(alpha,df$matCount,df$Good_depth-df$matCount+1))
  df$MAF_High = ifelse(df$matCount==df$Good_depth,1,qbeta(1-alpha,df$matCount+1,df$Good_depth-df$matCount))
  df$pValue = pbinom(df$matCount-1,df$Good_depth,0.05,lower.tail=FALSE)
  df$qValue = p.adjust(df$pValue,method='BH')
  df$Compartment = gsub('___.*','',df$Cluster)
  df$ClusterID = gsub('.*___','',df$Cluster)
  df$Sample = PDlab
  df$Type='Subs'
  sdf = df
  #Add in the CN if it exists
  fnom=file.path(fp_genotyping,sprintf('%s.segments.RDS',PDid))
  if(file.exists(fnom)){
    dat = readRDS(fnom)
    #Exclude non LOH regions.  Unless that's all we've got.
    if(!(PDlab %in% c('PapRCC','Wilms3'))){
      dat$cnts = dat$cnts[which(dat$cnts$SegmentID %in% which(gsegs[[PDlab]]$minor==0)),]
      dat$cCnts = dat$cCnts[which(dat$cCnts$SegmentID %in% which(gsegs[[PDlab]]$minor==0)),]
      dat$clCnts = dat$clCnts[which(dat$clCnts$SegmentID %in% which(gsegs[[PDlab]]$minor==0)),]
      #This one is a bit harder, need to find if gene is in approved region.  Use adjusted gCnts to define allowed geneIDs
      dat$sgCnts = dat$sgCnts[dat$sgCnts$geneID %in% dat$gCnts$geneID[which(dat$gCnts$SegmentID %in% which(gsegs[[PDlab]]$minor==0))],]
    }
    #Merge the chromosomes
    df = dat$cnts
    #Use the most up to date cluster definitions
    df$Cluster = clustTab$Cluster[match(df$DropletID,clustTab$DropletID)]
    df$ClusterID = gsub('.*___','',df$Cluster)
    df$Experiment = gsub('___.*','',df$Cluster)
    df = split(df,df$Cluster)
    df = lapply(df,function(e)
                data.frame(Cluster=unique(e$Cluster),
                           Good_depth=sum(e$Good_depth),
                           matCount = sum(e$matCount),
                           patCount = sum(e$patCount),
                           numSNPs = length(unique(paste0(e$Chr,':',e$Pos))),
                           numCells = length(unique(e$DropletID))
                           ))
    df = do.call(rbind,df)
    df$MAF = df$matCount/df$Good_depth
    df$MAF_Low = ifelse(df$matCount==0,0,qbeta(alpha,df$matCount,df$Good_depth-df$matCount+1))
    df$MAF_High = ifelse(df$matCount==df$Good_depth,1,qbeta(1-alpha,df$matCount+1,df$Good_depth-df$matCount))
    df$pValue = pbinom(df$matCount-1,df$Good_depth,0.5,lower.tail=FALSE)
    df$qValue = p.adjust(df$pValue,method='BH')
    df$Compartment = gsub('___.*','',df$Cluster)
    df$ClusterID = gsub('.*___','',df$Cluster)
    df$Sample = PDlab
    df$Type='CN'
    df = rbind(sdf,df)
  }
  bigdf[[PDlab]]=df
}
#Make the plot
bigdf = do.call(rbind,bigdf)
bigdf$ClusterID = paste0(compNames$clusterName[match(bigdf$Compartment,gsub('.*_','',compNames$dirtyName))],as.character(bigdf$ClusterID))
tmp = unique(bigdf$ClusterID) 
bigdf$ClusterID = factor(bigdf$ClusterID,levels=tmp[order(gsub('[0-9]+','',tmp),as.numeric(gsub('[A-Z]+','',tmp)))])
#Adjust the facet labels
bigdf$CompartmentLabs = compNames$prettyName[match(bigdf$Compartment,gsub('.*_','',compNames$dirtyName))]
#Filter out the rubbish
bigdf = dplyr::filter(bigdf,CompartmentLabs!='Indistinct')
gg = ggplot(bigdf,aes(ClusterID,MAF,group=Type)) +
  geom_point(aes(colour=qValue<.01,shape=Type),position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=MAF_Low,ymax=MAF_High,shape=Type),position=position_dodge(width=0.5)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,size=6))+
  facet_grid(Sample ~ CompartmentLabs,scales='free_x',labeller=label_wrap_gen())
pdf(file.path(plotDir,'genotypeOverview.pdf'),width=14,height=14)
plot(gg)
dev.off()

####################################
# cnGenotype and cnGenotypeControl #
####################################
#CN genotyping by cell for all compartments with information
#tgts = c('T','IT','N','IN','EN','PT','R')
tgts = c('T','IT')
figNoms = c('cnGenotype.pdf','cnGenotypeControl.pdf')
for(i in seq_along(tgts)){
  pdf(file.path(plotDir,figNoms[i]),width=14,height=14)
  m = match(tgts[i],compNames$clusterName)
  srcDir = compNames$dirtyName[m]
  print(srcDir)
  srcLab = compNames$prettyName[m]
  mDat = readRDS(file.path(baseDir,srcDir,'batchCorrected','metadata.RDS'))$meta.data
  rownames(mDat) = paste0(mDat$Sanger_study_ID,'___',gsub('^[0-9]+_','',rownames(mDat)))
  #Label at the level needed to find PDid
  mDat$PD_Label = paste0(mDat$Experiment,ifelse(mDat$TissueDiseaseState=='NephrogenicRest','NR',''))
  df = mDat
  df$Ncnts=NA
  df$Mcnts=NA
  #Get the Ncnts and Mcnts for each droplet
  for(PDlab in unique(df$PD_Label)){
    print(PDlab)
    #These ones we have no useful anything.
    if(PDlab=='Wilms1NR' | !(PDlab %in% names(PDids)))
      next
    #Get the PDid to load
    PDid = PDids[PDlab]
    #Now open the CN data for this sample
    CNdat = readRDS(file.path(fp_genotyping,sprintf('%s.segments.RDS',PDid)))
    #And extract the information about the cells we care about this time
    dat = CNdat$cnts[CNdat$cnts$DropletID %in% rownames(df[df$PD_Label==PDlab,]),]
    #Exclude non LOH regions.  Unless that's all we've got.
    if(!(PDlab %in% c('PapRCC','Wilms3')))
      dat = dat[which(dat$SegmentID %in% which(gsegs[[PDlab]]$minor==0)),]
    #Plot the MAF sorted by Ncnts for each Label
    Ncnts = sapply(split(dat$Good_depth,dat$DropletID),sum)
    Mcnts = sapply(split(dat$matCount,dat$DropletID),sum)
    df$Ncnts[match(names(Ncnts),rownames(df))]=Ncnts
    df$Mcnts[match(names(Mcnts),rownames(df))]=Mcnts
  }
  df$Mcnts[is.na(df$Mcnts)]=0
  df$Ncnts[is.na(df$Ncnts)]=0
  df$MAF = df$Mcnts/df$Ncnts
  df$MAF_Low = ifelse(df$Mcnts==0,0,qbeta(alpha,df$Mcnts,df$Ncnts-df$Mcnts+1))
  df$MAF_High = ifelse(df$Mcnts==df$Ncnts,1,qbeta(1-alpha,df$Mcnts+1,df$Ncnts-df$Mcnts))
  df = df[order(df$PD_Label,df$Ncnts),]
  df$Idx = seq(nrow(df))
  df$pValue = pbinom(df$Mcnts-1,df$Ncnts,0.5,lower.tail=FALSE)
  df$qValue = p.adjust(df$pValue,method='BH')
  df$sig = ifelse(df$pValue<1e-2,-2,log10(df$pValue))
  df$sig = ifelse(df$MAF<0.5,0.5,df$MAF)
  #Now make a plot showing the MAF for each point
  gg = ggplot(df,aes(tSNE1,tSNE2,colour=PD_Label,shape=pValue<.05)) +
    geom_point(data=df[df$Ncnts==0,],size=0.1) +
    geom_point(aes(fill=sig,size=log2(Ncnts)),alpha=1/2) +
    scale_shape_manual(breaks=c(TRUE,FALSE),values=c(21,22))+
    #scale_shape_manual(values=diseaseShapes,breaks=names(diseaseShapes))+
    scale_fill_gradient(high='black',low='#f2f2f2')+
    labs(shape='Significant?',
         fill='MAF',
         size='log2(#Counts)',
         colour='Experiment')+
    ggtitle(sprintf("CN Genotyping %s",srcLab))
  plot(gg)
  #Facet by experiment to show each sample more clearly
  #pdf(file.path(plotDir,sprintf('tumourGenotypingFacetedCN_%s.pdf',srcLab)),width=14,height=14)
  #plot(gg + facet_wrap(~PD_Label))
  #dev.off()
  dev.off()
}

######################################
# mutGenotype and mutGenotypeControl #
######################################
# Sub genotyping 
#tgts = c('T','IT','N','IN','EN','PT','R')
tgts = c('T','IT')
figNoms = c('mutGenotype.pdf','mutGenotypeControl.pdf')
for(i in seq_along(tgts)){
  pdf(file.path(plotDir,figNoms[i]),width=14,height=14)
  m = match(tgts[i],compNames$clusterName)
  srcDir = compNames$dirtyName[m]
  print(srcDir)
  srcLab = compNames$prettyName[m]
  mDat = readRDS(file.path(baseDir,srcDir,'batchCorrected','metadata.RDS'))$meta.data
  rownames(mDat) = paste0(mDat$Sanger_study_ID,'___',gsub('^[0-9]+_','',rownames(mDat)))
  #Label at the level needed to find PDid
  mDat$PD_Label = paste0(mDat$Experiment,ifelse(mDat$TissueDiseaseState=='NephrogenicRest','NR',''))
  df = mDat
  df$Ncnts=NA
  df$Mcnts=NA
  #Get the Ncnts and Mcnts for each droplet
  for(PDlab in unique(df$PD_Label)){
    print(PDlab)
    #These ones we have no useful anything.
    if(!(PDlab %in% names(PDids)))
      next
    #Get the PDid to load
    PDid = PDids[PDlab]
    #Now open the CN data for this sample
    subs = readRDS(file.path(fp_genotyping,sprintf('%s.subs.RDS',PDid)))
    #And extract the information about the cells we care about this time
    dat = subs$cnts[subs$cnts$DropletID %in% rownames(df[df$PD_Label==PDlab,]),]
    #And only include the clonal ones
    dat = dat[dat$MUT_AF>0.2,]
    if(nrow(dat)==0)
      next
    #Plot the MAF sorted by Ncnts for each Label
    Ncnts = sapply(split(dat$Good_depth,dat$DropletID),sum)
    Mcnts = sapply(split(dat$AltCount,dat$DropletID),sum)
    df$Ncnts[match(names(Ncnts),rownames(df))]=Ncnts
    df$Mcnts[match(names(Mcnts),rownames(df))]=Mcnts
  }
  df$Mcnts[is.na(df$Mcnts)]=0
  df$Ncnts[is.na(df$Ncnts)]=0
  df$MAF = df$Mcnts/df$Ncnts
  df$MAF_Low = ifelse(df$Mcnts==0,0,qbeta(alpha,df$Mcnts,df$Ncnts-df$Mcnts+1))
  df$MAF_High = ifelse(df$Mcnts==df$Ncnts,1,qbeta(1-alpha,df$Mcnts+1,df$Ncnts-df$Mcnts))
  df = df[order(df$PD_Label,df$Ncnts),]
  df$Idx = seq(nrow(df))
  df$pValue = pbinom(df$Mcnts-1,df$Ncnts,0.01,lower.tail=FALSE)
  df$qValue = p.adjust(df$pValue,method='BH')
  df$sig = ifelse(df$pValue<1e-2,-2,log10(df$pValue))
  df$sig = ifelse(df$MAF<0.5,0.5,df$MAF)
  #Now make a plot showing the MAF for each point
  gg = ggplot(df,aes(tSNE1,tSNE2,colour=PD_Label,shape=pValue<.05)) +
    geom_point(data=df[df$Ncnts==0,],size=0.1) +
    geom_point(aes(fill=MAF,size=log2(Ncnts)),alpha=1/2) +
    scale_shape_manual(breaks=c(TRUE,FALSE),values=c(21,22))+
    #scale_shape_manual(values=diseaseShapes,breaks=names(diseaseShapes))+
    scale_fill_gradient(high='black',low='#f2f2f2')+
    labs(shape='Significant?',
         fill='MAF',
         size='log2(#Counts)',
         colour='Experiment')+
    ggtitle(sprintf("Sub Genotyping %s",srcLab))
  #pdf(file.path(plotDir,sprintf('tumourGenotypingSubs_%s.pdf',srcLab)),width=14,height=14)
  plot(gg)
  ##Facet by experiment to show each sample more clearly
  #pdf(file.path(plotDir,sprintf('tumourGenotypingFacetedSubs_%s.pdf',srcLab)),width=14,height=14)
  #plot(gg + facet_wrap(~PD_Label))
  #dev.off()
  dev.off()
}


###############
# algoMarkers #
###############
#Aggregate all marker gene information into one table.  Also need to run Seurat DE for each I think.
tgts = c('T','IT','N','IN','F','IF','EN','PT')
bigTab=list()
for(i in seq_along(tgts)){
  m = match(tgts[i],compNames$clusterName)
  srcDir = compNames$dirtyName[m]
  #Cluster summaries
  cSums = grep('summaryOfCluster[0-9]+\\.RDS',list.files(file.path(baseDir,srcDir,'batchCorrected')),value=TRUE)
  #Load each in turn
  for(cSum in cSums){
    dat = readRDS(file.path(baseDir,srcDir,'batchCorrected',cSum))
    tmp = dat$markerGeneDat
    tmp$Cluster = paste0(tgts[i],gsub('summaryOfCluster([0-9]+)\\.RDS','\\1',cSum))
    tmp$Compartment = compNames$prettyComputerNames[m]
    tmp$EnsGeneID = gsub('.*_','',rownames(tmp))
    tmp$Symbol = gsub('_ENSG[0-9]+','',rownames(tmp))
    bigTab[[tmp$Cluster[1]]] = tmp
  }
}
bigTab = do.call(rbind,bigTab)
bigTab = dplyr::transmute(bigTab,
                   Cluster,
                   Symbol,
                   geneFrequencyInsideCluster=geneFrequency,
                   geneFrequencyOutsideCluster,
                   qval,
                   tfidf,
                   idf,
                   geneExpressionInsideCluster=geneExpression,
                   geneExpressionOutsideCluster,
                   geneFrequencyGlobal,
                   geneExpressionGlobal,
                   isKeyGene,
                   EnsGeneID)
write.table(bigTab,file.path(plotDir,'algoMarkers.tsv'),sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)

##################
# sampleManifest #
##################
#Clean up the manifest and save it out
tmp = mani
tmp = dplyr::transmute(tmp,
                SangerID=Sanger_study_ID,
                Channel10X=channel10X,
                Label=label,
                Experiment,
                TissueDiseaseState,
                Organ,
                Location1,
                Location2,
                BiologicalRepNo,
                TechnicalRepNo,
                Sort,
                AgeInMonthsPostConception=Age_donor.MonthsPostConception.,
                PatientDiseaseState,
                BulkID=Bulk_ID
                )
#Swap PDids for Wilms 2 and 3, they're the wrong way around
w = grep('PD3727[26]',tmp$BulkID)
mask = tmp$BulkID[w]
w1 = grep('PD37272',mask)
w2 = grep('PD37276',mask)
mask[w1] = gsub('PD37272','PD37276',mask[w1])
mask[w2] = gsub('PD37276','PD37272',mask[w2])
tmp$BulkID[w] = mask
#Add in the conversion from SangerID to EGA
iDat = read.delim('~/Projects/KidneySC/IrodsMetaData/Kidney_Single_Cell_Study_cram_manifest_INFO_from_iRODS.txt',sep='\t',header=TRUE)
iDat = rbind(iDat,read.delim('~/Projects/KidneySC/IrodsMetaData/Pilot_Fetal_Cell_Atlas_RNAseq_cram_manifest_INFO_from_iRODS.txt',sep='\t',header=TRUE))
#Now re-shape the columns to keep the ones that we want to keep
iDat = unique(dplyr::select(iDat,sample,sample_accession_number,study_accession_number,study_title))
iDat = dplyr::rename(iDat,
              SangerID=sample,
              EGA_ID=sample_accession_number,
              EGA_StudyID=study_accession_number,
              EGA_StudyTitle=study_title)
tmp = dplyr::left_join(tmp,iDat)
#Filter out the irrelevant samples
tmp = tmp[!grepl('^F15_',tmp$Label),]
tmp = tmp[!grepl('Ski|Dec|Liv',tmp$Label),]
tmp = tmp[!grepl('AV',tmp$Label),]
tmp = tmp[!grepl('Lym',tmp$Label),]
tmp = tmp[!grepl('_H_',tmp$Label),]
extraIDs=c(
`4602STDY7018923`='EGAN00001647933',
`4602STDY7018924`='EGAN00001647934',
`4602STDY7018925`='EGAN00001647935',
`4602STDY7018920`='EGAN00001647930',
`4602STDY7018921`='EGAN00001647931',
`4602STDY7018922`='EGAN00001647932',
`4602STDY7018927`='EGAN00001647937',
`4602STDY7090428`='EGAN00001647941',
`4602STDY7090429`='EGAN00001647942',
`4602STDY7090430`='EGAN00001647943',
`4602STDY7090425`='EGAN00001647938',
`4602STDY7090426`='EGAN00001647939',
`4602STDY7090427`='EGAN00001647940',
`4602STDY7090431`='EGAN00001647944',
`4602STDY6949193`='EGAN00001647951',
`4602STDY6949196`='EGAN00001647947',
`4602STDY6949198`='EGAN00001647949',
`4602STDY6949192`='EGAN00001647950',
`4602STDY6949195`='EGAN00001647953',
`4602STDY6949191`='EGAN00001647946',
`4602STDY6949194`='EGAN00001647952',
`4602STDY6949197`='EGAN00001647948',
`4602STDY6976426`='EGAN00001647958',
`4602STDY6976427`='EGAN00001647959',
`4602STDY6976428`='EGAN00001647960',
`4602STDY6976422`='EGAN00001647954',
`4602STDY6976423`='EGAN00001647955',
`4602STDY6976424`='EGAN00001647956',
`4602STDY6976425`='EGAN00001647957'
)
#Fix these ones up
w=is.na(tmp$EGA_ID)
tmp$EGA_ID[w] = extraIDs[tmp$SangerID[w]]
#Add in droplet counts for each
tmp$NoCellsBeforeQC = table(gsub('___.*','',colnames(rtoc)))[tmp$SangerID]
write.table(tmp,file.path(plotDir,'sampleManifest.tsv'),sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)


################
# cellManifest #
################
#Individual cell summary
#Work out what the bad channels are
cMani = read.table(file.path(plotDir,'sampleManifest.tsv'),sep='\t',header=TRUE)
#Start with Cluster table
cellTable = clustTab
cellTable = dplyr::filter(cellTable,Group!='notKidney')
rownames(cellTable) = cellTable$DropletID
#Clean up names
cellTable = dplyr::mutate(cellTable,channel=basename(channel),SangerID=gsub('___.*','',DropletID))
cellTable = dplyr::select(cellTable,barcode,channel,SangerID,DropletID,Group,Cluster,ClusterID,Compartment)
#And discard the shit channels that are not kidney, even though we haven't labelled them as such
cellTable = cellTable[!is.na(match(cellTable$SangerID,cMani$SangerID)),]
#Get the MT frac, gene count and UMI count for each cell
cellTable$nUMI = Matrix::colSums(rtoc[,cellTable$DropletID])
cellTable$nGenes = Matrix::colSums(rtoc[,cellTable$DropletID]>0)
cellTable$MTfrac = Matrix::colSums(rtoc[grep('^MT-',rownames(rtoc)),cellTable$DropletID])/cellTable$nUMI
#Print breakdown of cell numbers
if(TRUE){
  tmp = cellTable
  #Create a by biopsy category
  tmp$biopsy = gsub('_[0-9]$','',cMani$Label[match(tmp$SangerID,cMani$SangerID)])
  bLabs = unique(tmp$biopsy)
  cntsByBiopsy = list()
  cntsByBiopsy[['Start']] = ltable(tmp$biopsy,levels=bLabs)
  message(sprintf("To start with, their are %d cells.",nrow(tmp)))
  cnt = nrow(tmp)
  tmp = tmp[tmp$nGenes>200,]
  cntsByBiopsy[['GeneDrop']] = ltable(tmp$biopsy,levels=bLabs)
  message(sprintf("After dropping %d cells with fewer than 200 genes, there are %d cells.",cnt-nrow(tmp),nrow(tmp)))
  cnt = nrow(tmp)
  tmp = tmp[tmp$MTfrac<=0.2,]
  cntsByBiopsy[['MTDrop']] = ltable(tmp$biopsy,levels=bLabs)
  message(sprintf("After dropping %d cells with more than 20%% MT expression, there are %d cells.",cnt-nrow(tmp),nrow(tmp)))
  cnt = nrow(tmp)
  tmp = tmp[!(tmp$Group %in% c('unassigned','rubbish')),]
  cntsByBiopsy[['DoubletDrop']] = ltable(tmp$biopsy,levels=bLabs)
  message(sprintf("After dropping %d cells with ambiguous profiles, there are %d cells.",cnt-nrow(tmp),nrow(tmp)))
  message(sprintf("These are split into the following categories:"))
  w = grepl('Immune',tmp$Group)
  message(sprintf("Immune: %d",sum(w)))
  cntsByBiopsy[['Immune']] = ltable(tmp$biopsy[w],levels=bLabs)
  message(sprintf("Which is further split into:"))
  w = (tmp$Group=='foetalImmuneV3')
  cntsByBiopsy[['FoetalImmune']] = ltable(tmp$biopsy[w],levels=bLabs)
  message(sprintf(" - FoetalImmune: %d",sum(w)))
  w = (tmp$Group=='normalImmune')
  cntsByBiopsy[['NormalImmune']] = ltable(tmp$biopsy[w],levels=bLabs)
  message(sprintf(" - NormalImmune: %d",sum(w)))
  w =(tmp$Group=='tumourImmune')
  cntsByBiopsy[['TumourImmune']] = ltable(tmp$biopsy[w],levels=bLabs)
  message(sprintf(" - TumourImmune: %d",sum(w)))
  w = !grepl('Immune',tmp$Group)
  cntsByBiopsy[['Non-Immune']] = ltable(tmp$biopsy[w],levels=bLabs)
  message(sprintf("Non-Immune: %d",sum(w)))
  message(sprintf("Which is further split into:"))
  w = (tmp$Group=="foetalEpitheliumAndVascularV3")
  cntsByBiopsy[['Foetal']] = ltable(tmp$biopsy[w],levels=bLabs)
  message(sprintf(" - Foetal: %d",sum(w)))
  w = (tmp$Group=="tumourEpitheliumAndVascular")
  cntsByBiopsy[['Tumour']] = ltable(tmp$biopsy[w],levels=bLabs)
  message(sprintf(" - Tumour: %d",sum(w)))
  w = grepl('roximalTubular',tmp$Group)
  cntsByBiopsy[['Normal']] = ltable(tmp$biopsy[w],levels=bLabs)
  message(sprintf(" - Normal: %d",sum(w)))
  message(sprintf("This category is finally sub-divided into:"))
  w = (tmp$Group=='proximalTubularV2')
  cntsByBiopsy[['PT']] = ltable(tmp$biopsy[w],levels=bLabs)
  message(sprintf("   + Non-specific PT: %d",sum(w)))
  w = (tmp$Group=='normalEpitheliumAndVascularWithoutProximalTubularV2')
  cntsByBiopsy[['Non-PT']] = ltable(tmp$biopsy[w],levels=bLabs)
  message(sprintf("   + All other normal: %d",sum(w)))
  #Finalise the counts table
  cntsByBiopsy = do.call(rbind,cntsByBiopsy)
  write.table(cntsByBiopsy,file.path(plotDir,'cntsByBiopsy.tsv'),sep='\t',quote=FALSE,col.names=TRUE,row.names=TRUE)
  #Make a barplot of relative drops
  df = cntsByBiopsy[1:4,]
  df = t(t(df)/as.matrix(df)[1,])
  df = melt(df)
  #Order 
  ord = strsplit(levels(df$Var2),'_')
  df$Var2 = factor(as.character(df$Var2),levels=levels(df$Var2)[order(sapply(ord,`[`,1),sapply(ord,`[`,3),decreasing=TRUE)])
  gg = ggplot(df,aes(Var2,value)) +
    geom_col() +
    theme(axis.text.x = element_text(angle=90,hjust=1)) +
    facet_grid(~Var1) 
  pdf(file.path(plotDir,'FlowFraction.pdf'),width=28,height=14)
  plot(gg)
  dev.off()
}
#Change the rubbish to be all rubbish
cellTable$Group[cellTable$Group=='unassigned'] = 'rubbish'
cellTable = dplyr::mutate(cellTable,QCpass = MTfrac<=0.2 & nGenes>200 & Group != 'rubbish' )
#Do some final fixing up and re-ordering
cellTable = dplyr::mutate(cellTable,Compartment=compNames$prettyComputerNames[match(Group,compNames$shortDirtyNames)],ClusterID = ifelse(is.na(ClusterID),NA,paste0(compNames$clusterName,ClusterID)))
cellTable = dplyr::select(cellTable,-Group,-Cluster,-channel)
#Add in some basic information about the 10X channel.
cellTable = dplyr::inner_join(cellTable,dplyr::transmute(mani,SangerID=Sanger_study_ID,Source=label))
write.table(cellTable,file.path(plotDir,'cellManifest.tsv'),sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)

###############
# dnaManifest #
###############
# Table converting SangerID to EGA_ID
iDat = read.delim('~/Projects/KidneySC/IrodsMetaData/Orphan_Tumour_Study_-_familial_neuroblastoma_cram_manifest_INFO_from_iRODS.txt',sep='\t',header=TRUE)
iDat = rbind(iDat,read.delim('~/Projects/KidneySC/IrodsMetaData/Kidney_tumour_DNA_cram_manifest_INFO_from_iRODS.txt',sep='\t',header=TRUE))
#Get rid of the columns that are CRAM specific, we just want sample level info.
iDat = unique(dplyr::select(iDat,sample,sample_donor_id,sample_accession_number,study_accession_number,study_title))
iDat = dplyr::rename(iDat,
              SangerID=sample,
              PD_ID=sample_donor_id,
              EGA_ID=sample_accession_number,
              EGA_StudyID=study_accession_number,
              StudyTitle=study_title)
#Keep only the relevant ones
bulkIDs = c('PD37104a',
'PD37104b',
'PD37228c',
'PD37228f',
'PD35918h',
'PD35918g',
'PD36793a',
'PD36793c',
'PD36165d',
'PD36165b',
'PD36165e',
'PD37272a',
'PD37272g',
'PD37276a',
'PD37276g')
iDat = iDat[iDat$PD_ID %in% bulkIDs,]
#Add in the experiment name this corresponds to
cMani = read.table(file.path(plotDir,'sampleManifest.tsv'),sep='\t',header=TRUE)
iDat$Experiment = sapply(lapply(bulkIDs,function(e) cMani$Experiment[grep(gsub('[a-z]*','',e),cMani$BulkID)]),unique)[match(iDat$PD_ID,bulkIDs)]
write.table(iDat,file.path(plotDir,'dnaManifest.tsv'),sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)
#These are the studies that contain the sampleIDs we care about
#Wilms bulk DNA sequencing
#EGAS00001002171 - "Orphan Tumour Study - familial neuroblastoma"
#RCC bulk DNA sequencing
#EGAS00001002486 - "Kidney tumour_DNA"
#Post birth single cell data
#EGAS00001002325 - "Kidney Single Cell Study"
#Foetal single cell data
#EGAS00001002553 - "Pilot Fetal Cell Atlas_RNAseq"

##############
# cnSegments #
##############
#CN changes in bulk DNA
cnTab = GRangesList(gsegs)
cnTab = unlist(cnTab)
cnTab$Experiment = names(cnTab)
cnTab = data.frame(cnTab)
cnTab = dplyr::transmute(cnTab,
                  Chr=seqnames,
                  Start=start,
                  End=end,
                  TotalCN=total,
                  MinorCN=minor,
                  Experiment)
cnTab$PD_ID = basename(dirname(snpFiles))[match(cnTab$Experiment,names(snpFiles))]
write.table(cnTab,file.path(plotDir,'cnSegments.tsv'),sep='\t',quote=FALSE,col.names=TRUE,row.names=FALSE)

#################
# tableOfCounts #
#################
#Write out the final table of counts
toc = rtoc[,cellTable$DropletID]
#Save in the convenient format
saveRDS(toc,file.path(plotDir,'tableOfCounts.RDS'))
#And in matrix mart for readability/future proofing
writeMM(toc,file.path(plotDir,'tableOfCounts.mtx'))
#And save row and column name tables with it
df = data.frame(RowNumber = seq_along(rownames(toc)),
                GeneLabel = rownames(toc),
                Symbol = gsub('_.*','',rownames(toc)),
                EnsemblID = gsub('.*_','',rownames(toc))
                )
write.table(df,file.path(plotDir,'tableOfCounts_rowLabels.tsv'),sep='\t',quote=FALSE,col.names=TRUE,row.names=FALSE)
df = data.frame(ColNumber = seq_along(colnames(toc)),
                DropletID = colnames(toc),
                SangerID = gsub('___.*','',colnames(toc)),
                Barcode = gsub('.*___','',colnames(toc))
                )
write.table(df,file.path(plotDir,'tableOfCounts_colLabels.tsv'),sep='\t',quote=FALSE,col.names=TRUE,row.names=FALSE)

##########################################
# Convert symbolic names to Numeric ones #
##########################################
numConv = read.delim(fp_num,sep='\t',header=TRUE)
numConv$Target = with(numConv,paste0(Label,SubLabel,SubSubLabel))
#List all files
files = list.files(plotDir)
extensions = gsub('^[^\\.]*\\.','',files)
for(i in seq_len(nrow(numConv))){
  #Check if it's scripted and can be copied
  if(is.na(numConv$PurelyScripted[i]))
    next
  if(numConv$PurelyScripted[i]!='Yes')
    next
  #Find the file, if it exists
  g = grep(sprintf('^%s\\.',numConv$ShortName[i]),files)
  if(length(g)==0)
    next
  #Work out what to rename it to
  tgt = numConv$Target[i]
  if(any(startsWith(tgt,c('R','S',1:9)))){
    out = 'Figure'
  }else if(startsWith(tgt,'T')){
    out = 'Table'
    tgt = gsub('T','S',tgt)
  }else if(startsWith(tgt,'D')){
    out = 'Data'
    tgt = gsub('D','S',tgt)
  }else{
    stop("Unknown type")
  }
  out = paste0(out,tgt)
  file.copy(file.path(plotDir,files[g]),file.path(plotDir,paste0(out,'.',extensions[g])),overwrite=FALSE)
  print(out)
}
