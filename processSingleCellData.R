#################################
# Load libraries and parameters #
#################################
 
import(core,all=TRUE)
import(genome,all=TRUE)
library(grid)
library(Seurat)
library(Matrix)
library(GenomicFiles)
library(GenomicAlignments)
library(plot3D)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(Rtsne)
library(GenomicRanges)
library(org.Hs.eg.db)
library(gridExtra)
library(pheatmap)

###########################
# Define common functions #
###########################

#' Calculates the tf-idf for a set of target cells against a background of the cells given in "universe"
#'
#' @param data The data matrix to use.
#' @param target Columns that are the target
#' @param universe Columns that we should consider (target must be a subset)
tfidf = function(data,target,universe){
  if(!all(target %in% universe))
    stop('Target must be a subset of universe')
  nObs = Matrix::rowSums(data[,target,drop=FALSE]>0)
  nTot = Matrix::rowSums(data[,universe,drop=FALSE]>0)
  tf = nObs/length(target)
  idf = log(length(universe)/nTot)
  score = tf*idf
  #Calculate p-value for significance based on using a hypergeometric distribution to simulate the results of infinite random sampling
  pvals = phyper(nObs-1,nTot,length(universe)-nTot,length(target),lower.tail=FALSE)
  qvals = p.adjust(pvals,method='BH')
  ntf = (exp(-idf)*length(universe)-tf*length(target))/(length(universe)-length(target))
  return(data.frame(geneFrequency=tf,
                    geneFrequencyOutsideCluster=ntf,
                    geneFrequencyGlobal=exp(-idf),
                    geneExpression=Matrix::rowMeans(data[,target,drop=FALSE]),
                    geneExpressionOutsideCluster = Matrix::rowMeans(data[,universe[!(universe%in%target)],drop=FALSE]),
                    geneExpressionGlobal = Matrix::rowMeans(data),
                    idf=idf,
                    tfidf=score,
                    qval=qvals)[order(score,decreasing=TRUE),])
}
getAnnotations = function(genes,gene_ids='ensembl_gene_id',attributes=c("ensembl_gene_id","hgnc_symbol", "chromosome_name","gene_biotype", "start_position", "end_position","description"),dataset='hsapiens_gene_ensembl',biomart='ENSEMBL_MART_ENSEMBL',host='www.ensembl.org'){
  if(is.null(host)){
    bmart = biomaRt::useMart(biomart=biomart,dataset=dataset)
  }else{
    bmart = biomaRt::useMart(biomart=biomart,dataset=dataset,host=host)
  }
  dat = biomaRt::getBM(attributes=attributes,filters=gene_ids,values=genes,mart=bmart)
  #Now make sure that we have a line for every gene we asked for...
  bad = genes[!(genes%in%dat[,gene_ids])]
  #Add these ones back in 
  dat[bad,gene_ids]=bad
  return(dat)
}
extract_field  = function (string, field = 1, delim = "_")
{
    fields = as.numeric(unlist(strsplit(as.character(field),
        ",")))
    if (length(fields) == 1)
        return(strsplit(string, delim)[[1]][field])
    return(paste(strsplit(string, delim)[[1]][fields], collapse = delim))
}
Read10X=function(data.dir = NULL,use.symbols=FALSE){
  full_data <- list()
  for(i in seq_along(data.dir)){
    run <- data.dir[i]
    if (!dir.exists(run)){
      stop("Directory provided does not exist")
    }
    
    if(!grepl("\\/$", run)){
      run <- paste(run, "/", sep = "")
    }
    
    barcode.loc <- paste(run, "barcodes.tsv", sep ="")
    gene.loc <- paste(run, "genes.tsv", sep ="")
    matrix.loc <- paste(run, "matrix.mtx", sep ="")
    
    if (!file.exists(barcode.loc)){
      stop("Barcode file missing")
    }
    if (!file.exists(gene.loc)){
      stop("Gene name file missing")
    }
    if (!file.exists(matrix.loc)){
      stop("Expression matrix file missing")
    }
    
    data <- readMM(matrix.loc)
    cell.names <- readLines(barcode.loc)
    gene.names <- readLines(gene.loc)
    if(all(grepl("\\-1$", cell.names)) == TRUE) {
      cell.names <- as.vector(as.character(sapply(cell.names, extract_field, 1, delim = "-")))
    }
    rownames(data) <- make.unique(as.character(sapply(gene.names, extract_field, 2-as.numeric(!use.symbols), delim = "\\t"))) 
    
    if(is.null(names(data.dir))){
      if(i < 2){
        colnames(data) <- cell.names
      }
      else {
        colnames(data) <- paste0(i, "_", cell.names, sep = "") 
      }
    } else {
      colnames(data) <- paste0(names(data.dir)[i],"_",cell.names) 
    }
    full_data <- append(full_data, data)
  }
  full_data <- do.call(cbind, full_data)
  return(full_data)
}
buildManifest = function(srcDirs,manifestPaths = unique(file.path(dirname(gsub('/+$','',srcDirs)),'manifest.csv'))){
  #This is pulled directly from singleCellProcessV2.R and should be updated from there
  srcDirs = gsub('/+$','',srcDirs)
  #Get the channels from each source directory
  tmp = strsplit(srcDirs,'/')
  channelDat = data.frame(path = normalizePath(srcDirs),
                          channel = basename(srcDirs),
                          project = sapply(lapply(tmp,rev),`[`,2)
                          )
  #Load any extra data we can from the manifests
  for(manifestPath in manifestPaths){
    #Check if it exists
    if(file.exists(manifestPath)){
      #Read it in, manifests must load with this command and the first column must match srcDirs
      tmp = read.table(manifestPath,header=TRUE,row.names=NULL,sep=',')
      #Find path column
      if('path' %in% colnames(tmp)) pathCol=match('path',colnames(tmp)) else pathCol=1
      m = match(normalizePath(tmp[,pathCol]),normalizePath(srcDirs))
      if(all(is.na(m))){
        warnings(sprintf('Manifest %s contains no information about the data in srcDirs.  Skipping.',manifestPath))
        next
      }
      #Keep only the ones relavent to our data
      tmp = tmp[!is.na(m),]
      m = m[!is.na(m)]
      #Now add data from each column sequentially.  If duplicates occur, throw an error
      for(newCol in colnames(tmp)[-1]){
        #Check it's not a protected collumn name
        if(newCol %in% c('path','channel','project')){
          warning(sprintf('Column %s in manifest %s will be ignored.  Manifest cannot overwrite values of "path", "channel" or "project".',newCol,manifestPath))
          next
        }
        #Get the non-empty entries for this column
        w=which(!is.na(tmp[,newCol]))
        #Add it to the channelDat data.frame
        if(newCol %in% colnames(channelDat)){
          #Check that it's not over-writing an existing value
          if(any(!is.na(channelDat[m[w],newCol]))){
            stop("Duplicate manifest entries detected.  Please fix and re-run")
          }
        }
        #All fine, store the results
        channelDat[m[w],newCol]=tmp[w,newCol]
      }
    }
  }
  #Create default labels if no manifest is present
  defaultLabels = paste0(channelDat$project,'_',channelDat$channel)
  if('label' %in% colnames(channelDat)){
    w = which(is.na(channelDat$label))
    channelDat$label[w] = defaultLabels[w]
  }else{
    channelDat$label = defaultLabels
  }
  return(channelDat)
}

##################
# Pre-processing #
##################

#' pre-process data
#' 
#' This should operate on one channel only.  The output will be a count matrix with only cells included and any filtering of the soup done (or other pre-processing).  Can also save other summarising information such as the soupFractions/counts, parameters and code used to run the processing, etc.
#'
#' @param path The top level folder of the 10x data
#' @param soupBelow Any droplet with fewer than this many UMIs is considered pure background.
#' @param cellAbove Any droplet with more than this many UMIs is considered a cell.
#' @param minGenesForCell Any droplet without this any genes is flagged as probably not a cell (but not removed).
#' @param maxMT Any droplet with more than this MT fraction is flagged as probably not a cell (but not removed).
#' @param FDRcut The FDR threshold to use for selecting droplets using Poisson method.
#' @param verbose Be verbose?
#' @param doPlot Plot things?
#' @param BPPARAM Parallel execution of Poisson method.
preProcessChannel = function(path,soupBelow=100,cellAbove=1e4,minGenesForCell=200,maxMT=0.2,FDRcut=0.05,verbose=TRUE,doPlot=FALSE,BPPARAM=SerialParam()) {
  #Load the raw data
  if(verbose)
    message(header(sprintf('Processing channel %s',path)))
  srat = Read10X(file.path(path,'raw_gene_bc_matrices/GRCh38'))
  #Load the gene symbols
  symbolTable = read.table(file.path(path,'raw_gene_bc_matrices/GRCh38/genes.tsv'),sep='\t',header=FALSE)
  #Rename rows to include both
  rownames(srat) = paste0(symbolTable[match(rownames(srat),symbolTable[,1]),2],'_',rownames(srat))
  #Get the cellranger filtered barcodes
  bcodeFile = file.path(path,'filtered_gene_bc_matrices/GRCh38/barcodes.tsv')
  cRangerCells = readLines(bcodeFile)
  cRangerCells = gsub('-1$','',cRangerCells)
  #Summarise by UMIs and genes
  geneCounts = Matrix::colSums(srat>0)
  umiCounts = Matrix::colSums(srat)
  #And get the MT fraction
  mtFrac =Matrix::colSums(srat[grep('^MT-.*_',rownames(srat)),])/Matrix::colSums(srat)
  mtFrac[is.na(mtFrac)]=0
  ######################
  # Find soup and cell #
  ######################
  #Find the peak of the UMI distribution
  peak = findPeak(umiCounts,plot=FALSE)
  #Get the knee from the cumulative distribution
  if(sum(umiCounts>=soupBelow)<5){
    knee = 0
  }else{
    knee = EmptyDrops::findKneePoint(srat,lower=soupBelow)
  }
  #Get an enstimate of the knee from the CellRanger filtered cells
  kneeCR = min(umiCounts[match(cRangerCells,colnames(srat))],na.rm=TRUE)-1
  if(doPlot){
    #Plot the UMI distribution with some annotation
    y = sort(umiCounts,decreasing=TRUE)
    plot(log10(seq_along(y)),log10(y),'l',xlab='log10(Barcode Index)',ylab='log10(UMIs)',main=sprintf('UMI distribution'))
    abline(h=log10(c(kneeCR,knee)),col=c('green','red'),lty=1)
  }
  #Store which droplets pass the CR and Knee cut-offs
  wCellKnee = which(umiCounts > knee) 
  wCellCR = which(colnames(srat) %in% cRangerCells)
  #Now do the more liberal, Poisson method to find cells
  if(sum(umiCounts>=soupBelow)>5){
    tmp = EmptyDrops::testBarcodeAmbience(srat,wSoup=NA,lower=soupBelow,BPPARAM=BPPARAM)
    #Do FDR correction and identify true cells
    qval = tmp$PValue
    qval[tmp$Total > knee * 5] = 0
    qval = p.adjust(qval,method='BH')
  }else{
    qval = rep(1,length(wCellCR))
  }
  #Store the cells
  wCellPois = which(qval<FDRcut)
  #Make a plot of 
  if(doPlot){
    #Plot a diagnostic plot for all the cells
    pMask = rep(1,length(geneCounts))
    pMask[wCellPois] = 3
    cMask = mtFrac
    mtBreaks = c(0,.05,.1,.2,.5,1)
    cMask = cut(cMask,mtBreaks,include.lowest=TRUE)
    #Make a multi-panel plot showing the density
    layout(cbind(c(1,1,1,2)))
    plot(log10(umiCounts),log10(umiCounts/geneCounts),cex=ifelse(pMask==3,1,.1),col=as.numeric(cMask),main=sprintf('Cell overview'),pch=pMask)
    abline(v=log10(c(floor(peak/2),peak*2,knee,kneeCR)),col=c('black','black','red','green'),lty=c(2,2,1,1))
    abline(-log10(100),1,lty=3)
    abline(-log10(500),1,lty=4)
    abline(-log10(1000),1,lty=5)
    legend('topleft',levels(cMask),pch=1,col=seq_len(length(levels(cMask))),title='MT-Fraction')
    legend('left',c('Cell','Empty'),pch=c(3,1),title='Poisson Method',pt.cex=c(1,.1))
    legend('bottomright',c('soup boundaries','cumulative cut-off','cellranger cut-off','100 Genes','500 Genes','1000 Genes'),lty=c(2,1,1,3,4,5),col=c('black','red','green','black','black','black'),title='Bondaries')
    null = findPeak(umiCounts,plot=TRUE)
    #Restore default layout
    layout(1)
  }
  #Make sure the unambiguously highly expressed ones are included
  wCellHigh = which(umiCounts > cellAbove)
  #Get the big list of the ones to keep
  wCell = unique(c(wCellPois,wCellKnee,wCellCR,wCellHigh))
  #And keep a record of which ones look like crap
  wCellLow = which(geneCounts > minGenesForCell)
  #And those that have excess MT
  wCellMT = which(mtFrac < maxMT)
  goodBarcodes = data.frame(barCode = colnames(srat)[wCell],
                                    passPois = wCell %in% wCellPois,
                                    passKnee = wCell %in% wCellKnee,
                                    passCR = wCell %in% wCellCR,
                                    passHigh = wCell %in% wCellHigh,
                                    passNGenes = wCell %in% wCellLow,
                                    passMT = wCell %in% wCellMT)
  #########################
  # Any further filtering #
  #########################
  #Now save the output
  return(list(srat = srat[,wCell],contMask=contMask,goodBarcodes=goodBarcodes,mtFrac=mtFrac[wCell],umiCounts=umiCounts[wCell],geneCounts=geneCounts[wCell],knee=knee,kneeCR=kneeCR))
}

#################
# Main pipeline #
#################

#' Reload the project from a saved RDS
resumeProject = function(tgt,...){
  return(readRDS(tgt))
}

#' Loads pre-processed data and make Seurat object
#'
#' @param channelDat The channel meta data.
#' @param projectLabel Name for project.
#' @param plotDir Directory where output should live.
#' @param savedRDS The name of the saved RDS.  Will live in plotDir  
#' @param targetRDS Either the name of a column in \code{channelDat} that contains the path to an RDS to load (or save to if not processed) or a common suffix to add to channelDat$path.
#' @param decontaminate If we have a decontamination mask, should we use it to zero out the counts?
#' @param mtCut Any cell with MT expression over this threshold is dropped.
#' @param geneCut Any cell with fewer than this many genes detected is dropped.
#' @param keepPois Should we keep droplets that Poisson method thinks are cells?
#' @param filterFile The name of a file that will be found in the 10X directory that contains the barcodes to filter (in a "Barcode" column) and the reason for filtering them (in a "reason" column).  If not found, nothing is filtered beyond the defaults.
#' @param keepFile A file with a column named "barcode" giving a white list of barcodes to keep for each channel.  If specified, build will fail unless it is found.
#' @param processParams If the pre-processed RDS doesn't exist, it will be calculated by calling \code{preProcessChannel}.  If you want to pass any other arguments to this function, specify them in a list here.
#' @param dontLoad Don't load Seurat object from a saved RDS if present.
#' @param ImFeelingLucky Don't save the output to file, just fly by the seat of your pants.
#' @return Isn't it obvious?
buildSeurat = function(channelDat,projectLabel,plotDir,saveRDS='sratV2.RDS',targetRDS='preProcessedChannel.RDS',decontaminate=FALSE,mtCut=1.0,geneCut=0,keepPois=TRUE,filterFile=NULL,keepFile=NULL,processParams=list(),forceReprocess=FALSE,dontLoad=FALSE,ImFeelingLucky=FALSE){
  if(!dontLoad & file.exists(file.path(plotDir,saveRDS))){
    return(resumeProject(file.path(plotDir,saveRDS)))
  }
  if(!(targetRDS %in% colnames(channelDat))){
    targetRDS = file.path(channelDat$path,targetRDS)
    channelDat$targetRDS = targetRDS
    targetRDS = 'targetRDS'
  }
  toc = list()
  mDat = list()
  fDat = list()
  for(i in seq_len(nrow(channelDat))){
    tgt = channelDat[i,targetRDS]
    if(!file.exists(tgt) | forceReprocess){
      pdat = do.call(preProcessChannel,c(list(channelDat$path[i]),processParams))
      saveRDS(pdat,tgt)
    }else{
      message(sprintf('Loading table of counts for %s',tgt))
      pdat = readRDS(tgt)
    }
    #Now decide what to keep
    #Apply filters first
    dat = pdat$srat
    #Drop the ones with excess MT
    noms = names(which(pdat$mtFrac <= mtCut))
    dat = dat[,colnames(dat)%in% noms,drop=FALSE]
    #Drop the ones with too few genes
    nom = names(which(pdat$geneCounts > geneCut))
    dat = dat[,colnames(dat)%in%nom,drop=FALSE]
    #Throw the ones that are just here from Poisson filter
    if(!keepPois)
      dat = dat[,colnames(dat) %in% pdat$goodBarcodes$barCode[!(pdat$goodBarcodes$passPois & rowSums(pdat$goodBarcodes[,-match('barCode',colnames(pdat$goodBarcodes))])==1)]]
    #Try and find filter files
    filts = NULL
    if(!is.null(filterFile)){
      ffile = file.path(channelDat$path[i],filterFile)
      #If we do have one, drop the crap things
      if(file.exists(ffile)){
        filts = read.delim(ffile,sep='\t',header=TRUE)
        dat = dat[,!(colnames(dat) %in% filts$barcode),drop=FALSE]
      }
    }
    #Look for a white-list file and throw anything not on it
    if(!is.null(keepFile)){
      kfile = file.path(channelDat$path[i],keepFile)
      if(!file.exists(kfile))
        stop(sprintf("White list barcode filter file '%s' not found for directory '%s'.",keepFile,channelDat$path[i]))
      tmp = read.table(kfile,sep='\t',header=TRUE)
      dat = dat[,colnames(dat) %in% tmp$barcode,drop=FALSE]
    }
    #Any other criteria go here
    #Decontaminate?  Only possible if we have a contamination mask
    if(decontaminate & !is.null(pdat$contMask)){
      #Force both to same format
      pdat$srat = as(pdat$srat,'dgTMatrix')
      pdat$contMask = as(pdat$contMask,'dgTMatrix')
      #Copy counts over, but only for the unmasked entries
      dat = pdat$contMask
      m= match(paste0(dat@i,'_',dat@j),paste0(pdat$srat@i,'_',pdat$srat@j))
      dat@x = pdat$srat@x[m]
    }
    #Fix up the column names
    if(i>1){
      if(ncol(dat)>0)
        colnames(dat) = paste0(i,'_',colnames(dat))
      pdat$goodBarcodes$barCode = paste0(i,'_',pdat$goodBarcode$barCode)
    }
    #Add channel level meta.data
    channelDat[i,'knee']=pdat$knee
    channelDat[i,'kneeCR']=pdat$kneeCR
    #Store the output
    toc[[i]]=dat
    #And the meta-data that we have about each barcode
    mDat[[i]] = pdat$goodBarcodes[match(colnames(dat),pdat$goodBarcodes$barCode),]
    #Keep a record of what's been filtered
    fDat[[i]] = filts
  }
  #Now merge and create the Seurat object
  toc = do.call(cbind,toc)
  mDat = do.call(rbind,mDat)
  #Extract the ENSEMBL IDs and Symbols
  ensIDs = gsub('.*_','',rownames(toc))
  symbs = gsub('_.*','',rownames(toc))
  geneDat = getAnnotations(ensIDs)
  #Keep just the first match on ENS ID
  m = match(ensIDs,geneDat$ensembl_gene_id)
  geneDat = geneDat[m,]
  #In case there are NAs, store at least the little info we have
  geneDat$ensembl_gene_id = ensIDs
  #Add in the gene symbol from 10x table
  geneDat$symb = symbs
  rownames(geneDat) = rownames(toc)
  #Now do the construction and normalisation
  srat = CreateSeuratObject(toc,project=projectLabel,min.cells=0,min.genes=0)
   #Initialise the associated meta-data slot
  # geneDat - Meta data table about genes.  Rows of geneDat should correspond exactly to rows of srat@data.
  # channelDat - Meta data table about channels.
  # satDat - Saturation analysis matricies.
  # progressLevel - Marker indicating how far this object has progressed though the pipeline.
  # params - The params that the function was called with when this function was called.
  # maxPC - Resolved maxPC value to use.  That is, this will be set to either the passed value or a calculated one if NA is passed when the time is right.
  # plotDir - Output directory.
  # targetRDS - Name of the target RDS file to store progress.
  srat@misc = list(geneDat = geneDat,channelDat = channelDat, satDat = NULL, progressLevel='init',params=NULL,maxPC=NULL,plotDir=plotDir,targetRDS = file.path(plotDir,saveRDS),filterData=fDat,ImFeelingLucky=ImFeelingLucky,keepFile=keepFile,batchCorrected=FALSE)
  #Add basic metadata
  channelDatIdx = suppressWarnings(as.integer(gsub('_.*','',srat@cell.names)))
  channelDatIdx[is.na(channelDatIdx)] = 1
  tmp = channelDat[channelDatIdx,]
  rownames(tmp) = srat@cell.names
  srat = AddMetaData(srat,metadata=tmp)
  #Add the metadata gathered from pre-processing
  tmp = mDat
  rownames(tmp) = tmp$barCode
  tmp = tmp[srat@cell.names,]
  rownames(tmp) = tmp$barCode
  tmp = tmp[,-match('barCode',colnames(tmp))]
  srat = AddMetaData(srat,metadata=tmp)
  srat = NormalizeData(srat,scale.factor=1e4)
  if(!srat@misc$ImFeelingLucky)
    saveRDS(srat,srat@misc$targetRDS)
  return(srat)
}


#' Replace srat@data with batch regressed data.
#'
#' @param srat The Seurat object.
#' @param batch A vector of length ncol(srat@data) giving the batch to use for each cell.  Can also be a character vector of length 1 in which case the batch is set equal to this meta.data column.
#' @param minCells Don't apply batch correction to genes present in fewer than this many cells.  Helps to reduce the size of the dense matrix that needs to be allocated.
#' @param zeroCut Post batch corrected values smaller than this are considered to be equivalent to zero.
#' @param preserveZeros By default, the batch correction can take a gene that had zero counts and give it a non-zero value.  This will prevent this from happening by forcing these to remain zero.
#' @param reNormalise Rescale data for each cell post batch correction so that it can still be interpreted in the same way as the un-corrected data.
#' @param dontLoad By default will try and load a saved version from the target directory.
#' @param ... Extra parameters passed to ComBat.
batchRegression = function(srat,batch,minCells=3,zeroCut=-Inf,preserveZeros=FALSE,reNormalise=FALSE,regressionMethod=c('combat','mnn'),dontLoad=FALSE,...){
  regressionMethod = match.arg(regressionMethod)
  #Load it and return if it exists
  if(!dontLoad & file.exists(file.path(srat@misc$plotDir,'batchCorrected',basename(srat@misc$targetRDS)))){
    return(resumeProject(file.path(plotDir,'batchCorrected',basename(srat@misc$targetRDS))))
  }
  #Check the batch vector is correct
  if(length(batch)==1 & batch %in% colnames(srat@meta.data))
    batch = srat@meta.data[match(colnames(srat@data),rownames(srat@meta.data)),batch]
  if(length(batch)!=ncol(srat@data))
    stop("Batch must match matrix size")
  #There's no batch regression to do, so don't do any
  if(length(unique(batch))==1){
    srat@misc$batchCorrected=FALSE
    warning("Only one batch given, no regression to do.")
    return(srat)
  }
  #First get the smallest matrix we can
  tmp = srat@data
  gCnts = Matrix::rowSums(tmp>0)
  tmp = as.matrix(tmp[gCnts>minCells,])
  #Keep track of which ones were zero
  w = which(tmp==0)
  #Do the batch regression
  message(sprintf("Performing batch regression with %s",regressionMethod))
  if(regressionMethod=='combat'){
    tmp = sva::ComBat(dat=tmp,batch=batch,...)
  }else{
    #Do correction
    tmp = do.call(scran::mnnCorrect,c(lapply(unique(batch),function(e) tmp[,batch==e,drop=FALSE]),list(...)))
    mnnMatrix = tmp$num.mnn
    #Reconstruct the original matrix
    tmp = do.call(cbind,tmp$corrected)
    tmp = tmp[,colnames(srat@data)]
  }
  #Force them to remain zero post batch correction
  if(preserveZeros)
    tmp[w]=0
  #Set ones close to zero to zero
  tmp[tmp<zeroCut]=0
  #Re-normalise so that we can still interpret the data as log(1+nUMIs) normalised to 10,000 per cell
  if(reNormalise){
    message("Re-normalising data post batch correction.")
    #Get the scale factor
    scaleFac = mean(Matrix::colSums(exp(srat@data)-1))
    #Need to find the correction factor for each cell
    pb = txtProgressBar(min=0,max=ncol(tmp),style=3)
    normFacts=list()
    for(i in seq_len(ncol(tmp))){
      normFacts[[i]] = uniroot(function(e) (sum(exp(tmp[,i]*e)-1)-scaleFac),c(0,10))
      setTxtProgressBar(pb,i)
    }
    close(pb)
    tmp = t(t(tmp)*sapply(normFacts,function(e) e$root))
    #Should now satisfy Matrix::colSums(exp(tmp))==tgt for all cells
  }
  #Adjust the row names to add back in the removed genes
  w = which(tmp>0,arr.ind=TRUE)
  #Now Create a new spares matrix with this as the input
  tmp = sparseMatrix(i=which(gCnts>3)[w[,1]],j=w[,2],x=tmp[w],dims=dim(srat@data),dimnames=list(rownames(srat@data),colnames(srat@data)))
  #Save it back in Seurat object, destroy the existing srat@data, we can always re-make it if needed
  srat@data = tmp
  srat@misc$batchCorrected=TRUE
  #Modify output directory and create extra folder (if it doesn't exist)
  srat@misc$plotDir = file.path(srat@misc$plotDir,'batchCorrected')
  srat@misc$targetRDS = file.path(srat@misc$plotDir,basename(srat@misc$targetRDS))
  dir.create(srat@misc$plotDir,showWarnings=FALSE)
  #Now save the new object 
  if(!srat@misc$ImFeelingLucky)
    saveRDS(srat,srat@misc$targetRDS)
  return(srat)
}


 

calcMetadata = function(srat,mitoCutHigh=0.2,cytoPath='~/scratch/Common/cytoBand.txt',skipCellCycle=TRUE) {
  #Check if we need to do this.
  if(srat@misc$progressLevel!='init')
    return(srat)
  #Get MT expression
  mito.genes = which(srat@misc$geneDat[rownames(srat@raw.data),'chromosome_name']=='MT')
  b=Matrix::colSums(srat@raw.data)
  frac.mito = Matrix::colSums(srat@raw.data[mito.genes,])/b
  df = data.frame(frac.mito,frac.mito>mitoCutHigh,row.names=srat@cell.names)
  colnames(df) = c('frac.mito',sprintf('percent.mito.gt%d',as.integer(mitoCutHigh*100)))
  srat = AddMetaData(srat,metadata=df)
  #Save how far we've processed
  srat@misc$progressLevel='metadata'
  if(!srat@misc$ImFeelingLucky)
    saveRDS(srat,srat@misc$targetRDS)
  return(srat)
}

#' @param markerGenes A list of genes that we will force to be in the variable gene list.  Given as gene symbols.
selectFeatures = function(srat,markerGenes=c(),maxPC=25,storePC=maxPC){
  if(srat@misc$progressLevel!='metadata')
    return(srat)
  message('Feature selection and PCA...')
  pdf(file.path(srat@misc$plotDir,'meanVar.pdf'))
  srat = FindVariableGenes(srat,x.low.cutoff=0.0125,x.high.cutoff=3,y.cutoff=0.5,do.contour=FALSE)
  dev.off()
  if(length(markerGenes)>0){
    #Find the gene name
    m = match(markerGenes,gsub('_.*','',rownames(srat@data)))
    extras = rownames(srat@data)[m[!is.na(m)]]
    #Discard anything with zero expression
    extras = extras[Matrix::rowSums(srat@data[extras,,drop=FALSE])>0]
    #Force them to be in the variable genes
    srat@var.genes = unique(c(extras,srat@var.genes))
  }
  #Scale the data.  If we were to regress anything out, then the model would matter, but as far as I can tell it otherwise just centers and scales the rows to be center 0, sd 1 (the usual PCA pre-processing)
  srat = ScaleData(srat)
  #PCA
  srat = RunPCA(srat,pcs.compute=storePC,weight.by.var=FALSE,do.print=FALSE,maxit=1000)
  pdf(file.path(srat@misc$plotDir,'PCA_lanes.pdf'))
  PCAPlot(srat,1,2,group.by='label')
  dev.off()
  #How many principle components to use
  pdf(file.path(srat@misc$plotDir,'PCA_Elbow.pdf'))
  PCElbowPlot(srat)
  dev.off()
  if(is.na(maxPC)){
    #The jack straw plot takes aggggggggges
    srat = JackStraw(srat,num.pc=storePC,num.replicate=100,do.print=TRUE)
    pdf(file.path(srat@misc$plotDir,'sigPCA.pdf'),width=20,height=14)
    JackStrawPlot(srat,PCs=seq_len(storePC))
    dev.off()
    pAll = srat@jackStraw.empP
    score.thresh=1e-5
    pvals = rep(NA,100)
    for(i in seq(storePC)){
      pvals[i]=suppressWarnings(prop.test(c(length(which(pAll[, i] <= score.thresh)), floor(nrow(pAll) * score.thresh)), c(nrow(pAll), nrow(pAll)))$p.val)
    }
    if(any(pvals<.01)){
      maxPC=max(which(pvals<.01))
    }else{
      #Nothing is significant
      maxPC=0
    }
    #Set bounds on numbers of PCs to use
    if(maxPC<8)
      maxPC=8
    if(maxPC>storePC)
      maxPC=storePC
  }
  #Save maxPC for use on reload
  srat@misc$maxPC=maxPC
  srat@misc$progressLevel='featureSelection'
  if(!srat@misc$ImFeelingLucky)
    saveRDS(srat,srat@misc$targetRDS)
  return(srat)
}


reduceDimensions = function(srat,...){
  if(srat@misc$progressLevel!='featureSelection')
    return(srat)
  message('Running tSNE...')
  #Dimension reduction
  srat = RunTSNE(srat,dims.use=seq_len(srat@misc$maxPC),do.fast=TRUE,check_duplicates=FALSE,...)
  srat@misc$progressLevel='dimReduce'
  if(!srat@misc$ImFeelingLucky)
    saveRDS(srat,srat@misc$targetRDS)
  return(srat)
}

clusterCells = function(srat){
  if(srat@misc$progressLevel!='dimReduce')
    return(srat)
  message('Finding clusters...')
  #pdf(file.path(srat@misc$plotDir,'SNN_plot.pdf'),width=20,height=14)
  srat = FindClusters(srat,dims.use=seq_len(srat@misc$maxPC),resolution=1,print.output=0,save.SNN=TRUE,plot.SNN=FALSE)
  #srat = FindClusters(srat,dims.use=seq_len(srat@misc$maxPC),resolution=10^seq(-1,1,length.out=3),print.output=0,save.SNN=TRUE,plot.SNN=FALSE)
  #dev.off()
  srat@misc$progressLevel='findClusters'
  srat = SetAllIdent(srat,'res.1')
  if(!srat@misc$ImFeelingLucky)
    saveRDS(srat,srat@misc$targetRDS)
  return(srat)
}
  
callMarkerGenes = function(srat){
  if(srat@misc$progressLevel!='findClusters')
    return(srat)
  message("Identifying marker genes...")
  srat@misc$progressLevel='findMarkers'
  srat@misc$markers = tryCatch(FindAllMarkers(srat,test.use='roc',only.pos=TRUE,return.thresh=0.5,min.cells=0),error = function(e) {warning("Marker identification failed.");return(NULL)})
  if(!srat@misc$ImFeelingLucky)
    saveRDS(srat,srat@misc$targetRDS)
  return(srat)
}

makeStandardPlots = function(srat,
                             markerTable,
                             markerCellAnnotVars,
                             markerGeneAnnotVars,
                             markerCellAnnotLabs = markerCellAnnotVars,
                             markerGeneAnnotLabs = markerGeneAnnotVars,
                             markerGenes=unique(markerTable$Gene),
                             keyGenes = c('PTPRC','VCAM1','PECAM1','EPCAM','MKI67','SIX1','SIX2'),
                             maxMarkers = 100
                             ){
  if(!(srat@misc$progressLevel %in% c('findMarkers','findClusters')))
    return(srat)
  #Plot some standard markers
  message('Making standard plots...')
  for(markerGene in markerGenes){
    mark = rownames(srat@misc$geneDat)[match(markerGene,srat@misc$geneDat$symb)]
    if(is.na(mark) | !(mark %in% rownames(srat@data))){
      warning(sprintf("Marker gene %s not present in data, usually this is because it's not expressed.",markerGene))
      next
    }
    #Have to do this twice because Seurat plotting functions suck but I don't care enough to write my own
    pdf(file.path(srat@misc$plotDir,sprintf('expression_%s.pdf',markerGene)))
    #Apply a cut-off at the 99th quantile which truncates everything to the 99th
    df = data.frame(marker=srat@data[mark,],
                    tSNE1=srat@dr$tsne@cell.embeddings[,1],
                    tSNE2=srat@dr$tsne@cell.embeddings[,2])
    df$marker[df$marker>3]=3
    gg = ggplot(df,aes(tSNE1,tSNE2,colour=marker)) +
      geom_point(size=0.5) +
      labs(colour=markerGene) +
      scale_colour_gradientn(colors=c('grey','blue'),limits=c(0,3)) +
      theme(legend.title=element_blank()) +
      ggtitle(mark)
    plot(gg)
    dev.off()
    png(file.path(srat@misc$plotDir,sprintf('expression_%s.png',markerGene)),width=960,height=960)
    plot(gg)
    dev.off()
  }
  ###########
  # Heatmap #
  ###########
  #Check we can actually find all the markers
  w = markerTable$Gene %in% srat@misc$geneDat$symb
  if(!all(w)){
    warning(sprintf("Could not find marker genes (invalid IDs?): %s",paste0('"',markerTable$Gene[!w],'"',collapse=', ')))
  }
  mTable = markerTable[w,]
  #Work out which ones we have data for
  mTable$mark = rownames(srat@misc$geneDat)[match(mTable$Gene,srat@misc$geneDat$symb)]
  m = match(mTable$mark,rownames(srat@data))
  df = srat@data[m[!is.na(m)],,drop=FALSE]
  mTabIdx = which(!is.na(m))
  rownames(df) = seq_len(nrow(df))
  df = as.matrix(df)
  #Normalise each row to the maximum in it so we're not dominated by one highly expressed gene
  #df = df / rowMax(df)
  df[df>2] = 2
  #Add in the zeros
  df = rbind(df,matrix(0,nrow=sum(is.na(m)),ncol=ncol(srat@data)))#,dimnames=list(mTable$mark[is.na(m)],colnames(srat@data))))
  mTabIdx = c(mTabIdx,which(is.na(m)))
  rownames(df) = rownames(mTable)[mTabIdx]
  colnames(df) = colnames(srat@data)
  #Plot using pheatmap
  annotCell = srat@meta.data[colnames(df),markerCellAnnotVars,drop=FALSE]
  colnames(annotCell) = markerCellAnnotLabs
  for(nom in colnames(annotCell))
    annotCell[,nom] = factor(annotCell[,nom])
  annotGene = mTable[rownames(df),markerGeneAnnotVars,drop=FALSE]
  for(nom in colnames(annotGene))
    annotGene[,nom] = factor(annotGene[,nom])
  colnames(annotGene) = markerGeneAnnotLabs
  #Re-order cells by cluster
  o = do.call(order,annotCell[,markerCellAnnotLabs,drop=FALSE])
  df = df[,o]
  #And by gene groupings
  #Create another sorter which is the mean expression
  o = do.call(order,cbind(annotGene[,markerGeneAnnotLabs,drop=FALSE],-rowMeans(df)))
  df = df[o,]
  #Get default ggplot2 discrete colours (taken from here https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette)
  ggplotDefaultCols = function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[seq_len(n)]
  }
  #Make colour scheme for annotation
  annotCols = list()
  if('Cluster' %in% markerCellAnnotLabs)
    annotCols[['Cluster']] = setNames(ggplotDefaultCols(length(levels(annotCell[,'Cluster']))),levels(annotCell[,'Cluster']))
  #Set all colour schemes manually
  #for(nom in markerCellAnnotLabs){
  #  annotCols[[nom]] = setNames(ggplotDefaultCols(length(levels(annotCell[,nom]))),levels(annotCell[,nom]))
  #}
  #for(nom in markerGeneAnnotLabs){
  #  annotCols[[nom]] = setNames(ggplotDefaultCols(length(levels(annotGene[,nom]))),levels(annotGene[,nom]))
  #}
  colScheme = colorRampPalette(c('grey','blue'))(100)
  pheatmap(df,
                color = colScheme,
                cluster_rows=FALSE,
                cluster_cols=FALSE,
                legend=FALSE,
                annotation_row = annotGene,
                annotation_col = annotCell,
                annotation_legend = TRUE,
                annotation_names_row = TRUE,
                annotation_names_col = TRUE,
                drop_levels = FALSE,
                show_rownames = TRUE,
                show_colnames = FALSE,
                main = "Cell type overview",
                display_numbers = FALSE,
                labels_row = mTable[rownames(df),'Gene'],
                gaps_row = cumsum(Rle(annotGene[rownames(df),markerGeneAnnotLabs[1]])@lengths),
                gaps_col = cumsum(Rle(annotCell[colnames(df),markerCellAnnotLabs[1]])@lengths),
                annotation_colors = annotCols,
                filename=file.path(srat@misc$plotDir,'markerHeatmap.pdf'),
                width=21,
                height=21)
  #Make a second plot where we drop the low variance genes
  tmp = annotCell[colnames(df),markerCellAnnotLabs[1]]
  fracOn = apply(df,1,function(e) sapply(split(e,tmp),function(ee) sum(ee>0)/length(ee)))
  maxFracOn = apply(fracOn,2,max)
  rdf = df[maxFracOn>0.1,]
  pheatmap(rdf,
                color = colScheme,
                cluster_rows=FALSE,
                cluster_cols=FALSE,
                legend=FALSE,
                annotation_row = annotGene,
                annotation_col = annotCell,
                annotation_legend = TRUE,
                annotation_names_row = TRUE,
                annotation_names_col = TRUE,
                drop_levels = FALSE,
                show_rownames = TRUE,
                show_colnames = FALSE,
                main = "Cell type overview",
                display_numbers = FALSE,
                labels_row = mTable[rownames(rdf),'Gene'],
                gaps_row = cumsum(Rle(annotGene[rownames(rdf),markerGeneAnnotLabs[1]])@lengths),
                gaps_col = cumsum(Rle(annotCell[colnames(rdf),markerCellAnnotLabs[1]])@lengths),
                annotation_colors = annotCols,
                filename=file.path(srat@misc$plotDir,'markerHeatmapSmall.pdf'),
                width=21,
                height=21)
  ####Add in the zeros
  ###df = rbind(df,matrix(0,nrow=sum(is.na(m)),ncol=ncol(srat@data),dimnames=list(mTable$mark[is.na(m)],colnames(srat@data))))
  ####Melt it
  ###df = melt(df)
  ###noms = c('Gene','Cell','Expression')
  ###colnames(df) = noms
  ###m = match(df$Gene,mTable$mark)
  ####Add in anything extra in the mTable, but discard Gene, Cell and Expression
  ###df = cbind(df,mTable[m,!(colnames(mTable)%in%noms)])
  ####And the clusters
  ###df$Cluster = srat@ident[df$Cell]
  ####Convert gene names back to mTable values
  ###df$Gene = gsub('_.*','',df$Gene)
  #####Order Cells by cluster
  ####tmp = unique(as.character(df$Cell))
  ####df$Cell = factor(df$Cell,levels = tmp[order(srat@ident[tmp],tmp)])
  #####Order genes by annotation
  ####tmp = unique(as.character(df$Gene))
  ####df$Gene = factor(df$Gene,levels = tmp[order(mTable$Organ[match(tmp,mTable$Gene)])])
  ####Create the main plot
  ###markerGroupVar='Organ'
  ###cellAnnotVar='label'
  ###df$cellAnnot = srat@meta.data[df$Cell,cellAnnotVar]
  ###gg = ggplot(df,aes(Cell,Gene,fill=Expression)) + 
  ###  geom_tile() + 
  ###  facet_grid(sprintf('%s ~ Cluster',markerGroupVar),space='free',scales='free') +
  ###  guides(fill=FALSE) +
  ###  theme(panel.spacing.x=unit(0.1, "lines"),
  ###        panel.spacing.y=unit(0.1, "lines"),
  ###        axis.ticks.x=element_blank(),
  ###        axis.text.x=element_blank(),
  ###        axis.title.x=element_blank(),
  ###        strip.background=element_blank(),
  ###        strip.text.y=element_text(angle=0))
  ####Add extra annotation for cells
  ###df$cellAnnotLabel = cellAnnotVar
  ###gg_x = ggplot(df,aes(Cell,cellAnnotLabel,fill=cellAnnot))+
  ###  geom_tile() +
  ###  facet_grid(. ~ Cluster, space = 'free',scales='free') +
  ###  theme(panel.spacing.x=unit(0.1, "lines"),
  ###        panel.spacing.y=unit(0.1, "lines"),
  ###        strip.background = element_blank(),
  ###        strip.text=element_blank(),
  ###        axis.ticks.x=element_blank(),
  ###        axis.text.x=element_blank(),
  ###        axis.title.y=element_blank(),
  ###        legend.position='bottom')
  ####And align them
  ###g1 = ggplotGrob(gg)
  ###g2 = ggplotGrob(gg_x)
  ###colnames(g1) <- paste0(seq_len(ncol(g1)))
  ###colnames(g2) <- paste0(seq_len(ncol(g2)))
  ###pdf(file.path(srat@misc$plotDir,'markerHeatmap.pdf'))
  ###grid.draw(combine(g1, g2, along=2))
  ###dev.off()
  ###############
  #Plot some QC things
  #nUMI - log10
  srat@meta.data$log10_nUMI = log10(srat@meta.data$nUMI)
  pdf(file.path(srat@misc$plotDir,'log10_nUMI.pdf'))
  FeaturePlot(srat,'log10_nUMI',no.legend=FALSE,pt.size=0.5,cols.use=c('grey','blue'))
  dev.off()
  #Number of genes
  pdf(file.path(srat@misc$plotDir,'nGene.pdf'))
  FeaturePlot(srat,'nGene',no.legend=FALSE,pt.size=0.5,cols.use=c('grey','blue'))
  dev.off()
  #Number of genes - log10
  srat@meta.data$log10_nGene = log10(srat@meta.data$nGene)
  pdf(file.path(srat@misc$plotDir,'log10_nGene.pdf'))
  FeaturePlot(srat,'log10_nGene',no.legend=FALSE,pt.size=0.5,cols.use=c('grey','blue'))
  dev.off()
  #Pass based on CR cut-off
  pdf(file.path(srat@misc$plotDir,'passCR.pdf'))
  TSNEPlot(srat,group.by='passCR',pt.size=0.5)
  dev.off()
  #Cell cycle
  if('cell.cycle' %in% colnames(srat@meta.data)){
    pdf(file.path(srat@misc$plotDir,'cellCycle.pdf'))
    TSNEPlot(srat,group.by='cell.cycle',pt.size=0.5)
    dev.off()
  }
  ##The resolution = 1 clustering.
  pdf(file.path(srat@misc$plotDir,'clusters.res.1.pdf'))
  TSNEPlot(srat,group.by='res.1',do.label=TRUE,pt.size=0.5,no.legend=TRUE)
  dev.off()
  ##And the channels
  pdf(file.path(srat@misc$plotDir,'channels.pdf'),width=10.5,height=7)
  TSNEPlot(srat,group.by='channel',do.label=FALSE,pt.size=0.5,no.legend=FALSE)
  dev.off()
  #And labels
  pdf(file.path(srat@misc$plotDir,'labels.pdf'),width=10.5,height=7)
  TSNEPlot(srat,group.by='label',do.label=FALSE,pt.size=0.5,no.legend=FALSE)
  dev.off()
  #Do each of the annotation colums as a different plot
  for(var in markerCellAnnotVars){
    pdf(file.path(srat@misc$plotDir,sprintf("cellLabels_%s.pdf",markerCellAnnotLabs[match(var,markerCellAnnotVars)])))
    TSNEPlot(srat,group.by=var,do.label=FALSE,pt.size=0.5,no.legend=FALSE)
    dev.off()
  }
  #Mitochondrial expression
  pdf(file.path(srat@misc$plotDir,'mitoFrac.pdf'))
  FeaturePlot(srat,'frac.mito',no.legend=FALSE,pt.size=0.5,cols.use=c('grey','blue'))
  dev.off()
  #Boxplot showing mitochondrial expression fraction by label
  pdf(file.path(srat@misc$plotDir,'mitoFracByLabel.pdf'),width=14)
  gg = ggplot(srat@meta.data,aes(label,frac.mito)) +
    geom_boxplot() +
    stat_summary(fun.data = function(e) c(y=-0.02,label=length(e)),geom='text',fun.y=median) +
    stat_summary(fun.data = function(e) c(y=1.02,label=round(mean(e),2)),geom='text',fun.y=mean) +
    coord_flip() +
    xlab('Fraction of expression due to MT') +
    ylab('Label')+
    ggtitle('Mitochondrial expression (left#=cellNo,right#=medianMT)')
  plot(gg)
  dev.off()
  #Boxplot showing mitochondrial expression fraction by cluster
  pdf(file.path(srat@misc$plotDir,'mitoFracByCluster.pdf'),width=10)
  gg = ggplot(srat@meta.data,aes(res.1,frac.mito)) +
    geom_boxplot() +
    stat_summary(fun.data = function(e) c(y=-0.02,label=length(e)),geom='text',fun.y=median) +
    stat_summary(fun.data = function(e) c(y=1.02,label=round(mean(e),2)),geom='text',fun.y=mean) +
    coord_flip() +
    xlab('Fraction of expression due to MT') +
    ylab('Cluster#')+
    ggtitle('Mitochondrial expression (left#=cellNo,right#=medianMT)')
  plot(gg)
  dev.off()
  #Show the saturation curve for each cluster
  if(!is.null(srat@misc$satDat)){
    xx = sapply(split(names(srat@ident),srat@ident),function(e) apply(srat@misc$satDat$detectedReads[,e],1,median))
    yy = sapply(split(names(srat@ident),srat@ident),function(e) apply(srat@misc$satDat$detectedUMIs[,e],1,median))
    xx = melt(xx)
    yy = melt(yy)
    if(any(xx$Var1!=yy$Var1) | any(xx$Var2!=yy$Var2)){
      warning('Something went wrong with saturation curve building.  Skipping...')
    }else{
      df = data.frame(nReads=xx$value,nUMIs=yy$value,subSampleRate = srat@misc$satDat$subSampleRates[xx$Var1],cluster=factor(xx$Var2))
      #Plot the summary curve as a function of rates
      pdf(file.path(srat@misc$plotDir,'clusterSaturationByRate.pdf'))
      gg = ggplot(df,aes(subSampleRate,nUMIs,colour=cluster)) +
        geom_point() +
        geom_line() +
        ggtitle(sprintf('Sequencing saturation for each cluster')) +
        xlab('Fraction of sub-sampled reads') +
        ylab('Median number of detected genes in cluster')
      plot(gg)
      dev.off()
      #Plot the summary curve as a function of reads
      pdf(file.path(srat@misc$plotDir,'clusterSaturationByReads.pdf'))
      gg = ggplot(df,aes(nReads,nUMIs,colour=cluster)) +
        geom_point() +
        geom_line() +
        scale_x_log10() +
        scale_y_log10() +
        ggtitle(sprintf('Sequencing saturation for each cluster')) +
        xlab('Median number of reads in cluster') +
        ylab('Median number of detected genes in cluster')
      plot(gg)
      dev.off()
    }
  }
  #Summarise each cluster
  message('Creating cluster summaries...')
  srat@misc$clusterSummaries=list()
  for(cl in levels(srat@ident)){
    tmp = summariseCluster(srat,cl,sortGenes=keyGenes,outFile=file.path(srat@misc$plotDir,sprintf('summaryOfCluster%s.txt',cl)),maxWriteMarkers=maxMarkers,markerGenes=markerGenes)
    srat@misc$clusterSummaries[[cl]] = tmp
    saveRDS(tmp,file.path(srat@misc$plotDir,sprintf('summaryOfCluster%s.RDS',cl)))
  }
  #Make a summary of the markers and their expression in each cluster
  sharedMarkers = unique(unlist(lapply(srat@misc$clusterSummaries,function(e) rownames(e$markerGeneDat)[seq(min(nrow(e$markerGeneDat),50))]),use.names=FALSE))
  sharedMarkerTable = matrix(NA,nrow=length(sharedMarkers),ncol=length(levels(srat@ident)),dimnames=list(sharedMarkers,levels(srat@ident)))
  for(cl in levels(srat@ident)){
    cluster = names(srat@ident)[srat@ident %in% cl]
    sharedMarkerTable[,cl] = Matrix::rowSums(srat@raw.data[sharedMarkers,cluster,drop=FALSE]>0)/length(cluster)
  }
  #Add in the marker table
  tmp = sapply(srat@misc$clusterSummaries,function(e) e$markerGeneSummary$fracExpr)
  rownames(tmp) = rownames(srat@misc$clusterSummaries[[1]]$markerGeneSummary)
  #Record how many shared we have
  nShared=nrow(sharedMarkerTable)
  sharedMarkerTable = rbind(sharedMarkerTable,tmp[!(rownames(tmp)%in%rownames(sharedMarkerTable)),colnames(sharedMarkerTable)])
  #Add a column to record why it's here
  sharedMarkerTable = cbind(sharedMarkerTable,isSpecialMarker=rep(c(FALSE,TRUE),c(nShared,nrow(sharedMarkerTable)-nShared)))
  #Save it to file and to the R object
  srat@misc$sharedMarkerTable=sharedMarkerTable
  write.table(sharedMarkerTable,file.path(srat@misc$plotDir,'sharedMarkerTable.tsv'),quote=FALSE,sep='\t')
  #Finally, save a bunch of things
  saveRDS(srat,srat@misc$targetRDS)
  saveRDS(srat@data,file.path(srat@misc$plotDir,'tableOfCounts.RDS'))
  #Add tSNE to meta.data before saving
  srat@meta.data$tSNE1 = srat@dr$tsne@cell.embeddings[rownames(srat@meta.data),'tSNE_1']
  srat@meta.data$tSNE2 = srat@dr$tsne@cell.embeddings[rownames(srat@meta.data),'tSNE_2']
  saveRDS(list(meta.data=srat@meta.data,
               geneDat=srat@misc$geneDat,
               channelDat=srat@misc$channelDat,
               passedBarcodes = srat@misc$passedBarcodes,
               maxPC = srat@misc$maxPC,
               filterData =srat@misc$filterData,
               sharedMarkerTable = srat@misc$sharedMarkerTable,
               clusterSummaries = srat@misc$clusterSummaries),
               file.path(srat@misc$plotDir,'metadata.RDS'))
  saveRDS(srat@dr,file.path(srat@misc$plotDir,'dimensionReduction.RDS'))
  #All done!
  return(srat)
}
