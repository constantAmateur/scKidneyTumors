#' A module providing functions needed to genotype 10X cells from bulk sequencing data.
#' Requires AlleleCount (https://github.com/cancerit/alleleCount) to be installed
import(genome)
fix_chr = genome$fix_chr
import(sanger)
load_VCF = sanger$load_VCF
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
 
#' Runs allele Counter in parallel over a set of input and output files
parallelAlleleCounter = function(inputs,outputs,bams,nParallel=8,f=0,F=0,x=TRUE,d=TRUE,m=20,q=200,r='~/canpipe/ref/human/1000Genomes_hs37d5/genome.fa',bin='~/bin/alleleCounter',skipIfExists=FALSE){
  #Construct the base string
  baseCommand = sprintf('%s -l %%s -b %%s -o %%s -r %s -f %d -F %d -m %d -q %d',bin,r,f,F,m,q)
  if(x)
    baseCommand = paste0(baseCommand,' -x')
  if(d)
    baseCommand = paste0(baseCommand,' -d')
  if(length(inputs)==1)
    inputs=rep(inputs,length(bams))
  if(length(inputs)!=length(outputs) | length(inputs)!=length(bams))
    stop("inputs,outputs, and bams must all have the same length.")
  toRun = sprintf(baseCommand,inputs,bams,outputs)
  #Filter out anything where the output file exists
  if(skipIfExists)
    toRun = toRun[!file.exists(outputs)]
  #Run the commands...
  system(sprintf("parallel -j %d -- %s",nParallel,paste0('"',toRun,'"',collapse=' ')))
}

#' Run allele counters over every position in \code{grange} for the input files given and load the result.
#' 
#' @param outString The output filenames will be constructed using sprintf applied to this string.  The first entry is replaced by \code{inputFile} and the second by the bam names.
alleleCounterFromGRange = function(grange,bams,inputFile,outString='%s.%s.tsv',...){
  #Write out the positions to calculate
  write.table(cbind(as.character(seqnames(grange)),start(grange)),inputFile,sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
  #Get BAM names or make them
  if(is.null(names(bams))){
    noms = paste0('BAM',seq_along(bams))
  }else{
    noms = names(bams)
  }
  #Run allele counter on the BAMs
  outFiles = sprintf(outString,inputFile,noms)
  parallelAlleleCounter(inputFile,outFiles,bams,...)
  #Load each result
  out=list()
  for(i in seq_along(noms)){
    x = read.delim(outFiles[i],sep='\t',header=TRUE)
    #Get rid of suffix on barcodes
    x$Barcode = gsub('-[0-9]+','',x$Barcode)
    #Make globally unique reference
    x$DropletID = paste(rep(noms[i],nrow(x)),x$Barcode,sep='___')
    x$SampleIdx = rep(i,nrow(x))
    out[[i]]=x
  }
  out = do.call(rbind,out)
  #Store the names of the BAM files used
  out$BAM = bams[out$SampleIdx]
  #Fix column names
  noms = colnames(out)
  noms[noms=='X.CHR']='Chr'
  noms[noms=='POS']='Pos'
  colnames(out)=noms
  #Geth extra meta-data from grange if present
  m = match(paste0(out$Chr,':',out$Pos),paste0(as.character(seqnames(grange)),':',start(grange)))
  x = which(!(colnames(mcols(grange)) %in% colnames(out)))
  if(length(x)>0)
    out = cbind(out,as.data.frame(mcols(grange)[m,x]))
  #If we have REF and ALT, produce counts for these
  if(all(c('REF','ALT') %in% colnames(out))){
    out$RefCount = as.numeric(out[cbind(nrow=seq_len(nrow(out)),ncol=match(sprintf('Count_%s',out$REF),colnames(out)))])
    out$AltCount = as.numeric(out[cbind(nrow=seq_len(nrow(out)),ncol=match(sprintf('Count_%s',out$ALT),colnames(out)))])
    out$otherCount = out$Good_depth-out$RefCount-out$AltCount
  }
  return(out)
}
 
#' General summary function
#' 
#' @param aggFunctions  A list of lists.  Each sub-list should have two entries, the first giving the columns to aggregate in this way and the second giving the aggregation function.
summariseByMark = function(out,marks,dropThresh=1,aggFunctions=list()){
  gmark = apply(out[,marks,drop=FALSE],1,paste,collapse=':',sep='')
  tCnts = data.frame(Mark = unique(gmark))
  for(nom in colnames(out)){
    #Test if this column should be specially treated
    processed=FALSE
    for(i in seq_along(aggFunctions)){
      if(nom %in% aggFunctions[[i]][[1]]){
        tmp  = sapply(split(out[,nom],gmark),aggFunctions[[i]][[2]])
        processed=TRUE
        break
      }
    }
    #Decide what to do with those we don't have a special processing strategy for
    if(!processed){
      tmp  = split(out[,nom],gmark)
      nEntries = lengths(lapply(tmp,unique))
      if(any(nEntries>dropThresh))
        next
      #Aggregate and keep
      if(all(nEntries==1)){
        tmp = sapply(tmp,unique)
      }else{
        tmp =sapply(tmp,function(e) paste(unique(e),collapse=', '))
      }
    }
    #Store the summarised result
    if(length(tmp)>0)
      tCnts[match(names(tmp),tCnts$Mark),nom]=tmp
  }
  rownames(tCnts)=NULL
  #Return the result
  return(tCnts)
}

#' Special summary function for CN data.  Estimates aggregated Major Allele Frequencies and confidence limits.  Also calculates p-values for difference from allelic balance.
#' @param conf.level What confidence range to create.
#' @param minCoverage Don't test any location with coverage less than this value.
summariseCN = function(out,marks,conf.level=0.95,minCoverage=10){
  aggFuns = list(sum=list(c('Good_depth','patCount','matCount','RefCount','AltCount','otherCount'),sum),
                 median = list(c('normAF','mutAF'),median),
                 string = list(c('Chr','ClusterID','regionType'),function(e) paste(unique(e),collapse=', '))
                 )
  if(nrow(out)==0)
    return(NULL)
  #Stuff for confidence intervals
  p.L = function(x,n,alpha){
    ifelse(x==0,0,qbeta(alpha,x,n-x+1))
  }
  p.U = function(x,n,alpha){
    ifelse(x==n,1,qbeta(1-alpha,x+1,n-x))
  }
  alpha = (1-conf.level)/2
  #Do the main summarisation
  tCnts = summariseByMark(out,marks,dropThresh=1,aggFuns)
  #Make special extra columns
  tmarks = apply(out[,marks,drop=FALSE],1,paste,collapse=':',sep='')
  tCnts$numSNPs = sapply(split(paste(out$DropletID,out$Chr,out$Pos,sep=':'),tmarks),function(e) length(unique(e)))[tCnts$Mark]
  tCnts$numCells = sapply(split(out$DropletID,tmarks),function(e) length(unique(e)))[tCnts$Mark]
  rownames(tCnts)=tCnts$Mark
  tCnts$MAF = tCnts$matCount/tCnts$Good_depth
  tCnts$MAF_Low = p.L(tCnts$matCount,tCnts$Good_depth,alpha)
  tCnts$MAF_High = p.U(tCnts$matCount,tCnts$Good_depth,alpha)
  tCnts$pval = pmin(pbinom(tCnts$matCount-1,tCnts$Good_depth,0.5,lower.tail=FALSE),pbinom(tCnts$patCount-1,tCnts$Good_depth,0.5,lower.tail=FALSE))
  #Do restricted multiple hypothesis test
  w = which(tCnts$Good_depth>minCoverage)
  tCnts$qval = 1
  tCnts$qval[w] = p.adjust(tCnts$pval[w],method='BH')
  #Re-order by p-value
  tCnts = tCnts[order(tCnts$pval),]
  return(tCnts)
}

#' Add information about which gene and what type of region a SNP falls in
annotateByRegion = function(snps){
  gns = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
  gns = renameSeqlevels(gns,setNames(fix_chr(seqlevels(gns)),seqlevels(gns)))
  exns = exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
  exns = reduce(GRanges(fix_chr(seqnames(exns)),ranges(exns)))
  intrns = intronsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene)
  intrns = unlist(intrns)
  intrns = reduce(GRanges(fix_chr(seqnames(intrns)),ranges(intrns)))
  #Get the type of each SNP
  snps$regionType = 'InterGene'
  o = suppressWarnings(findOverlaps(snps,gns,ignore.strand=TRUE))
  snps$regionType[queryHits(o)]='Gene'
  #Also record which gene this is
  oo = split(subjectHits(o),queryHits(o))
  oo  = sapply(oo,function(e) paste0(gns$gene_id[e],collapse=','))
  snps$geneID = ''
  snps$geneID[as.numeric(names(oo))] = as.character(oo)
  o = suppressWarnings(findOverlaps(snps,exns))
  snps$regionType[queryHits(o)]='Exon'
  o = suppressWarnings(findOverlaps(snps,intrns))
  snps$regionType[queryHits(o)]='Intron'
  return(snps)
}

#' Gets the SNPs that are informative for the CN in the specified sample.
#'
#' Given a set of segments where a CN change with unbalanced alleles is known to occur, finds those SNPs that are heterozygous in the germline and so will be informative of the CN state.  These SNPs are phased by the BAF in the tumour.
#'
#' @param fn_snps Caveman SNP file for the individual with the CN change.
#' @param segments ASCAT style summary output (with first column removed) entered as as string.  This will define the segments where the CN change occurs and what that CN change is.
#' @param HetCut In order to be considered heterozygous in the normal the AF in the normal must be within this value of 0.5.
#' @param qCut Adjusted p-value cut-off for selecting SNPs that are have BAF in the tumour above/below 0.5 and so should be phased together.
#' @return A GRanges object with the informative SNPs together with some useful annotations.
getCN_SNPs = function(fn_snps,segments,HetCut=0.1,qCut=0.1){
  snps = load_VCF(fn_snps,verbose=FALSE)
  segments = as.data.frame(t(sapply(lapply(strsplit(segments,'\n')[[1]],strsplit,','),`[[`,1)))
  segments = GRanges(segments[,1],IRanges(as.numeric(segments[,2]),as.numeric(segments[,3])),total=as.integer(segments[,6]),minor=as.integer(segments[,7]))
  segments$MAF = (segments$total-segments$minor)/segments$total
  snps = subsetByOverlaps(snps,segments)
  #Get those that are HET in the normal
  snps = snps[which(abs(snps$NORM_AF-0.5) < HetCut)]
  #And that we can un-ambiguously call as being enriched/depleted
  p = pmin(pbinom(snps$MUT_ALT,snps$MUT,0.5),pbinom(snps$MUT_ALT-1,snps$MUT,0.5,lower.tail=FALSE))
  q = p.adjust(p,method='BH')
  good = snps[q<qCut]
  good$maternalAllele = ifelse(good$MUT_AF>0.5,good$ALT,good$REF)
  good$paternalAllele = ifelse(good$MUT_AF<0.5,good$ALT,good$REF)
  #Annotate the SNPs
  good = annotateByRegion(good)
  o = findOverlaps(good,segments)
  good$SegmentID[queryHits(o)] = subjectHits(o)
  good$BulkMAF[queryHits(o)] = segments$MAF[subjectHits(o)]
  return(good)
}

#' Look for presence of mutations given by subFile VCF in bams.
#' @param clonalAF Splits mutations into sub-clonal and clonal based on which side of this Allele frequency they lie.
genoTypeBySubs = function(subFile,bams,inputFile,metaData=NULL,clusterIDs=NULL,clonalAF=0.2,...){
  subs = load_VCF(subFile)
  subs = subs[subs$FILTER=='PASS']
  #Get the counts at the subs
  cnts = alleleCounterFromGRange(subs,bams,inputFile,d=FALSE,...)
  #Now add in extra meta-data
  x = which(!(colnames(metaData) %in% colnames(out)))
  if(length(x)>0)
    cnts = cbind(cnts,metaData[cnts$SampleIdx,x])
  if(!is.null(clusterIDs))
    cnts$Cluster = clusterIDs[cnts$DropletID]
  cnts$ClusterID = gsub('.*___','',cnts$Cluster)
  cnts$Experiment = gsub('___.*','',cnts$Cluster)
  #cnts$Experiment = gsub('.*_','',cnts$Experiment)
  #Designate a Sub as clonal or not
  cnts$Clonal = cnts$MUT_AF> clonalAF
  #Split by expriment
  aggFuns = list(sum=list(c('RefCount','AltCount'),sum),
                 median=list(c('MUT_AF'),median))
  cCnts = summariseByMark(cnts,c('Experiment','ClusterID','Clonal'),aggFunctions=aggFuns)
  sCnts = summariseByMark(cnts,c('Chr','Pos'),aggFunctions=aggFuns)
  #Add in meta-data
  if(!is.null(metaData)){
    cCnts = cbind(cCnts,metaData[cCnts$SampleIdx,])
    sCnts = cbind(sCnts,metaData[sCnts$SampleIdx,])
  }
  return(list(cnts=cnts,cCnts=cCnts,sCnts=sCnts,bams=bams))
}


genoTypeByCN = function(snps,bams,inputFile,clusterIDs=NULL,excludeExperiments=c(),excludeRegions=c('InterGene'),pvalCut=0.05,MAFcut=0.2,...){
  out = alleleCounterFromGRange(snps,bams,inputFile,...)
  out$matCount = as.numeric(out[cbind(nrow=seq_len(nrow(out)),ncol=match(sprintf('Count_%s',out$maternalAllele),colnames(out)))])
  out$patCount = as.numeric(out[cbind(nrow=seq_len(nrow(out)),ncol=match(sprintf('Count_%s',out$paternalAllele),colnames(out)))])
  #Populate the cluster ID column if we have that information
  if(!is.null(clusterIDs)){
    out$Cluster=clusterIDs[out$DropletID]
  }
  #Get experiment and clusterID
  out$ClusterID = gsub('.*___','',out$Cluster)
  out$Experiment = gsub('___.*','',out$Cluster)
  #out$Experiment = gsub('.*_','',out$Experiment)
  #Filter out any regions we want to ignore
  out = out[!(out$regionType %in% excludeRegions),]
  #Aggregate at SNPs across cells we expect to be mostly CN neutral.  Anything that has a bias to one allele should be excluded as it must either be germline homozygous or have strong ASE.
  badSNPs = summariseCN(out[!(out$Experiment %in% excludeExperiments),],c('Chr','Pos'))
  #Discard only those that are pretty bad
  out = out[!(paste0(out$Chr,':',out$Pos) %in% rownames(badSNPs)[badSNPs$qval<pvalCut & abs(badSNPs$MAF-0.5)>MAFcut]),]
  #Aggregate across genes to find dodgy genes
  badGenes = summariseCN(out[!(out$Experiment %in% excludeExperiments),],'geneID')
  #Drop the bad ones again
  out = out[!(out$geneID %in% badGenes$geneID[badGenes$qval<pvalCut & abs(badGenes$MAF-0.5)>MAFcut]),]
  #Create different summaries 
  gCnts = summariseCN(out,c('DropletID','Chr','SegmentID','geneID'))
  cCnts = summariseCN(out,c('DropletID','Chr','SegmentID'))
  sgCnts = summariseCN(out,c('geneID','Experiment'))
  clCnts = summariseCN(out,c('Cluster','Chr','SegmentID'))
  return(list(cnts=out,
              gCnts=gCnts,
              cCnts=cCnts,
              clCnts=clCnts,
              sgCnts=sgCnts,
              badGenes = badGenes,
              badSNPs = badSNPs))
}

plotGenotypedCN = function(cCnts,clCnts,sgCnts,baseName){
  #How large a gap to leave between things grouped by chrosomoe
  gap = 0.4
  ##########################################################
  # Plot individual cells by cluster.  Split by experiment
  gg_cells=list()
  for(experiment in unique(cCnts$Experiment)){
    if(is.na(experiment)){
      tmp = cCnts[is.na(cCnts$Experiment),]
      lab = 'NA'
    }else{
      tmp = cCnts[cCnts$Experiment==experiment & !is.na(cCnts$Experiment),]
      lab = experiment
    }
    #Order by average p-value across the chromosomes
    x = sapply(split(tmp$pval,tmp$DropletID),mean)
    x = x[order(tmp$ClusterID[match(names(x),tmp$DropletID)],x)]
    x = match(tmp$DropletID,names(x))
    #Add offset for chromosome
    offset = seq(gap/2,1-gap/2,length.out=length(unique(tmp$Chr)))[as.numeric(factor(tmp$Chr))]-0.5
    tmp$x = x+offset
    tmp$ClusterID = factor(tmp$ClusterID,levels=unique(tmp$ClusterID)[order(suppressWarnings(as.numeric(unique(tmp$ClusterID))))])
    gg_cells[[lab]] = ggplot(tmp,aes(x,MAF,colour=factor(Chr),shape=qval<.1)) +
      geom_point() + 
      geom_errorbar(aes(ymin=MAF_Low,ymax=MAF_High),alpha=0.5) +
      #geom_line(aes(x,BulkMAF),inherit.aes=FALSE,colour='black',linetype=2) +
      facet_wrap(~ClusterID,scales='free_x') +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
            #panel.grid.minor.x = element_line(colour='black',size=0.1)) +
      scale_x_continuous(breaks=sort(unique(x)),minor_breaks = c(0,sort(unique(x)))+0.5) +
      ylim(0,1) + 
      geom_hline(yintercept=0.5) +
      ggtitle(sprintf("CN over-view at cell level for experiment %s",experiment))
  }
  ####################################
  # Plot cells grouped into clusters
  #Make custom x-coord to put different chromosomes next to one another
  o = order(suppressWarnings(as.numeric(clCnts$ClusterID)))
  x = match(clCnts$ClusterID,unique(clCnts$ClusterID[o]))
  labs = unique(clCnts$ClusterID[o])
  offset = seq(0+gap/2,1-gap/2,length.out=length(unique(clCnts$Chr)))[as.numeric(factor(clCnts$Chr))]-0.5
  if(length(unique(clCnts$Chr))==1)
    offset = rep(0,length(offset))
  clCnts$x = x+offset
  gg_cl = ggplot(clCnts,aes(x,MAF,colour=factor(Chr),shape=qval<.1)) +
    geom_point() +
    #geom_line(aes(x,BulkMAF),inherit.aes=FALSE,colour='black',linetype=2) +
    geom_errorbar(aes(ymin=MAF_Low,ymax=MAF_High),alpha=1.0) + 
    scale_x_continuous(breaks=seq(length(unique(clCnts$ClusterID))),labels=labs,minor_breaks = seq(0,length(unique(clCnts$ClusterID)))+0.5) +
    facet_wrap(~Experiment,scales='free_x') +
    ylim(0,1) +
    geom_hline(yintercept=0.5) +
    theme(axis.text.x = element_text(angle=90,hjust=1),panel.grid.minor.x = element_line(colour='black',size=0.1)) +
    ggtitle("CN over-view at cluster level")
  ###########################################################
  # Plot genes aggregated across all cells in an experiment
  sgCnts$x = seq_len(nrow(sgCnts))
  gg_genes = ggplot(sgCnts,aes(x,MAF,colour=factor(Chr),shape=qval<.1)) +
    geom_point() +
    facet_wrap(~Experiment) +
    geom_errorbar(aes(ymin=MAF_Low,ymax=MAF_High),alpha=0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
    ylim(0,1) +
    geom_hline(yintercept=0.5) +
    ggtitle("CN over-view at gene level")
  pdf(sprintf('%s_cells.pdf',baseName),width=28,height=35)
  for(gg in gg_cells){
    plot(gg)
  }
  dev.off()
  pdf(sprintf('%s_clusters.pdf',baseName),width=28,height=35)
  plot(gg_cl)
  dev.off()
  pdf(sprintf('%s_genes.pdf',baseName),width=28,height=35)
  plot(gg_genes)
  dev.off()
  return(list(cells=gg_cells,clusters=gg_cl,genes=gg_genes))
}
