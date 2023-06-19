# Raphael Mourad
# 28/08/2019
# Adapted by Vincent ROCHER for Snakemake and trans interaction
# Create a big square matrix with all chtromosomes
.libPaths(c(.libPaths(),"/usr/local/lib/R/site-library"))
library(Matrix)
library(MASS)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)
library(plyranges)
library(keras)
#Functions
#Return a boolean
# Test if the first chromosome (x1) > to the second
# used to order column in our data frame
getOrder <- function(x1,x2){
  x1 <-x1 %>% str_remove("chr") %>% dplyr::recode(X="23",Y="24") %>% as.numeric()
  x2 <-x2 %>% str_remove("chr") %>% dplyr::recode(X="23",Y="24") %>% as.numeric()
  x1>x2
}
#Same but for our specific order with DIvA
#CHECK THE CHROMOSOME ORDER ON JUICER 
getorder_for_our_Run <- function(x1,x2){
  my_order <- paste0("chr",c(1:7,"X",8:18,20,"Y",19,22,21))
  grep(paste("^",x1,"$", sep=""),my_order) > grep(paste("^",x2,"$", sep=""),my_order)
}
#Return a 0 score if the file is empty (to fill the matrix with 0)
processingNULL <- function(x1,x2,SeqInfo){
  message(str_c("WARNING :",x1,"and",x2,"combination can't be loaded, generation of 0 score instead."),sep=" ")
  if(getorder_for_our_Run(x1,x2)){
    c(x1,x2) %<-% c(x2,x1)
  }
  indices_ranges.1 <- my.ranges %>% filter(seqnames == x1) %>% as_tibble() %>% pull(idx)
  indices_ranges.2 <- my.ranges %>% filter(seqnames == x2) %>% as_tibble() %>% pull(idx)
  
  tibble(
    idx1 = c(1,max(indices_ranges.1)),
    idx2 = c(1,max(indices_ranges.2)),
    V3 = c(0,0)
  )
}
#Return a data frame of index1 index2 score interaction (main code)
processingMATRIX <- function(x1,x2,datai,my.ranges,whichorder=getorder_for_our_Run){
  #For each chromosome, get chrom size and compute the last bin
  chrend1=seqlengths(SeqInfo[x1]) 
  binend1=ceiling(chrend1/binSize)
  chrend2=seqlengths(SeqInfo[x2])
  binend2=ceiling(chrend2/binSize)
  #Remove NA value (no interaction)
  datai=datai[!is.na(datai[,3]),]
  if(whichorder(x1,x2)){
    colnames.datai <- colnames(datai)
    datai <- datai[,c(2,1,3)]
    colnames(datai) <- colnames.datai
  }
  #Add the last bin (used to compute the SparseMatrix)
  datai=rbind(datai,c((binend1-1)*binSize,(binend2-1)*binSize,0))
  if(obserOE=="observed"){
    datai[,3]=round(datai[,3])
  }else{
    datai[,3]=round(datai[,3],2)
  }
  #Add 1 : conversion from base 0 to base 1
  datai[,1] <- datai[,1] + 1
  datai[,2] <- datai[,2] + 1
  # Get the index of the first chromosome
  indices_ranges <- my.ranges %>% filter(seqnames == x1) %>% as_tibble() %>% dplyr::select(start,idx)
  #Left join index and position
  datai <- datai %>% as_tibble() %>% left_join(indices_ranges,by = c("V1"="start")) %>%
    dplyr::rename(idx1 = idx) 
  # Get the index of the second chromosome
  indices_ranges <- my.ranges %>% filter(seqnames == x2) %>% as_tibble() %>% dplyr::select(start,idx)
  #Left join index and position
  datai <- datai %>% left_join(indices_ranges,by = c("V2"="start")) %>%
    dplyr::rename(idx2 = idx) 
  #Return the two index and the score
  datai %>% dplyr::select(idx1,idx2,V3)
}
# Resolution
binSize=as.numeric(snakemake@params[["my_res"]])
binSizeTxt=paste0(binSize/1e3,"kb")
# LOAD DATA AND RESULTS ----------------------------------------------------------
dumpDir <- snakemake@params[["dumpdir"]]
# Exp, normalization type 
obserOE= snakemake@params[["obserOE"]] # "observed" "OE"
Experiment = snakemake@params[["exp"]]
my_norm = snakemake@params[["norm"]]


# Chromosomes
# Maybe try to get the chromosome list from snakemake directly
Chr.V <- paste0("chr",c(1:22,"X"))
# Compute all combination
Chr.V.combn <- crossing(col1=Chr.V,col2=Chr.V)
# Get the chrom size
SeqInfo=seqinfo(BSgenome.Hsapiens.UCSC.hg19)
# Build the bins of all genome
my.ranges <- tileGenome(SeqInfo[Chr.V],tilewidth=binSize, cut.last.tile.in.chrom=TRUE)
# my.ranges <- my.ranges %>% as_tibble() %>% group_by(seqnames) %>% filter(end != max(end)) %>% as_granges()
# Create an index
my.ranges <- my.ranges %>% mutate(idx = 1:length(my.ranges))


#For each combination of chromosome, do
res <- apply(Chr.V.combn,1,function(x){
  x1 <- x[1]
  x2 <- x[2]
  message(paste(x1,x2,sep="vs"))
  
  
  
  #Load the file
  file <- paste(c(
    dumpDir,
    str_c(c("dump",obserOE,my_norm,Experiment,str_c(c(x1,x2),collapse="_"),binSizeTxt),
          collapse="_")
  ),collapse="/") %>% paste("txt.gz",sep=".")
  #using fread from data.table
  datai=as.matrix(fread(cmd=paste0("zcat ",file),sep='\t',header=F))
  if(nrow(datai) == 0){
    processingNULL(x1,x2,my.ranges)
  }else{
    if(snakemake@params[["whichorder"]] == "custom"){
      processingMATRIX(x1,x2,datai,my.ranges)
    }else{
      processingMATRIX(x1,x2,datai,my.ranges,whichorder=getOrder)
    }
    
  }
  
}) %>% bind_rows() %>% as.matrix()

#Create the sparse matrix from the data frame pos1 pos2 score
data.Mati=sparseMatrix(i=res[,1],j=res[,2],x=res[,3])
#Create a symmetrical matrix with only the first half 
# Fill the other half with the transpose matrix
data.MatSymi=data.Mati+t(data.Mati)
# And divide by 2 the diag
diag(data.MatSymi)=diag(data.MatSymi)/2
#Return the range and the matrix sparse
res <- list(
  "range"=my.ranges,
  "matrix"=data.MatSymi
)
saveRDS(object = res,file = snakemake@output[[1]])
