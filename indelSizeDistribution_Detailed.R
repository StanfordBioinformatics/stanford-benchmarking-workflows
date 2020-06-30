#!/usr/bin/env Rscript

##USAGE: Rscript indelSizeDistribution_Detailed.R TruthVCF_bcftoolsStatsOut.txt TPonly.vcf.gz FPonly.vcf.gz FNonly.vcf.gz IndelsizeDistribution.txt IndelsizeDistributionPlot.pdf


#args = commandArgs(trailingOnly=TRUE)
args = commandArgs(TRUE)
if (length(args)==0 | length(args)!=6 )
  stop("First argument is text input file (bcftools stats on Truth.vcf.gz), Second: *_TPonly.vcf.gz , Third: *_FPonly.vcf.gz, Fourth: *_FNonly.vcf.gz, Fifth: Output file(txt) Sixth: Output indel plot file (PDF)", call.=FALSE)
  
inFile <- read.table(file=args[1],stringsAsFactors=FALSE, sep="\t",skip=2,header=FALSE) ## first two rows are header info so skip
## count indel of size 1 or -1
range1Sum <- sum(inFile$V4[which(inFile$V3==-1 | inFile$V3==1)])

## count insertion or deletion of size 2 to 5
range2a=seq(-5,-2,by=1)
range2b=seq(2,5,by=1)
range2aSum <- sum(unlist(lapply(range2a, function(x) { inFile$V4[which(inFile$V3==x)] }) ))
range2bSum <- sum(unlist(lapply(range2b, function(x) { inFile$V4[which(inFile$V3==x)] }) ))
range2Sum=range2aSum+range2bSum

## count insertion or deletion of size 6 to 10
range3a=seq(-10,-6,by=1)
range3b=seq(6,10,by=1)
range3aSum <- sum(unlist(lapply(range3a, function(x) { inFile$V4[which(inFile$V3==x)] }) ))
range3bSum <- sum(unlist(lapply(range3b, function(x) { inFile$V4[which(inFile$V3==x)] }) ))
range3Sum=range3aSum+range3bSum

## count insertion or deletion of size 11 to 20
range4a=seq(-20,-11,by=1)
range4b=seq(11,20,by=1)
range4aSum <- sum(unlist(lapply(range4a, function(x) { inFile$V4[which(inFile$V3==x)] }) ))
range4bSum <- sum(unlist(lapply(range4b, function(x) { inFile$V4[which(inFile$V3==x)] }) ))
range4Sum=range4aSum+range4bSum

## count insertion or deletion of size 21 to 50
range5a=seq(-50,-21,by=1)
range5b=seq(21,50,by=1)
range5aSum <- sum(unlist(lapply(range5a, function(x) { inFile$V4[which(inFile$V3==x)] }) ))
range5bSum <- sum(unlist(lapply(range5b, function(x) { inFile$V4[which(inFile$V3==x)] }) ))
range5Sum=range5aSum+range5bSum

## count insertion or deletion of size 51 to 100
range6a=seq(-100,-51,by=1)
range6b=seq(51,100,by=1)
range6aSum <- sum(unlist(lapply(range6a, function(x) { inFile$V4[which(inFile$V3==x)] }) ))
range6bSum <- sum(unlist(lapply(range6b, function(x) { inFile$V4[which(inFile$V3==x)] }) ))
range6Sum=range6aSum+range6bSum

## Obtain the line number where the vcf header ends
lineNo<- as.numeric(system(paste("zgrep -n '#CHROM'",args[2],"|cut -d ':' -f1",sep=" "), intern=TRUE))

## Skip the headers and read the REF and ALT columns in TP, FP and FN INDEL annotation files [ Avoids using bcftools view -H TP.vcf.gz |cut -f4,5 as extra workflow step]
##Use text= instead of file= in read.table to directly read system call output
TPcounts <- read.table(text=system(paste("zcat ", args[2], "|tail -n +", lineNo, "|cut -f4,5",sep=""),intern=TRUE) ,sep="\t",header=TRUE,stringsAsFactors=FALSE)
FPcounts <- read.table(text=system(paste("zcat ", args[3], "|tail -n +", lineNo, "|cut -f4,5",sep=""),intern=TRUE) ,sep="\t",header=TRUE,stringsAsFactors=FALSE)
FNcounts <- read.table(text=system(paste("zcat ", args[4], "|tail -n +", lineNo, "|cut -f4,5",sep=""),intern=TRUE) ,sep="\t",header=TRUE,stringsAsFactors=FALSE)

##Report Indel size as number of bases in REF if the ALT column is <DEL>
TPcounts$INDELsize <- as.numeric( apply(TPcounts,1, function(x) { ifelse(x[2]=="<DEL>",nchar(x[1]),0) }) )
FPcounts$INDELsize <- as.numeric( apply(FPcounts,1, function(x) { ifelse(x[2]=="<DEL>",nchar(x[1]),0) }) )
FNcounts$INDELsize <- as.numeric( apply(FNcounts,1, function(x) { ifelse(x[2]=="<DEL>",nchar(x[1]),0) }) )

##Report Indel size as (No of bases in REF - No of bases in ALT) (one anchor base is common in ALT). Check if ALT has multiple entries separated by commas

##Obtain the rows where ALT has multiple entries  which(TP_multipleALTs!=0)  
TP_multipleALTs <- unlist(lapply(TPcounts$ALT, function(x) length(grep(x,pattern=",",value=TRUE)) ))
FP_multipleALTs <- unlist(lapply(FPcounts$ALT, function(x) length(grep(x,pattern=",",value=TRUE)) ))
FN_multipleALTs <- unlist(lapply(FNcounts$ALT, function(x) length(grep(x,pattern=",",value=TRUE)) ))

## For ALT with multiple values, INDEL size = nchar(REF) - nchar(first entry in ALT)  
TPcounts$INDELsize[which(TP_multipleALTs!=0)] <- abs( nchar(TPcounts[,1][which(TP_multipleALTs!=0)]) - nchar(sapply(strsplit(as.character(TPcounts$ALT[which(TP_multipleALTs!=0)]), ","), "[[", 1))  )                                        
FPcounts$INDELsize[which(FP_multipleALTs!=0)] <- abs( nchar(FPcounts[,1][which(FP_multipleALTs!=0)]) - nchar(sapply(strsplit(as.character(FPcounts$ALT[which(FP_multipleALTs!=0)]), ","), "[[", 1))  )                                        
FNcounts$INDELsize[which(FN_multipleALTs!=0)] <- abs( nchar(FNcounts[,1][which(FN_multipleALTs!=0)]) - nchar(sapply(strsplit(as.character(FNcounts$ALT[which(FN_multipleALTs!=0)]), ","), "[[", 1))  )                                        

## Rest of ALT which is not <DEL> and not multiple values, change indel size from zero to relevant ELSE retain INDEL size as is
## as.numeric conversion of INDEL size column is very important even though indelsize is numeric before apply is used on it
TPcounts$INDELsize <- as.numeric( apply(TPcounts,1, function(x) { ifelse( as.numeric(x[3])==0, abs( nchar(x[1]) - nchar(x[2]) ), x[3]) }) )
FPcounts$INDELsize <- as.numeric( apply(FPcounts,1, function(x) { ifelse( as.numeric(x[3])==0, abs( nchar(x[1]) - nchar(x[2]) ), x[3]) }) ) 
FNcounts$INDELsize <- as.numeric( apply(FNcounts,1, function(x) { ifelse( as.numeric(x[3])==0, abs( nchar(x[1]) - nchar(x[2]) ), x[3]) }) )

## Compute the number of TPs, FPs and FNs in each of the INDEL size ranges
## Compute the correspoding precision and recall statistics for each indel size bin 
## Recall = TP/(TP+FN)
## Precision = TP/(TP+FP)

IndelSize1_TPcounts <- length(which(TPcounts$INDELsize==1))
IndelSize1_FPcounts <- length(which(FPcounts$INDELsize==1))
IndelSize1_FNcounts <- length(which(FNcounts$INDELsize==1))
IndelSize1_Precision <- round(100*(IndelSize1_TPcounts/(IndelSize1_TPcounts + IndelSize1_FPcounts)),digits=2)
IndelSize1_Recall <- round(100*(IndelSize1_TPcounts/(IndelSize1_TPcounts + IndelSize1_FNcounts)),digits=2)

IndelSize2to5_TPcounts <- length(which(TPcounts$INDELsize>=2 & TPcounts$INDELsize<=5))
IndelSize2to5_FPcounts <- length(which(FPcounts$INDELsize>=2 & FPcounts$INDELsize<=5))
IndelSize2to5_FNcounts <- length(which(FNcounts$INDELsize>=2 & FNcounts$INDELsize<=5))
IndelSize2to5_Precision <- round(100*(IndelSize2to5_TPcounts/(IndelSize2to5_TPcounts + IndelSize2to5_FPcounts)),digits=2)
IndelSize2to5_Recall <- round(100*(IndelSize2to5_TPcounts/(IndelSize2to5_TPcounts + IndelSize2to5_FNcounts)),digits=2)

IndelSize6to10_TPcounts <- length(which(TPcounts$INDELsize>=6 & TPcounts$INDELsize<=10))
IndelSize6to10_FPcounts <- length(which(FPcounts$INDELsize>=6 & FPcounts$INDELsize<=10))
IndelSize6to10_FNcounts <- length(which(FNcounts$INDELsize>=6 & FNcounts$INDELsize<=10))
IndelSize6to10_Precision <- round(100*(IndelSize6to10_TPcounts/(IndelSize6to10_TPcounts + IndelSize6to10_FPcounts)),digits=2)
IndelSize6to10_Recall <- round(100*(IndelSize6to10_TPcounts/(IndelSize6to10_TPcounts + IndelSize6to10_FNcounts)),digits=2)

IndelSize11to20_TPcounts <- length(which(TPcounts$INDELsize>=11 & TPcounts$INDELsize<=20))
IndelSize11to20_FPcounts <- length(which(FPcounts$INDELsize>=11 & FPcounts$INDELsize<=20))
IndelSize11to20_FNcounts <- length(which(FNcounts$INDELsize>=11 & FNcounts$INDELsize<=20))
IndelSize11to20_Precision <- round(100*(IndelSize11to20_TPcounts/(IndelSize11to20_TPcounts + IndelSize11to20_FPcounts)),digits=2)
IndelSize11to20_Recall <- round(100*(IndelSize11to20_TPcounts/(IndelSize11to20_TPcounts + IndelSize11to20_FNcounts)),digits=2)

IndelSize21to50_TPcounts <- length(which(TPcounts$INDELsize>=21 & TPcounts$INDELsize<=50))
IndelSize21to50_FPcounts <- length(which(FPcounts$INDELsize>=21 & FPcounts$INDELsize<=50))
IndelSize21to50_FNcounts <- length(which(FNcounts$INDELsize>=21 & FNcounts$INDELsize<=50))
IndelSize21to50_Precision <- round(100*(IndelSize21to50_TPcounts/(IndelSize21to50_TPcounts + IndelSize21to50_FPcounts)),digits=2)
IndelSize21to50_Recall <- round(100*(IndelSize21to50_TPcounts/(IndelSize21to50_TPcounts + IndelSize21to50_FNcounts)),digits=2)

IndelSize51greater_TPcounts <- length(which(TPcounts$INDELsize>=51))
IndelSize51greater_FPcounts <- length(which(FPcounts$INDELsize>=51))
IndelSize51greater_FNcounts <- length(which(FNcounts$INDELsize>=51))
IndelSize51greater_Precision <- round(100*(IndelSize51greater_TPcounts/(IndelSize51greater_TPcounts + IndelSize51greater_FPcounts)),digits=2)
IndelSize51greater_Recall <- round(100*(IndelSize51greater_TPcounts/(IndelSize51greater_TPcounts + IndelSize51greater_FNcounts)),digits=2)

#Truth total/frequency is based on TruthVCF
write(x="Indel Size distribution", file=args[5],append=TRUE)
write(x="INDEL SIZE\t FREQUENCY\t #TP\t #FP\t #FN\t PRECISION(%)\t RECALL(%)", file=args[5],append=TRUE)
write(x=paste("1",range1Sum,IndelSize1_TPcounts,IndelSize1_FPcounts,IndelSize1_FNcounts,IndelSize1_Precision,IndelSize1_Recall,sep="\t"), file=args[5],append=TRUE)
write(x=paste("2 - 5", range2Sum,IndelSize2to5_TPcounts,IndelSize2to5_FPcounts,IndelSize2to5_FNcounts,IndelSize2to5_Precision,IndelSize2to5_Recall,sep="\t"), file=args[5],append=TRUE)
write(x=paste("6 - 10", range3Sum,IndelSize6to10_TPcounts,IndelSize6to10_FPcounts,IndelSize6to10_FNcounts,IndelSize6to10_Precision,IndelSize6to10_Recall,sep="\t"), file=args[5],append=TRUE)
write(x=paste("11 - 20", range4Sum,IndelSize11to20_TPcounts, IndelSize11to20_FPcounts,IndelSize11to20_FNcounts,IndelSize11to20_Precision,IndelSize11to20_Recall,sep="\t"), file=args[5],append=TRUE)
write(x=paste("21 - 50", range5Sum,IndelSize21to50_TPcounts,IndelSize21to50_FPcounts,IndelSize21to50_FNcounts,IndelSize21to50_Precision,IndelSize21to50_Recall,sep="\t"), file=args[5],append=TRUE)
write(x=paste("51 and greater",range6Sum,IndelSize51greater_TPcounts,IndelSize51greater_FPcounts,IndelSize51greater_FNcounts,IndelSize51greater_Precision,IndelSize51greater_Recall,sep="\t"), file=args[5],append=TRUE)    
    
indelPlot <- cbind(c(IndelSize1_TPcounts, IndelSize2to5_TPcounts,IndelSize6to10_TPcounts,IndelSize11to20_TPcounts,IndelSize21to50_TPcounts,IndelSize51greater_TPcounts),
                   c(IndelSize1_FPcounts, IndelSize2to5_FPcounts,IndelSize6to10_FPcounts,IndelSize11to20_FPcounts,IndelSize21to50_FPcounts,IndelSize51greater_FPcounts),
				   c(IndelSize1_FNcounts, IndelSize2to5_FNcounts,IndelSize6to10_FNcounts,IndelSize11to20_FNcounts,IndelSize21to50_FNcounts,IndelSize51greater_FNcounts))

##Convert each element to numeric else matrix will change a single vector
## apply transforms the matrix as well
indelPlot <- apply(indelPlot, 1, function(x) x <- as.numeric(x))
colnames(indelPlot) <- c("1", "2-5","6-10","11-20","21-50","51 or greater")
rownames(indelPlot) <- c("TP","FP","FN")

##Generate plot and store as pdf
pdf(args[6],width=10,height=10)
barplot(indelPlot,col=c("darkblue","orange","black"),main="Indel size distribution", xlab="Indel size", ylab= "Frequency",names.arg=colnames(indelPlot),legend=c("True Positives","False Positives", "False Negatives"),beside=TRUE)
dev.off()



