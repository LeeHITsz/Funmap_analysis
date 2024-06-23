library(CARMA)
library(data.table)
library(magrittr)
library(dplyr)
library(devtools)
library(R.utils)

index <- 'Timing'
gene <- 'TBX3'  # 'GENE8', 'GENE5', 'TBX3', 'C6orf211', 'MYC'
label <- '9'

L0 <- 10
lambda.list<-list()
lambda.list[[1]]<-L0

ld = fread(file = paste(index,"/",gene,"/data/ld.txt", sep=""), sep = " ", header = F)
ld.list<-list()
ld.list[[1]]<-as.matrix(ld)

print(paste('CARMA_anno', label, ' :', sep=""))
cat('\n')

for(i in 251:300){
    sumstat<-fread(file = paste(index,"/",gene,"/data/C",i,".txt", sep=""), sep = "\t", header = F)
    anno=fread(file = paste(index,"/",gene,"/data/C",i,"_anno.txt", sep=""), sep="\t", header = T, check.names = F, data.table = F, stringsAsFactors = F)
    z.list<-list()
    anno.list<-list()
    z.list[[1]]<-sumstat[,2]
    anno.list[[1]]<-anno[,-1]

    start_time <- Sys.time()
    CARMA.results<-CARMA(z.list,ld.list,lambda.list=lambda.list,w.list=anno.list,outlier.switch=F,rho.index=0.95)
    end_time <- Sys.time()
    cat('\n')
    print(end_time - start_time)
    cat('\n')
}
















