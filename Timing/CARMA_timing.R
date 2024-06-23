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

print(paste('CARMA', label, ' :', sep=""))
cat('\n')

for(i in 251:300){
    sumstat<-fread(file = paste(index,"/",gene,"/data/C",i,".txt", sep=""), sep = "\t", header = F)
    z.list<-list()
    z.list[[1]]<-sumstat[,2]

    start_time <- Sys.time()
    CARMA.results<-CARMA(z.list, ld.list, lambda.list=lambda.list, outlier.switch=F,rho.index=0.95)
    end_time <- Sys.time()
    cat('\n')
    print(end_time - start_time)
    cat('\n')
}
