library(susieR)
library(data.table)
library(magrittr)
library(dplyr)
library(devtools)
library(R.utils)

options (warn = -1)

index <- 'Timing'
gene <- 'MYC'  # 'GENE8', 'GENE5', 'TBX3', 'C6orf211', 'MYC'
label <- '5'

ld = fread(file = paste(index,"/",gene,"/data/ld.txt", sep=""), sep = " ", header = F)

print(paste('SuSiE', label, ' :', sep=""))
cat('\n')

for(i in 1:50){
    sumstat<-fread(file = paste(index,"/",gene,"/data/C",i,".txt", sep=""), sep = "\t", header = F)
    start_time <- Sys.time()
    fitted_rss1 <- susie_rss(z=c(sumstat[,2])$V2, n=50000, R=as.matrix(ld), L=10, estimate_residual_variance=TRUE)
    end_time <- Sys.time()
    print(end_time - start_time)
}
