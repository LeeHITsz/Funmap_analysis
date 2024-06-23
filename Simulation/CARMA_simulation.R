library(CARMA)
library(data.table)
library(magrittr)
library(dplyr)
library(devtools)
library(R.utils)

# Setting Prefix
index = "25"

L0 <- 10
lambda.list<-list()
lambda.list[[1]]<-L0

gene_name = fread(file = "Simulation loci.csv", sep = ",", header = T)$gene.names

for(gene in 1:10){
    gene = 10
    ld = fread(file = paste(index,"/",gene_name[gene],"/data/ld.txt", sep=""), sep = " ", header = F)
    ld.list<-list()
    ld.list[[1]]<-as.matrix(ld)
    for(i in 1:50){
        i <- 50
        sumstat<-fread(file = paste(index,"/",gene_name[gene],"/data/C",i,".txt", sep=""), sep = "\t", header = F)
        anno=fread(file = paste(index,"/",gene_name[gene],"/data/C",i,"_anno.txt", sep=""), sep="\t", header = T, check.names = F, data.table = F, stringsAsFactors = F)
        z.list<-list()
        anno.list<-list()
        z.list[[1]]<-sumstat[,2]
        anno.list[[1]]<-anno[,-1]

        CARMA.results<-CARMA(z.list, ld.list, lambda.list=lambda.list, outlier.switch=F,rho.index=0.95)
        result = data.frame(
            "SNPs" = sumstat[,1],
            "PIP" = CARMA.results[[1]]$PIPs
        )
        fwrite(x = result, file = paste(index,"/",gene_name[gene],"/CARMA_result_no/C",i,".txt.csv", sep=""), sep = ",", quote = F, na = "NA", row.names = F, col.names = F)
        
        sets = data.frame(
            "SNPs" = sumstat[,1],
            "CS" = 0
        )
        if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
            for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){
                sets$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
            }
        }
        fwrite(x = sets, file = paste(index,"/",gene_name[gene],"/CARMA_sets/C",i,".sets", sep=""), sep = ",", quote = F, na = "NA", row.names = F, col.names = F)
        print(paste(index,gene_name[gene],i,"no","completed"))
        
        CARMA.results<-CARMA(z.list,ld.list,lambda.list=lambda.list,w.list=anno.list,outlier.switch=F,rho.index=0.95)
        result = data.frame(
            "SNPs" = sumstat[,1],
            "PIP" = CARMA.results[[1]]$PIPs
        )
        fwrite(x = result, file = paste(index,"/",gene_name[gene],"/CARMA_result_anno/C",i,".txt.csv", sep=""), sep = ",", quote = F, na = "NA", row.names = F, col.names = F)
        
        sets = data.frame(
            "SNPs" = sumstat[,1],
            "CS" = 0
        )
        if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
            for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){
                sets$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
            }
        }
        fwrite(x = sets, file = paste(index,"/",gene_name[gene],"/CARMA_anno_sets/C",i,".sets", sep=""), sep = ",", quote = F, na = "NA", row.names = F, col.names = F)
        print(paste(index,gene_name[gene],i,"anno","completed"))
    }
}
