library(CARMA)
library(data.table)
library(magrittr)
library(dplyr)
library(devtools)
library(R.utils)

Chr = "21"

ld_list = fread(file = paste('chr',Chr,'/LD/LD_list_TG.txt', sep=""), header = F)$V1

L0 <- 10
lambda.list<-list()
lambda.list[[1]]<-L0

for(Pos in ld_list){

    ld = fread(file = paste('chr',Chr,'/data/chr',Chr,'_',Pos*1000000+1,'_',(Pos+3)*1000000+1,'.ld', sep=""), sep = " ", header = F)

    anno = fread(file = paste('chr',Chr,'/data/chr',Chr,'_',Pos*1000000+1,'_',(Pos+3)*1000000+1,'.annotations', sep=""), sep=" ", header = T, check.names = F, data.table = F, stringsAsFactors = F)
    sd_vec <- apply(anno, 2, sd)
    if(any(sd_vec == 0)){
        invariant_var_idx <- which(sd_vec == 0) 
        other_var_idx <- setdiff(seq_len(ncol(anno)), invariant_var_idx)
        other_vars <- anno[, other_var_idx]
    } else {
        other_vars <- anno
    }

    ld.list<-list()
    anno.list<-list()
    ld.list[[1]]<-as.matrix(ld)
    anno.list[[1]]<-other_vars
    for(trait in c('HDL','LDL','Cho')){
    # trait = 'TG'
        if (!file.exists(paste('chr',Chr,'/data/chr',Chr,'_',Pos*1000000+1,'_',(Pos+3)*1000000+1,'_',trait,'.txt', sep=""))){
            next
        }

        sumstat<-fread(file = paste('chr',Chr,'/data/chr',Chr,'_',Pos*1000000+1,'_',(Pos+3)*1000000+1,'_',trait,'.txt', sep=""), sep = "\t", header = F)
        z.list<-list()
        z.list[[1]]<-sumstat[,2]

        CARMA.results<-CARMA(z.list,ld.list,lambda.list=lambda.list,w.list=anno.list,outlier.switch=F,rho.index=0.95,Max.Model.Dim=2000,all.iter=4,all.inner.iter=4)
        
        result = data.frame(
            "SNPs" = sumstat[,1],
            "PIP" = CARMA.results[[1]]$PIPs
        )
        fwrite(x = result, file = paste('chr',Chr,'/CARMA_result_anno/','chr',Chr,'_',Pos*1000000+1,'_',(Pos+3)*1000000+1,'_',trait,'.csv', sep=""), sep = ",", quote = F, na = "NA", row.names = F, col.names = F)
        
        sets = data.frame(
            "SNPs" = sumstat[,1],
            "CS" = 0
        )
        if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
            for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){
                sets$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
            }
        }
        fwrite(x = sets, file = paste('chr',Chr,'/CARMA_anno_sets/','chr',Chr,'_',Pos*1000000+1,'_',(Pos+3)*1000000+1,'_',trait,'.sets', sep=""), sep = ",", quote = F, na = "NA", row.names = F, col.names = F)
    }
}
