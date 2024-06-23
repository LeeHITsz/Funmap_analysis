library(XPASS)
library(susieR)
library(data.table)

# Setting Prefix
index <- "26"

gene_name = fread(file = "Simulation loci.csv", sep = ",", header = T)$gene.names

n <- 50000
m <- 20
m_anno <- 10
wvar <- 0.01

K_true <- 2
sb <- 0.0075
se <- 1 - sb

maxLD <- 1

for(gene in 1:10){

  gene <- gene_name[gene]
  
  geno <- XPASS::read_data(paste0("dataX/",gene,"/data/",gene))
  snps_indp <- read.table(paste0("indpSNP/",gene,".prune.in"),header=F)$V1
  
  X <- geno$X
  X <- X[sample(nrow(X),n,replace=F),]
  X <- scale(X)
  p <- ncol(X)
  LD <- cor(X)
  
  out.LD = data.frame(LD)
  fwrite(x = out.LD, file = paste0(index,"/",gene,"/data/ld.txt"), sep = " ", quote = F, na = "NA", row.names = F, col.names = F)
  
  idx_prune <- match(snps_indp,geno$snps)
  
  LD_max <- apply(LD,2,function(x){sort(abs(x),decreasing = T)[2]})
  idx_candidate <- which(LD_max<=maxLD)
  idx_candidate <- idx_candidate[idx_candidate %in% idx_prune]
  
  nrep <- 50

  for (irep in 1:nrep) {
    
    set.seed(irep)
    
    A <- matrix(rnorm(m*p,0,1),p,m)
    out.anno = cbind(data.frame(geno$snps), data.frame(A))
    colnames(out.anno)[1] <- "SNP"
    fwrite(x = out.anno, file = paste(index,"/",gene,"/data/C",irep,"_anno.txt", sep=""), sep = "\t", quote = F, na = "NA", row.names = F, col.names = T)
    
    w <- rnorm(m,0,sqrt(wvar))
    w[sample(1:m, m_anno)] <- 0
    pi_unnormalized <- A[idx_candidate,] %*% w
    pi_unnormalized <- pi_unnormalized-min(pi_unnormalized)
    pi_normalized <- exp(pi_unnormalized)/sum(exp(pi_unnormalized))
    
    idx_causal <- sample(idx_candidate,K_true,replace = F,prob = pi_normalized)
    idx_causal_out <- rep(0, p)
    idx_causal_out[idx_causal] = 1
    out.trueSNP = data.frame(
      "SNP" = geno$snps,
      "label" = idx_causal_out
    )
    fwrite(x = out.trueSNP, file = paste(index,"/",gene,"/true_SNP/C",irep,"_true.txt", sep=""), sep = "\t", quote = F, na = "NA", row.names = F, col.names = F)
    
    b <- rep(0, p)
    b[idx_causal] <- rnorm(K_true, 0, sqrt(sb / K_true))
    y0 <- X %*% b
    y <- y0 + rnorm(n, 0, sqrt((1 - sb) * var(y0) / sb))
    
    sumstats <- univariate_regression(X, y)
    z_scores <- sumstats$betahat / sumstats$sebetahat
    
    out.z_scores = data.frame(
      "SNP" = geno$snps,
      "z" = z_scores
    )
    fwrite(x = out.z_scores, file = paste(index,"/",gene,"/data/C",irep,".txt", sep=""), sep = "\t", quote = F, na = "NA", row.names = F, col.names = F)
    print(paste0(index,"-",gene,": C",irep," completed"))
  }
}

