import os
import numpy as np
import pandas as pd
from Funmap import SUSIE
from Funmap import FUNMAP

np.set_printoptions(threshold=1e8)

# Setting Prefix
index = '26'

gene_name = pd.read_csv("Simulation loci.csv")["gene.names"].values


for gene in gene_name[0:10]:

    R = pd.read_csv(index + "/" + gene + '/data/ld.txt', sep='\s+', header=None)
    R = R.values

    n = 50000

    for k in range(0, 50):

        A = pd.read_csv(index + "/" + gene + '/data/C' + str(k + 1) + "_anno.txt", sep='\t')
        A = A.values[:, 1:].astype(float)
        p = A.shape[0]
        m = A.shape[1]

        z = pd.read_csv(index + "/" + gene + '/data/C' + str(k + 1) + ".txt", header=None, sep='\t')
        snp_name = z.values[:, 0]
        z = z.values[:, 1]

        result = SUSIE(z, R, n, L=10)
        s_result = pd.DataFrame({'SNP': snp_name, 'PIP': result.pip})
        s_result.to_csv(index + "/" + gene + "/susie_result_no/C" + str(k + 1) + ".txt.csv",
                        header=None, index=False)
        file = open(index + "/" + gene + "/SuSiE_sets/C" + str(k + 1) + ".sets", "w")
        file.write(str(result.sets))
        file.close()

        result = FUNMAP(z, R, A, n, L=10)
        s_result = pd.DataFrame({'SNP': snp_name, 'PIP': result.pip})
        s_result.to_csv(index + "/" + gene + "/Funmap_result/C" + str(k + 1) + ".txt.csv",
                        header=None, index=False)
        file = open(index + "/" + gene + "/Funmap_sets/C" + str(k + 1) + ".sets", "w")
        file.write(str(result.sets))
        file.close()
