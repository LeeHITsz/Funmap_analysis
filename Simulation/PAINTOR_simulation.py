import os
import numpy as np
import pandas as pd


loci = pd.read_csv("Simulation loci.csv")
gene_name = loci["gene.names"].values

# Setting Prefix
index = '26'    

m = 20
anno = ["X" + str(j) for j in range(1, m+1)]
# anno = ["M" + str(j) for j in range(1, m+1)]

def compute_paintor_credible_sets(posteriors, ld_matrix, threshold=0.95):

    n_snps = posteriors.shape[0]
    credible_sets = np.zeros(n_snps, dtype=int)
    credible_set_id = 0

    while np.sum(posteriors) > threshold:

        # Find SNP with max posterior
        max_idx = np.argmax(posteriors)

        # Get indices of SNPs in LD above threshold
        ld_indexes = np.where(ld_matrix[max_idx] > 0.5)[0]

        if ld_indexes.size != 0 and np.sum(posteriors[ld_indexes]) > threshold:
            # Subset posteriors and LD
            sub_posteriors = posteriors[ld_indexes]
            sub_posteriors /= sub_posteriors.sum()

            # Find SNPs until cumulative sum reaches threshold
            credible_indices = np.argsort(sub_posteriors)[
                               -(np.where(np.cumsum(np.sort(sub_posteriors)[::-1]) < threshold)[0].shape[0] + 1):]

            credible_set_id += 1
            credible_sets[ld_indexes[credible_indices]] = credible_set_id

        # Remove subsetted SNPs
        posteriors[ld_indexes] = 0

    return credible_sets

for gene in gene_name[0:10]:

    LD = pd.read_csv(index + "/" + gene + '/data/ld.txt', sep='\s+', header=None)
    np.savetxt(index+"/"+gene+"/data/LD.txt", LD, fmt='%.3f')

    cmd1 = "PAINTOR_V3.0/PAINTOR " \
           "-input " + index+"/"+gene+"/data/input.files " \
           "-in " + index+"/"+gene+"/data/ " \
           "-out " + index+"/"+gene+"/PAINTOR_result_no/ " \
           "-Zhead Zscore " \
           "-LDname ld " \
           "-mcmc " \
           "-annotations Coding"

    cmd2 = "PAINTOR_V3.0/PAINTOR " \
           "-input " + index+"/"+gene+"/data/input.files " \
           "-in " + index+"/"+gene+"/data/ " \
           "-out " + index+"/"+gene+"/PAINTOR_result_anno/ " \
           "-Zhead Zscore " \
           "-LDname ld " \
           "-mcmc " \
           "-annotations " + ",".join(anno)

    for i in range(50):

        np.savetxt(index + "/" + gene + "/data/input.files", np.array(['C'+str(i+1)]),
                   fmt='% s', delimiter=' ', newline='\n', header='')

        data_z = pd.read_csv(index + "/" + gene + "/data" + "/" + "C" + str(i + 1) + ".txt", header=None, sep="\t")
        np.savetxt(index + "/" + gene + "/data" + "/" + "C" + str(i + 1), data_z.values[:, 1],
                   fmt='%.6f', delimiter=' ', newline='\n', header='Zscore', comments='')

        os.system("cp " + index + "/" + gene + "/data" + "/LD.txt " + index + "/" + gene + "/data" + "/" + "C" + str(
            i + 1) + ".ld")

        data_anno = pd.read_csv(index + "/" + gene + '/data/C'+str(i+1)+'_anno.txt', sep='\t')
        data_anno.rename(columns={'SNP': 'Coding'}, inplace=True)
        data_anno['Coding'] = 1
        data_anno.to_csv(index + "/" + gene + "/data" + "/" + "C" + str(i + 1) + ".annotations", index=False, sep=" ")

        os.system(cmd1)
        print(gene+": " + "C" + str(i+1) + " no completed")
        os.system(cmd2)
        print(gene+": " + "C" + str(i+1) + " anno completed")

        os.system("rm "+index+"/"+gene+"/data"+"/input.files")

    for i in range(50):

        os.system("rm "+index+"/"+gene+"/data"+"/"+"C"+str(i+1))
        os.system("rm "+index+"/"+gene+"/data"+"/"+"C"+str(i+1)+".ld")
        os.system("rm "+index+"/"+gene+"/data"+"/"+"C"+str(i+1)+".annotations")
    
    os.system("rm "+index+"/"+gene+"/data/LD.txt")

    LD = pd.read_csv(index+'/'+gene+'/data/ld.txt', header=None, sep=' ').values
    for i in range(50):
        paintor_pip = pd.read_csv(index + '/' + gene + '/PAINTOR_result_no/C' + str(i + 1) + '.results',
                                  sep=' ').values[:, -1]
        np.savetxt(index + '/' + gene + '/PAINTOR_sets/C' + str(i + 1) + '.sets',
                   compute_paintor_credible_sets(paintor_pip, LD), fmt='%d')

        paintor_pip = pd.read_csv(index + '/' + gene + '/PAINTOR_result_anno/C' + str(i + 1) + '.results',
                                  sep=' ').values[:, -1]
        np.savetxt(index + '/' + gene + '/PAINTOR_anno_sets/C' + str(i + 1) + '.sets',
                   compute_paintor_credible_sets(paintor_pip, LD), fmt='%d')
