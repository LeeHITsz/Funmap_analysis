import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

gene_name = pd.read_csv("Simulation loci.csv")["gene.names"].values

for index in [20, 21, 22, 23, 24, 25]:

    index = str(index+1)

    s_total_index = np.array([])
    s_total_index_pre1 = np.array([])
    s_total_index_pre2 = np.array([])
    s_total_index_pre3 = np.array([])
    s_total_index_pre4 = np.array([])
    s_total_index_pre5 = np.array([])
    s_total_index_pre6 = np.array([])

    for gene in gene_name[:10]:

        for i in range(50):

            s_index = pd.read_csv(index + "/" + gene + "/true_SNP/C" + str(i + 1) + "_true.txt", header=None,
                                sep='\t').values[:, -1]

            s_index_pre1 = pd.read_csv(index + "/" + gene + "/susie_result_no/C" + str(i + 1) + ".txt.csv",
                                    header=None).values[:, -1]
            s_index_pre2 = pd.read_csv(index + "/" + gene + "/Funmap_result/C" + str(i + 1) + ".txt.csv",
                                    header=None).values[:, -1]
            s_index_pre4 = pd.read_csv(index+"/"+gene+"/CARMA_result_anno/C" + str(i + 1) + ".txt.csv", header=None).values[:, -1]
            s_index_pre6 = pd.read_csv(index+"/"+gene+"/PAINTOR_result_anno/C" + str(i + 1) + ".results", sep=' ').values[:, -1]

            s_total_index = np.append(s_total_index, s_index)
            s_total_index_pre1 = np.append(s_total_index_pre1, s_index_pre1)
            s_total_index_pre2 = np.append(s_total_index_pre2, s_index_pre2)
            s_total_index_pre4 = np.append(s_total_index_pre4, s_index_pre4)
            s_total_index_pre6 = np.append(s_total_index_pre6, s_index_pre6)

    idx_casual = np.where(s_total_index == 1)
    idx_noncasual = np.where(s_total_index == 0)

    _, ax = plt.subplots(1, 3, figsize=(15, 5))
    ax[0].set_xlabel('Funmap PIP', fontsize=14)
    ax[0].set_ylabel('SuSiE PIP', rotation=90, fontsize=14)
    ax[0].set_xlim([-0.05, 1.05])
    ax[0].set_ylim([-0.05, 1.05])
    ax[0].scatter(s_total_index_pre2[idx_casual], s_total_index_pre1[idx_casual], color='r', label="casual", s=5)
    ax[0].scatter(s_total_index_pre2[idx_noncasual], s_total_index_pre1[idx_noncasual], color='grey', label="non-casual", s=5, alpha=0.4, linewidths=0)
    ax[0].plot([0, 1], [0, 1], color='m', linestyle='--')
    ax[0].legend(loc='lower right')

    ax[1].set_xlabel('Funmap PIP', fontsize=14)
    ax[1].set_ylabel('CARMA with annotation PIP', rotation=90, fontsize=14)
    ax[1].set_xlim([-0.05, 1.05])
    ax[1].set_ylim([-0.05, 1.05])
    ax[1].scatter(s_total_index_pre2[idx_casual], s_total_index_pre4[idx_casual], color='r', label="casual", s=5)
    ax[1].scatter(s_total_index_pre2[idx_noncasual], s_total_index_pre4[idx_noncasual], color='grey', label="non-casual", s=5, alpha=0.4, linewidths=0)
    ax[1].plot([0, 1], [0, 1], color='m', linestyle='--')
    ax[1].legend(loc='lower right')

    ax[2].set_xlabel('Funmap PIP', fontsize=14)
    ax[2].set_ylabel('PAINTOR with annotation PIP', rotation=90, fontsize=14)
    ax[2].set_xlim([-0.05, 1.05])
    ax[2].set_ylim([-0.05, 1.05])
    ax[2].scatter(s_total_index_pre2[idx_casual], s_total_index_pre6[idx_casual], color='r', label="casual", s=5)
    ax[2].scatter(s_total_index_pre2[idx_noncasual], s_total_index_pre6[idx_noncasual], color='grey', label="non-casual", s=5, alpha=0.4, linewidths=0)
    ax[2].plot([0, 1], [0, 1], color='m', linestyle='--')
    ax[2].legend(loc='lower right')

    plt.subplots_adjust(left=0.05, right=0.98, bottom=0.12, top=0.95, wspace=0.28)
    plt.savefig(index+"/"+"scat"+index+".png", dpi=300)
