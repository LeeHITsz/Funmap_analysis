import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

gene_name = pd.read_csv("Simulation loci.csv")["gene.names"].values

for index in [17]:

    index = str(index+1)

    threshold = np.arange(0.05, 1, 0.1)

    def empirical_FDR(target, true):

        return 1-np.intersect1d(target, true).size/target.size

    def threshold_to_index(target):

        global_fdr = np.cumsum(np.sort(1 - target)) / np.arange(1, target.size + 1)
        casual_num = np.array([global_fdr[global_fdr < i].size for i in threshold])
        return [(np.argpartition(target, -casual_num[i])[-casual_num[i]:]) for i in range(threshold.size)]

    def get_empirical_FDR(target, true):

        return [empirical_FDR(threshold_to_index(target)[i], true) for i in range(threshold.size)]


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
            s_index_pre3 = pd.read_csv(index+"/"+gene+"/CARMA_result_no/C" + str(i + 1) + ".txt.csv", header=None).values[:, -1]
            s_index_pre4 = pd.read_csv(index+"/"+gene+"/CARMA_result_anno/C" + str(i + 1) + ".txt.csv", header=None).values[:, -1]
            s_index_pre5 = pd.read_csv(index+"/"+gene+"/PAINTOR_result_no/C" + str(i + 1) + ".results", sep=' ').values[:, -1]
            s_index_pre6 = pd.read_csv(index+"/"+gene+"/PAINTOR_result_anno/C" + str(i + 1) + ".results", sep=' ').values[:, -1]

            s_total_index = np.append(s_total_index, s_index)
            s_total_index_pre1 = np.append(s_total_index_pre1, s_index_pre1)
            s_total_index_pre2 = np.append(s_total_index_pre2, s_index_pre2)
            s_total_index_pre3 = np.append(s_total_index_pre3, s_index_pre3)
            s_total_index_pre4 = np.append(s_total_index_pre4, s_index_pre4)
            s_total_index_pre5 = np.append(s_total_index_pre5, s_index_pre5)
            s_total_index_pre6 = np.append(s_total_index_pre6, s_index_pre6)

    s_total_index = np.nonzero(s_total_index)
    s_total_eFDR1 = get_empirical_FDR(s_total_index_pre1, s_total_index)
    s_total_eFDR2 = get_empirical_FDR(s_total_index_pre2, s_total_index)
    s_total_eFDR3 = get_empirical_FDR(s_total_index_pre3, s_total_index)
    s_total_eFDR4 = get_empirical_FDR(s_total_index_pre4, s_total_index)
    s_total_eFDR5 = get_empirical_FDR(s_total_index_pre5, s_total_index)
    s_total_eFDR6 = get_empirical_FDR(s_total_index_pre6, s_total_index)

    plt.figure(figsize=(22, 4))
    plt.subplot(164)
    plt.xlabel('FDR level', fontsize=14)
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.scatter(threshold, s_total_eFDR1, color='#315CB5', label="SuSiE")
    plt.plot([0, 1], [0, 1], color='m', linestyle='--')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=11, loc='lower right')

    plt.subplot(161)
    plt.xlabel('FDR level', fontsize=14)
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.ylabel('Empirical FDR', fontsize=14)
    plt.scatter(threshold, s_total_eFDR2, color='#9B2E2B', label="Funmap")
    plt.plot([0, 1], [0, 1], color='m', linestyle='--')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=11, loc='lower right')

    plt.subplot(165)
    plt.xlabel('FDR level', fontsize=14)
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.scatter(threshold, s_total_eFDR3, color='#459CD7', label="CARMA")
    plt.plot([0, 1], [0, 1], color='m', linestyle='--')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=11, loc='lower right')

    plt.subplot(162)
    plt.xlabel('FDR level', fontsize=14)
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.scatter(threshold, s_total_eFDR4, color='#E2533D', label="CARMA with annotation")
    plt.plot([0, 1], [0, 1], color='m', linestyle='--')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=11, loc='lower right')

    plt.subplot(166)
    plt.xlabel('FDR level', fontsize=14)
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.scatter(threshold, s_total_eFDR5, color='c', label='PAINTOR')
    plt.plot([0, 1], [0, 1], color='m', linestyle='--')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=11, loc='lower right')

    plt.subplot(163)
    plt.xlabel('FDR level', fontsize=14)
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.scatter(threshold, s_total_eFDR6, color='gold', label='PAINTOR with annotation')
    plt.plot([0, 1], [0, 1], color='m', linestyle='--')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=11, loc='lower right')

    plt.subplots_adjust(left=0.03, right=0.98, bottom=0.14, top=0.97)
    # plt.savefig(index+"/"+"global FDR.pdf")
    # plt.savefig(index+"/"+"global FDR"+index+".pdf", dpi=300)
    plt.savefig(index+"/"+"global FDR"+index+".png", dpi=300)
    # plt.show()
