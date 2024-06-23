import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


gene_name = pd.read_csv("Simulation loci.csv")["gene.names"].values

num_casual = 2*50*10

# for index in [8, 9, 10, 14, 15, 16, 20, 21, 22]:
# for index in [11, 12, 13, 17, 18, 19, 23, 24, 25]:
for index in [17]:

    index = str(index+1)
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
                                    header=None).values[:, -1] * s_index
            s_index_pre2 = pd.read_csv(index + "/" + gene + "/Funmap_result/C" + str(i + 1) + ".txt.csv",
                                    header=None).values[:, -1] * s_index
            s_index_pre3 = pd.read_csv(index + "/" + gene + "/CARMA_result_no/C" + str(i + 1) + ".txt.csv",
                                    header=None).values[:, -1] * s_index
            s_index_pre4 = pd.read_csv(index + "/" + gene + "/CARMA_result_anno/C" + str(i + 1) + ".txt.csv",
                                    header=None).values[:, -1] * s_index
            s_index_pre5 = pd.read_csv(index + "/" + gene + "/PAINTOR_result_no/C" + str(i + 1) + ".results",
                                    sep=' ').values[:, -1] * s_index
            s_index_pre6 = pd.read_csv(index + "/" + gene + "/PAINTOR_result_anno/C" + str(i + 1) + ".results",
                                    sep=' ').values[:, -1] * s_index

            s_total_index_pre1 = np.append(s_total_index_pre1, s_index_pre1)
            s_total_index_pre2 = np.append(s_total_index_pre2, s_index_pre2)
            s_total_index_pre3 = np.append(s_total_index_pre3, s_index_pre3)
            s_total_index_pre4 = np.append(s_total_index_pre4, s_index_pre4)
            s_total_index_pre5 = np.append(s_total_index_pre5, s_index_pre5)
            s_total_index_pre6 = np.append(s_total_index_pre6, s_index_pre6)

    # height_80_susie = [np.sum(s_total_index_pre1 > 0.80)/num_casual, np.sum(s_total_index_pre2 > 0.80)/num_casual]
    # height_80_carma = [np.sum(s_total_index_pre3 > 0.80)/num_casual, np.sum(s_total_index_pre4 > 0.80)/num_casual]
    # height_80_paintor = [np.sum(s_total_index_pre5 > 0.80)/num_casual, np.sum(s_total_index_pre6 > 0.80)/num_casual]
    height_90_susie = [np.sum(s_total_index_pre1 > 0.90)/num_casual, np.sum(s_total_index_pre2 > 0.90)/num_casual]
    height_90_carma = [np.sum(s_total_index_pre3 > 0.90)/num_casual, np.sum(s_total_index_pre4 > 0.90)/num_casual]
    height_90_paintor = [np.sum(s_total_index_pre5 > 0.90)/num_casual, np.sum(s_total_index_pre6 > 0.90)/num_casual]
    height_95_susie = [np.sum(s_total_index_pre1 > 0.95)/num_casual, np.sum(s_total_index_pre2 > 0.95)/num_casual]
    height_95_carma = [np.sum(s_total_index_pre3 > 0.95)/num_casual, np.sum(s_total_index_pre4 > 0.95)/num_casual]
    height_95_paintor = [np.sum(s_total_index_pre5 > 0.95)/num_casual, np.sum(s_total_index_pre6 > 0.95)/num_casual]
    height_99_susie = [np.sum(s_total_index_pre1 > 0.99)/num_casual, np.sum(s_total_index_pre2 > 0.99)/num_casual]
    height_99_carma = [np.sum(s_total_index_pre3 > 0.99)/num_casual, np.sum(s_total_index_pre4 > 0.99)/num_casual]
    height_99_paintor = [np.sum(s_total_index_pre5 > 0.99)/num_casual, np.sum(s_total_index_pre6 > 0.99)/num_casual]

    barWidth = 0.4
    x_position = [1, 2.4]
    x_position_fmt = ['Including annotations', 'No annotations']
    bbox = dict(facecolor='grey', edgecolor='k', boxstyle='square', alpha=0.2)

    plt.figure(figsize=(10, 3.5))
    # plt.subplot(141)
    # plt.bar(1, height_80_susie[1], color='#9B2E2B', label='Funmap', width=barWidth)
    # plt.bar(1 + barWidth, height_80_carma[1], color='#E2533D', label='CARMA+anno', width=barWidth)
    # plt.bar(1 + 2 * barWidth, height_80_paintor[1], color='#F9E4A9', label='PAINTOR+anno', width=barWidth)
    # plt.bar(2.4, height_80_susie[0], color='#315CB5', label='SuSiE', width=barWidth)
    # plt.bar(2.4 + barWidth, height_80_carma[0], color='#459CD7', label='CARMA', width=barWidth)
    # plt.bar(2.4 + 2 * barWidth, height_80_paintor[0], color='#B2DFE3', label='PAINTOR', width=barWidth)
    # plt.text(1 + 0 * barWidth, height_80_susie[1], '%.2f' % height_80_susie[1], ha='center', va='bottom')
    # plt.text(1 + 1 * barWidth, height_80_carma[1], '%.2f' % height_80_carma[1], ha='center', va='bottom')
    # plt.text(1 + 2 * barWidth, height_80_paintor[1], '%.2f' % height_80_paintor[1], ha='center', va='bottom')
    # plt.text(2.4 + 0 * barWidth, height_80_susie[0], '%.2f' % height_80_susie[0], ha='center', va='bottom')
    # plt.text(2.4 + 1 * barWidth, height_80_carma[0], '%.2f' % height_80_carma[0], ha='center', va='bottom')
    # plt.text(2.4 + 2 * barWidth, height_80_paintor[0], '%.2f' % height_80_paintor[0], ha='center', va='bottom')
    # plt.xticks([i + 0.4 for i in x_position], x_position_fmt)
    # plt.ylim([0, 1])
    # plt.ylabel('Power', rotation=90)
    # plt.legend(bbox_to_anchor=(3.8, -0.08), ncol=6)
    # plt.title("PIP > 0.80", bbox=bbox)

    plt.subplot(131)
    plt.bar(1, height_90_susie[1], color='#9B2E2B', label='Funmap', width=barWidth)
    plt.bar(1 + barWidth, height_90_carma[1], color='#E2533D', label='CARMA+anno', width=barWidth)
    plt.bar(1 + 2 * barWidth, height_90_paintor[1], color='#F9E4A9', label='PAINTOR+anno', width=barWidth)
    plt.bar(2.4, height_90_susie[0], color='#315CB5', label='SuSiE', width=barWidth)
    plt.bar(2.4 + barWidth, height_90_carma[0], color='#459CD7', label='CARMA', width=barWidth)
    plt.bar(2.4 + 2 * barWidth, height_90_paintor[0], color='#B2DFE3', label='PAINTOR', width=barWidth)
    plt.text(1 + 0 * barWidth, height_90_susie[1], '%.2f' % height_90_susie[1], ha='center', va='bottom')
    plt.text(1 + 1 * barWidth, height_90_carma[1], '%.2f' % height_90_carma[1], ha='center', va='bottom')
    plt.text(1 + 2 * barWidth, height_90_paintor[1], '%.2f' % height_90_paintor[1], ha='center', va='bottom')
    plt.text(2.4 + 0 * barWidth, height_90_susie[0], '%.2f' % height_90_susie[0], ha='center', va='bottom')
    plt.text(2.4 + 1 * barWidth, height_90_carma[0], '%.2f' % height_90_carma[0], ha='center', va='bottom')
    plt.text(2.4 + 2 * barWidth, height_90_paintor[0], '%.2f' % height_90_paintor[0], ha='center', va='bottom')
    plt.xticks([i + 0.4 for i in x_position], x_position_fmt)
    plt.ylim([0, 1])
    plt.ylabel('Power', rotation=90)
    plt.title("PIP > 0.90", bbox=bbox)
    plt.legend(bbox_to_anchor=(3.1, -0.08), ncol=6)

    plt.subplot(132)
    plt.bar(1, height_95_susie[1], color='#9B2E2B', label='Funmap', width=barWidth)
    plt.bar(1 + barWidth, height_95_carma[1], color='#E2533D', label='CARMA+anno', width=barWidth)
    plt.bar(1 + 2 * barWidth, height_95_paintor[1], color='#F9E4A9', label='PAINTOR+anno', width=barWidth)
    plt.bar(2.4, height_95_susie[0], color='#315CB5', label='SuSiE', width=barWidth)
    plt.bar(2.4 + barWidth, height_95_carma[0], color='#459CD7', label='CARMA', width=barWidth)
    plt.bar(2.4 + 2 * barWidth, height_95_paintor[0], color='#B2DFE3', label='PAINTOR', width=barWidth)
    plt.text(1 + 0 * barWidth, height_95_susie[1], '%.2f' % height_95_susie[1], ha='center', va='bottom')
    plt.text(1 + 1 * barWidth, height_95_carma[1], '%.2f' % height_95_carma[1], ha='center', va='bottom')
    plt.text(1 + 2 * barWidth, height_95_paintor[1], '%.2f' % height_95_paintor[1], ha='center', va='bottom')
    plt.text(2.4 + 0 * barWidth, height_95_susie[0], '%.2f' % height_95_susie[0], ha='center', va='bottom')
    plt.text(2.4 + 1 * barWidth, height_95_carma[0], '%.2f' % height_95_carma[0], ha='center', va='bottom')
    plt.text(2.4 + 2 * barWidth, height_95_paintor[0], '%.2f' % height_95_paintor[0], ha='center', va='bottom')
    plt.xticks([i + 0.4 for i in x_position], x_position_fmt)
    plt.ylim([0, 1])
    plt.title("PIP > 0.95", bbox=bbox)

    plt.subplot(133)
    plt.bar(1, height_99_susie[1], color='#9B2E2B', label='Funmap', width=barWidth)
    plt.bar(1 + barWidth, height_99_carma[1], color='#E2533D', label='CARMA+anno', width=barWidth)
    plt.bar(1 + 2 * barWidth, height_99_paintor[1], color='#F9E4A9', label='PAINTOR+anno', width=barWidth)
    plt.bar(2.4, height_99_susie[0], color='#315CB5', label='SuSiE', width=barWidth)
    plt.bar(2.4 + barWidth, height_99_carma[0], color='#459CD7', label='CARMA', width=barWidth)
    plt.bar(2.4 + 2 * barWidth, height_99_paintor[0], color='#B2DFE3', label='PAINTOR', width=barWidth)
    plt.text(1 + 0 * barWidth, height_99_susie[1], '%.2f' % height_99_susie[1], ha='center', va='bottom')
    plt.text(1 + 1 * barWidth, height_99_carma[1], '%.2f' % height_99_carma[1], ha='center', va='bottom')
    plt.text(1 + 2 * barWidth, height_99_paintor[1], '%.2f' % height_99_paintor[1], ha='center', va='bottom')
    plt.text(2.4 + 0 * barWidth, height_99_susie[0], '%.2f' % height_99_susie[0], ha='center', va='bottom')
    plt.text(2.4 + 1 * barWidth, height_99_carma[0], '%.2f' % height_99_carma[0], ha='center', va='bottom')
    plt.text(2.4 + 2 * barWidth, height_99_paintor[0], '%.2f' % height_99_paintor[0], ha='center', va='bottom')
    plt.xticks([i + 0.4 for i in x_position], x_position_fmt)
    plt.ylim([0, 1])
    plt.title("PIP > 0.99", bbox=bbox)

    plt.subplots_adjust(left=0.06, right=0.98, top=0.90, bottom=0.15)
    # plt.savefig(index+"/"+"double bar"+index+".pdf")
    plt.savefig(index+"/"+"double bar"+index+".svg")
    # plt.savefig(index+"/"+"double bar.png", dpi=300)
