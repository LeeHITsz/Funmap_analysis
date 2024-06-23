import numpy as np
import pandas as pd
from matplotlib import rcParams
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, average_precision_score

gene_name = pd.read_csv("Simulation loci.csv")["gene.names"].values

index = '26'

s_total_index = np.array([])
s_total_index_pre1 = np.array([])
s_total_index_pre2 = np.array([])
s_total_index_pre3 = np.array([])
s_total_index_pre4 = np.array([])
s_total_index_pre5 = np.array([])
s_total_index_pre6 = np.array([])

for gene in gene_name[:10]:

    for i in range(50):

        s_index = pd.read_csv(index+"/"+gene+"/true_SNP/C"+str(i+1)+"_true.txt", header=None, sep='\t').values[:, -1]

        s_index_pre1 = pd.read_csv(index+"/"+gene+"/susie_result_no/C" + str(i + 1) + ".txt.csv", header=None).values[:, -1]
        s_index_pre2 = pd.read_csv(index+"/"+gene+"/Funmap_result/C" + str(i + 1) + ".txt.csv", header=None).values[:, -1]
        s_index_pre3 = pd.read_csv(index+"/"+gene+"/CARMA_result_no/C" + str(i + 1) + ".txt.csv", header=None).values[:, -1]
        s_index_pre4 = pd.read_csv(index+"/"+gene+"/CARMA_result_anno/C" + str(i + 1) + ".txt.csv", header=None).values[:, -1]
        s_index_pre5 = pd.read_csv(index+"/"+gene+"/PAINTOR_result_no/C" + str(i + 1) + ".results", sep=' ').values[:, -1]
        s_index_pre6 = pd.read_csv(index+"/"+gene+"/PAINTOR_result_anno/C" + str(i + 1) + ".results", sep=' ').values[:, -1]
 
        if np.isnan(s_index_pre2.astype(float)).any():
            continue

        s_total_index = np.append(s_total_index, s_index)
        s_total_index_pre1 = np.append(s_total_index_pre1, s_index_pre1)
        s_total_index_pre2 = np.append(s_total_index_pre2, s_index_pre2)
        s_total_index_pre3 = np.append(s_total_index_pre3, s_index_pre3)
        s_total_index_pre4 = np.append(s_total_index_pre4, s_index_pre4)
        s_total_index_pre5 = np.append(s_total_index_pre5, s_index_pre5)
        s_total_index_pre6 = np.append(s_total_index_pre6, s_index_pre6)

s_total_index = s_total_index.astype(int)
s_total_index_pre1 = s_total_index_pre1.astype(float)
s_total_index_pre2 = s_total_index_pre2.astype(float)
s_total_index_pre3 = s_total_index_pre3.astype(float)
s_total_index_pre4 = s_total_index_pre4.astype(float)
s_total_index_pre5 = s_total_index_pre5.astype(float)
s_total_index_pre6 = s_total_index_pre6.astype(float)

precision1, recall1, thresholds1 = precision_recall_curve(s_total_index, s_total_index_pre1)
PRC1 = average_precision_score(s_total_index, s_total_index_pre1)
precision2, recall2, thresholds2 = precision_recall_curve(s_total_index, s_total_index_pre2)
PRC2 = average_precision_score(s_total_index, s_total_index_pre2)
precision3, recall3, thresholds3 = precision_recall_curve(s_total_index, s_total_index_pre3)
PRC3 = average_precision_score(s_total_index, s_total_index_pre3)
precision4, recall4, thresholds4 = precision_recall_curve(s_total_index, s_total_index_pre4)
PRC4 = average_precision_score(s_total_index, s_total_index_pre4)
precision5, recall5, thresholds5 = precision_recall_curve(s_total_index, s_total_index_pre5)
PRC5 = average_precision_score(s_total_index, s_total_index_pre5)
precision6, recall6, thresholds6 = precision_recall_curve(s_total_index, s_total_index_pre6)
PRC6 = average_precision_score(s_total_index, s_total_index_pre6)

thresholds1 = np.append(thresholds1, 1.01)
thresholds2 = np.append(thresholds2, 1.01)
thresholds3 = np.append(thresholds3, 1.01)
thresholds4 = np.append(thresholds4, 1.01)
thresholds5 = np.append(thresholds5, 1.01)
thresholds6 = np.append(thresholds6, 1.01)

index_095_1 = np.argmin(np.abs(thresholds1 - 0.95))
index_095_2 = np.argmin(np.abs(thresholds2 - 0.95))
index_095_3 = np.argmin(np.abs(thresholds3 - 0.95))
index_095_4 = np.argmin(np.abs(thresholds4 - 0.95))
index_095_5 = np.argmin(np.abs(thresholds5 - 0.95))
index_095_6 = np.argmin(np.abs(thresholds6 - 0.95))

plt.figure(figsize=(15, 5))
plt.subplot(131)
plt.xlabel('FDR', fontsize=14)
plt.ylabel('power', fontsize=14)
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.step(1-precision2, recall2, color='#9B2E2B', label='Funmap')
plt.step(1-precision4, recall4, color='#E2533D', label='CARMA+anno')
plt.step(1-precision6, recall6, color='#F9E4A9', label='PAINTOR+anno')
plt.step(1-precision1, recall1, color='#315CB5', label='SuSiE')
plt.step(1-precision3, recall3, color='#459CD7', label='CARMA')
plt.step(1-precision5, recall5, color='#B2DFE3', label='PAINTOR')

plt.scatter(1-precision1[index_095_1], recall1[index_095_1], color='k', facecolors='none')
plt.text(1-precision1[index_095_1]+0.05, recall1[index_095_1]-0.05, '0.95', color='#315CB5')
plt.scatter(1-precision2[index_095_2], recall2[index_095_2], color='k', facecolors='none')
plt.text(1-precision2[index_095_2]+0.05, recall2[index_095_2]-0.05, '0.95', color='#9B2E2B')
plt.scatter(1-precision3[index_095_3], recall3[index_095_3], color='k', facecolors='none')
plt.text(1-precision3[index_095_3]+0.05, recall3[index_095_3]-0.05, '0.95', color='#459CD7')
plt.scatter(1-precision4[index_095_4], recall4[index_095_4], color='k', facecolors='none')
plt.text(1-precision4[index_095_4]+0.05, recall4[index_095_4]-0.05, '0.95', color='#E2533D')
plt.scatter(1-precision5[index_095_5], recall5[index_095_5], color='k', facecolors='none')
plt.text(1-precision5[index_095_5]+0.05, recall5[index_095_5]-0.05, '0.95', color='#B2DFE3')
plt.scatter(1-precision6[index_095_6], recall6[index_095_6], color='k', facecolors='none')
plt.text(1-precision6[index_095_6]+0.05, recall6[index_095_6]-0.05, '0.95', color='#F9E4A9')
plt.plot([0, 1], [0, 1], color='m', linestyle='--')
plt.legend(bbox_to_anchor=(2.6, -0.12), ncol=6)
# plt.legend(fontsize=8)

plt.subplot(132)
plt.xlabel('thresholds', fontsize=14)
plt.ylabel('FDR', fontsize=14)
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.step(thresholds2, 1-precision2, color='#9B2E2B', label='Funmap')
plt.step(thresholds4, 1-precision4, color='#E2533D', label='CARMA+anno')
plt.step(thresholds6, 1-precision6, color='#F9E4A9', label='PAINTOR+anno')
plt.step(thresholds1, 1-precision1, color='#315CB5', label='SuSiE')
plt.step(thresholds3, 1-precision3, color='#459CD7', label='CARMA')
plt.step(thresholds5, 1-precision5, color='#B2DFE3', label='PAINTOR')
# plt.legend(fontsize=8)

plt.subplot(133)
plt.xlabel('thresholds', fontsize=14)
plt.ylabel('power', fontsize=14)
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.step(thresholds2, recall2, color='#9B2E2B', label='Funmap')
plt.step(thresholds4, recall4, color='#E2533D', label='CARMA+anno')
plt.step(thresholds6, recall6, color='#F9E4A9', label='PAINTOR+anno')
plt.step(thresholds1, recall1, color='#315CB5', label='SuSiE')
plt.step(thresholds3, recall3, color='#459CD7', label='CARMA')
plt.step(thresholds5, recall5, color='#B2DFE3', label='PAINTOR')
# plt.legend(fontsize=8)

plt.subplots_adjust(left=0.05, right=0.98, top=0.95, bottom=0.18)
plt.savefig(index+"/"+"PR curves"+index+".pdf", dpi=300)
# plt.show()
