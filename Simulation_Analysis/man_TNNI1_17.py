import os
import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy import stats
from numpy import array


index = '18'
gene = 'TNNI1'
Ci = '17'

appear_width = 13
mask_width = 3
true_snps = np.where(pd.read_csv(index+'/'+gene+'/true_SNP/C'+Ci+'_true.txt', sep='\t', header=None).values[:,1]==1)[0]


def get_carma_sets(data):

    n_sets = np.max(data)
    carma_sets = {}

    if n_sets:
        for l in range(1, n_sets+1):
            carma_sets['L'+str(l)] = np.where(data == l)[0]
        return carma_sets
    else:
        return None


pip_susie = pd.read_csv(index+'/'+gene+'/susie_result_no/C'+Ci+'.txt.csv', header=None).values[:, 1].astype(float)
pip_funmap = pd.read_csv(index+'/'+gene+'/Funmap_result/C'+Ci+'.txt.csv', header=None).values[:, 1].astype(float)
pip_carma = pd.read_csv(index+'/'+gene+'/CARMA_result_no/C'+Ci+'.txt.csv', header=None).values[:, 1].astype(float)
pip_carma_anno = pd.read_csv(index+'/'+gene+'/CARMA_result_anno/C'+Ci+'.txt.csv', header=None).values[:, 1].astype(float)
pip_paintor = pd.read_csv(index+'/'+gene+'/PAINTOR_result_no/C'+Ci+'.results', sep=' ').values[:, 1].astype(float)
pip_paintor_anno = pd.read_csv(index+'/'+gene+'/PAINTOR_result_anno/C'+Ci+'.results', sep=' ').values[:, 1].astype(float)
zscore = pd.read_csv(index+'/'+gene+'/PAINTOR_result_anno/C'+Ci+'.results', sep=' ').values[:, 0].astype(float)
pval = stats.norm.sf(np.abs(zscore)) * 2


with open(index+'/'+gene+'/SuSiE_sets/C'+Ci+'.sets', 'r', encoding='utf-8') as f:
    file = f.read()

if re.search(r"'cs':\s*(\{[^}]+\})", file) is None:
    sets_susie = None
else:
    sets_susie = eval(re.search(r"'cs':\s*(\{[^}]+\})", file).group(1))

with open(index+'/'+gene+'/Funmap_sets/C'+Ci+'.sets', 'r', encoding='utf-8') as f:
    file = f.read()

if re.search(r"'cs':\s*(\{[^}]+\})", file) is None:
    sets_funmap = None
else:
    sets_funmap = eval(re.search(r"'cs':\s*(\{[^}]+\})", file).group(1))

sets_CAMRA_data = pd.read_csv(index+'/'+gene+'/CARMA_sets/C'+Ci+'.sets', header=None).values[:, -1]
sets_CAMRA = get_carma_sets(sets_CAMRA_data)

sets_CAMRA_anno_data = pd.read_csv(index+'/'+gene+'/CARMA_anno_sets/C'+Ci+'.sets', header=None).values[:, -1]
sets_CAMRA_anno = get_carma_sets(sets_CAMRA_anno_data)

sets_PAINTOR_data = np.loadtxt(index+'/'+gene+'/PAINTOR_sets/C'+Ci+'.sets', dtype=int)
sets_PAINTOR = get_carma_sets(sets_PAINTOR_data)

sets_PAINTOR_anno_data = np.loadtxt(index+'/'+gene+'/PAINTOR_anno_sets/C'+Ci+'.sets', dtype=int)
sets_PAINTOR_anno = get_carma_sets(sets_PAINTOR_anno_data)

R = pd.read_csv(index+'/'+gene+'/data/ld.txt', header=None, sep=' ').values
R1 = R[true_snps[0]-(appear_width-1)//2:true_snps[0]+(appear_width-1)//2+1,
       true_snps[0]-(appear_width-1)//2:true_snps[0]+(appear_width-1)//2+1]
R2 = R[true_snps[1]-(appear_width-1)//2:true_snps[1]+(appear_width-1)//2+1,
       true_snps[1]-(appear_width-1)//2:true_snps[1]+(appear_width-1)//2+1]

corr = np.zeros((appear_width+(mask_width+1)//2, appear_width+(mask_width+1)//2))
corr[np.tril_indices_from(corr, k=-(mask_width+1)//2)] = R1[np.tril_indices_from(R1)]
corr[np.triu_indices_from(corr, k=(mask_width+1)//2)] = R2[np.triu_indices_from(R2)]
corr = np.abs(corr)

bands_x = np.arange(true_snps[1]-(appear_width-1)//2, true_snps[1]+(appear_width-1)//2+1).astype(str)
bands_x = np.append(np.array([' ']*((mask_width+1)//2)), bands_x)
bands_y = np.arange(true_snps[0]-(appear_width-1)//2, true_snps[0]+(appear_width-1)//2+1).astype(str)
bands_y = np.append(np.array([' ']*((mask_width+1)//2)), bands_y)

mask = np.zeros_like(corr, dtype=np.bool_)
mask[np.tril_indices_from(mask, k=-(mask_width+1)//2)] = True
mask[np.triu_indices_from(mask, k=(mask_width+1)//2)] = True

cmap = 'Blues'
colors = ['g', 'c', 'm', 'y', 'pink', 'purple', 'orange', 'grey', 'darkblue', 'darkcyan',
          'darkgoldenrod', 'darkgray', 'darkgreen', 'darkkhaki', 'darkmagenta', 'darkolivegreen',
          'darkorange', 'darkorchid', 'darkred', 'darksalmon', 'darkseagreen', 'darkslateblue',
          'darkslategray', 'darkturquoise', 'darkviolet']

plt.figure(figsize=(22, 10))
plt.subplot(241)
ax = sns.heatmap(corr, mask=~mask, cmap=cmap, square=True, linewidths=0.5, vmin=0, vmax=1, cbar=True,
                 xticklabels=bands_x, yticklabels=bands_y)
ax.set_xlabel("Region of Casual SNP2", fontsize=14, labelpad=15)
ax.set_ylabel("Region of Casual SNP1", fontsize=14, labelpad=15)
ax.set_xticklabels(ax.get_xticklabels(), fontsize=6)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
ax.tick_params(top=False, left=False)
plt.xticks([])
plt.yticks([])
ax.add_patch(patches.Rectangle(((appear_width-1)//2+(mask_width+1)//2, (appear_width-1)//2), (appear_width+1)//2, 1,
                               linewidth=2, edgecolor='r', facecolor='none'))
ax.add_patch(patches.Rectangle((0, (appear_width-1)//2+(mask_width+1)//2), (appear_width+1)//2, 1,
                               linewidth=2, edgecolor='r', facecolor='none'))
ax.add_patch(patches.Rectangle(((appear_width-1)//2, (appear_width-1)//2+(mask_width+1)//2), 1, (appear_width+1)//2,
                               linewidth=2, edgecolor='r', facecolor='none'))
ax.add_patch(patches.Rectangle(((appear_width-1)//2+(mask_width+1)//2, 0), 1, (appear_width+1)//2,
                               linewidth=2, edgecolor='r', facecolor='none'))

plt.subplot(245)
plt.scatter(np.arange(pip_susie.shape[0]), -np.log10(pval), c='darkblue', s=48)
plt.scatter(np.arange(pip_susie.shape[0])[true_snps], -np.log10(pval)[true_snps], c='r', s=48)
plt.text(true_snps[0]+1, -np.log10(pval)[true_snps[0]]+1, 'SNP1', color='r', fontsize=14)
plt.text(true_snps[1]+1, -np.log10(pval)[true_snps[1]]+1, 'SNP2', color='r', fontsize=14)
plt.axhline(y=-np.log10(5e-8), color='r', linestyle='--', linewidth=1)
plt.xticks(np.arange(0, pip_susie.shape[0], 200), fontsize=14)
plt.xlabel('Variant', fontsize=14)
plt.ylabel(r'$-\log_{10}(p)$', fontsize=14)
plt.title('GWAS', fontsize=18, fontweight="bold")
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.yticks(fontsize=14)

plt.subplot(242)
plt.scatter(np.arange(pip_funmap.shape[0]), pip_funmap, c='darkblue', s=48)
if sets_funmap is not None:
    nonsets_index = np.concatenate(list(sets_funmap.values()))
    for i, group in enumerate(sets_funmap.values()):
        color = colors[i]
        for idx in group:
            plt.scatter(np.arange(pip_funmap.shape[0])[idx], pip_funmap[idx], c='darkblue', s=48,
                        edgecolors=color, linewidths=6)
    plt.scatter(np.arange(pip_funmap.shape[0])[nonsets_index], pip_funmap[nonsets_index],
                c='darkblue', s=48)
else:
    plt.scatter(np.arange(pip_funmap.shape[0]), pip_funmap, c='darkblue', s=48)
plt.scatter(np.arange(pip_funmap.shape[0])[true_snps], pip_funmap[true_snps], c='r', s=48)
plt.ylabel('PIP', fontsize=14)
plt.title('Funmap', fontsize=18, fontweight="bold")
plt.ylim([-0.05, 1.05])
plt.gca().axes.xaxis.set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.yticks(fontsize=14)

plt.subplot(246)
plt.scatter(np.arange(pip_susie.shape[0]), pip_susie, c='darkblue', s=48)
if sets_susie is not None:
    nonsets_index = np.concatenate(list(sets_susie.values()))
    for i, group in enumerate(sets_susie.values()):
        color = colors[i]
        for idx in group:
            plt.scatter(np.arange(pip_susie.shape[0])[idx], pip_susie[idx], c='darkblue', s=48,
                        edgecolors=color, linewidths=6)
    plt.scatter(np.arange(pip_susie.shape[0])[nonsets_index], pip_susie[nonsets_index],
                c='darkblue', s=48)
else:
    plt.scatter(np.arange(pip_susie.shape[0]), pip_susie, c='darkblue', s=48)
plt.scatter(np.arange(pip_susie.shape[0])[true_snps], pip_susie[true_snps], c='r', s=48)
plt.xticks(np.arange(0, pip_susie.shape[0], 200), fontsize=14)
plt.xlabel('Variant', fontsize=14)
plt.ylabel('PIP', fontsize=14)
plt.title('SuSiE', fontsize=18, fontweight="bold")
plt.ylim([-0.05, 1.05])
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.yticks(fontsize=14)

plt.subplot(243)
plt.scatter(np.arange(pip_carma_anno.shape[0]), pip_carma_anno, c='darkblue', s=48)
if sets_CAMRA_anno is not None:
    nonsets_index = np.concatenate(list(sets_CAMRA_anno.values()))
    for i, group in enumerate(sets_CAMRA_anno.values()):
        color = colors[i]
        for idx in group:
            plt.scatter(np.arange(pip_carma_anno.shape[0])[idx], pip_carma_anno[idx], c='darkblue', s=48,
                        edgecolors=color, linewidths=6)
    plt.scatter(np.arange(pip_carma_anno.shape[0])[nonsets_index], pip_carma_anno[nonsets_index],
                c='darkblue', s=48)
else:
    plt.scatter(np.arange(pip_carma_anno.shape[0]), pip_carma_anno, c='darkblue', s=48)
plt.scatter(np.arange(pip_carma_anno.shape[0])[true_snps], pip_carma_anno[true_snps], c='r', s=48)
plt.title('CARMA with annotation', fontsize=18, fontweight="bold")
plt.ylabel('PIP', fontsize=14)
plt.ylim([-0.05, 1.05])
plt.gca().axes.xaxis.set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.yticks(fontsize=14)

plt.subplot(247)
plt.scatter(np.arange(pip_carma.shape[0]), pip_carma, c='darkblue', s=48)
if sets_CAMRA is not None:
    nonsets_index = np.concatenate(list(sets_CAMRA.values()))
    for i, group in enumerate(sets_CAMRA.values()):
        color = colors[i]
        for idx in group:
            plt.scatter(np.arange(pip_carma.shape[0])[idx], pip_carma[idx], c='darkblue', s=48,
                        edgecolors=color, linewidths=6)
    plt.scatter(np.arange(pip_carma.shape[0])[nonsets_index], pip_carma[nonsets_index],
                c='darkblue', s=48)
else:
    plt.scatter(np.arange(pip_carma.shape[0]), pip_carma, c='darkblue', s=48)
plt.scatter(np.arange(pip_carma.shape[0])[true_snps], pip_carma[true_snps], c='r', s=48)
plt.xticks(np.arange(0, pip_susie.shape[0], 200), fontsize=14)
plt.xlabel('Variant', fontsize=14)
plt.title('CARMA', fontsize=18, fontweight="bold")
plt.ylabel('PIP', fontsize=14)
plt.ylim([-0.05, 1.05])
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.yticks(fontsize=14)

plt.subplot(244)
plt.scatter(np.arange(pip_paintor_anno.shape[0]), pip_paintor_anno, c='darkblue', s=48)
if sets_PAINTOR_anno is not None:
    nonsets_index = np.concatenate(list(sets_PAINTOR_anno.values()))
    for i, group in enumerate(sets_PAINTOR_anno.values()):
        color = colors[i]
        for idx in group:
            plt.scatter(np.arange(pip_paintor_anno.shape[0])[idx], pip_paintor_anno[idx], c='darkblue', s=48,
                        edgecolors=color, linewidths=6)
    plt.scatter(np.arange(pip_paintor_anno.shape[0])[nonsets_index], pip_paintor_anno[nonsets_index],
                c='darkblue', s=48)
else:
    plt.scatter(np.arange(pip_paintor_anno.shape[0]), pip_paintor_anno, c='darkblue', s=48)
plt.scatter(np.arange(pip_paintor_anno.shape[0])[true_snps], pip_paintor_anno[true_snps], c='r', s=48)
plt.title('PAINTOR with annotation', fontsize=18, fontweight="bold")
plt.ylabel('PIP', fontsize=14)
plt.ylim([-0.05, 1.05])
plt.gca().axes.xaxis.set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.yticks(fontsize=14)

plt.subplot(248)
plt.scatter(np.arange(pip_paintor.shape[0]), pip_paintor, c='darkblue', s=48)
if sets_PAINTOR is not None:
    nonsets_index = np.concatenate(list(sets_PAINTOR.values()))
    for i, group in enumerate(sets_PAINTOR.values()):
        color = colors[i]
        for idx in group:
            plt.scatter(np.arange(pip_paintor.shape[0])[idx], pip_paintor[idx], c='darkblue', s=48,
                        edgecolors=color, linewidths=6)
    plt.scatter(np.arange(pip_paintor.shape[0])[nonsets_index], pip_paintor[nonsets_index], c='darkblue', s=48)
else:
    plt.scatter(np.arange(pip_paintor.shape[0]), pip_paintor, c='darkblue', s=48)
plt.scatter(np.arange(pip_paintor.shape[0])[true_snps], pip_paintor[true_snps], c='r', s=48)
plt.xticks(np.arange(0, pip_susie.shape[0], 200), fontsize=14)
plt.xlabel('Variant', fontsize=14)
plt.title('PAINTOR', fontsize=18, fontweight="bold")
plt.ylabel('PIP', fontsize=14)
plt.ylim([-0.05, 1.05])
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.yticks(fontsize=14)

plt.subplots_adjust(left=0.05, right=0.98, bottom=0.08, top=0.95)
# plt.savefig(index+'/'+gene+'/man_C' + Ci + '.pdf')
plt.savefig(index+'/'+gene+'/man_C' + Ci + '.png', dpi=300)
# plt.show()
