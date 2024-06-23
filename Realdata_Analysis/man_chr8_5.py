import os
import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from numpy import array
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


Chr = '8'
Pos = 5
trait = 'Cho'


def get_ldname(chr, pos):

    return 'chr' + str(chr) + '_' + str(int(pos * 1e6 + 1)) + '_' + str(int((pos + 3) * 1e6 + 1))


def get_carma_sets(data):

    n_sets = np.max(data)
    carma_sets = {}

    if n_sets:
        for l in range(1, n_sets+1):
            carma_sets['L'+str(l)] = np.where(data == l)[0]
        return carma_sets
    else:
        return None

gene = pd.read_csv('geneChr8.csv')
gene = gene[['geneName', 'txStart', 'txEnd']].values

y_list = [5, 4, 3, 2, 1, 4.8, 3.8, 2.8, 1.8, 0.8, 4.6, 3.6, 2.6, 1.6, 0.6, 4.3, 3.3, 2.3, 1.3, 0.3, 4.0, 3.0, 2.0]

chr_data = pd.read_csv('chr' + Chr + '/zscore/Cholesterol.csv.gz').drop_duplicates(subset=['rsid'], keep='first')

pip_susie = pd.read_csv('chr' + Chr + '/SuSiE_result/' + get_ldname(Chr, Pos) + '_' + trait + '.csv',
                        header=None)
pip_funmap = pd.read_csv('chr' + Chr + '/Funmap_result/' + get_ldname(Chr, Pos) + '_' + trait + '.csv',
                         header=None)
pip_carma = pd.read_csv('chr' + Chr + '/CARMA_result/' + get_ldname(Chr, Pos) + '_' + trait + '.csv',
                        header=None)
pip_carma_anno = pd.read_csv(
    'chr' + Chr + '/CARMA_result_anno/' + get_ldname(Chr, Pos) + '_' + trait + '.csv', header=None)
pip_paintor = pd.read_csv(
    'chr' + Chr + '/PAINTOR_result/' + get_ldname(Chr, Pos) + '_' + trait + '.results', sep=' ')
pip_paintor_anno = pd.read_csv(
    'chr' + Chr + '/PAINTOR_result_anno/' + get_ldname(Chr, Pos) + '_' + trait + '.results', sep=' ')

pos = pd.merge(chr_data[['pos', 'rsid', 'pval']], pip_susie.iloc[:, 0].rename("rsid"), on='rsid')
pval = pos['pval'].values
rsid = pos['rsid'].values
pos = pos['pos'].values / 1e6

pip_susie = pip_susie.values[:, 1].astype(float)
pip_funmap = pip_funmap.values[:, 1].astype(float)
pip_carma = pip_carma.values[:, 1].astype(float)
pip_carma_anno = pip_carma_anno.values[:, 1].astype(float)
pip_paintor = pip_paintor.values[:, 1].astype(float)
pip_paintor_anno = pip_paintor_anno.values[:, 1].astype(float)

with open('chr' + Chr + '/SuSiE_sets/' + get_ldname(Chr, Pos) + '_' + trait + '.sets', 'r',
          encoding='utf-8') as f:
    file = f.read()

if re.search(r"'cs':\s*(\{[^}]+\})", file) is None:
    sets_susie = None
else:
    sets_susie = eval(re.search(r"'cs':\s*(\{[^}]+\})", file).group(1))

with open('chr' + Chr + '/Funmap_sets/' + get_ldname(Chr, Pos) + '_' + trait + '.sets', 'r',
          encoding='utf-8') as f:
    file = f.read()

if re.search(r"'cs':\s*(\{[^}]+\})", file) is None:
    sets_funmap = None
else:
    sets_funmap = eval(re.search(r"'cs':\s*(\{[^}]+\})", file).group(1))

if os.path.exists('chr' + Chr + '/CARMA_sets/' + get_ldname(Chr, Pos) + '_' + trait + '.sets'):

     sets_CARMA_data = pd.read_csv(
        'chr' + Chr + '/CARMA_sets/' + get_ldname(Chr, Pos) + '_' + trait + '.sets',
        header=None).values[:, -1]
     sets_CARMA = get_carma_sets( sets_CARMA_data)
else:
     sets_CARMA = None

if os.path.exists('chr' + Chr + '/CARMA_anno_sets/' + get_ldname(Chr, Pos) + '_' + trait + '.sets'):

     sets_CARMA_anno_data = pd.read_csv(
        'chr' + Chr + '/CARMA_anno_sets/' + get_ldname(Chr, Pos) + '_' + trait + '.sets',
        header=None).values[:, -1]
     sets_CARMA_anno = get_carma_sets( sets_CARMA_anno_data)
else:
     sets_CARMA_anno = None

if os.path.exists('chr' + Chr + '/PAINTOR_sets/' + get_ldname(Chr, Pos) + '_' + trait + '.sets'):

    sets_PAINTOR_data = np.loadtxt(
        'chr' + Chr + '/PAINTOR_sets/' + get_ldname(Chr, Pos) + '_' + trait + '.sets', dtype=int)
    sets_PAINTOR = get_carma_sets(sets_PAINTOR_data)
else:
    sets_PAINTOR = None

if os.path.exists('chr' + Chr + '/PAINTOR_anno_sets/' + get_ldname(Chr, Pos) + '_' + trait + '.sets'):

    sets_PAINTOR_anno_data = np.loadtxt(
        'chr' + Chr + '/PAINTOR_anno_sets/' + get_ldname(Chr, Pos) + '_' + trait + '.sets', dtype=int)
    sets_PAINTOR_anno = get_carma_sets(sets_PAINTOR_anno_data)
else:
    sets_PAINTOR_anno = None

colors = ['g', 'c', 'y', 'pink', 'm', 'purple', 'orange', 'grey', 'darkblue', 'darkcyan',
          'darkgoldenrod', 'darkgray', 'darkgreen', 'darkkhaki', 'darkmagenta', 'darkolivegreen',
          'darkorange', 'darkorchid', 'darkred', 'darksalmon', 'darkseagreen', 'darkslateblue',
          'darkslategray', 'darkturquoise', 'darkviolet']

"""
For the convenience of drawing, 
the coordinate point of rs2928617 
is temporarily moved out of the figure and drawn separately
"""

pip_susie_5266 = pip_susie[5266]
pip_funmap_5266 = pip_funmap[5266]
pip_carma_5266 = pip_carma[5266]
pip_carma_anno_5266 = pip_carma_anno[5266]
pip_paintor_5266 = pip_paintor[5266]
pip_paintor_anno_5266 = pip_paintor_anno[5266]
pval_5266 = pval[5266]

pip_susie[5266] = 5  
pip_funmap[5266] = 5
pip_carma[5266] = 5  
pip_carma_anno[5266] = 5  
pip_paintor[5266] = 5  
pip_paintor_anno[5266] = 5
pval[5266] = 1e-100  

R = pd.read_csv('chr' + Chr + '/data/' + get_ldname(Chr, Pos) + '.ld', header=None, sep=' ').values
corr = np.abs(R[5266-7:5267+7, 5266-7:5267+7])
bands_y = rsid[5266-7:5267+7]

length = 1000
r1 = 5266 - length//2 + np.where(R[5266, 5266-length//2:5267+length//2]**2>0.8)[0]
r2 = 5266 - length//2 + np.where((R[5266, 5266-length//2:5267+length//2]**2>0.6) & (R[5266, 5266-length//2:5267+length//2]**2<=0.8))[0]
r3 = 5266 - length//2 + np.where((R[5266, 5266-length//2:5267+length//2]**2>0.4) & (R[5266, 5266-length//2:5267+length//2]**2<=0.6))[0]
r4 = 5266 - length//2 + np.where((R[5266, 5266-length//2:5267+length//2]**2>0.2) & (R[5266, 5266-length//2:5267+length//2]**2<=0.4))[0]
# r5 = 5266 - length//2 + np.where(R[5266, 5266-length//2:5267+length//2]**2<=0.2)[0]

fig = plt.figure(figsize=(22, 13))
grid = plt.GridSpec(3, 4, figure=fig, height_ratios=[2, 2, 1], width_ratios=[1, 1, 1, 1])
plt.subplot(grid[0, 0])
ax = sns.heatmap(corr, cmap='Blues', square=True, linewidths=0.5, vmin=0, vmax=1, cbar=True, yticklabels=bands_y)
ax.set_xlabel("Neighbourhood of rs2928617", fontsize=14, labelpad=15)
ax.set_xticklabels(ax.get_xticklabels(), fontsize=6)
ax.xaxis.tick_top()
plt.xticks([])
ax.xaxis.set_label_position('top')
ax.tick_params(top=False, left=False)
ax.add_patch(patches.Rectangle((0, 7), 15, 1, linewidth=2, edgecolor='r', facecolor='none'))
ax.add_patch(patches.Rectangle((7, 0), 1, 15, linewidth=2, edgecolor='r', facecolor='none'))

plt.subplot(grid[1, 0])
plt.scatter(pos, -np.log10(pval), c='darkblue', s=12)
plt.scatter(pos[r1], -np.log10(pval)[r1], c='r', s=12)
plt.scatter(pos[r2], -np.log10(pval)[r2], c='orange', s=12)
plt.scatter(pos[r3], -np.log10(pval)[r3], c='g', s=12)
plt.scatter(pos[r4], -np.log10(pval)[r4], c='c', s=12)
# plt.scatter(pos[r5], -np.log10(pval)[r5], c='b', s=12)
plt.scatter(pos[5266], -np.log10(pval_5266), c='purple', s=24, marker='D')
plt.axhline(y=5, color='r', linestyle='--', linewidth=1)
plt.text(pos[5266]-0.35, -np.log10(pval_5266), 'rs2928617', color='k', fontsize=14)
# plt.xlabel('Position on Chromosome ' + Chr + ' (Mb)', fontsize=14)
plt.ylabel(r'$-\log_{10}(p)$', fontsize=14)
plt.xlim(6, 7)
plt.ylim((-0.05, 9))

axins = plt.gca().inset_axes((0.05, 0.58, 0.40, 0.25))
axins.scatter(pos, -np.log10(pval), c='darkblue', s=12)
axins.scatter(pos[r1], -np.log10(pval)[r1], c='r', s=12)
axins.scatter(pos[r2], -np.log10(pval)[r2], c='orange', s=12)
axins.scatter(pos[5266]+0.002, -np.log10(pval_5266), c='purple', s=24, marker='D')
axins.set_xlim(6.555, 6.572)
axins.set_ylim(7.5, 8.7)
axins.set_xticks([])
axins.set_yticks([])
mark_inset(plt.gca(), axins, loc1=2, loc2=4, fc="none", ec='k', lw=0.2)

plt.title('GWAS', fontsize=14, fontweight="bold")
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.xticks([])
plt.yticks(fontsize=12)
legend_elements = [Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.8, 1]$', markerfacecolor='r', markersize=12),
                   Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.6, 0.8]$', markerfacecolor='orange', markersize=12),
                   Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.4, 0.6]$', markerfacecolor='g', markersize=12),
                   Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.2, 0.4]$', markerfacecolor='c', markersize=12)]
plt.legend(handles=legend_elements, loc='upper right')

plt.subplot(grid[0, 1])
plt.scatter(pos, pip_funmap, c='darkblue', s=24)
plt.scatter(pos[5266], pip_funmap_5266, c='purple', s=64, edgecolors='g', linewidths=2, marker='D')
plt.text(pos[5266]-0.3, pip_funmap_5266, 'rs2928617', color='k', fontsize=14)
plt.ylabel('PIP', fontsize=14)
plt.title('Funmap', fontsize=14, fontweight="bold")
plt.xlim(6, 7)
plt.ylim([-0.05, 1.05])
plt.gca().axes.xaxis.set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.yticks(fontsize=12)
legend_elements = [Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.8, 1]$', markerfacecolor='r', markersize=12),
                   Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.6, 0.8]$', markerfacecolor='orange', markersize=12),
                   Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.4, 0.6]$', markerfacecolor='g', markersize=12),
                   Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.2, 0.4]$', markerfacecolor='c', markersize=12)]
plt.legend(handles=legend_elements, loc='upper right')

plt.subplot(grid[1, 1])
nonsets_index = np.concatenate(list(sets_susie.values()))
nonsets_index = np.delete(np.arange(8646), nonsets_index)
plt.scatter(pos[nonsets_index], pip_susie[nonsets_index], c='darkblue', s=24)
for i, group in enumerate(sets_susie.values()):
    color = colors[i]
    for idx in group:
        if R[5266, idx]**2>0.8:
            plt.scatter(pos[idx], pip_susie[idx], c='r', s=48, edgecolors=color, linewidths=2)
        elif (R[5266, idx]**2>0.6) and (R[5266, idx]**2<=0.8):
            plt.scatter(pos[idx], pip_susie[idx], c='orange', s=48, edgecolors=color, linewidths=2)
        elif (R[5266, idx]**2>0.4) and (R[5266, idx]**2<=0.6):
            plt.scatter(pos[idx], pip_susie[idx], c='g', s=48, edgecolors=color, linewidths=2)
        elif (R[5266, idx]**2>0.2) and (R[5266, idx]**2<=0.4):
            plt.scatter(pos[idx], pip_susie[idx], c='c', s=48, edgecolors=color, linewidths=2)
        elif R[5266, idx]**2<=0.2:
            plt.scatter(pos[idx], pip_susie[idx], c='darkblue', s=48, edgecolors=color, linewidths=2)
plt.scatter(pos[5266], pip_susie_5266, c='purple', s=64, edgecolors='g', linewidths=2, marker='D')
plt.text(pos[5266]-0.3, pip_susie_5266, 'rs2928617', color='k', fontsize=14)
plt.ylabel('PIP', fontsize=14)
plt.title('SuSiE', fontsize=14, fontweight="bold")
plt.xlim(6, 7)
plt.ylim([-0.05, 1.05])
# plt.xlabel('Position on Chromosome ' + Chr + ' (Mb)', fontsize=14)
plt.xticks([])
# plt.gca().axes.xaxis.set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.yticks(fontsize=12)
legend_elements = [Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.8, 1]$', markerfacecolor='r', markersize=12),
                   Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.6, 0.8]$', markerfacecolor='orange', markersize=12),
                   Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.4, 0.6]$', markerfacecolor='g', markersize=12),
                   Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.2, 0.4]$', markerfacecolor='c', markersize=12)]
plt.legend(handles=legend_elements, loc='upper right')

plt.subplot(grid[0, 2])
nonsets_index = np.concatenate(list(sets_CARMA_anno.values()))
nonsets_index = np.delete(np.arange(8646), nonsets_index)
plt.scatter(pos[nonsets_index], pip_carma_anno[nonsets_index], c='darkblue', s=24)
for i, group in enumerate(sets_CARMA_anno.values()):
    color = colors[i]
    for idx in group:
        if R[5266, idx]**2>0.8:
            plt.scatter(pos[idx], pip_carma_anno[idx], c='r', s=48, edgecolors=color, linewidths=2)
        elif (R[5266, idx]**2>0.6) and (R[5266, idx]**2<=0.8):
            plt.scatter(pos[idx], pip_carma_anno[idx], c='orange', s=48, edgecolors=color, linewidths=2)
        elif (R[5266, idx]**2>0.4) and (R[5266, idx]**2<=0.6):
            plt.scatter(pos[idx], pip_carma_anno[idx], c='g', s=48, edgecolors=color, linewidths=2)
        elif (R[5266, idx]**2>0.2) and (R[5266, idx]**2<=0.4):
            plt.scatter(pos[idx], pip_carma_anno[idx], c='c', s=48, edgecolors=color, linewidths=2)
        elif R[5266, idx]**2<=0.2:
            plt.scatter(pos[idx], pip_carma_anno[idx], c='darkblue', s=48, edgecolors=color, linewidths=2)
plt.scatter(pos[5266], pip_carma_5266, c='purple', s=64, edgecolors='c', linewidths=2, marker='D')
plt.text(pos[5266]+0.1, pip_carma_5266, 'rs2928617', color='k', fontsize=14)
plt.title('CARMA with annotation', fontsize=14, fontweight="bold")
plt.ylabel('PIP', fontsize=14)
plt.xlim(6, 7)
plt.ylim([-0.05, 1.05])
plt.gca().axes.xaxis.set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.yticks(fontsize=12)
legend_elements = [Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.8, 1]$', markerfacecolor='r', markersize=12),
                   Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.6, 0.8]$', markerfacecolor='orange', markersize=12),
                   Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.4, 0.6]$', markerfacecolor='g', markersize=12),
                   Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.2, 0.4]$', markerfacecolor='c', markersize=12)]
plt.legend(handles=legend_elements, loc='upper right')

plt.subplot(grid[1, 2])
nonsets_index = np.concatenate(list(sets_CARMA.values()))
nonsets_index = np.delete(np.arange(8646), nonsets_index)
plt.scatter(pos[nonsets_index], pip_carma[nonsets_index], c='darkblue', s=24)
for i, group in enumerate(sets_CARMA.values()):
    color = colors[i]
    for idx in group:
        if R[5266, idx]**2>0.8:
            plt.scatter(pos[idx], pip_carma[idx], c='r', s=48, edgecolors=color, linewidths=2)
        elif (R[5266, idx]**2>0.6) and (R[5266, idx]**2<=0.8):
            plt.scatter(pos[idx], pip_carma[idx], c='orange', s=48, edgecolors=color, linewidths=2)
        elif (R[5266, idx]**2>0.4) and (R[5266, idx]**2<=0.6):
            plt.scatter(pos[idx], pip_carma[idx], c='g', s=48, edgecolors=color, linewidths=2)
        elif (R[5266, idx]**2>0.2) and (R[5266, idx]**2<=0.4):
            plt.scatter(pos[idx], pip_carma[idx], c='c', s=48, edgecolors=color, linewidths=2)
        elif R[5266, idx]**2<=0.2:
            plt.scatter(pos[idx], pip_carma[idx], c='darkblue', s=48, edgecolors=color, linewidths=2)
plt.scatter(pos[5266], pip_carma_anno_5266, c='purple', s=64, edgecolors='g', linewidths=2, marker='D')
plt.text(pos[5266]+0.1, pip_carma_anno_5266, 'rs2928617', color='k', fontsize=14)
plt.title('CARMA', fontsize=14, fontweight="bold")
plt.ylabel('PIP', fontsize=14)
plt.xlim(6, 7)
plt.ylim([-0.05, 1.05])
# plt.xlabel('Position on Chromosome ' + Chr + ' (Mb)', fontsize=14)
plt.xticks([])
# plt.gca().axes.xaxis.set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.yticks(fontsize=12)
legend_elements = [Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.8, 1]$', markerfacecolor='r', markersize=12),
                   Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.6, 0.8]$', markerfacecolor='orange', markersize=12),
                   Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.4, 0.6]$', markerfacecolor='g', markersize=12),
                   Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.2, 0.4]$', markerfacecolor='c', markersize=12)]
plt.legend(handles=legend_elements, loc='upper right')

plt.subplot(grid[0, 3])
nonsets_index = np.concatenate(list(sets_PAINTOR_anno.values()))
nonsets_index = np.delete(np.arange(8646), nonsets_index)
plt.scatter(pos[nonsets_index], pip_paintor_anno[nonsets_index], c='darkblue', s=24)
for i, group in enumerate(sets_PAINTOR_anno.values()):
    color = colors[i]
    for idx in group:
        if R[5266, idx]**2>0.8:
            plt.scatter(pos[idx], pip_paintor_anno[idx], c='r', s=48, edgecolors=color, linewidths=2)
        elif (R[5266, idx]**2>0.6) and (R[5266, idx]**2<=0.8):
            plt.scatter(pos[idx], pip_paintor_anno[idx], c='orange', s=48, edgecolors=color, linewidths=2)
        elif (R[5266, idx]**2>0.4) and (R[5266, idx]**2<=0.6):
            plt.scatter(pos[idx], pip_paintor_anno[idx], c='g', s=48, edgecolors=color, linewidths=2)
        elif (R[5266, idx]**2>0.2) and (R[5266, idx]**2<=0.4):
            plt.scatter(pos[idx], pip_paintor_anno[idx], c='c', s=48, edgecolors=color, linewidths=2)
        elif R[5266, idx]**2<=0.2:
            plt.scatter(pos[idx], pip_paintor_anno[idx], c='darkblue', s=48, edgecolors=color, linewidths=2)
# plt.scatter(pos[5266], pip_paintor_anno[5266], c='r', s=24, edgecolors='g', linewidths=2)
plt.scatter(pos[5266], pip_paintor_5266, c='purple', s=64, marker='D')
plt.text(pos[5266]-0.3, pip_paintor_5266+0.02, 'rs2928617', color='k', fontsize=14)
plt.title('PAINTOR with annotation', fontsize=14, fontweight="bold")
plt.ylabel('PIP', fontsize=14)
plt.xlim(6, 7)
plt.ylim([-0.05, 1.05])
plt.gca().axes.xaxis.set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.yticks(fontsize=12)
legend_elements = [Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.8, 1]$', markerfacecolor='r', markersize=12),
                   Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.6, 0.8]$', markerfacecolor='orange', markersize=12),
                   Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.4, 0.6]$', markerfacecolor='g', markersize=12),
                   Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.2, 0.4]$', markerfacecolor='c', markersize=12)]
plt.legend(handles=legend_elements, loc='upper right')

plt.subplot(grid[1, 3])
nonsets_index = np.concatenate(list(sets_PAINTOR.values()))
nonsets_index = np.delete(np.arange(8646), nonsets_index)
plt.scatter(pos[nonsets_index], pip_paintor[nonsets_index], c='darkblue', s=24)
for i, group in enumerate(sets_PAINTOR.values()):
    color = colors[i]
    for idx in group:
        if R[5266, idx]**2>0.8:
            plt.scatter(pos[idx], pip_paintor[idx], c='r', s=48, edgecolors=color, linewidths=2)
        elif (R[5266, idx]**2>0.6) and (R[5266, idx]**2<=0.8):
            plt.scatter(pos[idx], pip_paintor[idx], c='orange', s=48, edgecolors=color, linewidths=2)
        elif (R[5266, idx]**2>0.4) and (R[5266, idx]**2<=0.6):
            plt.scatter(pos[idx], pip_paintor[idx], c='g', s=48, edgecolors=color, linewidths=2)
        elif (R[5266, idx]**2>0.2) and (R[5266, idx]**2<=0.4):
            plt.scatter(pos[idx], pip_paintor[idx], c='c', s=48, edgecolors=color, linewidths=2)
        elif R[5266, idx]**2<=0.2:
            plt.scatter(pos[idx], pip_paintor[idx], c='darkblue', s=48, edgecolors=color, linewidths=2)
# plt.scatter(pos[5266], pip_paintor_anno[5266], c='r', s=24, edgecolors='g', linewidths=2)
plt.scatter(pos[5266], pip_paintor_anno_5266, c='purple', s=64, marker='D')
plt.text(pos[5266]-0.3, pip_paintor_anno_5266+0.05, 'rs2928617', color='k', fontsize=14)
plt.title('PAINTOR', fontsize=14, fontweight="bold")
plt.ylabel('PIP', fontsize=14)
plt.xlim(6, 7)
plt.ylim([-0.05, 1.05])
# plt.xlabel('Position on Chromosome ' + Chr + ' (Mb)', fontsize=14)
plt.xticks([])
# plt.gca().axes.xaxis.set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.yticks(fontsize=12)
legend_elements = [Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.8, 1]$', markerfacecolor='r', markersize=12),
                   Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.6, 0.8]$', markerfacecolor='orange', markersize=12),
                   Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.4, 0.6]$', markerfacecolor='g', markersize=12),
                   Line2D([0], [0], marker='o', color='w', label='$r^2\in (0.2, 0.4]$', markerfacecolor='c', markersize=12)]
plt.legend(handles=legend_elements, loc='upper right')

plt.subplot(grid[2, 0])
for row, y in zip(gene, y_list):
    plt.hlines(y, xmin=row[1] / 1e6, xmax=row[2] / 1e6, color='grey', lw=8)
    if row[0] == 'AGPAT5':
        plt.text(row[1] / 1e6, y, row[0], fontsize=8,
                 verticalalignment="center", horizontalalignment="right", color='r')
    else:
        plt.text(row[1] / 1e6, y, row[0], fontsize=8, verticalalignment="center", horizontalalignment="right")
plt.xlim(6, 7)
plt.ylim(0, 5)
plt.yticks(fontsize=14)
plt.yticks([])
plt.xlabel('Position on Chromosome 8 (Mb)', fontsize=14)
plt.ylabel('Gene', fontsize=14)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

plt.subplot(grid[2, 1])
for row, y in zip(gene, y_list):
    plt.hlines(y, xmin=row[1] / 1e6, xmax=row[2] / 1e6, color='grey', lw=8)
    if row[0] == 'AGPAT5':
        plt.text(row[1] / 1e6, y, row[0], fontsize=8,
                 verticalalignment="center", horizontalalignment="right", color='r')
    else:
        plt.text(row[1] / 1e6, y, row[0], fontsize=8, verticalalignment="center", horizontalalignment="right")
plt.xlim(6, 7)
plt.ylim(0, 5)
plt.yticks(fontsize=14)
plt.yticks([])
plt.xlabel('Position on Chromosome 8 (Mb)', fontsize=14)
# plt.ylabel('Gene', fontsize=14)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

plt.subplot(grid[2, 2])
for row, y in zip(gene, y_list):
    plt.hlines(y, xmin=row[1] / 1e6, xmax=row[2] / 1e6, color='grey', lw=8)
    if row[0] == 'AGPAT5':
        plt.text(row[1] / 1e6, y, row[0], fontsize=8,
                 verticalalignment="center", horizontalalignment="right", color='r')
    else:
        plt.text(row[1] / 1e6, y, row[0], fontsize=8, verticalalignment="center", horizontalalignment="right")
plt.xlim(6, 7)
plt.ylim(0, 5)
plt.yticks(fontsize=14)
plt.yticks([])
plt.xlabel('Position on Chromosome 8 (Mb)', fontsize=14)
# plt.ylabel('Gene', fontsize=14)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

plt.subplot(grid[2, 3])
for row, y in zip(gene, y_list):
    plt.hlines(y, xmin=row[1] / 1e6, xmax=row[2] / 1e6, color='grey', lw=8)
    if row[0] == 'AGPAT5':
        plt.text(row[1] / 1e6, y, row[0], fontsize=8,
                 verticalalignment="center", horizontalalignment="right", color='r')
    else:
        plt.text(row[1] / 1e6, y, row[0], fontsize=8, verticalalignment="center", horizontalalignment="right")
plt.xlim(6, 7)
plt.ylim(0, 5)
plt.yticks(fontsize=14)
plt.yticks([])
plt.xlabel('Position on Chromosome 8 (Mb)', fontsize=14)
# plt.ylabel('Gene', fontsize=14)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

plt.subplots_adjust(left=0.05, right=0.98, bottom=0.08, top=0.95)
plt.savefig('man_chr8_5.png', dpi=900)
# plt.show()
