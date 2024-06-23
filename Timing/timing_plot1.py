import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


data = pd.read_csv('timing1.csv')
y1 = data.mean().values[[0, 6, 12, 18, 24]]
y2 = data.mean().values[[1, 7, 13, 19, 25]]
y3 = data.mean().values[[2, 8, 14, 20, 26]]
y4 = data.mean().values[[3, 9, 15, 21, 27]]
y5 = data.mean().values[[4, 10, 16, 22, 28]]
y6 = data.mean().values[[5, 11, 17, 23, 29]]

std1 = data.std().values[[0, 6, 12, 18, 24]]
std2 = data.std().values[[1, 7, 13, 19, 25]]
std3 = data.std().values[[2, 8, 14, 20, 26]]
std4 = data.std().values[[3, 9, 15, 21, 27]]
std5 = data.std().values[[4, 10, 16, 22, 28]]
std6 = data.std().values[[5, 11, 17, 23, 29]]

x = np.array([854, 1313, 1833, 2561, 4107])
# 创建图形
fig, ax = plt.subplots(figsize=(5, 4.5))

# 绘制折线图
ax.plot(x, y4, marker='o', color='#9B2E2B', label='Funmap')
ax.plot(x, y5, marker='o', color='#E2533D', label='CARMA+anno')
ax.plot(x, y6, marker='o', color='#F9E4A9', label='PAINTOR+anno')
ax.plot(x, y1, marker='o', color='#315CB5', label='SuSiE')
ax.plot(x, y2, marker='o', color='#459CD7', label='CARMA')
ax.plot(x, y3, marker='o', color='#B2DFE3', label='PAINTOR')
ax.plot(x, y4, marker='o', color='#9B2E2B')
# 添加误差棒
ax.errorbar(x, y4, yerr=0.5*std4, fmt='o', ecolor='#9B2E2B', capsize=3, color='#9B2E2B')
ax.errorbar(x, y5, yerr=0.5*std5, fmt='o', ecolor='#E2533D', capsize=3, color='#E2533D')
ax.errorbar(x, y6, yerr=0.5*std6, fmt='o', ecolor='#F9E4A9', capsize=3, color='#F9E4A9')
ax.errorbar(x, y1, yerr=0.5*std1, fmt='o', ecolor='#315CB5', capsize=3, color='#315CB5')
ax.errorbar(x, y2, yerr=0.5*std2, fmt='o', ecolor='#459CD7', capsize=3, color='#459CD7')
ax.errorbar(x, y3, yerr=0.5*std3, fmt='o', ecolor='#B2DFE3', capsize=3, color='#B2DFE3')
ax.errorbar(x, y4, yerr=0.5*std4, fmt='o', ecolor='#9B2E2B', capsize=3, color='#9B2E2B')

# 设置标题和坐标轴标签
ax.set_xlabel('Number of candidate variants (p)', fontsize=14)
ax.set_ylabel('Computational Time (Seconds)', fontsize=14)

# 显示图形
plt.legend(fontsize=12)
plt.show()
