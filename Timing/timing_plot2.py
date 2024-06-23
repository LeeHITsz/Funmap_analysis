import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


data = pd.read_csv('timing2.csv')
y4 = data.mean().values[[0, 3, 6, 9, 12]]
y5 = data.mean().values[[1, 4, 7, 10, 13]]
y6 = data.mean().values[[2, 5, 8, 11, 14]]

std4 = data.std().values[[0, 3, 6, 9, 12]]
std5 = data.std().values[[1, 4, 7, 10, 13]]
std6 = data.std().values[[2, 5, 8, 11, 14]]

x = np.array([10, 50, 100, 150, 200])
# 创建图形
fig, ax = plt.subplots(figsize=(5, 4.5))

# 绘制折线图
ax.plot(x, y4, marker='o', color='#9B2E2B', label='Funmap')
ax.plot(x, y5, marker='o', color='#E2533D', label='CARMA+anno')
ax.plot(x, y6, marker='o', color='#F9E4A9', label='PAINTOR+anno')

# 添加误差棒
ax.errorbar(x, y4, yerr=0.5*std4, fmt='o', ecolor='#9B2E2B', capsize=3, color='#9B2E2B')
ax.errorbar(x, y5, yerr=0.5*std5, fmt='o', ecolor='#E2533D', capsize=3, color='#E2533D')
ax.errorbar(x, y6, yerr=0.5*std6, fmt='o', ecolor='#F9E4A9', capsize=3, color='#F9E4A9')


# 设置标题和坐标轴标签
ax.set_xlabel('Number of functional annotations (m)', fontsize=14)
ax.set_ylabel('Computational Time (Seconds)', fontsize=14)
ax.set_xticks([0, 50, 100, 150, 200])

# 显示图形
plt.legend(fontsize=12)
plt.show()
