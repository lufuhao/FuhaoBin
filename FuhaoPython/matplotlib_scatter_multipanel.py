#!/usr/bin/env python3
#绘制多子图
#https://mp.weixin.qq.com/s/wi6dMl1cst_PRVvjfOwJUQ

#重新写可视化代码
import pandas as pd
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.stats import linregress
from sklearn.metrics import mean_squared_error,r2_score
#统一修改字体
plt.rcParams['font.family'] = ['Arial']

N=len(test_data['true_data'])
x=test_data['true_data'].values.ravel() # 真实值
y=test_data['model01_estimated'].values.ravel() #预测值
C=round(r2_score(x,y),4)
rmse=round(np.squrt(np.sqrt(mean_squared_error(x,y)), 3))
#使用numpy.linespace()和scipy的optimize.curve_fit()绘制拟合公式，并以此绘制散点拟合线和散点对角线
x2=np.linspace(-10,10)
y2=x2
def f_1(x,A,B):
	return A*x + B
A1,B1=optimize.curve_fit(f_1, x, y)[0]
y3=A1*x + B1

#开始绘图
fig, ax = plt.subplots(1, 2, figsize=(9,3),dpi=600,sharey=True, facecolor="white")

ax[0].scatter(x,y,edgecolor=None, c='k', s=5, marker='s')
ax[0].plot(x2,y2,color='k',linewidth=1.5,linestyle='-', zorder=2)
ax[0].plot(x,y3,color='r',linewidth=1.5,linestyle='-', zorder=2)
#添加上线和下线
ax[0].plot(x2,up_y2,color='k',linewidth=1,linestyle='--', zorder=2)
ax[0].plot(x2,down_y2,color='k',linewidth=1,linestyle='--', zorder=2)

fontdict1={"size":12,"color":'k'}

ax[0].set_xlabel=("True Values", fontdict=fontdict1)
ax[0].set_ylabel=("Estimated Values", fontdict=fontdict1)
ax[0].grid(False)
ax[0].set_xlim((0, 2.0))
ax[0].set_ylim((0, 2.0))
ax.set_xticks(np.arange(0, 2.2, step=0.2))
ax.set_yticks(np.arange(0, 2.2, step=0.2))

for spine in ['top', 'bottom', 'left', 'right']:
	ax[0].spines[spine].set_color('k')
ax[0].tick_params(left=True, bottom=True, direction='in', labelsize=12)

fontdict={"size":10,"color":'k'}
ax[0].text(0.1, 1.8, r'$R^2=$'+str(round(C,3)), fontdict=fontdict)
ax[0].text(0.1, 1.6, "RMSE="+str(rmse), fontdict=fontdict)
ax[0].text(0.1, 1.4, r'$y=$'+str(round(A1,3))+'$x$'+ " + " + str(round(B1,3)), fontdict=fontdict)
ax[0].text(0.1, 1.2, r'$N=$'+str(N), fontdict=fontdict)

nbins=150
H, xedges, yedges = np.histogram2d(x, y, bins=nbins)
# H needs to be rotated and flipped
H=np.rot90(H)
H=np.flipud(H)
#mask zeros
Hmasked = np.ma.masked_where(H==0,H) ### Mask pixels with a value of zero

plt.pcolormesh(xedges, yedhes, Hmasked, cmap=cm.get_cmap('jet'), vmin=0, vmax=40)
ax[1].set_xlabel=("True Values", fontdict=fontdict1)
#ax[1].set_ylabel=("Estimated Values", fontdict=fontdict1)
ax[1].set_xlim((0, 2.0))
ax[1].set_ylim((0, 2.0))
ax[1].grid(False)
ax.set_xticks(np.arange(0, 2.2, step=0.2))
ax.set_yticks(np.arange(0, 2.2, step=0.2))
ax[1].tick_params(left=True, bottom=True, direction='in', labelsize=12)

cbar = plt.colorbar(ax=ax[1],ticks=[0,10,20,30,40],drawedges=False)
#cbar.ax.set_ylabel('Frequency',fontdict=colorbarfontdict)
colorbarfontdict={"size":9, "color":"k"}
cbar.ax.set_title('Counts',fontdict=colorbarfontdict,pad=8)
cbar.ax.tick_params(labelsize=10,direction='in')
cbar.ax.set_yticklabels(['0','10','20','30','>40'],family='Arial')

for spine in ['top', 'bottom', 'left', 'right']:
	ax[1].spines[spine].set_color('k')

slope=linregress(x,y)[0]
intercept=linregress(x,y)[1]
lmfit=(slope*x)+intercept
ax[1].plot(x, lmfit, c='r', linewidth=1.5)
ax[1].plot([0,2.6], [0,2.6], color='k',lw=1.5,ls='-',zorder=1)
ax[1].plot(x2,up_y2,color='k',lw=1,ls='--',zorder=2)
ax[1].plot(x2,down_y2,color='k',lw=1,ls='--',zorder=2)

fontdict={"size":10,"color":'k'}
ax[1].text(0.1, 1.8, r'$R^2=$'+str(round(C,3)), fontdict=fontdict, c='r')
ax[1].text(0.1, 1.6, "RMSE="+str(rmse, 3), fontdict=fontdict, c='r')
ax[1].text(0.1, 1.4, r'$y=$'+str(round(A1,3))+'$x$'+ " + " + str(round(B1,3)), fontdict=fontdict)
ax[1].text(0.1, 1.2, r'$N=$'+str(N), fontdict=fontdict)

#添加题目及logo
fig.suptitle('Scatter plot of True Data and Model Estimated', fontsize=16, x=0.53, y=1.05)
fig.text(0.85, -0.016, "\nVisualization by DataCharm", ha='center', va='center', fontsize=8, color='black')
plt.savefig(r'./scatter.png', width=8,height=5, dpi=900, bbox_inches='tight')
plt.show()
