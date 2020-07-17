#!/usr/bin/env python3
#https://mp.weixin.qq.com/s/wi6dMl1cst_PRVvjfOwJUQ

#重新写可视化代码
import pandas as pd
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
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
fig, ax = plt.subplots(figsize=(7,5),dpi=600)
dian=plt.scatter(x,y,edgecolor=None, c='k',s=16, marker='s')
ax.plot(x2,y2,color='k',linewidth=1.5,linestyle='--')
ax.plot(x,y3,color='r',linewidth=2,linestyle='--')
fontdict1={"size":17,"color":'k',"family":'Arial'}
ax.set_xlabel=("True Values", fontdict=fontdict1)
ax.set_ylabel=("Estimated Values", fontdict=fontdict1)
ax.grid(False)
ax.set_xlim((0, 2.0))
ax.set_ylim((0, 2.0))
ax.set_xticks(np.arange(0, 2.2, step=0.2))
ax.set_yticks(np.arange(0, 2.2, step=0.2))
#设置刻度字体
labels=ax.get_xticklabls() + ax.get_yticklabels()
[label.set_fontname('Arial') for label in labels]

for spine in ['top', 'bottom', 'left', 'right']:
	ax.spines[spine].set_color('k')
ax.tick_params(left=True,bottom=True,direction='in',labelsize=14)
#添加题目
titlefontdict={"size":20,"color":'k',"family":'Arial'}
ax.set_title('Scatter plot of Ture Data and Model Estimated', titlefontdict, pad=20)
#ax.set_title()
fontdict={"size":16,"color":'k',"family":'Arial'}
ax.text(0.1, 1.8, r'$R^2=$'+str(round(C,3)), fontdict=fontdict)
ax.text(0.1, 1.6, "RMSE="+str(rmse), fontdict=fontdict)
ax.text(0.1, 1.4, r'$y=$'+str(round(A1,3))+'$x$'+ " + " + str(round(B1,3)), fontdict=fontdict)
ax.text(0.1, 1.2, r'$N=$'+str(N), fontdict=fontdict)

text_font={'family':'Arial', 'size':22, 'weight':'bold', 'color':'black'}
ax.text(.9, .9, "(a)", transform=ax.transAxes, fontdict=textfont, zorder=4)
ax.text(.8, .056, "\nVisualization by DataCharm", transform=ax.transAxes, ha='center', va='center', fontsize=10, color='black')

#计算密度范围并附上颜色
# Estimate the 2D histogram
nbins = 150
H, xedges, yedges = np.histogram2d(X, Y, bins=nbins)
# H needs to be rotated and flipped
H = np.rot90(H)
H = np.flipud(H)
# Mask zeros
Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
#开始绘图
plt.pcolormesh(xedges, yedges, Hmasked, cmap=cm.get_cmap('jet'), vmin=0, vmax=40)

#对colorbar进行定制化设置
colorbarfontdict={"size":9, "color":"k"}
cbar = plt.colorbar(ax=ax,ticks=[0,10,20,30,40],drawedges=False)
#cbar.ax.set_ylabel('Frequency',fontdict=colorbarfontdict)
cbar.ax.set_title('Counts',fontdict=colorbarfontdict,pad=8)
cbar.ax.tick_params(labelsize=12,direction='in')
cbar.ax.set_yticklabels(['0','10','20','30','>40'],family='Times New Roman')

#学术性相关性散点图还需添加拟合最佳上线(upper line)和下线(bottom line)，而两者的绘制依据为1：1 最佳线和误差 Δτ= ± (0.05+0.15 True data ,分别对应y=1.15×+0.05 (upper line) and y=0.85×−0.05 (bottom line)。
#用于绘制最佳拟合线
x2 = np.linspace(-10,10)
y2=x2
#绘制upper line
up_y2 = 1.15*x2 + 0.05
#绘制bottom line
down_y2 = 0.85*x2 - 0.05
#添加上线和下线
ax.plot(x2,up_y2,color='k',lw=1.5,ls='--',zorder=2)
ax.plot(x2,down_y2,color='k',lw=1.5,ls='--',zorder=2)


#上述结果是更改了matplotlib绘图风格，即在绘图之前添加如下代码：
plt.style.use('seaborn-darkgrid')

plt.savefig(r'./scatter.png', width=7,height=7, dpi=900, bbox_inches='tight')
plt.show()
