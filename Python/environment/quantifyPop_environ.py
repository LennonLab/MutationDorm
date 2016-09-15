from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
#import pylab
#import matplotlib as mpl
import os
import pandas as pd

mydir = os.path.expanduser("~/github/MutationDorm/")


content_list = []

for content in os.listdir(mydir + 'data/'):
    if content.endswith('.txt'):
        content_list.append(content)

fig = plt.figure()
ax1 = fig.add_subplot(111)

colors = ['green', 'orange', 'red']
seed_bank_size = ['Green = 1000 dormant cells', 'Orange = 100 dormant cells', 'Red = 10 dormant cells']
count = 0
for file_name in content_list:
    M_freq = file_name[:-4].split('_')
    if M_freq[1] == '0.5':
        print M_freq, colors[count]
        file_pd = pd.read_csv(mydir + 'data/'+ file_name, sep = '\t', index_col =0)
        file_pd_melt = pd.melt(file_pd)
        file_pd_melt.variable = file_pd_melt.variable.astype(int)
        ax1.scatter(file_pd_melt.ix[:,0].values, np.log10(file_pd_melt.ix[:,1].values),\
            color=colors[count],s=5, alpha=0.05, edgecolor='none', label= seed_bank_size[count])
        count += 1

plt.xlim([0,100])
plt.xlabel('Generations', fontsize=20)
plt.ylabel('Fitness, log10', fontsize=20)
plt.legend(loc='upper right')
fig.savefig(mydir + 'figs/test.png',  bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
