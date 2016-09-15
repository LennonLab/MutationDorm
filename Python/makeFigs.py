from __future__ import division
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os


mydir = os.path.expanduser("~/github/MutationDorm/")


def M_plot():
    Mms = [10, 100, 1000, 10000]
    pi = []
    x = []
    for M in Mms:
        IN = pd.read_csv(mydir + 'data/G10000_S100_N1000_M' + str(M) + '_c10.txt', \
            sep = ' ')
        L = len(IN.WT.values)
        pi.append(IN.WT.values)
        x.append( np.asarray([M] * L))
    fig, ax = plt.subplots()
    ax.set_xscale('log', basex=10)
    plt.scatter(x, pi)

    plt.savefig(mydir + 'figs/test.png')


def piMergeDF(M):
    df = pd.DataFrame()
    for i in range(0, 100):
        IN = pd.read_csv(mydir + 'data/sim_history/G10000_S100_N1000_M' + str(M) + '_c10_I' + str(i) + '.txt', \
            sep = ' ')
        df[i] = IN.PI
    df.to_csv(path_or_buf = mydir + 'data/sim_history/merged/G10000_S100_N1000_M' + str(M) + '_c10_Pi.txt', \
        sep = ' ', index_label = None, index = False, columns = None)

def piTimeSeries(M):
    IN = pd.read_csv(mydir + 'data/sim_history/merged/G10000_S100_N1000_M' + str(M) + '_c10_Pi.txt', \
        sep = ' ')
    x = np.arange(1,10001)
    fig, ax = plt.subplots()
    IN_mean = IN.mean(axis = 1)
    plt.plot(x, IN_mean, linewidth = 2, linestyle = '-', alpha=1.0, color = 'k')
    for column in IN:
        y = IN[column].values
        plt.plot(x, y, 'k-', alpha=0.1, color = 'b')
    plt.ylabel(r'$\pi$', fontsize=30, rotation=0, labelpad= 20)
    plt.xlabel('Generations', fontsize=20)
    plt.text(0.1, 0.9,'M = '+ str(M), \
     horizontalalignment='center', \
     verticalalignment='center', \
     transform = ax.transAxes, fontsize = 16)
    fig.savefig(mydir + '/figs/Simulation_M' + str(M)  + '.png', \
        bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def piTimeSeriesAll():
    Ms = [10, 100, 1000]
    Ex = [2.02, 2.2, 4]
    #Ms = [10]
    colors = ['#87CEEB', '#FFA500', '#FF6347']
    fig, ax = plt.subplots()
    x = np.arange(1,10001)
    for i, M in enumerate(Ms):
        IN = pd.read_csv(mydir + 'data/sim_history/merged/G10000_S100_N1000_M' + str(M) + '_c10_Pi.txt', \
            sep = ' ')
        IN_mean = IN.mean(axis = 1)
        IN_std = IN.std(axis = 1)
        ax.plot(x, IN_mean, lw=2, label='M = ' + str(M), color=colors[i], alpha = 0.9)
        ax.fill_between(x, IN_mean+IN_std, IN_mean-IN_std, facecolor=colors[i], alpha=0.5)
        ax.axhline(y = Ex[i], color=colors[i], ls = '--')
    ax.text(10, 1.8, r'$\mathbf{E}_{10}[\pi]$')
    ax.text(10, 2.35, r'$\mathbf{E}_{100}[\pi]$')
    ax.text(10 , 4.15, r'$\mathbf{E}_{1,000}[\pi]$')
    ax.legend(loc='upper left')
    ax.set_xlabel('Generations', fontsize=20)
    ax.set_ylabel(r'$\pi$', fontsize=30, rotation=0, labelpad= 20)
    ax.grid()
    ax.set_ylim([0,6])
    fig.savefig(mydir + '/figs/Simulation_All_Ms.png', \
        bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


#piMergeDF(10)
#piMergeDF(100)
#piMergeDF(1000)
#piTimeSeries(10)
#piTimeSeries(100)
#piTimeSeries(1000)

piTimeSeriesAll()
