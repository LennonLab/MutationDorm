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


def plot_fitness(gen = 60000):
    IN = pd.read_csv(mydir + 'data/fitness/T' + str(gen) + '.txt', sep = '\t')
    values = []
    print np.log10(np.mean(IN.fitness.values) / gen)
    in_zip = zip(IN.fitness.values.tolist(), IN.n_i.values.tolist())
    for i in in_zip:
        values.extend([round(i[0], 6)] * i[1])
    fig = plt.figure()
    plt.hist(values, bins=30, alpha = 0.8,  normed = False)
    plt.title('dist. of fitness')
    plt.xlabel('Fitness', fontsize=14)
    plt.ylabel('Probability', fontsize=14)
    #plt.xlim([0, max(120,  x_mean  + (x_mean * 2)) ])
    fig.tight_layout()
    out_plot = mydir + 'figs/T' + str(gen) + '.png'
    fig.savefig(out_plot, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def clean_v_dorm():
    OUT = open(mydir + "data/fitness/summary.txt", 'w')
    path = mydir + "data/fitness"
    print>> OUT, 'Pop', 'Time', 'Rep', 'Mean', 'Std'
    for file in os.listdir(path):
        if file.endswith(".txt"):
            file_path = os.path.join(path, file)
            if ('N_' in file_path or 'M_' in file_path) and '_T0' not in file_path:
                IN = pd.read_csv(file_path, sep = '\t')
                values = []
                file_split = file_path.split('_')
                if 'N_' in file_path:
                    in_zip = zip(IN.fitness.values.tolist(), IN.n_i.values.tolist())
                else:
                    in_zip = zip(IN.fitness.values.tolist(), IN.m_i.values.tolist())
                for i in in_zip:
                    values.extend([round(i[0], 6)] * i[1])
                mean = np.mean(values)
                std = np.std(values)
                #print file_split
                print>> OUT, file_split[0][-1], file_split[-2][1:], file_split[-1][:-4], str(mean), str(std)
    OUT.close()

def plot_v_dorm():
    IN = pd.read_csv(mydir + "data/fitness/summary.txt", sep = ' ')
    t1_M = IN.loc[(IN['Time'] == 80000) & (IN['Pop'] == 'M')]
    t2_M = IN.loc[(IN['Time'] == 90000) & (IN['Pop'] == 'M')]
    t1_N = IN.loc[(IN['Time'] == 80000) & (IN['Pop'] == 'N')]
    t2_N = IN.loc[(IN['Time'] == 90000) & (IN['Pop'] == 'N')]
    #t2 = IN.loc[IN['Time'] == 90000]
    #print t1_M
    print (np.mean(t2_M.Mean.values) - np.mean(t1_M.Mean.values))  / 10000
    print (np.mean(t2_N.Mean.values) - np.mean(t1_N.Mean.values))  / 10000

#piMergeDF(10)
#piMergeDF(100)
#piMergeDF(1000)
#piTimeSeries(10)
#piTimeSeries(100)
#piTimeSeries(1000)
#piTimeSeriesAll()
#plot_fitness()
#IN1 = pd.read_csv(mydir + 'data/fitness/T80000.txt', sep = '\t')
#IN2 = pd.read_csv(mydir + 'data/fitness/T90000.txt', sep = '\t')

#values1 = []
#in_zip1 = zip(IN1.fitness.values.tolist(), IN1.n_i.values.tolist())
#for i in in_zip1:
#    values1.extend([round(i[0], 6)] * i[1])
#values2 = []
#in_zip2 = zip(IN2.fitness.values.tolist(), IN2.n_i.values.tolist())
#for i in in_zip2:
#    values2.extend([round(i[0], 6)] * i[1])

#print np.mean(values2)
#print np.mean(values1)
#print np.log10((np.mean(values2) - np.mean(values1)) / 10000)
#clean_v_dorm()
plot_v_dorm()
