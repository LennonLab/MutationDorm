from __future__ import division
import generatePop as gp
import quantifyPop as qp
import numpy as np
import pandas as pd
import timeit
import os
from scipy.stats.stats import pearsonr
import  matplotlib.pyplot as plt
import matplotlib.patheffects as pe


mydir = os.path.expanduser("~/github/MutationDorm/")

mutation_rate = 0.0001
active_size = 1000
#dormant_size = 500
seq_length = 1000
generations = 10
c = 10

def get_history(mutation_rate, generations, seq_length, active_size, dormant_size, c):
    pop = gp.generate_pop(active_size, dormant_size, seq_length, random = False)
    history = gp.simulate_history(pop, mutation_rate, generations, seq_length, active_size, dormant_size, c)
    return history
    #qp.stacked_trajectory_plot(history, generations, \
    #    name = mydir + 'figs/trajectory_N_' + str(active_size)+ '_M_' + str(dormant_size) + '.png')

def get_history_custom(mutation_rate, generations, seq_length, active_size, dormant_size, c):
    history = []
    pop = gp.generate_WF_pop(mutation_rate, generations, seq_length, active_size, dormant_size)
    for gen in range(generations):
        print dormant_size, gen
        sim = gp.simulate(pop, mutation_rate, generations, seq_length, active_size, dormant_size, c)
        pop_merged = gp.merge_two_dicts(sim['Active'], sim['Dormant'])
        history.append(pop_merged)
    qp.stacked_trajectory_plot(history, generations, name = mydir + 'figs/trajectory_M_' + str(dormant_size) + '.png')

def simulateOUT(iterations, mutation_rate, generations, seq_length, active_size, dormant_size, c):
    OUT = open(mydir + 'data/G' + str(generations) + '_S' + str(seq_length) + '_N' + str(active_size) + \
        '_M' + str(dormant_size) +  '_c' + str(c) + '.txt', 'w+')
    print>> OUT, 'iteration', 'WT', 'PI', 'FT', 'TD', 'FD'
    for x in range(iterations):
        pop = gp.generate_pop(active_size, dormant_size, seq_length)
        sim = gp.simulate(pop, mutation_rate, generations, seq_length, active_size, dormant_size, c)
        print dormant_size, x
        pop_merged = gp.merge_two_dicts(sim['Active'], sim['Dormant'])
        FD = qp.fu_and_li_D(pop_merged, seq_length)
        TD = qp.tajimas_D(pop_merged, seq_length)
        WT = qp.wattersons_theta(pop_merged, seq_length)
        FT = qp.fu_and_li_theta(pop_merged, seq_length)
        PI = qp.tajimas_theta(pop_merged)
        print>> OUT, x, WT, PI, FT, TD, FD
    OUT.close()


def simulate_save_history(OUT, pop, mutation_rate, generations, seq_length, active_size, dormant_size, c):
    for g in range(generations):
        gp.time_step(pop, mutation_rate, seq_length, active_size, dormant_size, c)
        pop_merged = gp.merge_two_dicts(pop['Active'], pop['Dormant'])
        FD = qp.fu_and_li_D(pop_merged, seq_length)
        TD = qp.tajimas_D(pop_merged, seq_length)
        WT = qp.wattersons_theta(pop_merged, seq_length)
        FT = qp.fu_and_li_theta(pop_merged, seq_length)
        PI = qp.tajimas_theta(pop_merged)
        print>> OUT, g, WT, PI, FT, TD, FD
    return pop


def simulate_history_OUT(iterations, mutation_rate, generations, seq_length, active_size, dormant_size, c):
    for x in range(iterations):
        OUT = open(mydir + 'data/sim_history/Cs/G' + str(generations) + '_S' + str(seq_length) + '_N' + str(active_size) + \
            '_M' + str(dormant_size) +  '_c' + str(c) + '_I' + str(x) + '.txt', 'w+')
        print>> OUT, 'iteration', 'WT', 'PI', 'FT', 'TD', 'FD'
        pop = gp.generate_pop(active_size, dormant_size, seq_length)
        sim = simulate_save_history(OUT, pop, mutation_rate, generations, seq_length, active_size, dormant_size, c)
        print dormant_size, x
        OUT.close()

def correlationPlot(generations, active_size, dormant_size = 1000):
    #Cs = [10, 100, 1000]
    Cs = [1]
    r2s = []
    gens = range(generations+1)[:: int(generations/10)]
    r2_gens = [(gens[i], gens[i+1] ) for i in range(len(gens)-1) ]
    for C in Cs:
        r2_c = []
        for x in range(100):
            print C, x
            hist = get_history(mutation_rate, generations, seq_length, active_size, dormant_size, C)
            r2_c_x = []
            for r2_gen in r2_gens:
                hist_gen = pd.DataFrame([ hist[r2_gen[0]],   hist[r2_gen[1]]])
                hist_gen = hist_gen.fillna(0)
                #hist_r2_gen = hist.iloc[[[0], r2_gen[1]],:]
                hist_r2_gen_clean = hist_gen.loc[:, (hist_gen != 0).any(axis=0)]
                freq0 = hist_r2_gen_clean.values[0] / (active_size + dormant_size)
                freq1 = hist_r2_gen_clean.values[1] / (active_size + dormant_size)
                r2 = pearsonr(freq0, freq1)[0]
                r2_c_x.append(r2)
            #mean = sum(r2_c) / len(r2_c)
            r2_c.append(r2_c_x)
        df_T = pd.DataFrame(r2_c)
        #df_T.columns = ['10', '100', '1000']
        df_T.to_csv(mydir + 'data/corr_' + str(C) + '.txt', sep=' ', index = False)


def squared(x):
    # that, if x is a string,
    if type(x) is str:
        # just returns it untouched
        return x
    # but, if not, return it multiplied by 100
    elif x:
        return  x**2
    # and leave everything else
    else:
        return

def r2_time():
    df_1 = pd.read_csv(mydir + 'data/corr_1.txt', sep = ' ', header = 0)
    df_1 = df_1.applymap(squared)
    df_10 = pd.read_csv(mydir + 'data/corr_10.txt', sep = ' ', header = 0)
    df_10 = df_10.applymap(squared)
    df_100 = pd.read_csv(mydir + 'data/corr_100.txt', sep = ' ', header = 0)
    df_100 = df_100.applymap(squared)
    df_1000 = pd.read_csv(mydir + 'data/corr_1000.txt', sep = ' ', header = 0)
    df_1000 = df_1000.applymap(squared)
    std_1 = df_1.std(axis=0)
    mean_1 = df_1.mean(axis=0)
    std_10 = df_10.std(axis=0)
    mean_10 = df_10.mean(axis=0)
    std_100 = df_100.std(axis=0)
    mean_100 = df_100.mean(axis=0)
    std_1000 = df_1000.std(axis=0)
    mean_1000 = df_1000.mean(axis=0)

    x = [0,1,2,3,4,5,6,7,8,9]
    fig, ax = plt.subplots()
    ax.plot(x, mean_1, label='Time in seed bank = 1000', \
    path_effects=[pe.Stroke(linewidth=2.7, foreground='b'), pe.Normal()], lw = 2, color = '#9370DB')
    ax.plot(x, mean_10, label='Time in seed bank = 100', \
    path_effects=[pe.Stroke(linewidth=2.7, foreground='b'), pe.Normal()], lw = 2, color = '#87CEEB')
    ax.plot(x, mean_100, label='Time in seed bank = 10', \
        path_effects=[pe.Stroke(linewidth=2.7, foreground='b'), pe.Normal()], lw = 2, color = '#FFA500')
    ax.plot(x, mean_1000, label='Time in seed bank = 1', \
        path_effects=[pe.Stroke(linewidth=2.7, foreground='b'), pe.Normal()], lw = 2, color = '#FF6347')
    ax.fill_between(x, mean_1+std_1, mean_1-std_1, facecolor='#9370DB', alpha=0.5)
    ax.fill_between(x, mean_10+std_10, mean_10-std_10, facecolor='#87CEEB', alpha=0.5)
    ax.fill_between(x, mean_100+std_100, mean_100-std_100, facecolor='#FFA500', alpha=0.5)
    ax.fill_between(x, mean_1000+std_1000, mean_1000-std_1000, facecolor='#FF6347', alpha=0.5)

    x_ticks = ['0-100', '100-200', '200-300', '300-400', '400-500', '500-600', \
        '600-700', '700-800', '800-900', '900-1,000']
    plt.xticks(x, x_ticks, fontsize = 8)
    ax.set_ylim(0,1)
    ax.set_ylabel( r'$r^{2}_{t-1, t}$', fontsize=20)
    ax.set_xlabel('Time between sampling (generations)' , fontsize=20)
    ax.legend(loc='lower left', fontsize = 12)

    plt.savefig(mydir +'figs/r2_plot.png', bbox_inches='tight',  dpi = 600)
    plt.close()


#get_history(mutation_rate, 10, seq_length, 10, 0, 0)
#Ms = [10, 100, 1000, 10000]
#Ms = [10000]
#for M in Ms:
    #simulateOUT(100, mutation_rate, 10000, seq_length, active_size, M, c)
#    simulate_history_OUT(100, mutation_rate, 10000, seq_length, active_size, M, c, Mason = True)

    #get_history(mutation_rate, 10000, seq_length, active_size, M, c)

#Cs = [10, 100, 1000]
#for C in Cs:
#    simulate_history_OUT(100, mutation_rate, 10000, seq_length, active_size, 1000, C)
#1000
#correlationPlot(1000, active_size)
r2_time()
#correlationPlot(1000, active_size)
