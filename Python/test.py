from __future__ import division
import generatePop as gp
import quantifyPop as qp
import timeit
import os

mydir = os.path.expanduser("~/github/MutationDorm/")

mutation_rate = 0.00001
active_size = 1000
#dormant_size = 500
seq_length = 100
generations = 10
c = 10

def get_history(mutation_rate, generations, seq_length, active_size, dormant_size, c):
    pop = gp.generate_pop(active_size, dormant_size, seq_length, random = False)
    history = gp.simulate_history(pop, mutation_rate, generations, seq_length, active_size, dormant_size, c)
    qp.stacked_trajectory_plot(history, generations, \
        name = mydir + 'figs/trajectory_N_' + str(active_size)+ '_M_' + str(dormant_size) + '.png')

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


#get_history(mutation_rate, 10, seq_length, 10, 0, 0)
#Ms = [10, 100, 1000, 10000]
#Ms = [10000]
#for M in Ms:
    #simulateOUT(100, mutation_rate, 10000, seq_length, active_size, M, c)
#    simulate_history_OUT(100, mutation_rate, 10000, seq_length, active_size, M, c, Mason = True)

    #get_history(mutation_rate, 10000, seq_length, active_size, M, c)

Cs = [10, 100, 1000]
for C in Cs:
    simulate_history_OUT(100, mutation_rate, 10000, seq_length, active_size, 1000, C)
