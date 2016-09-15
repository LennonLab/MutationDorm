from __future__ import division
import generatePop as gp
import quantifyPop as qp
import timeit
import os

mutation_rate = 0.00001
active_size = 1000
dormant_size = 1000
seq_length = 100
generations = 10
c = 10

def simulate_history_OUT(iterations, mutation_rate, generations, seq_length, active_size, dormant_size, c):
    for x in range(iterations):
        OUT = open('/N/dc2/projects/Lennon_Sequences/2016_MutationDorm/data/sim_history/Cs/G' \
            + str(generations) + '_S' + str(seq_length) + '_N' + str(active_size) + \
            '_M' + str(dormant_size) +  '_c' + str(c) + '_I' + str(x) + '.txt', 'w+')
        print>> OUT, 'iteration', 'WT', 'PI', 'FT', 'TD', 'FD'
        pop = gp.generate_pop(active_size, dormant_size, seq_length)
        sim = simulate_save_history(OUT, pop, mutation_rate, generations, seq_length, active_size, dormant_size, c)
        print dormant_size, x
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

#Ms = [10000]
Cs = [10, 100, 1000]
for C in Cs:
    #simulateOUT(100, mutation_rate, 10000, seq_length, active_size, M, c)
    simulate_history_OUT(100, mutation_rate, 10000, seq_length, active_size, dormant_size, C)
