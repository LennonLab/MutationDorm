from __future__ import division
import numpy as np
import itertools
import math
#import matplotlib.pyplot as plt
#import pylab
#import matplotlib as mpl
import argparse
from collections import Counter
import random

'''
View individual fitness as a function of the number of mutations each individual has
Constant population size
'''

def generate_pop(pop_size):
    pop = {}
    pop[0] = pop_size
    return pop


def choose_by_weight(weights):
    weights = map(float, weights)
    rndm = random.random() * sum(weights)
    for i, j in enumerate(weights):
        rndm -= j
        if rndm < 0:
            return i

def get_mutation_count(U, pop_size):
    mean = U * pop_size
    return np.random.poisson(mean)


def get_random_indiv(pop):
    mutational_classes = pop.keys()
    size_step = sum(pop.values())
    frequencies = [x/float(size_step) for x in pop.values()]
    total = sum(frequencies)
    frequencies = [x / total for x in frequencies]
    return np.random.choice(mutational_classes, p=frequencies)

def mutation_event(pop):
    '''Returns individual of a specific mutational class'''
    indiv = get_random_indiv(pop)
    pop[indiv] -= 1
    if indiv + 1 in pop:
        pop[indiv + 1] += 1
    else:
        pop[indiv + 1] = 1

def mutation_drift_steps(pop, U):
    N = sum(pop.values())
    mutation_count = get_mutation_count(N, U)
    for i in range(mutation_count):
        mutation_event(pop)
    frequencies = [x/float(N) for x in pop.values()]
    return list(np.random.multinomial(N, frequencies))


def offspring_step(pop, U):
    counts = mutation_drift_steps(pop, U)
    for (mutationClass, count) in zip(pop.keys(), counts):
        if (count > 0):
            pop[mutationClass] = count
        else:
            del pop[mutationClass]
    #return pop

def good_bad_not_evil(generations, good_freq):
    bad_freq = 1 - good_freq
    good_bad_string = ''
    for i in range(0,generations):
        freqs = [good_freq, bad_freq]
        j = choose_by_weight(freqs)
        if j == 0:
            good_bad_string += '1'
        else:
            good_bad_string += '0'
    return good_bad_string


def time_step(pop, U):
    offspring_step(pop, U)


def dormancy_step(pop, seed_bank, c):
    N = sum(pop.values())
    M = sum(seed_bank.values())
    K = N / M
    #dorm_prob = c / N
    #resusc_prob = (c*K) / M
    c = int(round(c, 0))
    cK = int(round((K*c), 0))
    for i in range(0, c):
        active_i = get_random_indiv(pop)
        pop[active_i] -= 1
        if active_i in seed_bank:
            seed_bank[active_i] += 1
        else:
            seed_bank[active_i] = 1
    for i in range(0, cK):
        print seed_bank, c, cK
        dormant_i = get_random_indiv(seed_bank)
        seed_bank[dormant_i] -= 1
        if dormant_i in pop:
            pop[dormant_i] += 1
        else:
            pop[dormant_i] = 1
    #print seed_bank

def fitness(pop, s=0.01):
    pop_fitness = []
    N = sum(pop.values())
    for mutation_class, n in pop.iteritems():
        fitness_class = (1-s) ** mutation_class
        pop_fitness.append(fitness_class * n)
    return np.mean(pop_fitness)
    #print pop.values()

def simulate(N, U, delata_U,  generations, good_freq, c, M):
    pop = generate_pop(N)
    gen_list = good_bad_not_evil(generations, good_freq)
    if M != 0:
        seed_bank = generate_pop(M)
    for i in range(generations):
        environ_qual = int(gen_list[i])
        if (environ_qual == 0) and (M==0):
            time_step(pop, U + delata_U)
        elif (environ_qual == 0) and (M!=0):

            dormancy_step(pop, seed_bank, c)
            time_step(pop, U + delata_U)
        else:
            time_step(pop, U)
    # get mean and SD of W
    return fitness(pop, s=0.01)


def many_simulations(N, U, delata_U,  generations, good_freq, c, M):
    fitness_list = []
    for x in range(generations):
        fitness_x = simulate(N, U, delata_U,  generations, good_freq, c, M)
        #print fitness_x
        fitness_list.append(fitness_x)
    print np.mean(fitness_list)


U = 0.1
delata_U = 0.2
generations = 100
good_freq = 0
N = 100
c = 1
M = 50
many_simulations(N, U, delata_U,  generations, good_freq, c, M)


### So there's sampling error and the dormant pool decreases to zero
# need to fix this
