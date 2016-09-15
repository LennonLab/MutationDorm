from __future__ import division
import numpy as np
import itertools
import math
import argparse
from collections import Counter
import random
#from collections import defaultdict


def generate_base_haplotype(seq_length):
    base_haplotype = ''.join(["A" for i in range(seq_length)])
    return base_haplotype

def random_dna_sequence(seq_length):
    return ''.join(np.random.choice(list('ACTG')) for _ in range(seq_length))

def generate_pop(active_size, dormant_size, seq_length, random = False):
    pop = {}
    pop['Active'] = {}
    pop['Dormant'] = {}
    if random == True:
        for x in range(active_size):
            new_haplotype = random_dna_sequence(seq_length)
            if new_haplotype in pop['Active']:
                pop['Active'][new_haplotype] += 1
            else:
                pop['Active'][new_haplotype] = 1
        for x in range(dormant_size):
            new_haplotype = random_dna_sequence(seq_length)
            if new_haplotype in pop['Dormant']:
                pop['Dormant'][new_haplotype] += 1
            else:
                pop['Dormant'][new_haplotype] = 1

    else:
        base_haplotype = generate_base_haplotype(seq_length)
        pop['Active'][base_haplotype] = active_size
        pop['Dormant'][base_haplotype] = dormant_size
    return pop


def choose_by_weight(weights):
    weights = map(float, weights)
    rndm = random.random() * sum(weights)
    for i, j in enumerate(weights):
        rndm -= j
        if rndm < 0:
            return i

def good_bad_not_evil(generations, good_freq):
    bad_freq = 1 - good_freq
    good_bad_string = ''
    for i in range(0,generations):
        freqs = [good_freq, bad_freq]
        j = choose_by_weight(freqs)
        if j == 0:
            good_bad_string += 'G'
        else:
            good_bad_string += 'B'
    return good_bad_string

# mutation
def get_mutation_count(mutation_rate, active_size, seq_length):
    mean = mutation_rate * active_size * seq_length
    return np.random.poisson(mean)

def get_random_haplotype(pop, sub_pop):
    # Need random haplotype for active only
    if sub_pop == 'Active':
        #try:
        haplotypes = pop['Active'].keys()
        active_size_step = sum(pop['Active'].values())
        frequencies = [x/float(active_size_step) for x in pop['Active'].values()]
        total = sum(frequencies)
        frequencies = [x / total for x in frequencies]
        return np.random.choice(haplotypes, p=frequencies)
    elif sub_pop == 'Dormant':
        #except KeyError:
        haplotypes = pop['Dormant'].keys()
        dormant_size_step = sum(pop['Dormant'].values())
        frequencies = [x/float(dormant_size_step) for x in pop['Dormant'].values()]
        total = sum(frequencies)
        frequencies = [x / total for x in frequencies]
        return np.random.choice(haplotypes, p=frequencies)

def get_mutant(haplotype, seq_length):
    alphabet = ['A', 'T', 'G', 'C']
    site = np.random.randint(seq_length)
    possible_mutations = list(alphabet)
    possible_mutations.remove(haplotype[site])
    mutation = np.random.choice(possible_mutations)
    new_haplotype = haplotype[:site] + mutation + haplotype[site+1:]
    return new_haplotype

def mutation_event(pop, seq_length):
    sub_pop = 'Active'
    haplotype = get_random_haplotype(pop, sub_pop)
    if pop['Active'][haplotype] > 1:
        pop['Active'][haplotype] -= 1
        new_haplotype = get_mutant(haplotype, seq_length)
        if new_haplotype in pop['Active']:
            pop['Active'][new_haplotype] += 1
        else:
            pop['Active'][new_haplotype] = 1

def mutation_step(pop, mutation_rate, active_size, seq_length):
    mutation_count = get_mutation_count(mutation_rate, active_size, seq_length)
    for i in range(mutation_count):
        mutation_event(pop, seq_length)

# reproduce active pop, drift acts here
def get_offspring_counts(pop, active_size):
    haplotypes = pop['Active'].keys()
    frequencies = [x/float(active_size) for x in pop['Active'].values()]
    total = sum(frequencies)
    frequencies = [x / total for x in frequencies]
    return list(np.random.multinomial(active_size, frequencies))

def offspring_step(pop, active_size):
    counts = get_offspring_counts(pop, active_size)
    for (haplotype, count) in zip(pop['Active'].keys(), counts):
        if (count > 0):
            pop['Active'][haplotype] = count
        else:
            del pop['Active'][haplotype]


def dormancy_step(pop, active_size, dormant_size, c):
    if dormant_size <= 0:
        pass
    else:
        K = active_size / dormant_size
        for i in range(0, c):
            new_haplotype = get_random_haplotype(pop, 'Active')
            pop['Active'][new_haplotype] -= 1
            if new_haplotype in pop['Dormant']:
                pop['Dormant'][new_haplotype] += 1
            else:
                pop['Dormant'][new_haplotype] = 1
        for i in range(0, c):
            new_haplotype = get_random_haplotype(pop, 'Dormant')
            pop['Dormant'][new_haplotype] -= 1
            if new_haplotype in pop['Active']:
                pop['Active'][new_haplotype] += 1
            else:
                pop['Active'][new_haplotype] = 1


def time_step(pop, mutation_rate, seq_length, active_size, dormant_size, c):
    if mutation_rate != 0:
        mutation_step(pop, mutation_rate, active_size, seq_length)
    offspring_step(pop, active_size)
    dormancy_step(pop, active_size, dormant_size, c)



def merge_two_dicts(x, y):
    '''Given two dicts, merge them into a new dict as a shallow copy.'''
    x = Counter(x)
    y = Counter(y)
    z = x + y
    return z


def simulate(pop, mutation_rate, generations, seq_length, active_size, dormant_size, c):
    for i in range(generations):
        time_step(pop, mutation_rate, seq_length, active_size, dormant_size, c)
    return pop

def simulate_save_history(OUT, pop, mutation_rate, generations, seq_length, active_size, dormant_size, c):
    for i in range(generations):
        time_step(pop, mutation_rate, seq_length, active_size, dormant_size, c)
        return pop

def simulate_history(pop, mutation_rate, generations, seq_length, active_size, dormant_size, c):
    history = []
    history.append(merge_two_dicts(pop['Active'], pop['Dormant']))
    for i in range(generations):
        print active_size, dormant_size, c, i
        time_step(pop, mutation_rate, seq_length, active_size, dormant_size, c)
        pop_merged = merge_two_dicts(pop['Active'], pop['Dormant'])
        history.append(pop_merged)
    return history

def generate_WF_pop(mutation_rate, generations, seq_length, active_size, dormant_size):
    pop = generate_pop(active_size + dormant_size, 0, seq_length, random = False)
    simulate(pop, mutation_rate, generations, seq_length, active_size + dormant_size, 0, 0)
    countdown = dormant_size
    while countdown > 0:
        haplotype = random.sample(pop['Active'], 1)[0]
        print haplotype
        if pop['Active'][haplotype] > 0:
            countdown -= 1
            pop['Active'][haplotype] -= 1
            if haplotype in pop['Dormant']:
                pop['Dormant'][haplotype] += 1
            else:
                pop['Dormant'][haplotype] = 1
    return pop
