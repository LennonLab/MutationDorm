from __future__ import division
import os
import numpy as np
import scipy.stats as st
import scipy.special as sp
import pandas as pd
from collections import Counter

mydir = os.path.expanduser("~/github/MutationDorm/")


class dfe(st.rv_continuous):
    def _pdf(self, s, sigma, beta):
        term1 =  np.exp(- ((s/sigma)) ) / sigma
        term2 = sp.gamma(1 + (1/beta))
        return term1 / term2


def sim_adapt(U_b, N, gens = 100, rep = 1):
    # Normalized over its range, in this case [0,inf]
    my_cv = dfe(a=0, b='inf', name='dfe', shapes='sigma, beta')
    pop = {}
    pop[1.0] = N
    decimal_cutoff = int(np.log10(N)) + 1
    gens_save = range(gens)[::10000]
    gens_save.extend([gens])
    for gen in range(gens):
        N_gen = sum(pop.values())
        C = N / N_gen
        X_sum = 0
        for key, value in pop.iteritems():
            X_sum += key * value
        X_mean = X_sum / N_gen
        pop_items = pop.items()
        for pop_item in pop_items:
            fitness = pop_item[0]
            size = pop_item[1]
            # selection
            pois_mean = C * (1 + (fitness - X_mean)) * size
            lamba_i = np.random.poisson(lam=pois_mean)
            pop[fitness] = lamba_i
            # mutation
            muts = np.random.poisson(lam=U_b, size=lamba_i)
            index = np.argwhere(muts == 0)
            muts_keep = np.delete(muts, index)
            # fitness = mutation effects + background fitness
            fits = [round(fitness + sum(my_cv.rvs(size = x, sigma = 0.01, beta = 1)), decimal_cutoff) for x in muts_keep]
            pop[fitness] -= len(fits)
            for fit in fits:
                if fit in pop:
                    pop[fit] += 1
                else:
                    pop[fit] = 1
        pop = {k:v for k,v in pop.items() if v != 0}
        if gen in gens_save:
            print rep, gen
            df_gen = pd.DataFrame(pop.items())
            df_gen.columns = ['fitness', 'n_i']
            df_name = mydir + 'data/fitness/T' + str(gen) + '.txt'
            df_gen.to_csv(df_name, sep='\t', index = False)

def sim_adapt_dorm(U_b, N, M = 1000, c = 10, gens = 100, rep = 1):
    # Normalized over its range, in this case [0,inf]
    my_cv = dfe(a=0, b='inf', name='dfe', shapes='sigma, beta')
    pop_N = {}
    pop_M = {}
    pop_N[1.0] = N
    pop_M[1.0] = M
    decimal_cutoff = int(np.log10(N)) + 1
    gens_save = range(gens)[::10000]
    gens_save.extend([gens])
    print gens_save
    for gen in range(gens):
        N_gen = sum(pop_N.values())
        M_gen = sum(pop_M.values())
        #print N_gen, M_gen
        C = N / N_gen
        X_N_sum = 0
        for key, value in pop_N.iteritems():
            X_N_sum += key * value
        X_N_mean = X_N_sum / N_gen
        #print pop_M.values()
        #print np.asarray(pop_M.values())/ M_gen
        exit_dorm = np.random.choice(pop_M.keys(), p=np.asarray(pop_M.values()) / M_gen, size = c)
        enter_dorm = np.random.choice(pop_N.keys(), p=np.asarray(pop_N.values()) / N_gen, size = c)
        exit_dorm_fits = [round(x , decimal_cutoff) for x in exit_dorm]
        enter_dorm_fits = [round(x , decimal_cutoff) for x in enter_dorm]
        exit_dorm_fits_dict = Counter(exit_dorm_fits)
        enter_dorm_fits_dict = Counter(enter_dorm_fits)
        #print exit_dorm_fits_dict, enter_dorm_fits_dict
        for fit_exit in exit_dorm_fits_dict:
            pop_M[fit_exit] -= exit_dorm_fits_dict[fit_exit]
            if fit_exit in pop_N:
                pop_N[fit_exit] += exit_dorm_fits_dict[fit_exit]
            else:
                pop_N[fit_exit] = exit_dorm_fits_dict[fit_exit]

        for fit_enter in enter_dorm_fits_dict:
            pop_N[fit_exit] -= exit_dorm_fits_dict[fit_exit]
            if fit_enter in pop_M:
                pop_M[fit_enter] += enter_dorm_fits_dict[fit_enter]
            else:
                pop_M[fit_enter] = enter_dorm_fits_dict[fit_enter]
        pop_N = {k:v for k,v in pop_N.items() if v > 0}
        pop_M = {k:v for k,v in pop_M.items() if v > 0}
        pop_N_items = pop_N.items()
        for pop_N_item in pop_N_items:
            fitness = pop_N_item[0]
            size = pop_N_item[1]
            # selection
            pois_mean = C * (1 + (fitness - X_N_mean)) * size
            #if pois_mean< 0:
            #    print (fitness - X_N_mean), size, pois_mean
            lamba_i = np.random.poisson(lam=pois_mean)
            pop_N[fitness] = lamba_i
            # mutation
            muts = np.random.poisson(lam=U_b, size=lamba_i)
            index = np.argwhere(muts == 0)
            muts_keep = np.delete(muts, index)
            # fitness = mutation effects + background fitness
            fits = [round(fitness + sum(my_cv.rvs(size = x, sigma = 0.01, beta = 1)), decimal_cutoff) for x in muts_keep]
            pop_N[fitness] -= len(fits)
            for fit in fits:
                if fit in pop_N:
                    pop_N[fit] += 1
                else:
                    pop_N[fit] = 1
        if gen in gens_save:
            print rep, gen
            df_N_gen = pd.DataFrame(pop_N.items())
            df_M_gen = pd.DataFrame(pop_M.items())
            df_N_gen.columns = ['fitness', 'n_i']
            df_M_gen.columns = ['fitness', 'm_i']
            df_N_name = mydir + 'data/fitness/N_T' + str(gen) + '_R'+ str(rep) + '.txt'
            df_M_name = mydir + 'data/fitness/M_T' + str(gen) + '_R'+ str(rep) + '.txt'
            df_N_gen.to_csv(df_N_name, sep='\t', index = False)
            df_M_gen.to_csv(df_M_name, sep='\t', index = False)


#files = [10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000]
#for f in files:
#    f_name_1 = mydir + 'data/fitness/T' + str(f - 10000) + '.txt'
#    f_name_2 = mydir + 'data/fitness/T' + str(f) + '.txt'
#    df_1 = pd.read_csv(f_name_1, sep = '\t')
#    df_2 = pd.read_csv(f_name_2, sep = '\t')
#    P_1 = df_1.n_i.values / sum(df_1.n_i.values)
#    P_2 = df_2.n_i.values / sum(df_2.n_i.values)
#    delta_f = sum(P_2 * df_2.fitness.values) - sum(P_1 * df_1.fitness.values)
#    delta_t = 10000
#    print delta_f / delta_t
#    #mean = np.mean(df_2.Total.values) - np.mean(df_1.Total.values)
#    #print mean / 10000

#sim_adapt(U_b = 0.00001, N = 100000, gens =  100000)
for rep in range(1, 10):
    sim_adapt_dorm(U_b = 0.0001, N = 100000, M = 100000, gens =  100, c = 100, rep=rep)
