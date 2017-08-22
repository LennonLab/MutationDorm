from __future__ import division
import os
import numpy as np
import scipy.stats as st
import scipy.special as sp

mydir = os.path.expanduser("~/github/MutationDorm/")

#def dfe(beta, sigma):

#class dfe(st.rv_continuous):
#    def _pdf(self,s, sigma):

#        return np.exp(- (s/sigma) )  # Normalized over its range, in this case [0,1]

#my_cv = dfe(a=0, b='inf', name='my_pdf', shapes = 'sigma')
#print my_cv.rvs(size = 10, sigma = 2)

class dfe(st.rv_continuous):
    def _pdf(self, s, sigma, beta):
        #print type(st.gamma(1 + (1/beta)))
        #term1 = np.asarray(np.exp(- ((s/sigma) ** beta))) / np.asarray(st.gamma(1 + (1/beta)))
        #term1 = np.asarray(np.exp(- ((s/sigma) ** beta))) / np.asarray(st.gamma(1 + (1/beta)))
        #return (1/sigma) * term1 # Normalized over its range, in this case [0,inf]
        term1 =  np.exp(- ((s/sigma)) ) / sigma
        term2 = sp.gamma(1 + (1/beta))
        return term1 / term2

my_cv = dfe(a=0, b='inf', name='my_pdf', shapes='sigma, beta')
print my_cv.rvs(size = 10, sigma = 0.01, beta = 2)
#samples = my_cv.rvs(size = i, sigma = sigma)

#def get_mut(i, sigma, beta):
#    my_cv = dfe(a=0, b='inf', name='my_pdf', shapes='sigma')
#    samples = my_cv.rvs(size = i, sigma = sigma)
    #print samples
#    print  sp.gamma(1 + (1/beta))
#    samples=  np.asarray(list(samples))
#    return samples / 2.3

#print get_mut(i = 4, sigma = 0.01, beta = 4)

def sim_adapt(s, U_b, N):
    pop = np.zeros(N)
    pop[0] = N
    # selection THEN mutation

    print pop



#sim_adapt(0.1, 0.00001, 1000000)
