from __future__ import print_function

import numpy as np
from scipy.stats import norm
from scipy.stats import rv_continuous

import resource

def populate_triu(P, c_rv):
    '''
    Populates the upper triangle of P with connection probabilities
    sampled from c_rv. P_ij is the probability for a connection from
    node i to node j.  
    ----------
    P     :  matrix of connection probabiities
    c_rv  :  randomly returns connection probability
    '''

    rows, cols = np.triu_indices_from(P, k=1)
    
    for i,j in zip(rows,cols):
        P[i][j] = c_rv.rvs()

    assert(np.max(P)<=1)
    assert(np.min(P)>=0)
        
    return P
        

def populate_tril(P, c_rvy, params):
    '''
    Populates the lower triangle of P with connection probabilities
    sampled from c_rv, given the corr_method.
    ----------
    P           :  matrix of connection probabiities
    c_rvy      :  class for sampling correlated connection probabilities
    '''

    rows, cols = np.triu_indices_from(P, k=1)
    
    for i,j in zip(rows,cols):
        c_rvy.x = P[i][j]
        c_rvy.sigma = params["sigma"]
        c_rvy.norm_f = c_rvy.compute_norm_f()
        P[j][i] = c_rvy.rvs()
        
    assert(np.max(P)<=1)
    assert(np.min(P)>=0)

    return P
        

class C_rv_mult_norm(rv_continuous):
    '''
    normal modulated random variable for P_ji given P_ij=x
    and width sigma
    '''

    def __init__(self, c_rv, x, sigma):

        self.c_rv = c_rv
        self.x = x
        self.sigma = sigma

        self.norm_f = self.compute_norm_f()

        rv_continuous.__init__(self, a=0, b=1)

    def compute_norm_f(self):
        return integrate.quad(lambda y: self.c_rv.pdf(y)\
                              *norm.pdf(y, loc=self.x, scale=self.sigma),
                              0, 1)[0]
        
    def _pdf(self, y):
        return self.c_rv.pdf(y)*1/self.norm_f\
               * norm.pdf(y, loc=self.x, scale=self.sigma)
        
            

def connect_network(P):

    rel = np.random.uniform(size=np.shape(P))

    return (P>rel).astype(int)


def generate_network(N, c_rv, params):

    P = np.zeros((N,N))
    P = populate_triu(P, c_rv)
    c_rvy = C_rv_mult_norm(c_rv, 1, params["sigma"])
    P = populate_tril(P, c_rvy, params)
    G = connect_network(P)
    
    return G

def compute_overrep(G):

    N = len(G)
    
    U = (G+G.T)[np.triu_indices(N,1)]

    n_recip  = float(len(U[np.where(U==2)]))
    
    p_bar = np.sum(G)/float((N*(N-1)))                
    n_recip_bar = p_bar**2*N*(N-1)/2.

    return n_recip/n_recip_bar


  
if __name__ == '__main__':

    import os, pickle
    
    from gamma_functions import *

    a = 1.
    mu = 0.1
    
    c_rv = Trunc_gamma(a, mu)

    N = 250
    trials = 5

    df = []

    label = "a_1_sig_05_set5.p"
    with open("data/" + label, "wb") as pfile:
        pickle.dump(df, pfile)

    for sig in np.arange(0.05, 0.45, 0.05):
        params = {"sigma": sig}
        for i in range(trials):
            data = {}
            G = generate_network(N, c_rv, params)
            rho = compute_overrep(G)
            data["rho"] = rho
            data["sig"] = params["sigma"]
            data["mu"]  = float(np.sum(G))/(N*(N-1))
            data["N"] = len(G)
            data["alpha"] = c_rv.alph
            data["beta"] = c_rv.bet

            df.append(data)

            print( sig )
            print( "\t rho         : ", rho )
            print( "\t mu          : ", float(np.sum(G))/(N*(N-1)))
            
            os.rename("data/" + label, "tmp/" + label)
            with open("data/" + label, "wb") as pfile:
                pickle.dump(df, pfile)
