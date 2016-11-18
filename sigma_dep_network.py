
import numpy as np
from scipy.stats import norm
from scipy.stats import rv_continuous
from scipy.special import gamma
import scipy.integrate as integrate
from scipy.optimize import brentq


class Rv_Mult_Norm(rv_continuous):

    def __init__(self, rv, x, sigma):

        self.rv = rv
        self.x = x
        self.sigma = sigma

        self.norm_f = self.compute_norm_f()

        rv_continuous.__init__(self, a=0, b=1)

    def compute_norm_f(self):
        return integrate.quad(lambda y: self.rv.pdf(y)*norm.pdf(y, loc=self.x, scale=self.sigma), 0, 1)[0]
        
    def _pdf(self, y):
        return self.rv.pdf(y)*norm.pdf(y, loc=self.x, scale=self.sigma)/self.norm_f
        
            

def sample_rv_mult_norm(x, rv, sigma):

    rv_mult_norm = Rv_Mult_Norm(rv, x, sigma)

    return rv_mult_norm.rvs()



def K(a,b):

    int_f_0_1 = integrate.quad(
        lambda x: 1/(b**a*gamma(a))*x**(a-1)*np.exp(-x/b),0,1)
    K_ab = 1./(int_f_0_1[0])

    return K_ab


def fT(x,a,b):
    '''
    computes the probability density function $f^T_{\alpha, \beta}$ 
    of the truncated gamma distribution $\Gamma^{T}(\alpha,\beta)$
    '''

    if x < 0:
        fT_ab = 0
    elif x>1:
        fT_ab = 0
    else:
        K_ab = K(a,b)

        fT_ab = K_ab * 1/(b**a*gamma(a))*x**(a-1)*np.exp(-x/b)

    return fT_ab


class Trunc_Gamma_Rv(rv_continuous):

    def __init__(self, alph, bet, mu):

        self.alph = alph
        self.bet = bet
        self.mu = mu
        
        assert(abs(bet-self.get_b())<0.01)
        self.K_ab = self.get_K()

        rv_continuous.__init__(self, a=0, b=1)

    def get_b(self):
        root_f = lambda bet: integrate.quad(
            lambda x: x*fT(x,self.alph,bet),0,1)[0] - self.mu
        return brentq(root_f, 0.5*self.mu/self.alph, 5*self.mu/self.alph)

    def get_K(self):
        a,b = self.alph, self.bet
        int_f_0_1 = integrate.quad(
            lambda x: 1/(b**a*gamma(a))*x**(a-1)*np.exp(-x/b),0,1)
        K_ab = 1./(int_f_0_1[0])
        return K_ab


    def _pdf(self, x):
        if x < 0:
            fT_ab = 0
        elif x>1:
            fT_ab = 0
        else:
            fT_ab = self.K_ab * 1/(self.bet**self.alph*gamma(self.alph))*x**(self.alph-1)*np.exp(-x/self.bet)

        return fT_ab



   
if __name__ == '__main__':

    from params.sigma_dep_network_params import *
    import os, pickle

    data_frame = []
    with open("data/data_label.p", "wb") as pfile:
        pickle.dump(data_frame, pfile)

    for alph, bet in zip(alphas,betas):
                
        rv = Trunc_Gamma_Rv(alph, bet, mu)

        for sigma in np.arange(sig_low, sig_high, sig_step):

            for k in range(n_trials):

                print("a: ", alph, "\t sigma: ", sigma, "\t No.: ", k+1)
                
                data = {"alpha": alph, "beta": bet, "mu": mu, "sigma": sigma}
            
                # xs = rv.rvs(size=n_pairs)
                # ys = []
                # for x in xs:
                #     y = sample_rv_mult_norm(x, rv, sigma)
                #     ys.append(y)
                
                # data["xs"]=xs
                # data["ys"]=np.array(ys)

                data_frame.append(data)

                os.rename("data/data_label.p", "tmp/data_label.p")
                with open("data/data_label.p", "wb") as pfile:
                    pickle.dump(data_frame, pfile)
            

                
            



            
    # N = 500
    # trials = 3

    # params = {"sigma" : 0.1}

    # xs = rv.rvs(size=10000)
    # ys = []
    # for x in xs:
    #     y = f_dependent_mult_norm_crv(x, rv, params)
    #     ys.append(y)
    # #print np.mean(abs(xs-ys))
    # print np.mean(xs*ys)/(np.mean(np.concatenate((xs,ys)))**2)
