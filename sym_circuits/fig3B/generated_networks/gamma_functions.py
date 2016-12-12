
import numpy as np

from scipy.special import gamma
import scipy.integrate as integrate
from scipy.optimize import brentq

from scipy.stats import rv_continuous


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



class Trunc_gamma(rv_continuous):

    def __init__(self, alph, mu):

        self.alph = alph
        self.mu = mu
        self.bet = self.get_b()

        self.K_ab = K(self.alph,self.bet)

        rv_continuous.__init__(self, a=0, b=1)


    def get_b(self):

        root_f = lambda bet: integrate.quad(
            lambda x: x*fT(x,self.alph,bet),0,1)[0] - self.mu

        return brentq(root_f, 0.5*self.mu/self.alph, 5*self.mu/self.alph)


    def _pdf(self, x):
        if x < 0:
            fT_ab = 0
        elif x>1:
            fT_ab = 0
        else:
            fT_ab = self.K_ab * 1/(self.bet**self.alph*gamma(self.alph))*x**(self.alph-1)*np.exp(-x/self.bet)

        return fT_ab
