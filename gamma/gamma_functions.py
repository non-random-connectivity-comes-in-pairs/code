
from scipy.special import gamma
import scipy.integrate as integrate

import numpy as np


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




if __name__ == "__main__":

    import pylab as pl
    xs = np.arange(-0.2,1.2,0.01)
    pl.plot(xs,[f(x,5,0.025) for x in xs])
    pl.plot(xs,[f(x,1,0.1) for x in xs])
    pl.show()
