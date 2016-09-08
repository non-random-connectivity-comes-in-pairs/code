
import numpy as np

from gamma_functions import *

from scipy.optimize import brentq 
import scipy.integrate as integrate

from scipy.stats import gamma

# from gamma_figure_A
a = 0.248
b = 0.486852833106356

def get_rho(a,b):    
    numer = integrate.quad(lambda x:x**2*fT(x,a,b),0,1)[0]
    denom = (integrate.quad(lambda x: x*fT(x,a,b), 0,1)[0])**2    
    return numer/denom

assert abs(get_rho(a,b) - 4.) < 10e-3

p_cut = 0.01

# calculate 0.01 percentage
r = integrate.quad(lambda z: K(a,b)*gamma.pdf(z, a, scale=b), p_cut,1.)[0]

print "%.4f of pairs are connected with chance higher than %.2f" %(r,p_cut)

