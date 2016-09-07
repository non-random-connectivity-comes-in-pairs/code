
import numpy as np

from gamma_functions import *

import scipy.integrate as integrate


a_ns = np.load("data/a_ns.npy")
b_ns = np.load("data/b_ns.npy")

rhos = []

for a,b in zip(a_ns,b_ns):

    numer = integrate.quad(lambda x:x**2*fT(x,a,b),0,1)[0]
    denom = (integrate.quad(lambda x: x*fT(x,a,b), 0,1)[0])**2
    
    rho = numer/denom

    print a, "\t", b, "\t", rho

    rhos.append(rho)

np.save("data/rho_theoretical.npy", rhos)
