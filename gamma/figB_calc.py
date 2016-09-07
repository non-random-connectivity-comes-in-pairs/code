
import numpy as np

from gamma_functions import *

from scipy.optimize import brentq 
import scipy.integrate as integrate


# mu = \int_0^1 x f_ab(x) dx
# Given mu,a what is b?

def root_f(b):

    return integrate.quad(lambda x: x*fT(x,a,b),0,1)[0] - mu


mu = 0.1

a_ns = []
b_ns = []

for a in np.logspace(np.log10(0.2),np.log10(98),500):

    a_ns.append(a)

    b = brentq(root_f, 0.5*mu/a, 5*mu/a)

    b_ns.append(b)

    print a, " ", b


np.save("data/a_ns.npy", a_ns)
np.save("data/b_ns.npy", b_ns)
