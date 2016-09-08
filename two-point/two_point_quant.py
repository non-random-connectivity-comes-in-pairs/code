
import numpy as np

mu = 0.1
rho = 4.
x = 0.7

y = 1./(1/mu - x/mu**2) * (rho-x/mu)

print "y=", y

p = (mu-y)/(x-y)

print "p=", p


