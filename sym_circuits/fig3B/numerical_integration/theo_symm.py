
import numpy as np

from scipy import integrate
from scipy.optimize import brentq
from scipy.stats import norm

from gamma_functions import K, fT



class Gamma_network(object):

    def __init__(self, a, mu):

        self.a = a
        self.mu = mu

        self.b = self.get_b()

        self.sigma = 1.

        
    def get_b(self):

        root_f = lambda b: integrate.quad(
            lambda x: x*fT(x,self.a,b),0,1)[0] - self.mu

        return brentq(root_f, 0.5*self.mu/self.a, 5*self.mu/self.a)


    def xyfxfye2(self,x,y):
        numer =  x*y*fT(x,self.a,self.b)*fT(y,self.a,self.b)\
                 *norm.pdf(y,loc=x, scale=self.sigma) 
        denom =  integrate.quad(lambda y: fT(y,self.a,self.b)\
                                * norm.pdf(y,loc=x, scale=self.sigma),
                                0,1)[0]
        return numer/denom


    def yfxfye2(self,x,y):
        numer =  y*fT(x,self.a,self.b)*fT(y,self.a,self.b)\
                 *norm.pdf(y,loc=x, scale=self.sigma) 
        denom =  integrate.quad(lambda y: fT(y,self.a,self.b)\
                                *norm.pdf(y,loc=x, scale=self.sigma)
                                ,0,1)[0]
        return numer/denom
    

    def sim(self,sigma):

        self.sigma = sigma

        numer = integrate.dblquad(self.xyfxfye2, 0, 1,
                                  lambda x: 0, lambda x: 1)[0]
        denom_x = integrate.quad(lambda x: x*fT(x,self.a,self.b), 0,1)[0]
        denom_y = integrate.dblquad(self.yfxfye2, 0, 1,
                                    lambda x: 0, lambda x: 1)[0]
        return numer, denom_x, denom_y
        


if __name__=="__main__":

    import os, pickle
    
    label = "th_a1_mu01.p"

    gn = Gamma_network(1, 0.1)
    data_frame = []

    with open("data/" + label, "wb") as pfile:
        pickle.dump(data_frame, pfile)
    
    for sigma in np.arange(0.005,0.75005,0.005):

        numer, denom_x, denom_y  = gn.sim(sigma)
        print( sigma )
        print( "\t numer       : ", numer )
        print( "\t den_x, den_y: ", denom_x, ", ", denom_y )
        print( "\t rho         : ", numer/(denom_x*denom_y) )
        print( "\t alph, bet   : ", gn.a, ", ", gn.b )	

        data = {}
        data["alpha"] = gn.a
        data["beta"] = gn.b
        data["mu"] = gn.mu
        data["sigma"] = sigma
        data["numer"] = numer
        data["den_x"] = denom_x
        data["den_y"] = denom_y

        data_frame.append(data)
        
        os.rename("data/" + label, "tmp/" + label)
        with open("data/" + label, "wb") as pfile:
            pickle.dump(data_frame, pfile)

