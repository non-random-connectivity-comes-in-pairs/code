
import matplotlib as mpl
mpl.use('Agg')
import pylab as pl

import numpy as np
from gamma_functions import fT, K
from scipy import integrate

from matplotlib import rc

rc('text', usetex=True)
pl.rcParams['text.latex.preamble'] = [
    r'\usepackage{tgheros}',    # helvetica font
    r'\usepackage{sansmath}',   # math-font matching  helvetica
    r'\sansmath'                # actually tell tex to use it!
    r'\usepackage{siunitx}',    # micro symbols
    r'\sisetup{detect-all}',    # force siunitx to use the fonts
]  


fig, ax = pl.subplots(1,1)
fig.set_size_inches(7.5*0.5,2.3)

#fig.suptitle(r'$\mu = 0.1$', fontsize=15)

xs = np.arange(0,1.+0.001, 0.001)

linestyles = ['-',':','--','-.']

from scipy.stats import gamma
from scipy.stats import norm


a_ns =  [0.248, 1.0, 2.0, 15.0]
b_ns =  [0.486852833106356, 0.10004560945459162,
         0.0500000206118813, 0.006666666666666667] 

# idx = 0
# xloc = 0.15
# sig = 0.055

# idx =0
# xloc = 0.085
# sig = 0.035



idx =0
xloc = 0.15
sig = 0.025

rs = []
Nf = integrate.quad(lambda s: fT(s,a_ns[idx],b_ns[idx])*norm.pdf(s,loc=xloc, scale=sig),0,1)[0]
for x in xs:
    rs.append(fT(x, a_ns[idx], b_ns[idx]) *norm.pdf(x,loc=xloc,scale=sig)/Nf)


idx =0
xloc = 0.15
sig = 1.

us = []
Nf = integrate.quad(lambda s: fT(s,a_ns[idx],b_ns[idx])*norm.pdf(s,loc=xloc, scale=sig),0,1)[0]
for x in xs:
    us.append(fT(x, a_ns[idx], b_ns[idx]) *norm.pdf(x,loc=xloc,scale=sig)/Nf)

# ys = []
# for x in xs:
#     ys.append(fT(x, a_ns[idx], b_ns[idx]))

idx =0
xloc = 0.15
sig = 0.065
    
zs = []
Nf = integrate.quad(lambda s: fT(s,a_ns[idx],b_ns[idx])*norm.pdf(s,loc=xloc, scale=sig),0,1)[0]
for x in xs:
    zs.append(fT(x, a_ns[idx], b_ns[idx]) *norm.pdf(x,loc=xloc,scale=sig)/Nf)

zcs = []
for x in xs:
    zcs.append(fT(x, a_ns[idx], b_ns[idx]) *norm.pdf(x,loc=xloc,scale=sig))
    
ws = []
for x in xs:
    ws.append(norm.pdf(x,loc=xloc,scale=sig))
    
ax.plot(xs,rs, 'k', linestyle=':', label=r'$\sigma = 0.025$')
ax.plot(xs,zs, 'k', label=r'$\sigma=0.065$')
ax.plot(xs,us, 'k', linestyle='--', label=r'$\sigma=1$')
#ax.plot(xs,ys)
#ax.plot(xs,zcs)
#ax.plot(xs,ws)

ax.set_ylim(0,16)
ax.set_xlim(0,0.4)

pl.xticks([0.,0.15,0.3,0.4], ['0.0', r'$x=0.15$', '0.3','0.4'])
ax.set_title(r'$f_{P_{ij}}(y) = f_{\alpha,\beta}^T(y),\, \alpha=0.248$', size=13.)

ax.set_ylabel(r'$f_{P_{ji} | P_{ij}}(y \mid x)$')
ax.set_xlabel(r'connection probability $y$')
# ax_left.set_xlim(0,0.3)
# ax_left.set_ylim(0,16)
# ax_left.legend()

# label = ax_left.set_xlabel(r'connection probability $P_{ij}$')
# ax_left.xaxis.set_label_coords(0.8, -0.135)
# ax_left.set_ylabel(r'$f(P_{ij})$')


# for a,b,rho,ls in zip(a_ns, b_ns,rhos,linestyles):

#     ax_right.plot(x, gamma.pdf(x, a, scale=b), 'k', linestyle=ls)

# ax_right.set_xlim(0.3,1.0)
# ax_right.set_ylim(0,0.5)

# ax_right.set_xticks([0.3, 0.5, 0.8, 1.0])
# ax_right.set_yticks([0., 0.1, 0.2, 0.3,0.4, 0.5])

# ax_right.yaxis.tick_right()
pl.legend(prop={'size':12})

pl.savefig('fig3A.pdf', dpi=600, bbox_inches='tight')
