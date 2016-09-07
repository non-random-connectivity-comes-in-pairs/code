

import pylab as pl
import numpy as np

from scipy.stats import gamma

import matplotlib.ticker as ticker
from matplotlib import rc

rc('text', usetex=True)
pl.rcParams['text.latex.preamble'] = [
    r'\usepackage{tgheros}',    # helvetica font
    r'\usepackage{sansmath}',   # math-font matching  helvetica
    r'\sansmath'                # actually tell tex to use it!
    r'\usepackage{siunitx}',    # micro symbols
    r'\sisetup{detect-all}',    # force siunitx to use the fonts
]  


a_ns = np.load("data/a_ns.npy")
rhos = np.load("data/rho_theoretical.npy")

fig, ax = pl.subplots(1,1)
fig.set_size_inches(7.2*0.5,2.2)

ax.plot(a_ns,  rhos , color='k', label=r'$\varrho$')
ax.plot(a_ns, [1+1./a for a in a_ns], color='k',
        linestyle='dashed', label=r'$1+\frac{1}{\alpha}$')

ax.set_title(r'$\mu = 0.1$')

ax.set_xscale('log')
ax.set_xlabel(r'shape parameter $\alpha$')
ax.set_ylabel(r'relative occurrence $\varrho$')

pl.yticks(list(pl.yticks()[0]),
          ['%.1f' % tick for tick in pl.yticks()[0]])

pl.xticks(sorted(list(pl.xticks()[0]) + [0.2]),
          sorted(['0.2']+[str(int(x)) for x in list(pl.xticks()[0])]))

ax.set_xlim(0.2,100)

ax.legend()


pl.savefig('gamma_figC.pdf', dpi=600, bbox_inches='tight')


