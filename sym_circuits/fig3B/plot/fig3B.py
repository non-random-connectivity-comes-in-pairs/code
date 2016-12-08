
import matplotlib as mpl
mpl.use('Agg')
import pylab as pl

import numpy as np
from scipy.stats import sem
import pickle

from matplotlib import rc

rc('text', usetex=True)
pl.rcParams['text.latex.preamble'] = [
    r'\usepackage{tgheros}',    # helvetica font
    r'\usepackage{sansmath}',   # math-font matching helvetica
    r'\sansmath'                # actually tell tex to use it!
    r'\usepackage{siunitx}',    # micro symbols
    r'\sisetup{detect-all}',    # force siunitx to use the fonts
]  


# -----  data from numerical integration -----

# a=1
with open("data/th_a1_mu01.p", "rb") as pfile:
    df_a1 = pickle.load(pfile)

# data for sigma=0 from Figure 2    
t_sigs_a1 = [0]
t_rhos_a1 = [1.996]

t_sigs_a1 += [d["sigma"] for d in df_a1]
t_rhos_a1 += [d["numer"]/(1./4*(d["den_x"]+d["den_y"])**2) for d in df_a1]

# a=2
with open("data/th_a2_mu01.p", "rb") as pfile:
    df_a2 = pickle.load(pfile)
    
# data for sigma=0 from Figure 
t_sigs_a2 = [0]
t_rhos_a2 = [1.500]

t_sigs_a2 = [d["sigma"] for d in df_a2]
t_rhos_a2 = [d["numer"]/(1./4*(d["den_x"]+d["den_y"])**2) for d in df_a2]


# -----  data from sampling connection probabilities -----

with open("data/4_alpha_0-05_sigma.p", "rb") as pfile:
    df = pickle.load(pfile)

alph = 0.248
x = [d for d in df if d["alpha"]==alph]

sample_sigs = list(set([d["sigma"] for d in x]))
rhos_a0248 = []
rhos_a0248_sem = []

for sigma in sample_sigs:    
    df_sig = [d for d in x if d["sigma"]==sigma]
    rhos_sig = [np.mean(d["xs"]*d["ys"])\
                * 1/(np.mean(np.concatenate((d["xs"],d["ys"])))**2)\
                for d in df_sig]
    rhos_a0248.append(np.mean(rhos_sig))
    rhos_a0248_sem.append(sem(rhos_sig))


alph = 1.
x = [d for d in df if d["alpha"]==alph]

sample_sigs_a1 = list(set([d["sigma"] for d in x]))
rhos_a1 = []
rhos_a1_sem = []

for sigma in sample_sigs_a1:
    df_sig = [d for d in x if d["sigma"]==sigma]
    rhos_sig = [np.mean(d["xs"]*d["ys"])\
                * 1/(np.mean(np.concatenate((d["xs"],d["ys"])))**2)\
                for d in df_sig]
    rhos_a1.append(np.mean(rhos_sig))
    rhos_a1_sem.append(sem(rhos_sig))

    
# ----- fit for alpha=0.248 sampled data -----

xs = np.arange(0.,0.7,0.001)
ys = 1.086317 + (4.043159 - 1.086317)/(1 + (xs/0.2587529)**3.275628)


# ----- data from generated networks -----

with open("data/gn_a1_sig05_1.p", "rb") as pfile:
    df = pickle.load(pfile, encoding='latin1')

gn_sample_sigs = list(set([d["sig"] for d in df]))

gn_a1_rhos = []
gn_a1_rho_sems = []
for sig in gn_sample_sigs:
    df_sig  = [d for d in df if d["sig"]==sig]
    gn_a1_rhos.append(np.mean([d["rho"] for d in df_sig]))
    gn_a1_rho_sems.append(sem([d["rho"] for d in df_sig]))


    

fig, ax = pl.subplots(1,1)
fig.set_size_inches(7.5*0.5,2.3)

pl.plot(xs,ys, 'k', label=r'$\alpha=0.248$')
(_, caps, _) = pl.errorbar(sample_sigs, rhos_a0248, yerr=rhos_a0248_sem,
                           fmt=None, ecolor='k', elinewidth=1.5,)
for cap in caps:
    cap.set_markeredgewidth(0.8)

pl.plot(t_sigs_a1, t_rhos_a1, 'k', linestyle=':', label=r'$\alpha=1$')
pl.plot(t_sigs_a2, t_rhos_a2, 'k', linestyle='--', label=r'$\alpha=2$')

(_, caps, _) = pl.errorbar(gn_sample_sigs, gn_a1_rhos, yerr=gn_a1_rho_sems,
                           fmt=None, ecolor = 'r', elinewidth=1.5,)
for cap in caps:
    cap.set_markeredgewidth(0.8)

pl.xticks([0.,0.1,0.2,0.3,0.4,0.5,0.6])

ax.set_title(r'$f_{P_{ij}}(y) = f_{\alpha,\beta}^T(y)$', size=13.)
pl.xlim(0,0.6)
pl.ylim(1.,4.05)

pl.legend(prop={'size':12})

pl.xlabel(r'width $\sigma$')
pl.ylabel(r'relative occurrence $\varrho$')

pl.savefig('fig3B.pdf', dpi=600, bbox_inches='tight')
