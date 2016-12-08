
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


# with open("../data/data_f_mult_norm_1.p", "rb") as pfile:
#     data = pickle.load(pfile, encoding='latin1')

# sigmas, rho_means, rho_sems = [], [], []

# for sigma in data.keys():
#     sigmas.append(sigma)
#     rho_means.append(np.mean(data[sigma]))
#     rho_sems.append(sem(data[sigma]))


# -----  data from numerical integration -----

# a=1
with open("../data/th_a1_mu01.p", "rb") as pfile:
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

t_sigs_a1 = [d["sigma"] for d in df_a1]
t_rhos_a1 = [d["numer"]/(1./4*(d["den_x"]+d["den_y"])**2) for d in df_a1]


# -----  data from sampling connection probabilities -----

with open("../data/4_alpha_0-05_sigma.p", "rb") as pfile:
    df = pickle.load(pfile)

alph = 0.248

x = [d for d in df if d["alpha"]==alph]

sample_sigs = list(set([d["sigma"] for d in x]))
rhos = []
rhos_sem = []
mus = []
xavgs = []
yavgs = []

for sigma in sample_sigs:
    
    df_sig = [d for d in x if d["sigma"]==sigma]
    #assert(len(df_sig) == 5)
    rhos_sig = [np.mean(d["xs"]*d["ys"])/(np.mean(np.concatenate((d["xs"],d["ys"])))**2) for d in df_sig]
    rhos.append(np.mean(rhos_sig))
    rhos_sem.append(sem(rhos_sig))
    
    mus_sig = [np.mean(np.concatenate((d["xs"],d["ys"]))) for d in df_sig]
    mus.append(np.mean(mus_sig))
    
    xavg_sig = [np.mean(d["xs"]) for d in df_sig]
    xavgs.append(np.mean(xavg_sig))

    yavg_sig = [np.mean(d["ys"]) for d in df_sig]
    yavgs.append(np.mean(yavg_sig))


alph = 1.

x = [d for d in df if d["alpha"]==alph]

sample_sigs_a1 = list(set([d["sigma"] for d in x]))
rhos_a1 = []
rhos_a1_sem = []
mus = []
xavgs = []
yavgs = []

for sigma in sample_sigs_a1:
    
    df_sig = [d for d in x if d["sigma"]==sigma]
    #assert(len(df_sig) == 5)
    rhos_sig = [np.mean(d["xs"]*d["ys"])/(np.mean(np.concatenate((d["xs"],d["ys"])))**2) for d in df_sig]
    rhos_a1.append(np.mean(rhos_sig))
    rhos_a1_sem.append(sem(rhos_sig))


fig, ax = pl.subplots(1,1)
fig.set_size_inches(7.5*0.5,2.3)

xs = np.arange(0.,0.7,0.001)
ys = 1.086317 + (4.043159 - 1.086317)/(1 + (xs/0.2587529)**3.275628)

pl.plot(xs,ys, 'k', label=r'$\alpha=0.248$')
(_, caps, _) = pl.errorbar(sample_sigs, rhos, yerr=rhos_sem, fmt=None, ecolor='k', elinewidth=1.5,)
for cap in caps:
    cap.set_markeredgewidth(0.8)
# for i in range(len(rhos)):
#     print(sample_sigs[i], rhos[i])
pl.plot(t_sigmas_a1, t_rhos_a1, 'k', linestyle=':', label=r'$\alpha=1$')
pl.plot(t_sigmas_a2, t_rhos_a2, 'k', linestyle='--', label=r'$\alpha=2$')

pl.xticks([0.,0.1,0.2,0.3,0.4,0.5,0.6])

ax.set_title(r'$f_{P_{ij}}(y) = f_{\alpha,\beta}^T(y)$', size=13.)
pl.xlim(0,0.6)
pl.ylim(1.,4.05)

pl.legend(prop={'size':12})

pl.xlabel(r'width $\sigma$')
pl.ylabel(r'relative occurrence $\varrho$')

pl.savefig('fig3B.pdf', dpi=600, bbox_inches='tight')
