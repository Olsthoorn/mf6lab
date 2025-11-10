# Further analysis
# De vraag is dus, wat moet je doen in de verdere analyse?
# het blijkt uit de simulatie dat de mu van groot belang is om de curve te reproduceren.
# De sigma lijkt bijna niets te doen. Dat kunnen we ook simuleren.
# Vraag is vervolgens hoe dat zit met de overige simulaties.
# Kan je de resultaten van een simulatie op elkaar laten vallen door delen door de mean t.
# In elk geval bleek de mean t strikt lineair.
# Je kan ook de verblijftijd omrekenen naar een gemiddelde snelheid per deeltje.
# Dat betekent dat je door de plaats van de deeltjes op elk moment kan bepalen.
# Kan je analytisch de doorslagtijd omrekenen naar snelheid en plaats?
# Deeltjes blijven ook na omkering gesorteerd.
# De pluim is de inverse van de doorslagkromme. Beide bevatten dezelfde informatie.



# %% 
import os
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from itertools import cycle
from scipy.stats import lognorm
from scipy.optimize import curve_fit

from src.mf6tools import  Dirs


def get_Ihs(fnames):
    """Return de horizontal integration lengths burried in the file names.
    
    The result will be a sorted list of tuples [(Ih, basename), (.., ..), ...]
    
    """
    Ihs = []
    for fname in fnames:
        Ihs.append((get_Ih(fname), os.path.basename(fname)))
    return sorted(Ihs)

def get_Ih(fname):
    """Return the hoirizontal integration lenght buffied in the file name.
    
    In "fname_exp24.ext" the exp is the case name and the 24 the Ih in meters.
    """
    return int(fname[fname.find('exp')+3:fname.find('.')])


LENGTH_UNITS = 'meter'
TIME_UNITS = 'days'

# %% get the directories straight

HOME = '/Users/Theo/GRWMODELS/python/mf6lab/Projects/Dispersion/'
assert os.path.isdir(HOME), "Can't find the directory {}".format(HOME)

dirs = Dirs(HOME)

# Get casename and add the directories for working with this case
sim_name = 'wim2'
section_name = sim_name

dirs = dirs.add_case(sim_name)
os.chdir(dirs.case)

# %% Collect all data in dictionaries

# The data come from the Modpath mfpath file analysis and are all in separate nd.arrays
# with their dtypes. Look at the dtypes to see what data they contain.

# Get the file names in tuples (k, fname)
btt_tuples = get_Ihs([os.path.basename(fname) for fname in glob('./data/BTarrays_t*.npz')])
btx_tuples = get_Ihs([os.path.basename(fname) for fname in glob('./data/BTarrays_x*.npz')])
pst_tuples = get_Ihs([os.path.basename(fname) for fname in glob('./data/pass_times_*.npy')])
psx_tuples = get_Ihs([os.path.basename(fname) for fname in glob('./data/pass_xvals_*.npy')])
psi_tuples = get_Ihs([os.path.basename(fname) for fname in glob('./data/psi_*.npy')])

# Breakthrough times data
btt = dict()
for k, fname in btt_tuples:
    dta = np.load(os.path.join('./data/', fname))
    btt[k] = {'props': dta['bt_cdf_props_t'], 'times':dta['bt_times'], 'mode':dta['bt_mode']}

# x location break through data
btx = dict()
for k, fname in btx_tuples:
    dta = np.load(os.path.join('./data/', fname))
    btx[k] = {'props': dta['bt_cdf__xprops'], 'times':dta['bt_times'], 'mode':dta['bt_xmode'],
              'lnt':dta['lognorm_t_pars'], 'lnx':dta['lognorm_x_pars']}

# Pass times data
pst = dict()
for k, fname in pst_tuples:
    pst[k] = np.load(os.path.join('./data/', fname))

# Pass xvals data
psx = dict()
for k, fname in psx_tuples:
    psx[k] = np.load(os.path.join('./data/', fname))

# Psi data
psi = dict()
for k, fname in psi_tuples:
    psi[k] = np.load(os.path.join('./data/', fname))
   
# %% 
# Show the development of the properties of the breakthrough curves for each case
clrs = cycle("brgkmc")
for k in btt.keys():
    clr = next(clrs)
    p = btt[k]['props']
    a, b, c = np.polyfit(p['x'], p['mu'], 2)
    plt.scatter(p['x'], p['mu'],     color=clr, label=f'data Ih = {k}')
    # plt.plot(p['x'], a * p['x'] ** 2 + b * p['x'] + c, color=clr, label=f'Fit: y={a:.2f}x^2 + {b:.2f}x + {c:.2f}')

plt.legend()
plt.show()

# %% 
# Show the development of the properties of the breakthrough curves for each case
clrs = cycle("brgkmc")
x, y = [], []
for k in btt.keys():
    p = btt[k]['props']
    x += list(p['x'])
    y += list(p['mu'])

plt.scatter(x, y,     color=clr, label=f'data Ih = {k}')

a, b, c = np.polyfit(x, y, 2)

xp = np.logspace(-3, np.log10(2000), 100)

plt.plot(xp, a * xp ** 2 + b * xp + c , color=clr, label=f'Fit: y={a:.2f}x^2 + {b:.2f}x + {c:.2f}')

plt.legend()
plt.show()

# %% fit power function

def model_pow(x, a, b):
    return a * x  ** b
def model_sqrt(x, n, a, b):
    return n + a * np.sqrt(b * x)
def model_lim(x, n, a, b):
    return n + b * (1 - np.exp(-x / a))

x, y = [], []
for k in btt.keys():
    p = btt[k]['props']
    x += list(p['x'])
    y += list(p['mu'])
    plt.plot(p['x'], p['mu'], '.', label=f"Ih={k} m")
    
params, _ = curve_fit(model_lim, x, y, p0=(6, 5, 500.))
n, a, b = params

#params, _ = curve_fit(model_pow, x, y, p0=(1.5, 0.5))
#a, b = params

# params, _ = curve_fit(model_sqrt, x, y, p0=(6., 1.5, 0.5))
# n, a, b = params

# plt.scatter(x, y, label="x, mu")

xp = np.logspace(-3, np.log10(2000.), 100)
plt.plot(xp, model_lim(xp, n, a, b), label=fr"$y = {n:.2f} + {b:.2f} \left(1 - e^{{-\frac{{x}}{{{a:.2f}}}}}\right)$")

# plt.plot(xp, model_pow(xp, a, b), label=f"y = {a:.2f} * x * {b:.2f}")

# plt.plot(p['x'], model_sqrt(p['x'], n, a, b), label=f"y = {n:.2f} + {a:.2f} np.sqrt({b:.2f} * x)")

plt.legend()
plt.show()

# %% Sigma

# %% fit sigma

x, y = [], []
for k in btt.keys():
    p = btt[k]['props']
    x += list(p['x'])
    y += list(p['sigma'])
    plt.plot(p['x'], p['sigma'], '.-', label=f"Ih={k} m")
plt.legend()
plt.show()

# %% fit loc

x, y = [], []
for k in btt.keys():
    p = btt[k]['props']
    x += list(p['x'])
    y += list(p['loc'])
    plt.plot(p['x'], p['loc'], '.-', label=f"Ih={k} m")
plt.legend()
plt.show()

# %% fit t_mean

x, y = [], []
for k in btt.keys():
    p = btt[k]['props']
    x += list(p['t_mean'])
    y += list(p['t_mean'])
    plt.plot(p['x'], p['t_mean'], '.-', label=f"Ih={k} m")
plt.legend()
plt.show()

# %% fit t_med

x, y = [], []
for k in btt.keys():
    p = btt[k]['props']
    x += list(p['t_med'])
    y += list(p['t_med'])
    plt.plot(p['x'], p['t_med'], '.-', label=f"Ih={k} m")
plt.legend()
plt.show()

# %% fit t_mode

x, y = [], []
for k in btt.keys():
    p = btt[k]['props']
    x += list(p['t_mode'])
    y += list(p['t_mode'])
    plt.plot(p['x'], p['t_mode'], '.-', label=f"Ih={k} m")
plt.legend()
plt.show()

# %% Normal and lognormal pdfs

def norm(x, mu, sigma):
    return 1  /(sigma * np.sqrt(2 * np.pi)) * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))
def lnorm(x, mu, sigma):
    return 1  /(x * sigma * np.sqrt(2 * np.pi)) * np.exp(-(np.log(x) - mu) ** 2 / (2 * sigma ** 2))

x = np.linspace(0, 2000., 1000)

for k in btt.keys():
    p = btt[k]['props']
    mu = p['mu'][-1]
    sigma = p['sigma'][-1]
    plt.plot(x, lnorm(x, mu, sigma), '.-', label=f"Ih={k} m")
plt.legend()
plt.show()


# %% =========S P R E A D ==================================

xObs = pst[1]['xObs'][0]
std = np.zeros((len(pst), len(xObs)))
mu  = np.zeros_like(std)
med = np.zeros_like(std)
for ik, k in enumerate(pst.keys()):    
    for io, xo in enumerate(xObs):
        std[ik, io] = np.std(pst[k]['time'][:, io])
        mu [ik, io] = np.mean(pst[k]['time'][:, io])
        med[ik, io] =np.median(pst[k]['time'][:, io])
        

print("std")
print(np.round(std / 1000., 0))
print("mu")
print(np.round(mu / 1000., 0))        
print("med")
print(np.round(med / 1000., 0))

# %% The spread of the breakthough time data

fig, ax = plt.subplots(figsize=(10, 6))
ax.set(title="std", xlabel='xObs', ylabel='std')
for k, std_ in zip(pst.keys(), std):
    ax.plot(xObs, std_, label=f"Ih = {k} m")
ax.grid()
ax.legend()
plt.show()      


# %% Plot the lognormal distribution
mu, sigma = 10., 0.5
    
fig, ax = plt.subplots(figsize=(10, 6))
ax.set(title=f"cdfs computed, for mu={mu}, sigma={sigma}", xlabel='t [d]',  ylabel='pdf')
    
mus = np.array([8.63, 10.13, 10.39, 10.47, 10.78, 11.07, 11.02, 10.98, 10.94, 11.17])
sigmas = np.array([0.8, 0.44, 0.47, 0.56, 0.53, 0.45, 0.50, 0.57, 0.63, 0.58])
t0s = np.array([5000., 9800., 26130., 46440. , 57500., 67000., 93000., 119100., 144100.])
t = np.linspace(0, 300000., 2000)

clrs = cycle("brgkmc")
for mu, sigma, t0 in zip(mus,  sigmas, t0s):
    clr = next(clrs)
    ax.plot(t[t > t0], lognorm.cdf(t[t > t0], sigma, loc=t0, scale=np.exp(mu)),
            color=clr, label=f"mu={mu}, sigma={sigma}, t0={t0}d")
    mmu = np.mean(mus)
    ax.plot(t[t > t0], lognorm.cdf(t[t > t0], sigma, loc=t0, scale=np.exp(mmu)), '--',
            color=clr, label=f"mu={mmu}, sigma={sigma}, t0={t0}d")
    msig = np.mean(sigmas)
    ax.plot(t[t > t0], lognorm.cdf(t[t > t0], msig, loc=t0, scale=np.exp(mu)), ':',
            color=clr, label=f"mu={mu}, sigma={msig}, t0={t0}d")
    
ax.grid()
# ax.legend()

# %%
files = glob(os.path.join(dirs.data, '*t_exp*.npz'))
files
fff = []
for file in files:
    basename = os.path.basename(file)
    Ih = int(basename[basename.find("exp")+3:basename.find('.')])
    fff.append((Ih, basename))
fff = sorted(fff)
    
# %%
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 12))
ax1.set(title='mu sigma', xlabel="", ylabel='mu, sigma (log)')
ax2.set(title='times', xlabel=['x [m]'], ylabel='times [d]')

for Ih, fname in fff:
    print(f"{Ih:>10d}   {fname}")
    D = np.load('./data/' + fname)
    tpr = D['bt_cdf_props_t']
    ax1.plot(tpr['x'], tpr['mu'], '-', label=f'mu_{Ih}')
    ax1.plot(tpr['x'], tpr['sigma'], '--', label=f'sigma_{Ih}')
    ax1.plot(tpr['x'], np.log(tpr['t_mean'] - tpr['loc']),  ':', label=f'log(mean_{Ih} - loc_{Ih}')
    
    ax2.plot(tpr['x'], tpr['loc'] / tpr['t_mean'] , label=f'loc_{Ih}')
    #ax2.plot(tpr['x'], tpr['t_mode'] / tpr['t_mean'], label=f'mode_{Ih}')
    #ax2.plot(tpr['x'], tpr['t_med'] / tpr['t_mean'], label=f'med_{Ih}')
    #ax2.plot(tpr['x'], tpr['t_mean'], label=f'mean_{Ih}')
    
    
ax1.grid()
ax2.grid()
ax1.legend()
ax2.legend()

plt.show()
# %%

# fname = os.path.join(dirs.data, 'psi_exp*.npy')
# psifiles = glob(fname)
# Ihs = get_Ihs(psifiles)

# psi = []
# for (Ih, fname) in Ihs:
#     psi_ = np.max(np.load(os.path.join(dirs.data, fname)))
#     print(f"{Ih:>5d}, {fname:>15s}, {np.max(psi_):.4g}")
#     psi.append(psi_)
# psi = np.array(psi)
# psi[0] / psi

# %%
