# Analyting pass times

# %%
import os

os.chdir('/Users/Theo/GRWMODELS/python/mf6lab/Projects/Dispersion/cases/wimsaquif')

import numpy as np
import matplotlib.pyplot as plt
from settings import dirs, props
from scipy.interpolate import UnivariateSpline
from scipy import stats
from itertools import cycle

os.chdir(dirs.MP7)

def round_up_nice(x, m=1):
    return 10 ** (np.floor(np.log10(x)) - m) * np.ceil(
        np.round(10 ** (np.log10(x) + m - np.floor(np.log10(x))), 2))

# %%
k_field_str = props['k_field_pars']['k_field_str']

name = props['k_field_pars']['name']
pass_times = np.load('pass_times_' + name + '.npy')

cumulatives = dict()

for j, xo in enumerate(pass_times['xObs'][0]):
    idx = np.argsort(pass_times['time'][:, j])    
    t = pass_times['time'][:, j][idx]
    cumulatives[int(xo)] = t

last_times = cumulatives[list(cumulatives.keys())[-1]]
t_last = last_times[int(0.995 * len(last_times))]

xlim = (0, round_up_nice(t_last, m=2))

# dtype([('id', '<i4'), ('xObs', '<f8', (10,)), ('time', '<f8', (10,))])

# %%

fig, ax = plt.subplots(figsize=(10, 6))
ax.set(title='Breakthrough curves from Modpath with the folowing k-field parameters:\n' + k_field_str, xlabel='time [d]', ylabel='p[t < t[p]]')
p = pass_times['id'] / 1000.
for k in cumulatives:
    ax.plot(cumulatives[k], p, label=f'x = {k} m')
    
ax.set_xlim(-1000, 150000)
ax.grid()
ax.legend()

# %% == show log normal distribution

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, figsize=(10, 11))
ax1.set(title='cdfs with the following k-field parameters\n'
        + k_field_str, xlabel='t', ylabel='p[t < t(p)]') 
ax2.set(title='pdfs', xlabel='t', ylabel='d(cdf)/dt')

ax1.set_xlim(xlim)
ax2.set_xlim(xlim)

dtype=np.dtype([('x', '<i4'), ('mu', '<f8'), ('sigma', '<f8'), ('loc', '<f8'),
                ('mode_t', '<f8'), ('median', '<f8'), ('mean_t', '<f8'),
                ('vmode', '<f8'), ('vmean', '<f8')])
dtype2=np.dtype([('x', '<i4'), ('t05', '<f8'), ('t16', '<f8'), ('t50', '<f8'),
                ('t84', '<f8'), ('t95', '<f8')])

musigmaloc = np.zeros(len(cumulatives), dtype=dtype)
bt_times = np.zeros(len(cumulatives), dtype=dtype2) # Breakthrough times
clrs = cycle(["blue", "red", "green", "black", "gray", "magenta", "cyan", "orange", "brown", "darkgray"])

print("{:5s} {:10s} {:10s} {:10s}".format('x', 'mu', 'sigma', 'scale'))

for ix, x in enumerate(cumulatives):
    if x not in [200, 500, 1000]:
        continue
    clr = next(clrs)
    t = cumulatives[x]
    cdf_o = p

    params_lognorm = stats.lognorm.fit(t)    
    cdf_lognorm    = stats.lognorm.cdf(t, *params_lognorm)
    pdf_lognorm    = stats.lognorm.pdf(t, *params_lognorm)

    shape, loc, scale = params_lognorm

    dist = stats.lognorm(s=shape, loc=loc, scale=scale)
    t05, t16, t50, t84, t95 = dist.ppf(0.05), dist.ppf(0.1587), dist.ppf(0.5),    dist.ppf(0.8413), dist.ppf(0.95)
    
    bt_times[ix] = np.array([(x,
                              np.round(t05),
                              np.round(t16),
                              np.round(t50),
                              np.round(t84),
                              np.round(t95))], dtype=dtype2)

    mu = np.log(scale)
    sigma = shape 
    t0=loc
    mode_t = t0 + np.exp(mu - sigma ** 2)
    median_t = loc + np.exp(mu)
    mean_t = t0 + np.exp(mu + 0.5 * sigma ** 2)
    vmean = x / mean_t
    vmode = x / mode_t
    musigmaloc[ix] = np.array([(x, mu, sigma,
                                np.round(t0),
                                np.round(mode_t),
                                np.round(median_t),
                                np.round(mean_t),
                                np.round(vmode, 4),
                                np.round(vmean, 4)
                                )], dtype=dtype)
    
    print("{:5d} {:10.3g} {:10.3g} {:10.3g}".format(x, mu, sigma, loc))
    
    pdf_mode = stats.lognorm.pdf(mode_t, *params_lognorm)

    
    ax1.plot(t, cdf_o, '--', color=clr, label=f"cdf {x} m (original data)")
    ax1.plot(t,
        cdf_lognorm, color=clr, label=fr"x={x:4d} m, $\mu,\sigma,\,t_0$,loc=({mu:.2f}, {sigma:.2f}, {loc:.3e} d)")
    ax2.plot(t, pdf_lognorm, color=clr, label=fr"x={x:4d} m, $\mu,\sigma,\,t_0$=({mu:.2f}, {sigma:.2f}, {loc:.3e} d)")

    alpha = np.array([0.2, 0.5, 1.0, 2.0, 5.0])
    ta = t0 + (mode_t - t0) * alpha
    
    ax1.plot(ta, stats.lognorm.cdf(ta, *params_lognorm), 'o', mfc=clr, mec='k')
    ax2.plot(ta, stats.lognorm.pdf(ta, *params_lognorm), 'o', mfc=clr, mec='k')

    ax1.plot(median_t, stats.lognorm.cdf(median_t, *params_lognorm), 's', ms=7, mfc=clr, mec='k')
    ax2.plot(median_t, stats.lognorm.pdf(median_t, *params_lognorm), 's', ms=7, mfc=clr, mec='k')
    
    ax1.plot(mean_t, stats.lognorm.cdf(mean_t, *params_lognorm), '*', ms=10, mfc=clr, mec='k')
    ax2.plot(mean_t, stats.lognorm.pdf(mean_t, *params_lognorm), '*', ms=10, mfc=clr, mec='k')

    # ax1.vlines([t05, t16, t50, t84, t95], 0, ymax=0.95, colors=clr)
    # ax2.vlines([t05, t16, t50, t84, t95], 0, ymax=pdf_mode, colors=clr)    
    
    # ax1.vlines(ta, 0, ymax=0.95, colors=clr)
    # ax2.vlines(ta, 0, ymax=pdf_mode, colors=clr)
    
    ax1.text(0.6, 0.6, 
         r"$\blacksquare = t_{median},\ \bigstar = t_{mean}$", 
         transform=ax1.transAxes, 
         fontsize=14, 
         color='green',
         ha='left', va='top')
    ax2.text(0.6, 0.6,
         r"$\blacksquare = t_{median},\ \bigstar = t_{mean}$", 
         transform=ax2.transAxes, 
         fontsize=14, 
         color='green',
         ha='left', va='top')

for ax in [ax1, ax2]:
    ax.grid(True)
    ax.legend()

# %% == show normal distribution
ttl1 = 'square = medium, * = mean'
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, figsize=(10, 11))
ax1.set(title='ln(cdf)' + ttl1, xlabel='ln(t)', ylabel='p[ln(t) < ln(t(p))]') 
ax2.set(title='ln(pdf)' + ttl1, xlabel='ln(t)', ylabel='d(ln(cd)f)/dln(t)') 

print("{:5s} {:10s} {:10s}".format('x', 'mu', 'sigma'))
for x in cumulatives:
    t = cumulatives[x]
    cdf_o = p

    params_norm = stats.norm.fit(np.log(t))    
    cdf_norm    = stats.norm.cdf(np.log(t), *params_norm)
    pdf_norm    = stats.norm.pdf(np.log(t), *params_norm)

    mu, sigma = params_norm
    print("{:5d} {:10.3g} {:10.3g}".format(x, mu, sigma))
        
    ax1.plot(np.log(t), cdf_norm, label=fr"x={x} m, $\mu,\sigma$={mu:.2f}, {sigma:.2f}")
    ax2.plot(np.log(t), pdf_norm, label=fr"x={x} m, $\mu,\sigma$={mu:.2f},{sigma:.2f}")

for ax in [ax1, ax2]:
    ax.grid(True)
    ax.legend()


# %%
