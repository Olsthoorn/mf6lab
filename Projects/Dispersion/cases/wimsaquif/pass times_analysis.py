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

# dtype of the fitted lognormal cdf data and extracted key values
dtype1t=np.dtype([('x', '<i4'), ('mu', '<f8'), ('sigma', '<f8'), ('loc', '<f8'),
                ('t_mode', '<f8'), ('t_med', '<f8'), ('t_mean', '<f8'),
                ('v_mode', '<f8'), ('v_mean', '<f8')])
dtype1x=np.dtype([('t', '<i4'), ('mu', '<f8'), ('sigma', '<f8'), ('loc', '<f8'),
                ('x_mode', '<f8'), ('x_med', '<f8'), ('x_mean', '<f8'),
                ('v_mode', '<f8'), ('v_mean', '<f8')])


# dtype for the times the cdf reaches specific quantiles (eg. t95 = t of 95%  quantitle)
dtype2t=np.dtype([('x', '<i4'), ('t05', '<f8'), ('t16', '<f8'), ('t50', '<f8'),
                ('t84', '<f8'), ('t95', '<f8')])
dtype2x=np.dtype([('t', '<i4'), ('x05', '<f8'), ('x16', '<f8'), ('x50', '<f8'),
                ('x84', '<f8'), ('x95', '<f8')])

# dtype of the characteristics of the pdf's mode. tm50 for instance is the
# time at 50% between t0=loc and t_mode. cm50 is its cdf value at tm50 and
# pm50 is the pdf value at tm50.
dtype3t=np.dtype([('x', '<i4'),
        ('tm02', '<f8'), ('tm05', '<f8'), ('tm10', '<f8'), ('tm20', '<f8'), ('tm50', '<f8'),
        ('cm02', '<f8'), ('cm05', '<f8'), ('cm10', '<f8'), ('cm20', '<f8'), ('cm50', '<f8'),
        ('pm02', '<f8'), ('pm05', '<f8'), ('pm10', '<f8'), ('pm20', '<f8'), ('pm50', '<f8')])

dtype3x=np.dtype([('t', '<i4'),
        ('xm02', '<f8'), ('xm05', '<f8'), ('xm10', '<f8'), ('xm20', '<f8'), ('xm50', '<f8'),
        ('cm02', '<f8'), ('cm05', '<f8'), ('cm10', '<f8'), ('cm20', '<f8'), ('cm50', '<f8'),
        ('pm02', '<f8'), ('pm05', '<f8'), ('pm10', '<f8'), ('pm20', '<f8'), ('pm50', '<f8')])


# Often needed to mark figure and outputfiles
k_field_str = props['k_field_pars']['k_field_str']
case_name = props['k_field_pars']['name']

# Output data will be stored in dirs.data and are marked with the case name.
# Saved figure are stored in dirs.images and are marked with the case name.

# ===================== Let's go =========================================

# %% === Get the dat for times at different x-locations and show them ====

# pass_times dtype([('id', '<i4'), ('xObs', '<f8', (10,)), ('time', '<f8', (10,))])
# pass_times has the times of each particle (id) for a set of x-locations
pass_times = np.load(os.path.join(dirs.data, 'pass_times_' + case_name + '.npy'))

nxObs = len(pass_times['xObs'][0])

# Get just the x and times so we can sort the times (eliminate the ids)
dtype0t = np.dtype([('xObs', '<i4', (nxObs,)), ('time', '<f8', (nxObs,))])
cumul_t = np.zeros(len(pass_times), dtype=dtype0t)  # dict()
cumul_t['xObs'] = pass_times['xObs']

# Sort the times
for j in range(nxObs):
    idx = np.argsort(pass_times['time'][:, j])    
    t = pass_times['time'][:, j][idx]    
    cumul_t['time'][:, j] = t
    
# Shorten the xlim for better picture
last_times = cumul_t['time'][:, -1]
t_last = last_times[int(0.995 * len(last_times))]

xlim = (0, round_up_nice(t_last, m=2))


# %% === Show the raw BTC's breakthrough (p vs sorted times one curve per xObs)

# Raw break th
fig, ax = plt.subplots(figsize=(10, 6))
ax.set(title='BTCs from Modpath (case ' + case_name + ') with the folowing k-field parameters:\n' + k_field_str, xlabel='time [d]', ylabel='p[t < t[p]]')
p = pass_times['id'] / 1000.
for X, t in zip(cumul_t['xObs'].T, cumul_t['time'].T):
    ax.plot(t, p, label=f'x = {X[0]} m')
    
ax.set_xlim(-1000, 150000)
ax.grid()
ax.legend()

fig. savefig(os.path.join(dirs.images, 'BTcurves_' + case_name + '.png'))


# %% == Fit data to lognormal cdf's analyze and show them =========

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, figsize=(10, 11))
ax1.set(title='cdfs (case ' + case_name + ') with the following k-field parameters\n'
        + k_field_str, xlabel='t', ylabel='p[t < t(p)]') 
ax2.set(title='pdfs', xlabel='t', ylabel='d(cdf)/dt')

ax1.set_xlim(xlim)
ax2.set_xlim(xlim)


bt_cdf_props_t = np.zeros(cumul_t['xObs'].shape[-1], dtype=dtype1t)
bt_times       = np.zeros(cumul_t['xObs'].shape[-1], dtype=dtype2t) # Breakthrough times
bt_tmode      = np.zeros(cumul_t['xObs'].shape[-1], dtype=dtype3t) # Breakthrough times
clrs = cycle(["blue", "red", "green", "black", "gray", "magenta", "cyan", "orange", "brown", "darkgray"])

print("{:5s} {:10s} {:10s} {:10s}".format('x', 'mu', 'sigma', 'scale'))

for ix, (X, t) in enumerate(zip(cumul_t['xObs'].T, cumul_t['time'].T)):
    x = X[0].item()
    clr = next(clrs)    
    cdf_o = p

    par_log_t = stats.lognorm.fit(t)    
    cdf_log_t = stats.lognorm.cdf(t, *par_log_t)
    pdf_log_t = stats.lognorm.pdf(t, *par_log_t)

    shape_t, loc_t, scale_t = par_log_t
    mu_t, sig_t, t0 = np.log(scale_t), shape_t, loc_t

    print("x = {:5d}, mu = {:10.3g}, sigma = {:10.3g} loc = t0 = {:10.3g}".format(x, mu_t, sig_t, t0))

    # Get lognormal distribution properties.
    t_mode = t0 + np.exp(mu_t - sig_t ** 2)
    t_med  = t0 + np.exp(mu_t)
    t_mean = t0 + np.exp(mu_t + 0.5 * sig_t ** 2)
    v_mean = x / t_mean
    v_mode = x / t_mode
    
    # Store them in bt_cdf_props_t recarray.
    bt_cdf_props_t[ix] = np.array([(x, mu_t, sig_t, t0, t_mode, t_med, t_mean, v_mode, v_mean)], dtype=dtype1t)

    # Inversely get the times at which the cdf has given p values
    
    # First get the lognormal distribution with the fitted parameters
    dist = stats.lognorm(s=shape_t, loc=loc_t, scale=scale_t)
    
    # Use this distibution to get the times at which the cdf has these values
    # (useful for normal distribution less so for lognormal distribution)
    t05, t16, t50, t84, t95 = (dist.ppf(0.05), dist.ppf(0.1587), dist.ppf(0.5),
                               dist.ppf(0.8413), dist.ppf(0.95))
    
    # Store them, and for convenience round them in a recarray with dtype2t
    bt_times[ix] = np.array([(x, t05, t16, t50, t84, t95)], dtype=dtype2t)

        
    # Get the times at fraction alpha between t0 and t_mode    
    alpha = np.array([0.2, 0.5, 1.0, 2.0, 5.0]) # fraction of t_mode - t0
    tma = t0 + (t_mode - t0) * alpha # tma = that are at alpha between t0 and t_mode
    
    # Then compute the cdf and pdf at these points
    ca = stats.lognorm.cdf(tma, *par_log_t)
    pa = stats.lognorm.pdf(tma, *par_log_t)
    
    # Finally store the values in bt_mode with dtype3
    L = [x] + list(tma) + list(ca) + list(pa)
    nms = [nm[0] for nm in dtype3t.descr]
    for nm, v in zip(nms, L):
        bt_tmode[ix][nm] = v
    
    # Plot the results (not all to prevent clutter)
    if x in [200, 400, 600, 800, 1000]:
        ax1.plot(t, cdf_o, '--', color=clr, label=f"cdf {x} m (original data)")
        ax1.plot(t, cdf_log_t, color=clr,
                 label=fr"x={x:4d} m, $\mu,\sigma,\,t_0$,loc=({mu_t:.2f}, {sig_t:.2f}, {loc_t:.3e} d)")
        ax2.plot(t, pdf_log_t, color=clr, 
                 label=fr"x={x:4d} m, $\mu,\sigma,\,t_0$=({mu_t:.2f}, {sig_t:.2f}, {loc_t:.3e} d)")
  
        ax1.plot(tma, stats.lognorm.cdf(tma, *par_log_t), 'o', mfc=clr, mec='k')
        ax2.plot(tma, stats.lognorm.pdf(tma, *par_log_t), 'o', mfc=clr, mec='k')

        ax1.plot(t_med, stats.lognorm.cdf(t_med, *par_log_t), 's', ms=7, mfc=clr, mec='k')
        ax2.plot(t_med, stats.lognorm.pdf(t_med, *par_log_t), 's', ms=7, mfc=clr, mec='k')
        
        ax1.plot(t_mean, stats.lognorm.cdf(t_mean, *par_log_t), '*', ms=10, mfc=clr, mec='k')
        ax2.plot(t_mean, stats.lognorm.pdf(t_mean, *par_log_t), '*', ms=10, mfc=clr, mec='k')
        
        ax1.text(0.6, 0.75, 
            r"$\blacksquare = t_{median},\ \bigstar = t_{mean}$", 
            transform=ax1.transAxes, 
            fontsize=14, 
            color='green',
            ha='left', va='top')
        ax2.text(0.6, 0.55,
            r"$\blacksquare = t_{median},\ \bigstar = t_{mean}$", 
            transform=ax2.transAxes, 
            fontsize=14, 
            color='green',
            ha='left', va='top')

    for ax in [ax1, ax2]:
        ax.grid(True)
        ax.legend()

fig. savefig(os.path.join(dirs.images, 'BTanalysis_t_' + case_name + '.png'))

# Just show the data with less digits for better overview
with np.printoptions(precision=5, suppress=True):
    print("bt_cdf_props_t:")
    print(bt_cdf_props_t)
    print("bt_times:")
    print(bt_times)
    print("bt_mode")
    print(bt_tmode)

#  ==== Save the data arrays to dirs.data =====
# Saved are:  bt_cdf_props_t, bt_times and bt_mode
fname = os.path.join(dirs.data, "BTarrays_t_" + case_name + '.npz')
np.savez(fname, bt_cdf_props_t=bt_cdf_props_t, bt_times=bt_times, bt_mode=bt_tmode)    

# %% == show normal distribution

# The normal distribution is the lognormal distribution versus ln(t) instead of t.
# It's iteresting, but we don't actually need it as the cdf tells us all we
# need and is often clearer. Also, we included the pdf on the real time axis
# already aboven. Nevertheless, it's nice to see nice standard pdf's pop up
# at the ln(t) scale.

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 11))
ax1.set(title='ln(cdf) (case ' + case_name + ') with the following k-field parameters\n' + k_field_str, xlabel='ln(t)', ylabel='p[ln(t) < ln(t(p))]') 
ax2.set(title='ln(pdf)', xlabel='ln(t)', ylabel='d(ln(cd)f)/dln(t)') 

print("{:5s} {:10s} {:10s}".format('x', 'mu', 'sigma'))
for ix, (X, t) in enumerate(zip(cumul_t['xObs'].T, cumul_t['time'].T)):
    x = X[0]
    
    cdf_o = p

    par_norm_t = stats.norm.fit(np.log(t))    
    cdf_norm_t = stats.norm.cdf(np.log(t), *par_norm_t)
    pdf_norm_t = stats.norm.pdf(np.log(t), *par_norm_t)

    mu_t, sig_t = par_norm_t
    print("{:5d} {:10.3g} {:10.3g}".format(x, mu_t, sig_t))
        
    ax1.plot(np.log(t), cdf_norm_t, label=fr"x={x} m, $\mu,\sigma$=({mu_t:.2f}, {sig_t:.2f})")
    ax2.plot(np.log(t), pdf_norm_t, label=fr"x={x} m, $\mu,\sigma$=({mu_t:.2f}, {sig_t:.2f})")

for ax in [ax1, ax2]:
    ax.grid(True)
    ax.legend()

fig. savefig(os.path.join(dirs.images, 'BTanalysis_ln_t' + case_name + '.png'))



# ======================= Now the x values instead of time ========================


# %% === Compute and plot breakthrough (x locations at a given time) ===

pass_xvals = np.load(os.path.join(dirs.data, 'pass_xvals_' + case_name + '.npy'))

ntObs = len(pass_xvals['tObs'][0])

dtype0x = np.dtype([('tObs', '<f8', (ntObs,)), ('x', '<f8', (ntObs,))])
cumul_x = np.zeros(len(pass_xvals), dtype=dtype0x)  # dict()
cumul_x['tObs'] = pass_xvals['tObs']

for j in range(ntObs):
    idx = np.argsort(pass_xvals['x'][:, j])    
    x = pass_xvals['x'][:, j][idx]    
    cumul_x['x'][:, j] = x

last_xvals = cumul_x['x'][:, -1]
x_last = last_xvals[int(0.995 * len(last_xvals))]

xlim = (0, round_up_nice(x_last, m=2))

# %% === Ploting the cumulative curves for x-vals:
fig, ax = plt.subplots(figsize=(10, 6))
ax.set(title="Cumulative location curves (case " + case_name + ")\n" + k_field_str, xlabel="x [m]", ylabel="Frac")

tObs = pass_xvals['tObs'][0]
N = len(pass_xvals)
p = (np.arange(N) + 0.5) / N
for j, to in enumerate(tObs):
    idx = np.argsort(pass_xvals['x'][:, j])    
    x = pass_xvals['x'][:, j][idx]
    ax.plot(x, p, label=f"t = {tObs[j]} d")

x_last = x[int(0.995 * len(x))]
xlim = (0, round_up_nice(x_last, m=2))

ax.set_xlim(xlim)
ax.grid()
ax.legend()

fig.savefig(os.path.join(dirs.images, "cumul_x_" + case_name +".png"))

# %% === Fit a cdf through the data


fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, figsize=(10, 11))
ax1.set(title='cdfs for x (case ' + case_name + ') with the following k-field parameters\n'
        + k_field_str, xlabel='x', ylabel='p[x < x(p)]') 
ax2.set(title='pdfs for x', xlabel='x [m]', ylabel='d(cdf)/dx')

# ax1.set_xlim(xlim)
# ax2.set_xlim(xlim)

bt_cdf__xprops = np.zeros(cumul_x['tObs'].shape[-1], dtype=dtype1x)
bt_xvals     = np.zeros(cumul_x['tObs'].shape[-1], dtype=dtype2x) # Breakthrough x-values
bt_xmode      = np.zeros(cumul_x['tObs'].shape[-1], dtype=dtype3x) # Breakthrough x-values
clrs = cycle(["blue", "red", "green", "black", "gray", "magenta", "cyan", "orange", "brown", "darkgray"])

print("Data from pass_xvals, properties of fitted lognormal distribution for each tObs:")

for it, (T, x) in enumerate(zip(cumul_x['tObs'].T, cumul_x['x'].T)):
    t = T[0].item()
    clr = next(clrs)    
    cdf_o = p

    par_log_x = stats.lognorm.fit(x)    
    cdf_log_x = stats.lognorm.cdf(x, *par_log_x)
    pdf_log_x = stats.lognorm.pdf(x, *par_log_x)

    shape_x, loc_x, scale_x = par_log_x
    mu_x, sig_x, x0 = np.log(scale_x), shape_x, loc_x

    print("tObs = {:<.0f} d, mu_x = {:<.3g} m, sigma_x = {:<.3g} m, loc = x0 = {:<.3g} m".format(t, mu_x, sig_x, x0))

    # Get lognormal distribution properties.
    x_mode = x0 + np.exp(mu_x - sig_x ** 2)
    x_med =  x0 + np.exp(mu_x)
    x_mean = x0 + np.exp(mu_x + 0.5 * sig_x ** 2)
    v_mean = x_mean / t
    v_mode = x_mean / t
    
    # Store them in bt_cdf_props_x recarray.
    bt_cdf__xprops[it] = np.array([(t, mu_x, sig_x, x0, x_mode, x_med, x_mean, v_mode, v_mean)], dtype=dtype1x)

    # Inversely get the times at which the cdf has given p values
    
    # First get the lognormal distribution with the fitted parameters
    dist = stats.lognorm(s=shape_x, loc=loc_x, scale=scale_x)
    
    # Use this distibution to get the times at which the cdf has these values
    # (useful for normal distribution less so for lognormal distribution)    
    x05, x16, x50, x84, x95 = (dist.ppf(0.05), dist.ppf(0.1587), dist.ppf(0.5),
                               dist.ppf(0.8413), dist.ppf(0.95))
    
    xpvals = (x05, x16, x50, x84, x95)
    
    # Store them, and for convenience round them in a recarray with dtype2
    bt_xvals[it] = np.array([(t, x05, x16, x50, x84, x95)], dtype=dtype2x)

    # Get the times at fraction alpha between t0 and t_mode    
    alpha = np.array([0.2, 0.5, 1.0, 2.0, 5.0]) # fraction of x_mode - x0
    xma = x0 + (x_mode - x0) * alpha # xma = that are at alpha between x0 and x_mode
    
    # Then compute the cdf and pdf at these points
    ca = stats.lognorm.cdf(tma, *par_log_x)
    pa = stats.lognorm.pdf(tma, *par_log_x)
    
    # Finally store the values in bt_mode with dtype3
    L = [t] + list(xma) + list(ca) + list(pa)
    nms = [nm[0] for nm in dtype3x.descr]
    for nm, v in zip(nms, L):
        bt_xmode[it][nm] = v
    
    # Plot the results (not all to prevent clutter)
    if t in cumul_x['tObs'][0][::2]:
        ax1.plot(x, cdf_o, '--', color=clr, label=f"cdf {t} m (original data)")
        ax1.plot(x, cdf_log_x, color=clr,
                 label=fr"x={t:.0f} m, $\mu,\sigma,\,t_0$,loc=({mu_x:.2f}, {sig_x:.2f}, {loc_x:.3e} d)")
        ax2.plot(x, pdf_log_x, color=clr, 
                 label=fr"t={t:.0f} m, $\mu,\sigma,\,t_0$=({mu_x:.2f}, {sig_x:.2f}, {loc_x:.3e} d)")
  
        # ax1.plot(xma, stats.lognorm.cdf(xma, *par_log_x), 'o', mfc=clr, mec='k')
        # ax2.plot(xma, stats.lognorm.pdf(xma, *par_log_x), 'o', mfc=clr, mec='k')

        ax1.plot(xpvals, stats.lognorm.cdf(xpvals, *par_log_x), 'o', mfc=clr, mec='k')
        ax2.plot(xpvals, stats.lognorm.pdf(xpvals, *par_log_x), 'o', mfc=clr, mec='k')


        ax1.plot(x_med, stats.lognorm.cdf(x_med, *par_log_x), 's', ms=7, mfc=clr, mec='k')
        ax2.plot(x_med, stats.lognorm.pdf(x_med, *par_log_x), 's', ms=7, mfc=clr, mec='k')
        
        ax1.plot(x_mean, stats.lognorm.cdf(x_mean, *par_log_x), '*', ms=10, mfc=clr, mec='k')
        ax2.plot(x_mean, stats.lognorm.pdf(x_mean, *par_log_x), '*', ms=10, mfc=clr, mec='k')
        
        ax1.text(0.6, 0.75, 
            r"$\blacksquare = x_{median},\ \bigstar = x_{mean}$", 
            transform=ax1.transAxes, 
            fontsize=14, 
            color='green',
            ha='left', va='top')
        ax2.text(0.6, 0.55,
            r"$\blacksquare = x_{median},\ \bigstar = x_{mean}$", 
            transform=ax2.transAxes, 
            fontsize=14, 
            color='green',
            ha='left', va='top')

    for ax in [ax1, ax2]:
        ax.set_xlim(xlim)
        ax.grid(True)
        ax.legend()

fig. savefig(os.path.join(dirs.images, 'BTanalysis_x_' + case_name + '.png'))

#  ==== Save the data arrays to dirs.data =====
# Saved are:  bt_cdf_props_x, bt_times and bt_mode
fname = os.path.join(dirs.data, "BTarrays_x_" + case_name + '.npz')
np.savez(fname, bt_cdf__xprops=bt_cdf__xprops, bt_times=bt_xvals, bt_xmode=bt_xmode)  

# %% 
plt.show()
