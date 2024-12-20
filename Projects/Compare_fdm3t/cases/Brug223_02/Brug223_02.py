# Solution of Bruggeman 220_03, flow from cylinder of given R after sudden head change at t=0.
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j0 as J0, y0 as Y0, j1 as J1, y1 as Y1
from scipy.integrate import quad

import settings

def h(u=None, r=None, R=None, **_):
    rho = r / R
    return (J0(u) * Y0(u * rho) - Y0(u) * J0(u * rho)) / (J0(u) ** 2 + Y0(u) ** 2)

def dhdr(u=None, r=None, R=None, **_):
    """Return partial h/ partial r."""
    rho = r / R
    return u / R * ((Y0(u) * J1(rho * u) - J0(u) * Y1(rho * u))) / (J0(u) ** 2 + Y0(u) ** 2)

def argFi(z=None, r=None, R=None, beta=None, t=None, **_):
    u = np.exp(z)
    return h(u, r, R) * np.exp(- (u / (beta * R)) ** 2 * t)

def Fint(a, b, args=None):
    """Return phi(t, r).
    Parameters
    ----------
    a, b: limit of integration (should be -inf, inf)
    args: tuple of floats of remaining parameters.
        args = (r, R, beta, t)    
    """
    return 1. - 2. / np.pi * quad(argFi, a, b, args=args, limit=200)[0]

def argQ(z=None, r=None, R=None, beta=None, t=None, **kw):
    u = np.exp(z)
    return u * dhdr(u, r, R) * np.exp(- (u / (beta * R)) ** 2 * t)

def Qint(a, b, args=None):
    """Return Q(t, r).
    
    Parameters
    ----------
    a, b: floats, limits of the integration
    args: tuple of parameters
        args = (r, R, beta, t)    
    """
    return quad(argQ, a, b, args=args, limit=200)[0]

kD = np.sum(settings.props['kr'] * settings.props['D'])
S  = np.sum(settings.props['ss'] * settings.props['D'])
beta = np.sqrt(S / kD)

z = np.linspace(-7.5, 5., 1001)

R = settings.props['r'][1]
r = 12 * R
rs = np.array([1.0, 2.0, 3.0, 4.0]) * R

figsize = (10., 10.)

times = settings.props['t']

fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=figsize)
fig.suptitle(f"argFi, R={R:.4g} m,kD={kD:.4g} m2d, S={S:.4g}, beta={beta:.4g}")

for r, ax in zip(rs, axs.ravel()):
    ax.set_title(f"argFi, r={r:.4g} m, r/R={r/R:.4g}")
    ax.set_ylabel('argFi')
    ax.grid(True)
    
    for t in times[::10]:
        ax.plot(z, argFi(z=z, r=r, R=R, beta=beta, t=t), label=f't={t:4g} d')
    if np.isclose(r, 30):
        ax.legend(loc='upper left')

        
fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=figsize)
fig.suptitle(f"argFi, R={R:.4g} m,kD={kD:.4g} m2d, S={S:.4g}, beta={beta:.4g}")
for r, ax in zip(rs, axs.ravel()):
    ax.set_title(f"argQ, r={r:.4g} m, r/R={r/R:.4g}")
    ax.set_ylabel('argQ')
    ax.set_xlabel('z = np.exp(u)')
    ax.grid(True)

    for t in times[::10]:
        ax.plot(z, argQ( z=z, r=r, R=R, beta=beta, t=t), label=f't={t:.3g} d')
    if np.isclose(r, 30):
        ax.legend(loc='upper left')

_, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=figsize)
ax1.set_title(f"Bruggeman 220_03, head R={R:.4g} m, kD={kD:.4g} m2d, S={S:.4g}, beta={beta:.4g}")
ax1.set_xlabel('t [d]')
ax1.set_ylabel('phi [m]')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.grid(True)

ax2.set_title(f"Bruggeman 220_03 Q, R={R:.4g} m,kD={kD:.4g} m2d, S={S:.4g}, beta={beta:.4g}")
ax2.set_xlabel('t [d]')
ax2.set_ylabel('Q [m3/d]')
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.grid(True)

B = dict()
B['Q'] = dict()
B['F'] = dict()
B['t'] = dict()
F = dict()
a, b = -7.5, 5
for r in rs:
    ri = np.round(r)
    B['Q'][ri] = np.zeros_like(times)
    B['F'][ri] = np.zeros_like(times)
    B['t'][ri] = np.zeros_like(times)
    for it, t in enumerate(times):
        args = (r, R, beta, t)
        B['Q'][ri][it] = 4 * kD * settings.hb * Qint(a, b, args=args)
        B['F'][ri][it] = settings.hb * Fint(a, b, args=args)
        B['t'][ri][it] = t
    ax1.plot(times, B['F'][ri], label=f'Brug220_00 r={r:.4g} m')
    ax2.plot(times, B['Q'][ri], label=f'Brug220_00 r={r:.4g} m')

ax1.legend()
ax2.legend()

#plt.show()

