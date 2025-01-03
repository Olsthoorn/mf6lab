# Solution of Bruggeman 220_03, flow from cylinder of given R
# after sudden head change at t=0.
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j0 as J0, y0 as Y0, j1 as J1, y1 as Y1
from scipy.special import k0 as K0, k1 as K1
from scipy.integrate import quad # , quad_explain
from math import factorial as fac
import settings
from etc import color_cycler
from mf_adapt import dirs

# Integral solution
def h(u=None, r=None, R=None, **_):
    rho = r / R
    return ((J0(u) * Y0(u * rho) - Y0(u) * J0(u * rho)) /
            (J0(u) ** 2 + Y0(u) ** 2))

def dhdr(u=None, r=None, R=None, **_):
    """Return partial h/ partial r."""
    rho = r / R
    return (u / R * ((Y0(u) * J1(rho * u) - J0(u) * Y1(rho * u))) /
            (J0(u) ** 2 + Y0(u) ** 2))

def argFi_z(z=None, r=None, R=None, S=None, kD=None, t=None, **_):
    u = np.exp(z)
    beta = np.sqrt(S / kD)
    return h(u, r, R) * np.exp(- (u / (beta * R)) ** 2 * t)

def Fint_z(a, b, hb=None, r=None, R=None, S=None, kD=None, t=None):
    """Return phi(t, r).
    Parameters
    ----------
    a, b: limit of integration (should be -inf, inf)
    args: tuple of floats of remaining parameters.
        args = (hb, r, R, S, kD, t)    
    """
    args = (r, R, S, kD, t)
    return hb * (1. - 2. / np.pi * quad(
        argFi_z, a, b, args=args, limit=200)[0])

def Fint_z2(a, b, hb=None, r=None, R=None, S=None, kD=None, t=None):
    """Return phi(t, r).
    Parameters
    ----------
    a, b: limit of integration (should be -inf, inf)
    args: tuple of floats of remaining parameters.
        args = (hb, r, R, S, kD, t)    
    """
    args = (r, R, S, kD, t)
    left  = quad(argFi_z, -100, a, args=args, limit=200)[0]
    right = quad(argFi_z,    a, b, args=args, limit=200)[0]
    return hb * (1. - 2. / np.pi * (left + right))


def argQ_z(z=None, r=None, R=None, S=None, kD=None, t=None):
    u = np.exp(z)
    beta = np.sqrt(S / kD)
    return dhdr(u, r, R) * np.exp(- (u / (beta * R)) ** 2 * t)

def Qint_z(a, b, hb=None, r=None, R=None, S=None, kD=None, t=None):
    """Return Q(t, r).
    
    Parameters
    ----------
    a, b: floats, limits of the integration
    """
    return 4 * kD * hb * r * quad(
        argQ_z, a, b, args=(r, R, S, kD, t), limit=200)[0]
    
def Qint_z2(a, b, hb=None, r=None, R=None, S=None, kD=None, t=None):
    """Return Q(t, r).
    
    Parameters
    ----------
    a, b: floats, limits of the integration
    """
    left  = 4 * kD * hb * r * quad(argQ_z, -100, a, args=(r, R, S, kD, t), limit=200)[0]
    right = 4 * kD * hb * r * quad(argQ_z,    a, b, args=(r, R, S, kD, t), limit=200)[0]
    return left + right

def Qint_z_simpson(a, b, Np=None, hb=None, r=None, R=None, S=None, kD=None, t=None):
    """Return Q(t, r).
    
    Parameters
    ----------
    a, b: floats, limits of the integration
    """
    z = np.linspace(a, b, Np)
    F = argQ_z(z, r=r, R=R, S=S, kD=kD, t=t)
    integral = np.sum(np.diff(z) * 0.5 * (F[:-1] + F[1:]))
    return  4 * kD * hb * r * integral


# Stehfest for numerical back transformation of the Laplace transform
def sf_coefs(N):
    """Graver-Stehfest coefficients to numerically invert laplace transform.
    
    Paraneters:
    -----------
    N: int, must be even and rather <= 18 to prevent los of digits
        Stehfests number of coefficients
    """
    assert N % 2 == 0, "N must be even integer <= 18"
    z = np.zeros(N, dtype=float)
    for i in range(1, N + 1):                
        j1, j2 = int((i + 1) // 2), int(min(i, N // 2))        
        for j in range(j1, j2 + 1):
            z[i - 1] += (j ** (N // 2) * fac(2 * j) /
                (fac(N // 2 - j) * fac(j) * fac(j - 1) *
                 fac(i - j) * fac(2 * j - i))
            )                                        
        z[i - 1] *= (-1) ** (N // 2 + i)
    return z

def fhat(p, hb=None, r=None, R=None, S=None, kD=None):
    """Laplace transform of Phi Burgeman 223_02"""
    beta = np.sqrt(S / kD)
    return hb / p * K0(beta * r * np.sqrt(p) / K0(beta * R * np.sqrt(p)))

def qhat(p, hb=None, r=None, R=None, S=None, kD=None):
    """Lapace transform of Q of Bruggeman 223_02"""
    beta = np.sqrt(S / kD)
    return 2 * np.pi * r * np.sqrt(S * kD) * hb  / np.sqrt(p) * K1(beta * r * np.sqrt(p)) /K0(beta * R * np.sqrt(p))
    
                      
def Fback(lapl_func, times, N, args):
    """Stehfest back transformation."""
    zeta = sf_coefs(N)
    if np.isscalar(times):
        times = np.array([times])
    s = np.zeros_like(times)
    for it, t in enumerate(times):
        for k in range(1, N + 1):            
            p = k * np.log(2) / t
            s[it] += zeta[k - 1] * lapl_func(p, *args)
        s[it] *= np.log(2) / t
    return s

if __name__ == '__main__':
    #if True:
    kD = np.sum(settings.props['kr'] * settings.props['D'])
    S  = np.sum(settings.props['ss'] * settings.props['D'])
    hb = settings.hb
    u = np.logspace(-8, 4, 100000)
    z = np.log10(u)
    aL, (a, b) = -100., z[[0, -1]]
    R = settings.props['r'][1]
    rs = np.array([1.0, 2.0, 3.0, 4.0]) * R

    Nstehfest = 10
    
    figsize = (10., 10.)

    times = settings.props['t']
    times = np.logspace(-3, 9, 121)

    # Argument of phi    
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=figsize)
    fig.suptitle(f"argFi, R={R:.4g} m,kD={kD:.4g} m2d, S={S:.4g}")

    for r, ax in zip(rs, axs.ravel()):
        ax.set_title(f"argFi, r={r:.4g} m, r/R={r/R:.4g}")
        ax.set_ylabel('argFi')
        ax.set_xlabel('z = np.exp(u)')
        ax.grid(True)    
        for t in times[::10]:
            ax.plot(z, argFi_z(z=z, r=r, R=R, S=S, kD=kD, t=t), label=f't={t:4g} d')
        if np.isclose(r, 30):
            ax.legend(loc='upper left')

    # Argument of Q
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=figsize)
    fig.suptitle(f"argFi, R={R:.4g} m,kD={kD:.4g} m2d, S={S:.4g}")
    for r, ax in zip(rs, axs.ravel()):
        ax.set_title(f"argQ, r={r:.4g} m, r/R={r/R:.4g}")
        ax.set_ylabel('argQ')
        ax.set_xlabel('z = np.exp(u)')
        ax.grid(True)
        for t in times[::10]:
            ax.plot(z, argQ_z( z=z, r=r, R=R, S=S, kD=kD, t=t), label=f't={t:.3g} d')
        if np.isclose(r, 30):
            ax.legend(loc='upper left')

    title = f"""Bruggeman(1999) sol. 223_03.
Compare integration over different extents using scipy.integrate.quad and simpson.
h={hb:.4g}, R= {R:.4g}, S={S:.4g}, kD={kD:.4g}
"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
    fig.suptitle(title)
    ax1.set_title("Phi")
    ax1.set_ylabel("Phi")
    ax2.set_title("Q")    
    ax1.set_xlabel('t')
    ax2.set_ylabel('Q m3/d')
    ax1.grid(True)
    ax2.grid(True)
    ax1.set_xscale('log')
    ax2.set_xscale('log')
    
    Np = 100000
    clrs = color_cycler()
    for r in rs:
        clr = next(clrs)
        FiQuad   = np.zeros_like(times)
        FiQuad2  = np.zeros_like(times)
        Qsimpson = np.zeros_like(times)
        Qquad    = np.zeros_like(times)
        Qquad2   = np.zeros_like(times)
        for it, t in enumerate(times):
            args = (hb, r, R, S, kD, t)
            FiQuad[ it] = Fint_z( a, b, hb=hb, r=r, R=R, S=S, kD=kD, t=t)
            FiQuad2[it] = Fint_z2(a, b, hb=hb, r=r, R=R, S=S, kD=kD, t=t)
            
            Qsimpson[it] = Qint_z_simpson(a, b, Np=Np, hb=hb, r=r, R=R, S=S, kD=kD, t=t)
            Qquad[it]    = Qint_z( a, b, hb=hb, r=r, R=R, S=S, kD=kD, t=t)
            Qquad2[it]   = Qint_z2(a, b, hb=hb, r=r, R=R, S=S, kD=kD, t=t)
            
        ax1.plot(times, FiQuad,   '+-', color=clr, label=f"Fi quad,    r={r:.4g} (a, b)=({a:.4g}, {b:.4g})")
        ax1.plot(times, FiQuad2,  '.-', color=clr, mfc='none', label=f"Fi quad2,    r={r:.4g}, (a, b)=({aL:.4g}, {b:.4g})")
        ax2.plot(times, Qsimpson, 'x-', color=clr, label=f"Q simpson, r={r:.4g}, (a, b)=({aL:.4g}, {b:.4g}), Np={Np}")
        ax2.plot(times, Qquad,    '+-', color=clr, label=f"Q quad,    r={r:.4g}, (a, b)=({a:.4g}, {b:.4g})")
        ax2.plot(times, Qquad2,   '.-', color=clr, mfc='none', label=f"Q quad2,    r={r:.4g}, (a, b)=({aL:.4g}, {b:.4g})")
        ax2.plot(times, Fback(qhat, times, Nstehfest, (hb, r, R, S, kD)),
                 'o-', color=clr, mfc='none', label=f"Laplace, r={r:.4g} m")
    ax1.legend(loc="lower right", fontsize="xx-small")
    ax2.legend(loc='upper right', fontsize="xx-small") 
    
    plt.savefig(os.path.join(dirs.images, "CompareIntegration.png"))
    
    plt.show()

