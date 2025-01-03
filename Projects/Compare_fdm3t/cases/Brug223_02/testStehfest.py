

# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import exp1, erfc, k0 as K0, k1 as K1, kv as Kv, factorial as fac
import mpmath as mpm
from functools import partial
from inverselaplace import *

# %%

intro = """
Test Stehfest's Laplace back-transformation.

To test we need a function of which we know both the real solution and the solution
in the Laplace space. With these two functions, we can evaluate the Laplace solution
in the Stehfest routine to retain the approximate true solution numerically. This result
be plotted together with the real solution for comparison.

Two solutions are tested.
    1) Sudden head at x=0, t=0 in a half-infinite aquifer (h(x,t) = h erfc(kappa / (2 sqrt(t))),
        with kappa = x sqrt(kD / S))
    2) Theis well drawdown. s(r, t) = Q / (4 pi kD) exp1(r^2 S / (4 kD t))

Both examples show the match is good.

TO 20230306
"""

def stehfest_coefs(N=10):
    """Return the N Stehfest coefficients
    
    Note that Stehfest used k^(N/2 + 1) instead of k^(N/2). The latter
    as used by Hemker and Maas (1987) etc. is correct.
    
    @ TO 20230306  
    """
    v = np.zeros(N, dtype=float)
    for i in range(1, N + 1):
        j1, j2 = int((i + 1) // 2), int(min(i, N // 2))
        for k in range(j1, j2+1):
            v[i - 1] += k ** (N // 2) * fac(2 * k) / (
                fac(N // 2 - k) * fac(k) *
                fac(k - 1) * fac(i - k) *
                fac(2 * k - i)
            )
        v[i - 1] *= (-1) ** (i + N//2)
    return v

# Compute the coefficients once
vStehfest = stehfest_coefs(N=10)

# Laplace transform of sudden head change by 1 m at x=0, t=0
def FL_sudden(p, kappa=None):
    """Return the Laplace transform for sudden head change by 1 m at x=0, t=0.
    
    Parameters
    ----------
    p: float
        Laplace parameter
    kappa: float
        x sqrt(S/ kD)
    """
    return np.exp(-kappa * np.sqrt(p)) / p

def F_sudden(ts, kappa=None):
    """Return drawdown due to sudden head change by 1 m at x=0, t=0.
    
    Parameters
    ----------
    ts: np.ndarray
        times
    kappa: float
        x sqrt(S/ kD)
    """
    return erfc(kappa / (2 * np.sqrt(ts)))


def FL_theis(p=None, Q=None, kD=None, S=None, r=None, **kw):
    """Return Laplace-transformed Theis drawdown.
    
    Parameters
    ----------
    p: float [1/T]
        Laplace variable
    Q, kD, S, r: floats
        as usual
    kw: dict
        capture superfluous parameters
    """
    return Q / (2 * np.pi * kD * p) * K0(r * np.sqrt(p * S / kD))

def F_theis(Q=None, kD=None, S=None, r=None, ts=None):
    """Return Theis drawdown.
    
    Parameters
    ----------
    p: float [1/T]
        Laplace variable
    Q, kD, S, r: floats
        as usual
    ts: np.ndarray
        times
    """
    return Q / (4 * np.pi * kD) * exp1(r ** 2 * S /(4 * kD * ts))

def stehfest(lfunc, pars, ts):
    """Back-tranform lfunc(**pars) at times ts.
    
    lfunc: function
        Laplace transform as a function
    pars: dict
        parameters that enter lfunc
    ts: np.ndarray
        times for which the back-transform is requested
    """
    N = len(vStehfest)
    f = np.zeros_like(ts)
    for it, t in enumerate(ts):
        for v, i in zip(vStehfest, range(1, N + 1)):
            p = i * np.log(2) / t
            f[it] += v * lfunc(p, **pars)
        f[it] *= np.log(2) / t
    return f

# Laplace transforms of Bruggeman's solution 223_02
def fhat64(p, hb=None, r=None, R=None, S=None, kD=None):
    """Laplace transform of Phi Burgeman 223_02.
    
    using 64 bit float accuracy.
    """
    beta = np.sqrt(S / kD)
    return (hb / p * Kv(0, beta * r * np.sqrt(p)) /
                     Kv(0, beta * R * np.sqrt(p)))

def qhat64(p, hb=None, r=None, R=None, S=None, kD=None):
    """Lapace transform of Q of Bruggeman 223_02.
    
    using normal 64 bit float accuracy.
    """
    beta = np.sqrt(S / kD)
    return (2 * np.pi * r * np.sqrt(S * kD) * hb  / np.sqrt(p) *
            Kv(1, beta * r * np.sqrt(p)) /
            Kv(0, beta * R * np.sqrt(p)))

def fback_stehfest(lapl_func, times, N, args):
    """Stehfest back transformation.
    
    Simple straightforward back transformation by Stehfest.
    """
    zeta = stehfest_coefs(N)
    if np.isscalar(times):
        times = np.array([times])
    s = np.zeros_like(times)
    for it, t in enumerate(times):
        for k in range(1, N + 1):            
            p = k * np.log(2) / t
            s[it] += zeta[k - 1] * lapl_func(p, *args)
        s[it] *= np.log(2) / t
    return s

 
 # Laplace transforms of Bruggeman's solution 223_02
def fhat_mpm(p, hb=None, r=None, R=None, S=None, kD=None):
    """Laplace transform of Phi Burgeman 223_02"""
    beta = mpm.sqrt(S / kD)
    return (hb / p * mpm.besselk(0, beta * r * mpm.sqrt(p)) /
                     mpm.besselk(0, beta * R * mpm.sqrt(p)))

def qhat_mpm(p, hb=None, r=None, R=None, S=None, kD=None):
    """Lapace transform of Q of Bruggeman 223_02"""
    beta = mpm.sqrt(S / kD)
    return (2 * mpm.pi() * r * mpm.sqrt(S * kD) * hb  / mpm.sqrt(p) *
            mpm.besselk(1, beta * r * mpm.sqrt(p)) /
            mpm.besselk(0, beta * R * mpm.sqrt(p)))

 # %%
      
if __name__ == '__main__':  

    kD, S = 600, 0.001
    ts = np.logspace(-4, 3, 71)
    
    # Sudden change at x=0, t=0:
    x = 100.
    kappa = x * np.sqrt(S / kD)
    pars = {'kappa': kappa}
    f = stehfest(FL_sudden, pars, ts)
    
    fig, ax = plt.subplots()
    ax.set_title("Test Stehfest on sudden head change (erfc)")
    ax.set_xlabel('x [m]')
    ax.set_ylabel('s [m]')
    ax.set_xscale('log')
    ax.grid()
    ax.plot(ts, F_sudden(ts, **pars), label='Original')
    ax.plot(ts, f, '.', label='Stehfest')
    ax.set_ylim(ax.get_ylim()[::-1])
    ax.legend()
    
    
    # Theis drawdown
    pars_theis = {'Q': 600, 'kD': 600, 'S': 0.001, 'ts': ts, 'r': 100}
    f = stehfest(FL_theis, pars_theis, ts)
    fig, ax = plt.subplots()
    ax.set_title("Test Stehfest on Theis drawdown (exp1)")
    ax.set_xlabel('t [d]')
    ax.set_ylabel('s [m]')
    ax.set_xscale('log')
    ax.grid()
    ax.plot(ts, F_theis(**pars_theis), label='Original_theis')
    ax.plot(ts, f, '.', label='Stehfest')
    ax.set_ylim(ax.get_ylim()[::-1])
    ax.legend()
    
# %%
if True:
    method = 'stehfest'
    # method = 'talbot'
    method = 'dehoog'
    #method = 'cohen'
    print(f"Laplace inversion: {method}")
    
    ts = np.logspace(-4, 6, 36)
    
    hb, r, R, S, kD = 1., 30., 30., 0.1, 1000.
    rs = 30., 60., 90, 120
    
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 10))
    fig.suptitle(f"Bruggeman (1999) Solution 223_02\nNumerieke inversie volgens {method}\nhb={hb}, S={S}, kD={kD}")
    ax1.set_title("Head, different numerical backtransformations")
    ax2.set_title("Q, different numerical backtransformations")
    ax1.set_ylabel("Phi")
    ax1.set_ylim(-0.1, 1.1)
    ax2.set_ylabel("Q")
    ax2.set_xlabel("t [d]")
    ax2.set_yscale('log')
    ax2.set_xscale('log')
    ax2.set_ylim(1e-3, 1e4)
    ax1.grid(True)
    ax2.grid(True)


    invertors = {'FixedTalbot': FixedTalbot,
             'Stehfest': Stehfest,
             'DeHoog': DeHoog,
             'Cohen': Cohen}

    methods = ['Stehfest', 'FixedTalbot']    
    methods = ['FixedTalbot']    
        
    for r in rs:
        print(f"FP(r={r})")
        FP = partial(fhat_mpm, hb=hb, r=r, R=R, S=S, kD=kD)
        print(f"QP(r={r})")
        QP = partial(qhat_mpm, hb=hb, r=r, R=R, S=S, kD=kD)

        fi = np.zeros_like(ts)
        qi = np.zeros_like(ts)
        # This is using mpmath
        #for it, t in enumerate(ts):
        #    print(f"t = {t}")
        #    fi[it] = mpm.invertlaplace(FP, t, method=method)
        #    qi[it] = mpm.invertlaplace(QP, t, method=method)

        #ax1.plot(ts, fi, label=f"r={r}, R={R}")
        #ax2.plot(ts, qi, label=f"r={r}, R={R}")
        
        # ax1.plot(ts, fback_stehfest(fhat64, ts, N=12, args=(hb, r, R, S, kD)), '.-',
        #          label=f'Eigen stehfest, r={r}')
        # ax2.plot(ts, fback_stehfest(qhat64, ts, N=12, args=(hb, r, R, S, kD)), '.-', 
        #          label=f'Eigen stehfest, r={r}')   
        
        fpF = partial(fhat64, hb=hb, r=r, R=rs[0], S=S, kD=kD)
        fpQ = partial(qhat64, hb=hb, r=r, R=rs[0], S=S, kD=kD)

        for method in methods:
            invertor = invertors[method]()
        
            F_ = invertor.calc_time_domain_solution(fpF, ts, degree=10)
            Q_ = invertor.calc_time_domain_solution(fpQ, ts, degree=10)
        
            ax1.plot(ts, F_, 'o-', mfc='none', label=f'{method}, r={r}')
            ax2.plot(ts, Q_, 'o-', mfc='none', label=f'{method}, r={r}')
        
        print("F: ", F_[:5])
        print("Q: ", Q_[:5])
        print()
        
    print("Done")

    ax1.legend()
    ax2.legend()

# fig.savefig(f"{method}.png")


plt.show()
# %%
