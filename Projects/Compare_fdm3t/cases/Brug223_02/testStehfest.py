import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la
import scipy.special as sp

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
            v[i - 1] += k ** (N // 2) * sp.factorial(2 * k) / (
                sp.factorial(N // 2 - k) * sp.factorial(k) *
                sp.factorial(k - 1) * sp.factorial(i - k) *
                sp.factorial(2 * k - i)
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
    return sp.erfc(kappa / (2 * np.sqrt(ts)))


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
    return Q / (2 * np.pi * kD * p) * sp.k0(r * np.sqrt(p * S / kD))

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
    return Q / (4 * np.pi * kD) * sp.exp1(r ** 2 * S /(4 * kD * ts))

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
    plt.show()
    

