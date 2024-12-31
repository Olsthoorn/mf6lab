import numpy as np
import matplotlib.pyplot as plt
from scipy.special import k0 as K0, k1 as K1
from math import factorial as fac
from testStehfest import stehfest_coefs
from fdm.fdm3t import fdm3t, Grid, dtypeH
from etc import color_cycler

# Stehfest
def sf_coefs(N):
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

def fi_hat(p, args):
    """Laplace transformed head L{phi(r, t)}
    
    Note args: hb, rm, R, S, kD
    """
    h, r, R, S, kD = args
    beta = np.sqrt(S / kD)
    return h / p * K0(beta * r * np.sqrt(p) / K0(beta * R * np.sqrt(p)))

def q_hat(p, args):
    """Return Laplace transformed flow L{Q(r, t)}.
    
    Note args: hb, rm, R, S, kD
    """
    h, r, R, S, kD = args
    beta = np.sqrt(S  /kD)
    return 2 * np.pi * r * np.sqrt(S * kD) * h  / np.sqrt(p) * K1(beta * r * np.sqrt(p)) /K0(beta * R * np.sqrt(p))
    
class Brug223:

    def __init__(self, gr=None, t=None, hb=None, kD=None, S=None):        
        D = gr.DZ[0, 0, 0]     
        k = gr.const(kD / D)
        ss = gr.const(S / D)
        hi = gr.const(0.)
        idomain = gr.const(1, dtype=int)
        Ifh = gr.NOD[0:1, 0:1, 0].flatten()
        fh = np.zeros(len(Ifh), dtype=dtypeH)
        fh['I'], fh['h'] = Ifh, hb
        self.out = fdm3t(gr=gr, t=t,
                    k=(k, k, k), ss=ss, fh={0: fh}, hi=hi, idomain=idomain)
        return None
    
    @property
    def heads(self):
        return self.out['Phi']
    
    @property
    def flows(self):
        return self.out['Qx']

                      
def stehfback(fhat, t, N, args):
    """Stehfest back transformation for sequende of t values.
    
    Parameters:
    ----------
    fhat: function
        the Laplace transformed function that has to be back-tranformed
    t: sequence
        times
    N: int, must be even.
        number of Stehfest coefficients (must be even)
    """
    csf = sf_coefs(N)
    if np.isscalar(t):
        t = np.array([t])
    s = np.zeros_like(t)
    for it, t_ in enumerate(t):
        for k in range(1, N + 1):            
            p = k * np.log(2) / t_
            s[it] += csf[k - 1] * fhat(p, args)
        s[it] *= np.log(2) / t_
    return s
    
if __name__ == '__main__':
    
    N = 10  
    scn = sf_coefs(N)
    sc0 = stehfest_coefs(N)
    print(scn - sc0)

    hb, r0, R, S, kD, D = 1., 30., 30., 0.2, 600, 50.
    
    rs = np.hstack((R - 0.01, np.logspace(np.log10(R), 6, 60)))
    times = np.logspace(-2, 2, 41)
    
    gr = Grid(rs, None, [0, -D], axial=True)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 10))
    ax1.set_title(f"Heads, hb={hb}, R={R}, kD={kD}, S={S}")
    ax2.set_title(f"Flows, hb={hb}, R={R}, kD={kD}, S={S}")
    ax2.set_xlabel("times [d]")
    ax1.set_ylabel("head [m]")
    ax2.set_ylabel("flow [m3/d]")
    ax2.set_xscale("log")
    ax1.grid(True)
    ax2.grid(True)
    
    brug = Brug223(gr=gr, t=times, hb=hb, kD=kD, S=S)
    
    clrs = color_cycler()
    for ir in range(0, min(gr.nx, 30), 5):
        clr = next(clrs)
        heads = brug.heads[:, -1, -1, ir]
        flows = brug.flows[:, -1, -1, ir]
        rm = gr.xm[ir]
    
        fi = stehfback(fi_hat, times, N, (hb, rm, R, S, kD))
        Q  = stehfback(q_hat,  times, N, (hb, rm, R, S, kD))
        
        ax1.plot(times, heads[:], '.-', color=clr, mfc='none',
                 label=f"fdm3t, r={rm:.4g}")
        ax1.plot(times, fi, 's-', color=clr, mfc='none',
                 label=f"brug,  r={rm:.4g}")
        ax2.plot(times[1:], flows[:], '.-', color=clr, mfc='none',
                 label=f"fdm3t, r={rm:.4g}")
        ax2.plot(times, Q, 's-', color=clr, mfc='none', label=f"brug, r={rm:.4g}")
    ax1.legend()
    ax2.legend()
    
    plt.show()
    

    
   