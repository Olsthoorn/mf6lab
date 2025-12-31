"""Entry/exit resistance between ground and ditch is computed in this module.

This is done first by computing it numerically with the steady state
program fdm3.py under tools/fdm/src/fdm3d.py.

Thereafter we check it with the analytical Schwarz-Cristoffel conformal
transformation.

The idea is to derive a suitable analytial formula for the resistance in the
form of  w = 2 /(pi k) ln(D / Omega). The resistance is the extra drawdown
beyond the contact face between ditch and ground not the center of the ditch.

Numerically one can compute the drawdown at a given distance x due to a unit q
and compare it with the case  in which the ditch is fully penetrating.

The ditch bottom and ditch wall resistance is not included yet.
"""
# %%
import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from fdm.src.fdm3 import fdm3
from fdm.src.mfgrid import Grid
from matplotlib.patches import PathPatch, Path
from pathlib import Path as LPath
from scipy.interpolate import RegularGridInterpolator

# %%

#  De tabel coordinaten moeten alleen verfijnd worden naar 0 toe.

dpp_table = np.array([    
 [1.00E+06, 1.23E+00, 8.97E-01, 7.33E-01, 6.36E-01, 5.73E-01, 5.34E-01, 5.05E-01, 4.85E-01, 4.71E-01, 4.60E-01],
 [1.20E+00, 7.14E-01, 5.57E-01, 4.71E-01, 4.17E-01, 3.82E-01, 3.62E-01, 3.46E-01, 3.35E-01, 3.28E-01, 3.22E-01],
 [7.53E-01, 4.78E-01, 3.79E-01, 3.24E-01, 2.90E-01, 2.69E-01, 2.57E-01, 2.48E-01, 2.42E-01, 2.38E-01, 2.36E-01],
 [5.04E-01, 3.23E-01, 2.58E-01, 2.22E-01, 2.00E-01, 1.87E-01, 1.82E-01, 1.77E-01, 1.74E-01, 1.72E-01, 1.71E-01],
 [3.38E-01, 2.13E-01, 1.70E-01, 1.47E-01, 1.34E-01, 1.27E-01, 1.25E-01, 1.23E-01, 1.21E-01, 1.20E-01, 1.20E-01],
 [2.29E-01, 1.40E-01, 1.11E-01, 9.71E-02, 8.98E-02, 8.58E-02, 8.60E-02, 8.50E-02, 8.44E-02, 8.42E-02, 8.40E-02],
 [1.41E-01, 8.11E-02, 6.41E-02, 5.67E-02, 5.33E-02, 5.17E-02, 5.31E-02, 5.28E-02, 5.27E-02, 5.26E-02, 5.26E-02],
 [7.70E-02, 4.03E-02, 3.20E-02, 2.91E-02, 2.79E-02, 2.72E-02, 2.92E-02, 2.92E-02, 2.91E-02, 2.91E-02, 2.91E-02],
 [3.39E-02, 1.51E-02, 1.24E-02, 1.15E-02, 1.11E-02, 1.07E-02, 1.28E-02, 1.28E-02, 1.28E-02, 1.28E-02, 1.28E-02],
 [8.73E-03, 2.81E-03, 2.27E-03, 1.87E-03, 1.51E-03, 1.18E-03, 3.26E-03, 3.26E-03, 3.26E-03, 3.26E-03, 3.26E-03],
 [0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00]
 ])

D, dx, dy, nx, ny = 10., 0.1, 0.1, 10, 10
bs = np.linspace(0, D, int(D / (dx * nx)) + 1)
hs = np.linspace(0, D, int(D / (dy * ny)) + 1)
dpp_interp = RegularGridInterpolator((bs, hs), dpp_table, method='cubic')
dpp_interp((0, 3))


# %%
def dPP(out, k, L , D):
    Phi_end = out['Phi'][1:, 0, -1].mean()
    Q = out['Q'][1:, 0, -1].sum()
    return (Phi_end - Q*L / (k*D)) / Q

def dpp_contraction(h, D):
    """Rreturn dpp by contraction according to Verruijt(1970) p119."""
    return - 2 / np.pi * np.log(np.sin(np.pi * h / (2 * D)))

def strfun(gr, out):
    S = np.zeros((gr.nz + 1, gr.nx-1))
    S[:-1] = np.cumsum(out['Qx'][::-1, 0, :], axis=0)[::-1]
    return S

def compute_dpp_table(L=100, D=10, dx=0.1, dy=0.1, nx=10, ny=10):
    HUGE = 1000.
    L, D = 20, 10
    dx, dy = 0.1, 0.1

    irow = 0
    dpp = []    
    for h in np.arange(0, D + 0.01, ny * dy):
        dpp_line = []        
        for b in np.arange(0, D + 0.01, nx * dx):
            
            x = np.linspace(-b - dx, L + dx, int((L + b + 2*dx) / dx) + 1)
            z = np.linspace(dy, -D, int((D + dy) / dy + 1))
            gr = Grid(x, None, z)

            # xy = ((-b - dx, dy), (-b -dx, -h), (0, -h), (0, dy), (-b -dx, dy))
            # p = PathPatch(Path(xy), color='blue')

            DITCH = np.logical_and(gr.XM < 0, gr.ZM > -h)

            IBOUND = gr.const(1, dtype=int)
            IBOUND[:, :, 0] = 0
            IBOUND[0, :, :] = 0
            IBOUND[DITCH] = -1

            k = 1.
            K = gr.const(k)
            K[0, :, :] = 1e-8
            K[:, :, 0] = 1e-8
            K[DITCH] = HUGE
            K[1:, :, -1] = HUGE

            FQ = gr.const(0.)

            kdz = IBOUND[:, :, -1] * K[:, :, -1] * gr.DZ[: ,:, -1]
            FQ[:, :, -1] = -kdz / np.sum(kdz)

            HI = gr.const(0.)

            out= fdm3(gr, K=K, c=None, FQ=FQ, HI=HI, IBOUND=IBOUND, GHB=None)
            print('.', end="")
            
            dpp_line.append(dPP(out, k, L , D))
        irow += 1
        print(irow)
        dpp.append(dpp_line)
    return np.array(dpp)


# %%

def sim_one_ditch(L=100, D=10, dx=0.1, dy=0.1, b=5, h=3):
    HUGE = 1000.
    L, D = 20, 10
    dx, dy = 0.1, 0.1
    
    x = np.linspace(-b - dx, L + dx, int((L + b + 2*dx) / dx) + 1)
    z = np.linspace(dy, -D, int((D + dy) / dy + 1))
    gr = Grid(x, None, z)

    xy = ((-b - dx, dy), (-b -dx, -h), (0, -h), (0, dy), (-b -dx, dy))
    p = PathPatch(Path(xy), color='blue')

    DITCH = np.logical_and(gr.XM < 0, gr.ZM > -h)

    IBOUND = gr.const(1, dtype=int)
    IBOUND[:, :, 0] = 0
    IBOUND[0, :, :] = 0
    IBOUND[DITCH] = -1

    k = 1.
    K = gr.const(k)
    K[0, :, :] = 1e-8
    K[:, :, 0] = 1e-8
    K[DITCH] = HUGE
    K[1:, :, -1] = HUGE

    FQ = gr.const(0.)

    kdz = IBOUND[:, :, -1] * K[:, :, -1] * gr.DZ[: ,:, -1]
    FQ[:, :, -1] = -kdz / np.sum(kdz)

    HI = gr.const(0.)

    out= fdm3(gr, K=K, c=None, FQ=FQ, HI=HI, IBOUND=IBOUND, GHB=None)
    
    dpp_Q = dPP(out, k=k, L=L, D=D)
    
    dpp_Qa = 2 / (np.pi * k) * np.log(D / (2 * (b + h)))
    dpp_contr = dpp_contraction(h, D)
    
    S = strfun(gr, out)
    levels = np.linspace(S.min(), S.max(), 51)

    fig, ax = plt.subplots(figsize=(12, 7))
    fig.suptitle("Flow to ditch (no bottom resistance)")
    ax.set_title(f"h/D={h/D:.2f} m, b/D={b/D:.2f} m, D={D:.1f} m"
                + "\n"
                 + r"$\frac{d \phi}{Q}$"
                 + f"={dpp_Q:.3f} m/(m2/d_), "
                 + f"dpp_contr{dpp_contr:.3f}, "
                 + r"$\frac{2}{\pi} \ln(\frac{H}{\Omega})$"
                 + f"={dpp_Qa:.3f}")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    Cs = plt.contour(gr.xm, gr.zm, out['Phi'][:, 0, :], levels=-levels[::-1])
    Cs = plt.contour(gr.x[1:-1], gr.z, S, levels=levels)

    ax.add_patch(p)

    ax.set_aspect(1)
    ax.set_xlim(-b -dx, 8)
    
    parts = LPath(os.getcwd()).parts    
    home = os.path.join(*parts[:parts.index('GGOR') + 1], 'cases', 'AAN_GZK', )
    fig.savefig(os.path.join(home, "images", "ditch_isolines.png"))
    
    plt.show()

    return out

# %%
if __name__ == '__main__':


        ,     if True:
        table = compute_dpp_table(L=100, D=10, dx=0.1, dy=0.1, nx=10, ny=10)
        print(table)

        fname = os.path
        parts = LPath(os.getcwd()).parts    
        home = os.path.join(*parts[:parts.index('GGOR') + 1], 'cases', 'AAN_GZK', )
        pthnm = os.path.join(home, "data", "ditch_data.txt")

        with open(pthnm, 'w') as f:            
            f.write(str(table))
        print('Done', 'dpp table written to pthnm')

    # %%
    fig, ax = plt.subplots()

    C = ax.contour(table, levels=100)
    ax.clabel(C, levels=C.levels)


# %%
