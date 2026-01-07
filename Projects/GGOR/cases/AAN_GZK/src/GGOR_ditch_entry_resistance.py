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
import pickle
from scipy.optimize import least_squares
from itertools import cycle
from functools import lru_cache


# %%

def get_home_folder():
    """Return home folder"""
    parts = LPath(os.getcwd()).parts
    if 'GGOR' not in parts:  
        raise FileNotFoundError("'GGOR not in folder tree!")
    return os.path.join(*parts[:parts.index('GGOR') + 1], 'cases', 'AAN_GZK', )

def dPP(out, k):
    """Return Dphi/Q due to partial penetration."""
    phi0 = 0. # Head in the ditch
    
    # --- head at x=L
    phiL = out['Phi'][1:, 0, -1].mean()
    
    # --- Distance from ditch to extraction
    L = out['gr'].x[-2] - out['gr'].x[out['gr'].x <= 0][-1]
    
    # --- Total extraction
    Q = out['Q'][1:, 0, -1].sum()
    
    # --- Aquifer thickness
    D = out['gr'].DZ[1:, 0, -1].sum()
    
    # --- Dpp/Q
    dPhi_Q = k * (phiL - phi0) / Q - L / D    
    return dPhi_Q 

def dpp_contraction(h, D):
    """Rreturn dpp by contraction according to Verruijt(1970) p119."""
    return - 2 / np.pi * np.log(np.sin(np.pi * h / (2 * D)))

def strfun(out):
    """Return the stream function."""
    S = np.zeros((out['gr'].nz + 1, out['gr'].nx-1))
    S[:-1] = np.cumsum(out['Qx'][::-1, 0, :], axis=0)[::-1]
    return S

def compute_dpp_table(L=20, D=10, dx=0.1, dy=0.1, Q=-1, k=1, N=20):
    """Return dPhi_pp/Q for given grid.
    Parameters
    ----------
    L, D: scalars
        Length and Thickness of the aquifer.
    dx, dy: scalars
        Cell widths of the FDM mesh    
    Q: float
        total extraction
    k: float
        Hydraulic conductivity
    """
    HUGE = 1000.
    
    # --- for D is 10 and dx = dy = 0.1 this works to get a detailed table
    hs = D * np.unique(np.round(np.logspace(-2, 0, N + 1), 2))
    bs = D * np.unique(np.round(np.logspace(-2, np.log10(2), N + 1), 2))
    hs = np.hstack((0, hs))
    bs = np.hstack((0, bs))

    irow = 0
    dpp = []
    print(f"\nTable will have shape ({len(hs)}, {len(bs)})")
    for h in hs:
        print(f"{irow:>2}: ", end="") # --- Follow progress
        dpp_line = []        
        for b in bs:
            
            # --- x from -b - dx to L + dx
            x = np.linspace(-b - dx, L + dx, int(np.round((L + b + 2*dx) / dx)) + 1)
            
            # --- z from dy to -D
            z = np.linspace(dy, -D, int(np.round((D + dy) / dy)) + 1)
            z = np.unique(np.hstack((-h, z)))[::-1] # include h
            
            if (h == 0) and (b == 0):
                # --- Add a dummy value to prevent nan
                b = 0.5 * dx
                pb = 2 / (np.cosh(np.pi * b / D) - 1)
                dpp_ = b + np.log(pb) / np.pi
                dpp_line.append(dpp_)
                continue
            
            gr = Grid(x, None, z)

            DITCH = np.logical_and(gr.XM < 0, gr.ZM > -h)

            # --- Fixed heads and inactive cells
            IBOUND = gr.const(1, dtype=int)
            IBOUND[:, :, 0] = 0
            IBOUND[0, :, :] = 0
            IBOUND[DITCH] = -1
            
            # --- Conductivities in and outside the ditch
            K = gr.const(k)
            K[0, :, :] = 1e-8 # --- Inactive top row (allows zero depth ditch)
            K[:, :, 0] = 1e-8 # --- Inactive first col (allows zero width ditch)
            K[DITCH] = HUGE   # --- Ditch cells must not cause any resistance
            K[1:, :, -1] = HUGE # --- Right column below inact top row extractions

            # --- Extractions only in outer right column below inactive top row
            FQ = gr.const(0.)            
            # --- IBOUND used to keep FQ in inactive top row 0
            kdz = IBOUND[:, :, -1] * K[:, :, -1] * gr.DZ[: ,:, -1]
            FQ[:, :, -1] = Q * kdz / np.sum(kdz)

            # --- Initial heads
            HI = gr.const(0.)

            # --- Simulate and capture results in dict out
            out= fdm3(gr, K=K, c=None, FQ=FQ, HI=HI, IBOUND=IBOUND, GHB=None)
            out['gr'] = gr # --- Add gr to dict

            # --- Follow progress
            print('.', end="")

            dpp_line.append(dPP(out, k))
        irow += 1
        print()
        
        dpp.append(dpp_line)
    
    dPP_table = np.array(dpp)
    print("Done")
    
    # ---Round off all zeros for h==D
    dPP_table[-1] = 0
    
    # --- Top left corner (0,0) otherwise nan
    dPP_table[0, 0] = 0.5 * (dPP_table[0, 1] + dPP_table[1, 0])
    
    fdPP = RegularGridInterpolator((hs/D, bs/D), dPP_table, method='cubic')
    
    # --- Pickle the dPP_table as a dictionary
    pfile = os.path.join(get_home_folder(), 'data', 'dPP_table.pckl')
    with open(pfile,'wb') as f:
        dPP_table_dict = {'bs':bs/D, 'hs':hs/D, 'dPP_table':dPP_table}
        pickle.dump(dPP_table_dict, f)
    # Load the pickled file
    # with open(pfile,'rb') as f:        
    #     dPP_table_dict = pickle.load(f)
    
    return hs/D, bs/D, dPP_table, fdPP


def sim_one_ditch(L=100, Q=None, k=None, D=10, dx=0.1, dy=0.1, b=None, h=None):
    """Simulate and show the cross section for a single combination of b and h."""
    HUGE = 1000.
    L, D = 20, 10
    dx, dy = 0.1, 0.1
    
    hb = ((0, 5), (5, 0), (3, 6), (2, 2))
    
    fig, axs = plt.subplots(2, 2, figsize=(14, 13))
    fig.suptitle("Flow to ditch (no bottom resistance)\n"
                    f"D = {D} m, k={k} m/d, Q={Q} m2/d"
                    )

    for ia, ((h, b), ax) in enumerate(zip(hb, axs.flatten())):
    
        Q = -1.0 # m^2/d
       
        if ia in [2, 3]:  
            ax.set_xlabel("x [m]")
        if ia in [0, 2]:
            ax.set_ylabel("y [m]")

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
        FQ[:, :, -1] = Q * kdz / np.sum(kdz)

        HI = gr.const(0.)

        out= fdm3(gr, K=K, c=None, FQ=FQ, HI=HI, IBOUND=IBOUND, GHB=None)
        out['gr'] = gr
        
        dpp_Q = dPP(out, k=k)
        
        S = strfun(out)
        levels = np.linspace(S.min(), S.max(), 51)
    
        ax.set_title(f"h/D={h/D:.2f} m, b/D={b/D:.2f} m, "
                    fr"$d \phi={dpp_Q:.3f} \times Q$ m")
    
        ax.contour(gr.xm, gr.zm, out['Phi'][:, 0, :], levels=-levels[::-1])
        ax.contour(gr.x[1:-1], gr.z, S, levels=levels)

        ax.add_patch(p)

        ax.set_aspect(1)
        ax.set_xlim(-b -dx, 8)
    
    fig.savefig(os.path.join(get_home_folder(), "images", "ditch_isolines.png"))
    
    plt.show()

    return out

def show_dPP_table(bs=None, hs=None, dPP_table=None, fPP=None):
    """Contour the dPP_table.
    
    Parameters
    ----------
    bs: ditch widths
    hs: ditch depths
    dPP_table: np.ndarray | None
        Table with dPhi_pp/Q values
    dPP: RegularInterpolator | None
        Interpolator from which the dPP_table can be regenerated
    """
    
    assert not ((dPP_table is None) and (fPP is None)), (
        "dPP_table and fPP must not both be None."
    )

    # --- If dPP_table is None the interpolator must be given
    # ----and we generate the table
    if dPP_table is None:
        fPP = RegularGridInterpolator((hs, bs), dPP_table, method='cubic')
        dPP_table = np.array([fPP((h, b)) for h in hs for b in bs])
        dPP_table = dPP_table.reshape(len(hs), len(bs))
        dPP_table[0, 0] = np.nan
                         
    fig, ax = plt.subplots(figsize=(8, 8))    
    fig.suptitle("Partial penetration of a rectangular ditch in half infinite field")
    ax.set_title("dphi/Q (partial penetration (h/D, b/D) of rectangular ditch)")
    ax.set_xlabel('b/D')
    ax.set_ylabel('h/D)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    C = ax.contour(bs, hs, dPP_table, levels=25)
    ax.clabel(C, levels=C.levels)
    
    ax.set_xlim(0.1, 20)
    ax.set_ylim(0.1, 10)
    ax.invert_yaxis()
    ax.set_aspect(1)
    ax.grid(True, which='both')
    
    fig.savefig(os.path.join(get_home_folder(), 'images', 'ditch_pp.png'))
    plt.show()
    return ax

def fit_dPP_table():
    """Fit the dPP_table"""

    def f(params, bs, hs):
        """Return approx dPP/Q.
        
        Parameters
        ----------
        params: tuple
            parameters to calibrate
        bs: np.ndarray of real (0 < bs)
            Ditch width devided by aquifer thickness.
        hs: np.ndarray of real (0 < hs <= 1)
            Vertical gap width divied by aquifer thickness        
        """
        alpha, beta, gamma = params
        Rb = np.cosh(np.pi * bs[None, :]) - 1
        Rh = np.sin(np.pi * hs[:, None] / 2)**2
        p = 1 / (alpha * Rb + beta * Rh + gamma * Rb * Rh)        
        return bs + np.log(p) / np.pi
        
    def residuals(params, bs, hs):
        r = dPP_table - f(params, bs, hs)
        return r.flatten()[1:]   # drop the (hs≈0, bs≈0) corner

    # --- Pickle the dPP_table as a dictionary
    pfile = os.path.join(get_home_folder(), 'data', 'dPP_table.pckl')
    with open(pfile,'rb') as file:        
        dPP_table_dict = pickle.load(file)
        b = dPP_table_dict['bs']
        h = dPP_table_dict['hs']
        dPP_table = dPP_table_dict['dPP_table']
        D = h.max()
        bs, hs = b/D, h/D

    # Initial guess (important!)
    p0 = np.array([1.0, 1.0, 1.0])

    res = least_squares(
        residuals,
        p0,
        args = [bs, hs],
        bounds=([0, 0, 0], [np.inf, np.inf, np.inf])  # all must be positive
    )

    alpha, beta, gamma = res.x

    print("alpha, beta, gamma =", alpha, beta, gamma)
    print("RMS error =", np.sqrt(np.mean(res.fun**2)))

    fig, axis = plt.subplots()
    axis.plot(dPP_table.flatten(), f([alpha, beta, gamma], bs, hs).flatten(), '.')
    axis.set_xlabel("numerical")
    axis.set_ylabel("model")
    axis.axis('equal')
    axis.grid(True)
    
    clrs = cycle('brgkmc')
    fig, axis = plt.subplots()
    for ih in [0, 5, 10, 15, 20, 25]:
        clr = next(clrs)
        axis.plot(bs, dPP_table[ih], 'o', color=clr, label=f'ih={ih}')
        axis.plot(bs, f([alpha, beta, gamma], bs, hs)[ih], '.-', color=clr, label=f"ih={ih}")
    axis.set_title("Horizontal")
    axis.set_xlabel("b")
    axis.set_ylabel("dPP/Q")
    axis.grid(True)
    axis.legend()

    clrs = cycle('brgkmc')
    fig, axis = plt.subplots()
    for ib in [0, 5, 10, 15, 20]:
        clr = next(clrs)
        axis.plot(hs, dPP_table.T[ib], 'o', color=clr, label=f'ib={ib}')
        axis.plot(hs, f([alpha, beta, gamma], bs, hs).T[ib], '.-', color=clr, label=f"ib={ib}")
    axis.set_title("Vertical")
    axis.set_xlabel("h")
    axis.set_ylabel("dPP/Q")
    axis.grid(True)
    axis.legend()
    
    return alpha, beta

def show_num_vs_analytic_boundaries():
    """Compare boundary cases."""

    def f(params, bs, hs):
        """Return approx dPP/Q.
        
        Parameters
        ----------
        params: tuple
            parameters to calibrate
        bs: np.ndarray of real (0 < bs)
            Ditch width devided by aquifer thickness.
        hs: np.ndarray of real (0 < hs <= 1)
            Vertical gap width divied by aquifer thickness        
        """
        alpha, beta = params
        return bs[None, :] - 1 / np.pi * np.log(
            alpha * np.sin(np.pi * hs[:, None]/2)**2
            + beta*(np.cosh(np.pi * bs[None, :]) - 1)
        )

    # --- Pickle the dPP_table as a dictionary
    pfile = os.path.join(get_home_folder(), 'data', 'dPP_table.pckl')
    with open(pfile,'rb') as file:        
        dPP_table_dict = pickle.load(file)
        bs = dPP_table_dict['bs']
        hs = dPP_table_dict['hs']
        dPP_table = dPP_table_dict['dPP_table']

    fig, ax = plt.subplots()
    ax.set_title(r"$\Delta\Phi/Q$")
    
    ax.set(xlabel="b/D", ylabel="h/D")    
    ax.plot(bs, dPP_table[0], 'r.', label="Numeric hor.")
    ax.plot(hs, dPP_table[:, 0], 'b.', label="Numeric vert.")

    pb = 2 / (np.cosh(np.pi * bs) - 1)
    ph = 1 / np.sin(np.pi * hs / 2) ** 2
    
    ax.plot(bs, bs + np.log(pb) / np.pi, 'r-', label="Analytic hor.")
    ax.plot(hs, np.log(ph) / np.pi, label="Analytic vert.")
    
    # --- the analytic approximation in Huisman (1972, p57) is not very good
    # ax.plot(bs, np.log(1 / (2 * bs)), label='Huis hor')
    # ax.plot(hs, np.log(1 / (2 * hs)), label='Huis vert')
    
    ax.grid(True)
    ax.legend()


@lru_cache(maxsize=1)
def _build_dPhi_Q_interpolator():
    """Return interpolator for DeltaPhi/Q."""

    # --- Get the pickled dPP_table as a dictionary
    pfile = os.path.join(get_home_folder(), 'data', 'dPP_table.pckl')
    
    with open(pfile,'rb') as file:        
        dPP_table_dict = pickle.load(file)
        bs = dPP_table_dict['bs']
        hs = dPP_table_dict['hs']
        dPP_table = dPP_table_dict['dPP_table']
        
    return RegularGridInterpolator((hs, bs), dPP_table, method='cubic')


def dPhi_Q(hb):
    """Get ΔΦ/Q by interpolation.
    
    Parameters
    ----------
    hb : tuple or sequence of tuples (h_D, b_D)
     h_D: relative ditchdepth (h/D)
        ditch depth  0 <= h <= 1
     b_D: relative ditch width (b/D)
        ditch width b_D >= 0
        For values b_D > 2, b_D = 2 is used
    """
    hb = np.atleast_2d(hb)

    interp = _build_dPhi_Q_interpolator()

    # extract bounds from interpolator
    hmin, hmax = interp.grid[0][0], interp.grid[0][-1]
    bmin, bmax = interp.grid[1][0], interp.grid[1][-1]

    hb[:, 0] = np.clip(hb[:, 0], hmin, hmax)
    hb[:, 1] = np.clip(hb[:, 1], bmin, bmax)

    out = interp(hb)
    return out[0] if out.size == 1 else out
    

# %%
if __name__ == '__main__':

    if False:
        Q, k, D, L = -1, 1, 10, 20
        
        h, b = 3, 6 # Height-width for single case
            
        # out = sim_one_ditch(L=L, k=k, Q=Q, D=D, dx=0.1, dy=0.1, b=b, h=h)
        
        hs, bs, dPP_table, fpp = compute_dpp_table(L=L, D=D, dx=0.1, dy=0.1, Q=Q, k=k)
        
        show_dPP_table(bs=bs, hs=hs, dPP_table=dPP_table)
        plt.show()
        
        print("Done")
    
    # alpha, beta = fit_dPP_table()
    
    
    hs, bs, dPP_table, fDpp = compute_dpp_table(L=20, D=10, dx=0.1, dy=0.1, Q=-1, k=1)
    
    show_num_vs_analytic_boundaries()
    plt.show()
    
    print("Done")

