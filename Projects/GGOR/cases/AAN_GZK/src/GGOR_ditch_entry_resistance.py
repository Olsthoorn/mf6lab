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
from dataclasses import dataclass


# %%

def get_home_folder():
    """Return home folder"""
    parts = LPath(os.getcwd()).parts
    if 'GGOR' not in parts:  
        raise FileNotFoundError("'GGOR not in folder tree!")
    return os.path.join(*parts[:parts.index('GGOR') + 1])

def dPP(out, k):
    """Return Dphi/Q due to partial penetration from fdm model run."""
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

@dataclass
class Xsec:
    """X-section objext.
    
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
    N: int
        Number of points in creating 
    """
    L: float
    D: float
    dx: float
    dy: float
    Q: float
    k: float
    N: int
    
    def __post_init__(self):
        assert self.L > 0
        assert self.D > 0
        assert self.dx > 0
        assert self.dy >= 0
        assert self.k > 0
        assert self.N > 10
    
    @property
    def kD(self) -> float:
        """Transmissivity kD."""
        return self.k * self.D
      
    def summary(self):
        """Return current aquifer properties as a dict."""
        return dict(
            L=self.L,
            D=self.D,
            dx=self.dx,
            dy=self.dy,
            Q=self.Q,
            k=self.k,
            N=self.N,
        )
    def __str__(self):
        return f"L={self.L}, D={self.D}, dx={self.dx}, dy={self.dy}, Q={self.Q}, k={self.k}, N={self.N}"


class DPP():
    
    def __init__(self, xsec):
        """Instantiate DPP_table object.
        
        Parameters
        ----------
        xsec: dataclass object
            cross section data.
        """
        self.ditch = ditch
        

    def get_dpp_table(self, pkl_file='dPP_table.pckl'):
        xsec = self.xsec
        try:
            # --- Pickle the dPP_table as a dictionary
            pfile = os.path.join(get_home_folder(), 'data', pkl_file)
            with open(pfile,'rb') as file:        
                dPP_table_dict = pickle.load(file)
                bs = dPP_table_dict['bs']
                hs = dPP_table_dict['hs']
                dPP_table = dPP_table_dict['dPP_table']
                fdPP = RegularGridInterpolator((hs/xsec.D, bs/xsec.D), dPP_table, method='cubic')
                return hs, bs, dPP_table, fdPP                
        except Exception:
            return self.generate_dPP_table()

    def generate_dpp_table(self):
        """Generate dPhi/Q table by numerical simulation."""         

        xsec = self.xsec

        HUGE = 1000.
        
        # --- for D is 10 and dx = dy = 0.1 this works to get a detailed table
        hs = xsec.D * np.unique(np.round(np.logspace(-2, 0, xsec.N + 1), 2))
        bs = xsec.D * np.unique(np.round(np.logspace(-2, np.log10(2), xsec.N + 1), 2))
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
                x = np.linspace(-b - xsec.dx,
                                xsec.L + xsec.dx,
                                int(np.round((xsec.L + b + 2*xsec.dx) / xsec.dx)) + 1)
                
                # --- z from dy to -D
                z = np.linspace(xsec.dy, -xsec.D, int(np.round((xsec.D + xsec.dy) / xsec.dy)) + 1)
                z = np.unique(np.hstack((-h, z)))[::-1] # include h
                
                if (h == 0) and (b == 0):
                    # --- Add a dummy value to prevent nan
                    b = 0.5 * xsec.dx
                    pb = 2 / (np.cosh(np.pi * b / xsec.D) - 1)
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
                K = gr.const(xsec.k)
                K[0, :, :] = 1e-8 # --- Inactive top row (allows zero depth ditch)
                K[:, :, 0] = 1e-8 # --- Inactive first col (allows zero width ditch)
                K[DITCH] = HUGE   # --- Ditch cells must not cause any resistance
                K[1:, :, -1] = HUGE # --- Right column below inact top row extractions

                # --- Extractions only in outer right column below inactive top row
                FQ = gr.const(0.)            
                # --- IBOUND used to keep FQ in inactive top row 0
                kdz = IBOUND[:, :, -1] * K[:, :, -1] * gr.DZ[: ,:, -1]
                FQ[:, :, -1] = xsec.Q * kdz / np.sum(kdz)

                # --- Initial heads
                HI = gr.const(0.)

                # --- Simulate and capture results in dict out
                out= fdm3(gr, K=K, c=None, FQ=FQ, HI=HI, IBOUND=IBOUND, GHB=None)
                out['gr'] = gr # --- Add gr to dict

                # --- Follow progress
                print('.', end="")

                dpp_line.append(dPP(out, xsec.k))
            irow += 1
            print()
            
            dpp.append(dpp_line)
        
        dPP_table = np.array(dpp)
        print("Done")
        
        # ---Round off all zeros for h==D
        dPP_table[-1] = 0
        
        # --- Top left corner (0,0) otherwise nan
        dPP_table[0, 0] = 0.5 * (dPP_table[0, 1] + dPP_table[1, 0])
        
        fdPP = RegularGridInterpolator((hs/xsec.D, bs/xsec.D), dPP_table, method='cubic')
        
        # --- Pickle the dPP_table as a dictionary
        pfile = os.path.join(get_home_folder(), 'data', 'dPP_table.pckl')
        with open(pfile,'wb') as f:
            dPP_table_dict = {'bs':bs/xsec.D, 'hs':hs/xsec.D, 'dPP_table':dPP_table}
            pickle.dump(dPP_table_dict, f)
        
        return hs/xsec.D, bs/xsec.D, dPP_table, fdPP
    
    
    def show_dPP_table(self):
        
        hs, bs, dPP_table, _ = self.get_dpp_table()                     
        fig, ax = plt.subplots(figsize=(8, 8))    
        ax.set_title(r"$d\Phi/Q$ (partial penetration (h/D, b/D)")
        ax.set_xlabel('b/D')
        ax.set_ylabel('h/D)')
        
        C = ax.contour(bs, hs, dPP_table, levels=np.linspace(0, 2.5, 26))
        ax.clabel(C, levels=C.levels[::5])        
        ax.invert_yaxis()
        ax.set_aspect(1)
        ax.grid(True, which='both')
        
        fig.savefig(os.path.join(get_home_folder(), 'doc', 'images', 'ditch_pp.png'))
        plt.show()
        return ax

    @lru_cache(maxsize=1)
    def _build_dPhi_Q_interpolator(self):
        """Return interpolator for DeltaPhi/Q."""
        _, _, _, interp = self.get_dpp_table()        
        return interp

    def dPhi_Q(self, hb):
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
        interp = self._build_dPhi_Q_interpolator()
        
        hb = np.atleast_2d(hb)
        
        # extract bounds from interpolator
        hmin, hmax = interp.grid[0][0], interp.grid[0][-1]
        bmin, bmax = interp.grid[1][0], interp.grid[1][-1]

        hb[:, 0] = np.clip(hb[:, 0], hmin, hmax)
        hb[:, 1] = np.clip(hb[:, 1], bmin, bmax)

        out = interp(hb)
        return out[0] if out.size == 1 else out
    

    def sim_one_ditch(self, hb=None, L=100, D=10, Q=-1, k=1):
        """Simulate and show the cross section for a single combination of b and h."""
        xsec = self.xsec
        xsec.L, xsec.D, xsec.Q, xsec.k = L, D, Q, k

        HUGE = 1000.
        
        if hb is None:
            hb = np.array([(0, 5), (5, 0), (3, 6), (2, 2)])
        hb = np.atleast_2d(hb)
        assert np.all(hb.T[0] >=0) and np.all(hb.T[1] >=0) and np.all(hb.T[0] <= xsec.D), f"not all 0<=h<={xsec.D} and/or not all b>=0"
        
        fig, axs = plt.subplots(2, 2, figsize=(14, 13))
        fig.suptitle("Flow to ditch (no bottom resistance)\n"
                        f"D = {xsec.D} m, k={xsec.k} m/d, Q={xsec.Q} m2/d"
                        )

        for ia, ((h, b), ax) in enumerate(zip(hb, axs.flatten())):
        
            Q = -1.0 # m^2/d
        
            if ia in [2, 3]:  
                ax.set_xlabel("x [m]")
            if ia in [0, 2]:
                ax.set_ylabel("y [m]")

            x = np.linspace(-b - xsec.dx, L + xsec.dx, int((L + b + 2*xsec.dx) / xsec.dx) + 1)
            z = np.linspace(xsec.dy, -xsec.D, int((xsec.D + xsec.dy) / xsec.dy + 1))
            gr = Grid(x, None, z)

            xy = ((-b - xsec.dx, xsec.dy), (-b -xsec.dx, -h), (0, -h), (0, xsec.dy), (-b -xsec.dx, xsec.dy))
            p = PathPatch(Path(xy), color='blue')

            DITCH = np.logical_and(gr.XM < 0, gr.ZM > -h)

            IBOUND = gr.const(1, dtype=int)
            IBOUND[:, :, 0] = 0
            IBOUND[0, :, :] = 0
            IBOUND[DITCH] = -1

            k = xsec.k
            K = gr.const(k)
            K[0, :, :] = 1e-8
            K[:, :, 0] = 1e-8
            K[DITCH] = HUGE
            K[1:, :, -1] = HUGE

            FQ = gr.const(0.)

            kdz = IBOUND[:, :, -1] * K[:, :, -1] * gr.DZ[: ,:, -1]
            FQ[:, :, -1] = xsec.Q * kdz / np.sum(kdz)

            HI = gr.const(0.)

            out= fdm3(gr, K=K, c=None, FQ=FQ, HI=HI, IBOUND=IBOUND, GHB=None)
            out['gr'] = gr
            
            dpp_Q = dPP(out, k=xsec.k)
            
            S = strfun(out)
            levels = np.linspace(S.min(), S.max(), 51)
        
            ax.set_title(f"h/D={h/xsec.D:.2f} m, b/D={b/xsec.D:.2f} m, "
                        fr"$d \phi={dpp_Q:.3f} \times Q$ m")
        
            ax.contour(gr.xm, gr.zm, out['Phi'][:, 0, :], levels=-levels[::-1])
            ax.contour(gr.x[1:-1], gr.z, S, levels=levels)

            ax.add_patch(p)

            ax.set_aspect(1)
            ax.set_xlim(-b -xsec.dx, 8)
        
        fig.savefig(os.path.join(get_home_folder(), "doc", "images", "ditch_isolines.png"))
        
        plt.show()

        return out

   
def fit_dPP_table():
    """Fit the dPP_table.
    
    Try calibrating the numerically computed dPP_table
    by fitting some function that resembles the analytic
    solutions along the top and left edge of the X-section.

    Several function have been tried, guided by ChatGPT
    but none was satisfactory.
    
    I finally gave up trying and now use the
    interpolation in the numerically derived
    table.
    
    So, this function is obsolete.    
    """

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
    dpp = DPP()
    hs, bs, dPP_table, _ = dpp.get_dpp_table()
    
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

def show_num_exact_vs_huisman():
    """Compare boundary cases with Huisman (1972, p57).
    
    Huisman (1972, p57) provides an approximation for
    DeltaPhi/Q for shallow ditches with arbitrary
    cross section.
    
    This function plots the exact solutions for
    the zero-depth ditch along the top of the X-section
    and the contraction along the left edge of the X-section
    together with the expression given by Huisman.

    Huisman uses D/Omega, for the flat ditch and the contraction
    Omega = 2b or 2h respectively, so that is used and is show
    to fit well only for b, h < 0.1D
    """
    def dppHuisman(D_over_Omega):
        """Return dPhi/Q acc to Huisman(1972)."""
        return 2 / np.pi * np.log(D_over_Omega)
    
    dpp = DPP()
    hs, bs, dPP_table, _= dpp.get_dpp_table()

    fig, ax = plt.subplots()
    ax.set_title(r"$\Delta\Phi/Q$")        
    ax.set(xlabel="b/D, h/D", ylabel=r"$\Delta\Phi/Q$")
    
    # --- Use the numerically computed values directly from the dPP_table
    ax.plot(bs, dPP_table[0], 'r.', label="Numeric hor.")
    ax.plot(hs, dPP_table[:, 0], 'b.', label="Numeric vert.")

    # --- Using the exact analytical solution
    pb = 2 / (np.cosh(np.pi * bs) - 1)
    ph = 1 / np.sin(np.pi * hs / 2) ** 2
    
    ax.plot(bs, bs + np.log(pb) / np.pi, 'r-', label="Analytic hor.")
    ax.plot(hs, np.log(ph) / np.pi, label="Analytic vert.")
    
    # --- Use Huismans approximation (is not very good).
    # ax.plot(bs, bs + dppHuisman(1/(2 * bs), label='Huisman horizontal')
    ax.plot(hs, dppHuisman(1/(2 * hs)), label='Huisman')
    
    ax.grid(True)
    ax.legend()
    
    fig.savefig(os.path.join(get_home_folder(), 'doc', 'images',
                             'dPP_analytic_vs_huisman.png'))


    

# %%
if __name__ == '__main__':
    
    xsec = Xsec(L=20, D=10, dx=0.1, dy=0.1, Q=-1, k=1, N=20)
    
    dpp=DPP(xsec)

    if False:
        Q, k, D, L = -1, 1, 10, 20
        
        h, b = 3, 6 # Height-width for single case
            
        out = dpp.sim_one_ditch(hb=None, L=L, k=k, Q=Q, D=D)
        
        hs, bs, dPP_table, fpp = dpp.get_dpp_table()
        
        dpp.show_dPP_table(bs=bs, hs=hs, dPP_table=dPP_table)
        
        plt.show()
    if False:        
        alpha, beta = fit_dPP_table()
    if False:
        hs, bs, dPP_table, fDpp = dpp.get_dpp_table()
    if True:
        dpp.show_num_vs_analytic_boundaries()
        
    plt.show()
    
    print("Done")

