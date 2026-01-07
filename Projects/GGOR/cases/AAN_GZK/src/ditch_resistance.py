# Ditch resistance conformal transformation

# %%
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

parts = Path(os.getcwd()).parts
images = os.path.join(*parts[:parts.index('GGOR') + 1], 'doc', 'images')

# %%
@dataclass
class Aquifer:
    Q: float
    k: float
    D: float # thickness of aquifer. Y assumed between 0>=y>=-D
    DL: float # Most left poin of ditch (e.g. 5 + D * 1j)
    DR: float # Most right point of ditch (e.g. 10 + D * 1j)

    def __post_init(self):
        assert self.k > 0
        assert self.D > 0
        assert (self.DR.real > self.DL.real) or (self.DR.imag > self.DL.imag)
        
    def __str__(self):        
        return f"Aquifer(Q={self.Q}, k={self.k}, D={self.D}, DL={self.DL}, DR={self.DR}"
    
class Ditches:        
    def __init__(self, aq: Aquifer)->None:
        self.aq = aq
        
    def zGrid(self, x: np.ndarray, y: np.ndarray)-> np.ndarray:
        """Return a regular grid of z-coordinates bound by aq.
    
        Parameters
        ---------
        x: np.ndarray of float
            grid x coordinates
        y: np.ndarray of float
            grid.z coordinates
        """
        aq = self.aq
        x = np.unique(np.hstack((x, aq.DL.real, aq.DR.real)))
        y.sort(); y = y[::-1]
        y = y[np.logical_and(y >= 0, y <= aq.D)]
        assert len(y) > 0, "No y-coordinates between 0 and -D."
        
        X, Y = np.meshgrid(x, y)
        return X + 1j * Y
        
    def omGrid(self, phi: np.array, psi:np.array) ->np.array:
        """Return Omega grid Phi + 1j * Psi.
        
        It is assumed that 0= < psi <= Q
        Q = 0 will be at the bottom, Q=1 at the top.
        """
        aq = self.aq
        psi = psi[np.logical_and(psi >= 0, psi <= aq.Q)]
        psi = psi[::-1]
        assert len(psi) > 0, "psi not compliant with 0 < psi < 1"        
        Phi, Psi = np.meshgrid(phi, psi)
        return Phi + 1j * Psi
    
    @property
    def pq(self)->tuple:
        aq = self.aq
        zDL = aq.DL
        zDR = aq.DR
        ztaL = self.zeta1(zDL)
        ztaR = self.zeta1(zDR)
        p = 2  /(ztaL - ztaR)
        q = -(ztaL + ztaR) / (ztaL - ztaR)
        return p, q
    
    def zeta(self, z: complex|np.ndarray)->complex|np.ndarray:
        p, q, = self.pq
        zta1 = self.zeta1(z)
        zta = p * zta1 + q
        return zta if zta.size > 1 else zta.item()
    
    def plot(self, zta, figsize=(10, 8), ax=None, **kwargs):
        """Plot any complex grid."""
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        ax.plot(zta.real, zta.imag, **kwargs)
        ax.plot(zta.real.T, zta.imag.T, **kwargs)
        ax.grid(True)
        return ax
        
    def contour(self, Z, omega, levels=20, figsize=(10, 8), ax=None, **kwargs):
        """Contour stream and potential lines."""
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        ax.contour(Z.real, Z.imag, omega.real, levels=levels, **kwargs)
        ax.contour(Z.real, Z.imag, omega.imag, levels=levels, **kwargs)
        ax.set_aspect(1)
        return ax


class Ditch_exp(Ditches):        

    def zeta1(self, z:complex | np.ndarray)->complex | np.ndarray:
        aq = self. aq
        z = np.atleast_1d(z)
        zta1 = np.exp(np.pi / aq.D *z)
        return zta1 if zta1.size > 1 else zta1.item()
    
    def zeta(self, z: complex|np.ndarray)->complex|np.ndarray:
        p, q, = self.pq
        zta1 = self.zeta1(z)
        zta = p * zta1 + q
        return zta if zta.size > 1 else zta.item()

    def zeta3(self, z:complex|np.ndarray)->complex|np.ndarray:
        zta = self.zeta(z)
        return np.arcsin(zta) + np.pi/2
    
    def Omega_A(self, z:complex|np.ndarray)->np.ndarray:
        aq = self.aq
        zta3 = self.zeta3(z)
        return -1j * aq.Q / np.pi * zta3
    
    def Omega(self, z:complex|np.ndarray)->np.ndarray:
        aq = self.aq        
        return -1j * aq.Q / np.pi * (
            np.arcsin(p * np.exp(np.pi / aq.D * z) + q + 0j) + np.pi/2)
    
    def z_fr_om(self, Omega:complex|np.ndarray)->complex|np.ndarray:
        aq = self.aq
        p, q = self.pq
        arg = -1j * np.pi / aq.Q * Omega -np.pi/2
        z = aq.D  / np.pi * np.log((np.sin(arg + 0j) - q) / p)
        return z

class Ditch_sin(Ditches):        

    def zeta1(self, z:complex | np.ndarray)->complex | np.ndarray:
        aq = self. aq
        z = np.atleast_1d(z)
        zta1 = np.sin(1j * np.pi / aq.D *z + np.pi /2)
        return zta1 if zta1.size > 1 else zta1.item()
    
    def zeta3(self, z:complex|np.ndarray)->complex|np.ndarray:
        zta = self.zeta(z)
        return np.arcsin(zta)
    
    def Omega_A(self, z:complex|np.ndarray)->np.ndarray:
        aq = self.aq
        zta3 = self.zeta3(z)
        return -1j * aq.Q / np.pi * (zta3 - np.pi / 2)
    
    def Omega(self, z:complex|np.ndarray)->np.ndarray:
        aq = self.aq
        p, q = self.pq   
        return -1j * aq.Q / np.pi * (
            np.arcsin(p * np.sin(1j * np.pi / aq.D * z + np.pi / 2) + q) - np.pi/2)
    
    def z_fr_om(self, Omega:complex|np.ndarray)->complex|np.ndarray:
        aq = self.aq
        p, q = self.pq
        arg = -1j * np.pi / aq.Q * Omega + np.pi/2
        z = -1j * aq.D  / np.pi * (np.arcsin((np.sin(arg) - q) / p)-np.pi / 2)
        return z
    
    # --- Symptotic behavior
    def asymptote(self, z):
        aq = self.aq
        p, _ = self.pq
        return aq.Q/ aq.D * z+ (aq.Q / np.pi) * np.log(p)
    

def show_sin_based(DL=0+5j, DR=10 + 10j, D=10, Q=1, k=1,
                   N=50, Nlevels=20, plot_what={'omega'}):
    """Stream and contour lines in a half-infinite X-section.
    
    X section: 
    1.     0 <= x <= infinity
    2.     9 <= y <= D

    1. The ditch is a section of the aquifer contour between
    2. DR and DL where the head is fixed at zero.
    3. The X-section runs from closed boundary at x to inifinity.
    4. A flow equalto Q [L2/T] is directed to the right.
    
    The edges if the ditch given by points DL and DR can be
    anywhere along any of theoutside edges of the X-section.
    The only condition is, that when looked at
    the cross section and following the edges in a
    clock-wise fashion, DR must be to the right of DL.
    DR left of DL at bottom DR above DL along left edge
    and DR right of DL at bottom. But DR and DL do not
    have to be both on the same edge.
    
    --------DL---- -> ----DR-------
    |
    DR
    |
    ^
    |
    DL
    |
    --------DR---- <- ----DL---------
    
    Also possible:
    
    ---------DR-------     ---------DR----
    |                      |
    |                      DL
    |                      |
    ---DL-------------     ----------------
    """
    
    # --- Start data defining the cross section with the ditch.
    if D is None:
        # --- D is not provide, assume DR at top of aquifer
        D = DR.imag
    else:
        # --- assert DL and DR are compatible with provided D
        assert (DR.imag == D) or (DR.real == 0) or (DR.imag == 0), (
            f"DR not on any of the edges of the X-section (x=0, y={D} y=0).")
        assert (DL.imag == D) or (DL.real == 0) or (DL.imag == 0), (
            f"DL not on any of the edges of the X-section (x=0, y={D} y=0).")
        assert np.angle((DR - (D + 0.5 * D * 1j)) / (DL - (D + 0.5 * D * 1j))) < 0, (
            "DR not clockwise of DL")
        
    # --- z-grid
    # --- Exp case
    x = np.linspace(-D, 2* D, 3 * N + 1)
    
    # --- Sin case
    x = x[x >=0].clip(1e-3)
    
    y = np.linspace(0, D, N + 1).clip(1e-3, D - 1e-3)

    # --- Define aquifer
    aq = Aquifer(Q=Q, k=k, D=D, DL=DL, DR=DR)

    # --- Instantiate the ditch
    ditch = Ditch_sin(aq)
    
    # --- Generate the Z-grid 
    Z = ditch.zGrid(x=x, y=y)

    # --- Compute the 𝛀 for this Z-grid
    Om = ditch.Omega(Z)
    
    # --- Plotting intermedate and final planes
    all_ = {'zeta1', 'zeta', 'zeta3', 'omega', 'omega_cont'}
    assert len(plot_what.difference(all_)) == 0, (
        f'what = {plot_what} is not subset of {all_}'
    )    
    # --- Zeta 1(sin transform, not yet shifted)
    if 'zeta1' in plot_what:
        ax = ditch.plot(ditch.zeta1(Z))
        ditch.plot(ditch.zeta1([aq.DL, aq.DR]), ax=ax, marker='o', mfc='r')
        ax.set_title(r"$\zeta_1$, $Z$-lines in the $\zeta_1$-plane")
        ax.set(xlabel=r'$\Re(\zeta_1)$', ylabel=r'$\Im(\zeta_1)\times i$', aspect=1)

    # --- Zeta (sin-tranform shifted)
    if 'zeta' in plot_what:
        ax = ditch.plot(ditch.zeta(Z))
        ditch.plot(ditch.zeta([aq.DL, aq.DR]), ax=ax, marker='o', mfc='r')
        ax.set_title(r"$\zeta$, $Z$-line in the $\zeta$-plane.")
        ax.set(xlabel=r'$\Re(\zeta)$', ylabel=r'$\Im(\zeta)\times i$', aspect=1)

    # --- arsin applied on zeta, just before rotating back to 𝛀
    if 'zeta3' in plot_what:
        ax = ditch.plot(ditch.zeta3(Z))
        ditch.plot(ditch.zeta3([aq.DL, aq.DR]), ax=ax, marker='o', mfc='r')
        ax.set_title(r"$\zeta_3$, $Z$-lines in the $\zeta_3$ plane")
        ax.set(xlabel=r'$\Re(\zeta_3)$', ylabel=r'$\Im(\zeta_3)\times i$', aspect=1)
    
    # --- Final 𝛀 plane (plotting lines of constant 𝛀)
    if 'omega' in plot_what:
        ax = ditch.plot(ditch.Omega(Z))
        ax.set_title(r"$\Omega(Z)$: Distored $Z$-grid lines in the $\Omega$ plane")
        ax.set(xlabel=r'$\Phi$', ylabel=r'$\Psi\times i$', aspect=1)

    # --- Same thing, but instead of plotting,
    # --- contouring 𝛀 on the distorted z-grid
    if 'omega' in plot_what:
        ax = ditch.contour(Z=Z, omega=ditch.Omega(Z), levels=Nlevels)
        ax.set_title(r"$\Omega$-contoured on the regular $Z$ in the $Z$-plane")
        ax.set(xlabel=r'$x$', ylabel=r'$iy$', aspect=1)

    if 'omega_cont' in plot_what:
        # --- Starting from Omega
        # --- First generate regular 𝛀 field/grid
        phi = np.linspace(0, 2 * Q,  2 * Nlevels + 1)
        psi = np.linspace(0, Q, Nlevels + 1).clip(1e-3, aq.Q - 1e-3)
        Om = ditch.omGrid(phi, psi)

        # --- We can then contour 𝛀 given a distorted z-grid generated from 𝛀
        ax = ditch.contour(Z=ditch.z_fr_om(Om), omega=Om, levels=Nlevels)
        ax.set_title(r"$\Omega$ contoured on a distorted $Z$-grid")
        ax.set(xlabel=r'$x$', ylabel=r'$iy$', aspect=1)

    if 'omega' in plot_what:
        # --- Or we just plot the z generated from 𝛀 directly
        ax = ditch.plot(ditch.z_fr_om(Om))
        ax.set_title(r"$z$-lines from $\Omega$ grid directly plotted on the $Z$-plane")
        ax.set(xlabel=r'$x$', ylabel=r'$iy$', aspect=1)

# %%
def resistance(DL=0+5j, DR=10 + 10j, D=10, Q=1, k=1, N=50, ax=None):
    """Show the resistance (drawdown along top and bottom of the X-section.
    """
    
    # --- Start data defining the cross section with the ditch.
    if D is None:
        # --- D is not provide, assume DR at top of aquifer
        D = DR.imag
    else:
        # --- assert DL and DR are compatible with provided D
        assert (DR.imag == D) or (DR.real == 0) or (DR.imag == 0), (
            f"DR not on any of the edges of the X-section (x=0, y={D} y=0).")
        assert (DL.imag == D) or (DL.real == 0) or (DL.imag == 0), (
            f"DL not on any of the edges of the X-section (x=0, y={D} y=0).")
        assert np.angle((DR - (D + 0.5 * D * 1j)) / (DL - (D + 0.5 * D * 1j))) < 0, (
            "DR not clockwise of DL")

    # --- Define aquifer
    aq = Aquifer(Q=Q, k=k, D=D, DL=DL, DR=DR)
        
    # --- z-grid
    # --- Exp case
    x = np.linspace(0, 3* D, 3 * N + 1).clip(1e-3)    
    y = np.linspace(0, D, 3).clip(1e-3, D - 1e-3)[-1:]

    # --- Instantiate the ditch
    ditch = Ditch_sin(aq)
    Z = ditch.zGrid(x, y)
    Omega =ditch.Omega(Z)
    
    # --- Delta Phi op x=DR.real
    p, _ = ditch.pq
    b = max(aq.DL.real, aq.DR.real)
    h = aq.D - max(aq.DL.imag, aq.DR.imag)
    
    # --- Extra drawdown due to partial penetration of ditch (of zero depth)
    dPhi = aq.Q / aq.D * b + aq.Q / np.pi * np.log(p)
    
    # --- Extra length to travel to get the same dPhi
    dL = b  + aq.D / np.pi * np.log(p)
    
    om_asymp = ditch.asymptote(Z)
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 10))        
        title1 = fr"D={D} m ditch width (b), ditch depth (h)={h}, Q={aq.Q} m2/d, k={aq.k} m/d"
        ax.set_title("Potential along top middle and bottom of X-section\n"
                 + title1
                 )
    ax.set(xlabel='x', ylabel='Phi')
    ax.grid(True)

        
    for z, omega in zip(Z, Omega):
        ax.plot(z.real, omega.real, label=f'drawdown at y={z[0].imag:.1f}, b={b}, h={h}')
    ax.plot(Z[0].real, om_asymp[0].real, '-.', label='Asymptotic drawdown')
    ax.plot([b, b], [0, dPhi], 'o-', lw=2, label=r'$d \Phi_{pp}$')
    ax.plot([b - dL, b], [0, dPhi], 'x--', color='k', lw=2, label=r'extra length $\Delta L$')
    
    ax.legend(loc='lower left', fontsize='small')
    
    ax.figure.savefig(os.path.join(images, "w_entry_dphi-dL.png"))
    
    return ax
    

def dPhi_PP(Q=1, k=1, D=10, l_D=None, axis=0):
    """Return dPhi/Q along axis
    """
    l_D = l_D[l_D > 0]
    if axis==0:
        l_D = l_D[l_D <= 1]

    dPhi_Q = []
    for ld in l_D:
        if axis == 0:            
            DR = 0  + 1j * D
            DL = DR - 1j * D * ld
            bD = 0                        
        elif axis==1:
            DL = 0  + D * 1j
            DR = DL + ld * D
            bD = ld
        else:
            raise ValueError("axis must be 0 or 1")
        
        aq = Aquifer(Q=Q, k=k, D=D, DL=DL, DR=DR)
    
        # --- Instantiate the ditch
        ditch = Ditch_sin(aq)
    
        # --- Delta Phi op x=DR.real
        p, _ = ditch.pq
        print(f" {p.real:.4g}", end="")
                    
        # --- Extra drawdown due to partial penetration of ditch (of zero depth)
        dPhi_Q.append(np.real(bD + np.log(p) / np.pi))
        
    print(f", axis={axis}")
    
    return l_D, np.array(dPhi_Q)
    
def get_dphi_Q(Q=-1, k=1, D=10, l_D=None):
    bs_D, dPhi_Q_x = dPhi_PP(Q=Q, k=k, D=D, l_D=l_D, axis=1)
    hs_D, dPhi_Q_y = dPhi_PP(Q=Q, k=k, D=D, l_D=l_D, axis=0)

    fig, ax = plt.subplots()
    ax.plot(bs_D, dPhi_Q_x, label=r"dPP$_x=d\Phi/Q$ horizontal ditch")
    ax.plot(hs_D, dPhi_Q_y, label=r"dPP$_y=d\Phi/Q$ vertical ditch")
    ax.grid(True)
    ax.set_title(
        r"$\frac{\Delta \Phi}{Q}$ due to partial penetration: "        
        r"$\frac{\Delta\Phi}{Q}=\frac{k \, \Delta\phi}{Q}=\frac{\Delta L}{D}$"
        "\n"
        f"D={D}, Q={Q}, k={k}")
    ax.set(xlabel='h/D or b/d', ylabel=r"$\frac{d\Phi}{Q}=\frac{\Delta L}{D}$")
    ax.legend()
    ax.figure.savefig(os.path.join(images, "dPhi_pp_Q_hor_vert.png"))


# %%
if __name__ == '__main__':
    if False:
        # --- Examples of placing the disk at any position along the edges  
        show_sin_based(DL=0+5j, DR=10+10j, D=10, Nlevels=40, plot_what={'omega_cont'})
        show_sin_based(DL=5+10j, DR=10+10j, D=10, Nlevels=40, plot_what={'omega_cont'})
        show_sin_based(DL=0+7j, DR=0+8j, D=10, Nlevels=40, plot_what={'omega_cont'})
        show_sin_based(DL=2+0j, DR=0+2j, D=10, Nlevels=40, plot_what={'omega_cont'})
    if False:
        resistance(DL=0+10j, DR=10 + 10j, D=10, Q=-1, k=1, N=50)
        resistance(DL=0+5j, DR=3 + 10j, D=10, Q=-1, k=1, N=50)
        resistance(DL=0+6j, DR=0 + 8j, D=10, Q=-1, k=1, N=50)
        resistance(DL=3+0j, DR=0 + 2j, D=10, Q=-1, k=1, N=50)
    if True:
        D = 10
        DL = 0 + D * 1j
        ax = resistance(DL=DL, DR=1 + D * 1j, D=D, Q=-1, k=1, N=50)
        resistance(DL=DL, DR=2 + D * 1j, D=D, Q=-1, k=1, N=50, ax=ax)
        resistance(DL=DL, DR=5 + D * 1j, D=D, Q=-1, k=1, N=50, ax=ax)
        resistance(DL=DL, DR=10 + D * 1j, D=D, Q=-1, k=1, N=50, ax=ax)
        resistance(DL=DL, DR=15 + D * 1j, D=D, Q=-1, k=1, N=50, ax=ax)
        resistance(DL=DL, DR=20 + D * 1j, D=D, Q=-1, k=1, N=50, ax=ax)
    if True:        
        get_dphi_Q(Q=1, k=1, D=10, l_D=np.logspace(-2, np.log10(1.5), 50))
        
    plt.show()
    print("Done")
    # %%