# Ditch resistance conformal transformation

# %%
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

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
        assert (self.DR.real > self.DL.real) or (DR.imag > DL.imag)
        
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
        return -1j * aq.Q / np.pi * (
            np.arcsin(p * np.sin(1j * np.pi / aq.D * z + np.pi / 2) + q) - np.pi/2)
    
    def z_fr_om(self, Omega:complex|np.ndarray)->complex|np.ndarray:
        aq = self.aq
        p, q = self.pq
        arg = -1j * np.pi / aq.Q * Omega + np.pi/2
        z = -1j * aq.D  / np.pi * (np.arcsin((np.sin(arg) - q) / p)-np.pi / 2)
        return z
    
# %%
Q, D = 1., 10.
DL, DR =0 + 0.5 * D * 1j, 5 + D * 1j

x = np.linspace(-D, 2* D, 151)
x = x[x >=0].clip(1e-3)
y = np.linspace(0, D, 51).clip(1e-3, D - 1e-3)

aq = Aquifer(Q=Q, k=1, D=10, DL=DL, DR=DR)

ditch = Ditch_sin(aq)
Z = ditch.zGrid(x=x, y=y)

p, q = ditch.pq

Om = ditch.Omega(Z)

# %%
ax = ditch.plot(ditch.zeta1(Z))
ditch.plot(ditch.zeta1([aq.DL, aq.DR]), ax=ax, marker='o', mfc='r')
ax.set_title("zeta 1")

ax = ditch.plot(ditch.zeta(Z))
ditch.plot(ditch.zeta([aq.DL, aq.DR]), ax=ax, marker='o', mfc='r')
ax.set_title("zeta")

ax = ditch.plot(ditch.zeta3(Z))
ditch.plot(ditch.zeta3([aq.DL, aq.DR]), ax=ax, marker='o', mfc='r')
ax.set_title("zeta 3")

ax = ditch.plot(ditch.Omega(Z))
ax.set_title("Omega")

ax = ditch.contour(Z=Z, omega=ditch.Omega(Z), levels=20)

phi = np.linspace(0, 2 * Q,  21)
psi = np.linspace(0, Q, 11).clip(1e-3, aq.Q - 1e-3)
Om = ditch.omGrid(phi, psi)
ax = ditch.contour(Z=ditch.z_fr_om(Om), omega=Om, levels=20)
ax.set_title("Omega, contours")

plt.show()
print("Don")