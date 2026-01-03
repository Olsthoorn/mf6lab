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
    xL: float # assumed xL.imag = 0
    xR: float # assumed xR.imag = 0

    def __post_init(self):
        assert self.k > 0
        assert self.D > 0
        assert self.xR > self.xL
        
    @property
    def DL(self):
        """Right point of ditch."""
        return self.xL + self.D * 1j
    
    @property
    def DR(self):
        """Left point of ditch."""
        return self.xR + self.D * 1j
    
    def __str__(self):        
        return f"Aquifer(Q={self.Q}, k={self.k}, D={self.D}, xL={self.xL}, xR={self.xR}"
    
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
        x = np.unique(np.hstack((x, aq.xL, aq.xR)))
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
Q, D, xL, xR =1., 10, 0, 5

x = np.linspace(-D, 2* D, 151)
x = x[x >=0].clip(1e-3)
y = np.linspace(0, D, 51).clip(1e-3, D - 1e-3)

aq = Aquifer(Q=Q, k=1, D=10, xL=xL, xR=xR)

ditch = Ditch_sin(aq)
Z = ditch.zGrid(x=x, y=y)

p, q = ditch.pq

Om = ditch.Omega(Z)

# %%
fig, ax = plt.subplots()
zta1 = ditch.zeta1(Z)
zta1pnts = ditch.zeta1([aq.DL, aq.DR])
ax.set_title("zeta 1")
ax.plot(zta1.real, zta1.imag)
ax.plot(zta1.real.T, zta1.imag.T)
ax.plot(zta1pnts.real, zta1pnts.imag, 'ro')
ax.grid(True)

fig, ax = plt.subplots()
zta = ditch.zeta(Z)
ztapnts = ditch.zeta([aq.DL, aq.DR])
ax.set_title("zeta")
ax.plot(zta.real, zta.imag)
ax.plot(zta.real.T, zta.imag.T)
ax.plot(ztapnts.real, ztapnts.imag, 'ro')
ax.grid(True)

fig, ax = plt.subplots()
zta = ditch.zeta3(Z)
ztapnts = ditch.zeta3([aq.DL, aq.DR])
ax.set_title("zeta 3")
ax.plot(zta.real, zta.imag)
ax.plot(zta.real.T, zta.imag.T)
ax.plot(ztapnts.real, ztapnts.imag, 'ro')
ax.grid(True)

fig, ax = plt.subplots()
ax.set_title("Omega")
Om = ditch.Omega(Z)
ax.plot(Om.real, Om.imag)
ax.plot(Om.real.T, Om.imag.T)
ax.grid(True)

fig, ax = plt.subplots()
ax.set_title("Omega, contours")
Om = ditch.Omega(Z)

ax.contour(Z.real, Z.imag, Om.real, levels=20)
ax.contour(Z.real, Z.imag, Om.imag, levels=20)

fig, ax = plt.subplots()
ax.set_title("z from Omega, direct")
Om = ditch.Omega(Z)

phi = np.linspace(0, 2 * Q,  21)
psi = np.linspace(0, Q, 11).clip(1e-3, aq.Q - 1e-3)
Om = ditch.omGrid(phi, psi)
Z1 = ditch.z_fr_om(Om)
ax.plot(Z1.real, Z1.imag)
ax.plot(Z1.real.T, Z1.imag.T)
ax.grid(True)

plt.show()
print("Don")