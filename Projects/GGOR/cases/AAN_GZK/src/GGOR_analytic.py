"""This file implements analytic solution that are and can be used with GGOR.

The GGOR (desired groundwater and surface water regime) can be simulated
numerically with the GGOR-tool and analytically. The analytical
implementation of the simulators are in this module.

@TO 2025-12-19
"""
# %% Imports
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from abc import ABC, abstractmethod
import pandas as pd
from itertools import cycle
import ggor_meteo

# %%

@dataclass
class Aquifer:
    k: float
    D: float
    c: float
    mu: float
    b: float
    
    def __post_init__(self):
        assert self.k > 0
        assert self.D > 0
        assert self.c > 0
        assert self.mu > 0
        assert self.b > 0
    
    @property
    def kD(self) -> float:
        """Transmissivity kD."""
        return self.k * self.D
    
    @property
    def lam(self) -> float:
        """Spreading length λ = √(k D c)."""
        return np.sqrt(self.k * self.D * self.c)
    
    @property
    def T(self)->float:
        """Characteristic time of leaky top aquifer.groundwater system.
        
        T = mu c (b / lambda) * coth(b / lambda)        
        """
        L = self.lam
        return self.mu * self.c * self.b / L * __class__.coth(self.b / L)
    
    @property
    def G(self)->float:
        """Factor popping up in the steady case derivation.
        
        G = b / Lambda * coth( b / lambda) - 1        
        """
        L = self.lam
        return self.b / L * __class__.coth(self.b / L) - 1
    
    @classmethod
    def coth(cls, x:float | np.ndarray)->float | np.ndarray:
        """Return coth (not in numpy or scipy)"""
        return 1 / np.tanh(x)
        
    def summary(self):
        """Return current aquifer properties as a dict."""
        return dict(
            k=self.k,
            D=self.D,
            c=self.c,
            mu=self.mu,
            b=self.b,
            kD = self.kD,
            lam=self.lam,
            T=self.T,
            G=self.G,
        )

aq = Aquifer(k=10.0, D=5.0, c=0.2, mu=0.15, b=50)
print(aq)
print(aq.summary())

@dataclass
class Tdata:
    tdata: pd.DataFrame

    def __post_init__(self):
        assert pd.api.types.is_datetime64_any_dtype(self.tdata.index)
        assert 'RH' in self.tdata.columns
        assert 'EV24' in self.tdata.columns

        self.tdata['R'] = self.tdata['RH'] - self.tdata['EV24']
        self.tdata = self.tdata[['RH', 'EV24', 'R']]

# %%

class AnalyticalSolution(ABC):
    
    def __init__(self, aq: Aquifer) -> None:
        """Instantiate analytic simulator."""
        self.aq = aq
        
    @abstractmethod
    def steady(self, Nx: int, phi: float, hLR: float, R: float, **kwargs)->tuple:
        """Return steady-state solution (h, x)"""
        pass

    
    @abstractmethod
    def transient(self, time_props):
        """Return transient solufion()"""
        pass        

# %% Implementation base case

class Base_case(AnalyticalSolution):
    
    def __init__(self, aq: Aquifer)-> None:
        """Instantiate analytic simulator.
        
        Parameters
        ----------
        aq: Aquifer object
            Aquifer properties (see Aquifer)
        """
        self.aq = aq        


    def steady(self, x: int, phi: float, hLR: float, R: float, **kwargs)->tuple:
        """Return steady state solution of cross section x.
        
        x: float | np.ndarray [m] 
            x-coordinates  -b <= x <= b            
        phi: float [L]
            uniform head in regional aquifer
        hLR: float [L]
            water level in the ditch
        R: float [L/T]
            recharge rate
        """
        aq = self.aq
        x = np.atleast_1d(x)    
        h = (phi + R * aq.c) + (hLR - (phi + R * aq.c)) * (
            np.cosh(x/aq.lam) / np.cosh(aq.b/aq.lam)
        )
        return h if len(h) > 1 else h.item()
    
    def steady_avg(self, phi: float, hLR: float, R: float, **kwargs)->float:
        """Return the averages steady-state head.
        
        >>> aq = Aquifer()
        >>> x = np.linspace(-aq.b, +aq.b, 100)
        >>> bc = Base_case(aq)
        >>> h = bc.Steady(x=x, phi=0., hLR=0., R=0.001)
        >>> havg = bc.steady_avg(phi=0., hLR=0., R=0.001)
        >>> np.isclose(havg, np.mean(h))
        >>> True
        """
        aq = self.aq
        return (phi + R * aq.c) + (hLR - (phi + R * aq.c)) * (
            np.tanh(aq.b/aq.lam) / (aq.b/aq.lam)
        )

    def transient(self, rch, h_summer=None, h_winter=None, q=None):
        """Return result of dynamic simulation.

        Parameters
        ----------
        rch: pd.Series [m/d]
            recharge with index = pd.Timestamps
        h_summer, h_winter: floats
            summer and winter ditch levels
        q: float or series:
            upward seepage rate
        """        
        tdata = pd.DataFrame(rch, columns=['rch'])
        # --- Add column with ditch water level according to summer or winter
        
        summer = np.logical_and(tdata.index.month >= 4, tdata.index.month <= 9)
        tdata['hLR'] = h_winter
        tdata.loc[summer, 'hLR'] = h_summer
                        
        # --- initialize the column with the computed heads
        tdata['h'] = tdata['hLR']
        
        # --- Add q column
        tdata['q'] = q
                
        # --- end of the previous day
        t0 = tdata.index[0] - np.timedelta64(1, 'D')
        
        # --- Head at the start of the first day
        h0 = tdata['hLR'][0]
        
        aq = self.aq
                
        hcol = list(tdata.columns).index('h')
        
        for it,(t, R, hLR, q) in enumerate(tdata[['rch', 'hLR', 'q']].itertuples()):

            # --- allow timesteps to vary
            dt = (t - t0) / np.timedelta64(1, 'D')
            
            exp = np.exp(-dt / aq.T)
            havg = hLR +  (h0 - hLR) * exp +  (R + q) * aq.c * aq.G * (1 - exp)
            
            # --- havg is for the start of the next day, because it used the data of today.
            try:
                tdata.iloc[it + 1, hcol] = havg
            except Exception:
                print(f"Done, it={it}")
                break
            h0 = havg
            t0 = t
            
        tdata['phi'] = tdata['h'] + tdata['q'] * aq.c
        return tdata
    
    def transient_0(self, time=None, R=None, h0=0., hLR=0., q=None):
        """Return result of dynamic simulation for constant inputs.
        
        Inputs are constant, time are an np.ndarray of floats [T].
        This allows observing asympotic behavior.

        Parameters
        ----------
        time: np.ndata | float
            time or times
        h0: float
            initial head
        hLR: float
            ditch water level
        q: float
            upward seepage rate
        """
        time = np.atleast_1d(time)
        h = np.zeros_like(time) + h0
        
        t0 = time[0]
        
        aq = self.aq
                
        exp = np.exp(-(time - t0) / aq.T)
        h =hLR +  (h0 - hLR) * exp +  (R + q) * aq.c * aq.G * (1 - exp)   
        return h
    
    
    def asymptote(self, R:float, q: float, hLR: float)->float:
        """Return the transient solution for t at infinity (sready X-section average)
        
        havg = hLR + (R + q) c G
        
        """
        aq = self.aq        
        return hLR + (R + q) * aq.c * aq.G        
        
def base_case_steady():
    """Show the base case steady and compare with single layer."""
    
    from itertools import cycle
    
    # --- Input data
    aq = Aquifer(k=10, D=10, c=200, mu=0.2, b=50)
    phi, hLR, R = 0., 0., 0.001

    # --- Points along X-section
    x = np.linspace(-aq.b, aq.b, 101)
    
    # --- Setup figure
    fig, ax = plt.subplots(figsize=(10, 10))
    fig.suptitle("Base case analytical solution, steady-state")
    
    ax.set_title('Analytial base-case, steady state' + "\n"
                 + str(aq).replace(' c=200,', '')
                 )
    ax.set(xlabel='x [m]', ylabel='h [m]')
    
    # --- Simulate and show for different c-values.
    clrs = cycle('brgkmcy')
    for c in [10, 30, 100, 300, 1000, 3000, 10000]:
        
        # --- Manage graph color.
        clr = next(clrs)

        # --- Replace c, this gives new aquifer and new model.
        aq.c = c
        mdl = Base_case(aq=aq)
    
        # --- head along X-section and in center x=0
        hx = mdl.steady(x=x, phi=phi, hLR=hLR, R=R)
        h0 = mdl.steady(x=0, phi=phi, hLR=hLR, R=R)
        
        ax.plot(x, hx, color=clr, label=f'c={c} d')
        ax.plot(0, h0, 'o', ms=8, mec=clr, mfc='none')
    
    # --- Steady analytical solution single layer, with no leakage.
    ha = R * (mdl.aq.b **2 - x ** 2) / (2 * mdl.aq.kD)
    ax.plot(x, ha, '.', color='k', label='steady one-layer')    
    ax.grid(True)
    ax.legend(loc='upper right')
    plt.show()
    
def base_case_transient(b=50, R=0.001, h0=0, hLR=0, q=0):
    """Show head development for steady inputs together with asymptote
    
    Parameters
    ----------
    b: float [L]
        half-width of the X-section
    R: float [L/T]
        Recharge
    h0: float
        Initial head at t=0
    hLR: float
        Ditch water level.
    q: float [L/T]
        Upward seepage
    """
    # --- Input data
    aq = Aquifer(k=10, D=10, c=200, mu=0.2, b=b)
        
    title1 = str(aq).replace(", c = 200", "")
    title2 = f"R={R} m/d, h0={h0} m, hLR={hLR} m, q={q} m/d"
    
    time = np.logspace(0, 4, 81)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_title("Head development voor stationary inputs" + "\n" + 
                 title1 + "\n" + title2)
    ax.set(xlabel='time [d]', ylabel='head [m]', xscale='log')
    
    clrs = cycle('brgkmcy')
    for c in [10, 30, 100, 300, 1000, 3000]:
        clr = next(clrs)
        aq.c = c
        mdl = Base_case(aq=aq)
        
        h = mdl.transient_0(time=time, R=R, h0=h0, hLR=hLR, q=q)
        hinf = mdl.asymptote(R=R, hLR=hLR, q=q)

        ax.plot(time[1:], h[1:], color=clr, label=f'c={aq.c:7.0f} d, T={aq.T:8.3g} d')
        ax.plot(time[-1], hinf, 'o', mfc=clr)
    ax.grid(True)
    ax.legend(loc="lower right")    

def base_transient(rch=None, b=50, h0=0, h_summer=-0.9, h_winter=-1.1, q=0):
    """Show transient head development driven by meteo
    
    Parameters
    ----------
    tdata: pd DataFrame with fields RH (precip) and EV24 (evapotranspiraton)
        The meteo data
    b: float [L]
        half-width of the X-section
    h0: float
        Initial head at t=0
    h_winter: float
        Ditch water level during winter.
    h_summer: float
        Ditch water level during summer.
    q: float [L/T]
        Upward seepage
    """
    # --- Input data
    aq = Aquifer(k=10, D=10, c=200, mu=0.2, b=b)
        
    title1 = str(aq).replace(", c = 200", "")
    title2 = f"h0={h0} m, h_summer={h_summer} m, h_winter={h_winter} m, q={q} m/d"
    
    # --- Get meteo data
    if rch is None:
        meteo = ggor_meteo.Meteo()
        rch = meteo.recharge
        
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_title("Head driven by meteo data" + "\n" + 
                 title1 + "\n" + title2)
    ax.set(xlabel='time', ylabel='head [m]', xscale='log')
    
    clrs = cycle('brgkmcy')
    for c in [10, 1000]: # 30, 100, 300, 1000, 3000]:
        clr = next(clrs)
        aq.c = c
        mdl = Base_case(aq=aq)
        
        tdata = mdl.transient(rch, h_summer=h_summer, h_winter=h_winter, q=q)
        
        ax.plot(tdata.index, tdata['h'], color=clr,
                label=f'c={aq.c:7.0f} d, T={aq.T:8.3g} d')
    ax.grid(True)
    ax.legend(loc="lower right")    

if __name__ == "__main__":
    if False:
        base_case_steady()
    if False:
        base_case_transient(b=50, R=0.001, h0=0, hLR=0, q=0)
        base_case_transient(b=500, R=0.001, h0=0, hLR=0, q=0)
        base_case_transient(b=5000, R=0.001, h0=0, hLR=0, q=0)    
        base_case_transient(b=50, R=0.0, h0=3, hLR=0, q=0)
        base_case_transient(b=50, R=0.0, h0=0, hLR=3, q=0)
        base_case_transient(b=50, R=0.0, h0=0, hLR=0, q=0.001)
    if True:
        rch = ggor_meteo.Meteo().recharge
        h_summer, h_winter = -0.9, -1.1
        h0, q = 0, 0
        base_transient(rch, b=50, h0=3, h_summer=0, h_winter=0, q=q)
        
        
    
    plt.show()
    

