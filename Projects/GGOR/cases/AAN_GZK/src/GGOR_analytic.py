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
from pathlib import Path
from scipy.optimize import brentq
from scipy.signal import lfilter
# %%

@dataclass
class Aquifer:
    k: float
    D: float
    c: float
    w: float
    mu: float
    b: float
    
    def __post_init__(self):
        assert self.k > 0
        assert self.D > 0
        assert self.c > 0
        assert self.w >= 0
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
    def T1L(self):
        """Return T voor single layer (Dupuit) aquifer"""
        return self.b**2 / (3 * self.kD) + self.w * self.b / self.D 
        
    @property
    def T2L(self):
        """Return T for aquifer on top of regional aquifer (base-case)."""
        return self.c * ((self.w / self.c) * (self.b / self.D) + self.G)
        
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
            w = self.w,
            b=self.b,
            kD = self.kD,
            lam=self.lam,            
            G=self.G,
        )

# aq = Aquifer(k=10.0, D=5.0, c=0.2, mu=0.15, b=50, w=0)
# print(aq)
# print(aq.summary())
        
def root_x_tan_x(a, k):
    """Return k roots of x tan(x) = a.
    
    Required for manu solutions with Bessel functions and
    solutions with sums involving entry resistance.
    """
    def f(x):
        return x*np.tan(x) - a
    eps = 1e-12
    left  = k*np.pi + eps
    right = k*np.pi + np.pi/2 - eps
    return brentq(f, left, right)


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
    def steady_avg(self, Nx: int, phi: float, hLR: float, R: float, **kwargs)->tuple:
        """Return steady-state solution (h, x)"""
        pass

    @abstractmethod
    def transient_avg(self, time_props)-> pd.DataFrame:
        """Return X-section average transient head."""
        pass

    @abstractmethod
    def transient(self, time_props)-> pd.DataFrame:
        """Return transient head."""
        pass
    
    @abstractmethod
    def transient_pd(self, time_props)-> np.ndarray | float:
        """Return transient head with pandas input and output"""
        pass
    
    @abstractmethod
    def asymptote(self, time_props)->float:
        """Return the transient solution for t at infinity (sready X-section average)
        """
        pass
    
    def block_response(self, time: np.ndarray,
                       R: float=0,
                       dh:float=0,
                       q: float=0) -> np.ndarray:
        """Return the block response for R, dh or q."""
        L = (R, dh, q)
        assert all(x in (0, 1) for x in L) and sum(L) == 1, (
        """
        To get the correct block response make sure that:
            To get BR for recharge use R=1 and the rest zero
            To get BR for dh, dh=1 and the rest 0.
            To get BR for q: use q=1 and the rest 0.
        """
        )
        h = self.transient(time=time, R=R, dh=dh, q=q)
        h[1:] -= h[:-1]
        return h[h >= 0.001]

     
    def sim_by_lfilter(self, rch, h_summer=0, h_winter=0, q=0):
        """Return head simulated by convolution."""
        tdata = pd.DataFrame(rch, columns=['rch'])
        time = (tdata.index - tdata.index[0]) / np.timedelta64(1, 'D')
        
        summer = rch.index.month >3 and rch.index.month < 10
        hLR = np.zeros(len(rch)) + h_winter
        hLR[summer] = h_summer
        tdata['hLR'] = hLR
        tdata['q'] = q
        
        b_R  = self.block_response(time, R=1)
        b_dh = self.block_response(time, dh=1)        
        b_q  = self.block_response(time, q=1)
        
        h = (  lfilter(b_R,  1, tdata['rch'])
             + lfilter(b_dh, 1, tdata['hLR'] - h_winter)
             + lfilter(b_q,  1, tdata['q'])
             + h_winter
        )
        tdata['h'] = h
        return tdata
  

# %% Implementation base case

class Dupuit(AnalyticalSolution):
    """
    The Dupuit case is a single layer with constant kD and
    no regional aquifr below it, but entry resistance
    is included.
    
    The transient case is only for X-section-average head.
    """
    
    def steady(self, x: int, hLR: float, R: float)->tuple:
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
        h = hLR + R / (2 * aq.kD) * (aq.b**2 - x**2) + R * aq.w * aq.b / aq.D
        return h if len(h) > 1 else h.item()
    
    def steady_avg(self, hLR: float, R: float)->float:
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
        return hLR + R * aq.T1L

    def transient_pd(self, rch: pd.Series,
            h_summer: float=0, h_winter: float=0, q: float=0) -> pd.DataFrame:
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
                        
        # --- end of the previous day
        t0 = tdata.index[0] - np.timedelta64(1, 'D')
        
        # --- Head at the start of the first day
        h0 = tdata['hLR'][0]
        
        aq = self.aq
                
        hcol = list(tdata.columns).index('h')
        
        for it, (t, R, hLR) in enumerate(tdata[['rch', 'hLR']].itertuples()):

            # --- allow timesteps to vary
            dt = (t - t0) / np.timedelta64(1, 'D')
            
            T = aq.T1L
            
            exp = np.exp(-dt / (aq.mu * T))
            havg = hLR +  (h0 - hLR) * exp +  R * T * (1 - exp)
            
            # --- havg is for the start of the next day, because it used the data of today.
            try:
                tdata.iloc[it + 1, hcol] = havg
            except Exception:
                print(f"Done, it={it}")
                break
            h0 = havg
            t0 = t
            
        return tdata
    
    def transient_avg(self, time=None, R=None, h0=0, hLR=0):
        return self.transient(time=time, R=R, h0=h0, hLR=hLR)
    
    def transient(self, time=None, R=None, h0=0, hLR=0):
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
        """
        time = np.atleast_1d(time)
        h = np.zeros_like(time) + h0
        
        t0 = time[0]
        
        aq = self.aq
        
        T = aq.T1L
                
        exp = np.exp(-(time - t0) / (aq.mu * T))
        h =hLR +  (h0 - hLR) * exp +  R * T * (1 - exp)   
        return h
    
    def asymptote(self, R:float, hLR: float)->float:
        """Return the transient solution for t at infinity (sready X-section average)
        
        havg = hLR + R T
        
        """
        aq = self.aq        
        return hLR + R * aq.T1L
    
class Base_case(AnalyticalSolution):
    """The GGOR-tool base case is for a single layer  of constant kD
    above a regional aquifer with fixed head phi separated by
    a vertical resistance c. The current base case includes
    entry resistance.
    
    The transient case is only for X-section average head.
    """
    
    def steady_phi(self, x: float | np.ndarray,
                       phi: float, hLR: float, R: float)->tuple:
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
        chx = np.cosh(x / aq.lam)
        chb = np.cosh(aq.b / aq.lam)
        shb = np.sinh(aq.b / aq.lam)  
        h =hLR + (R  + (phi - hLR) / aq.c) * aq.c * (
            1 - chx /(aq.k * aq.w / aq.lam * shb + chb))
        return h if len(h) > 1 else h.item()
        
    def steady(self, x: float | np.ndarray, hLR: float,
                                    R: float, q: float)->tuple:
        """Return steady state solution of cross section x for given q.
        
        For this regional head phi has to be converted to seepage phi.
        This is done by first computing the X-section averaged water
        table, and then taking phi = havg + q c, after which the
        head-method is called.
        
        x: float | np.ndarray [m] 
            x-coordinates  -b <= x <= b            
        hLR: float [L]
            water level in the ditch
        R: float [L/T]
            recharge rate
        q: float [L/T]
            Upward seepage from regional aquifer
        """
        aq = self.aq
        havg = hLR + (R + q) * aq.T2L
        phi = q * aq.c + havg
        
        return self.steady(x, phi=phi, hLR=hLR, R=R)
    
    def steady_avg(self, hLR: float, R: float, q: float)->float:
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
        return hLR + (R + q) * aq.T2L        

    def transient_pd(self, rch: pd.Series,
                  h_summer: float, h_winter: float, q: float)->pd.DataFrame:
        """Return result of dynamic simulation.
        
        Notice: hLR is not an input parameters, it's taken form h_winter

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
            
            T = aq.T2L
            
            exp = np.exp(-dt / T)
            havg = hLR +  (h0 - hLR) * exp +  (R + q) * T * (1 - exp)
            
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
        
    def transient(self, time: float | np.ndarray,
                  R: float=0,
                  h0: float=0,
                  hLR: float=0,
                  q: float=0)->float | np.ndarray:
        """Return result of dynamic simulation for constant inputs."""
        return self.transient_avg(time=time, R=R, h0=h0, hLR=hLR, q=q)
    
    def transient_avg(self, time: float | np.ndarray,
                      R: float=0,
                      h0: float=0,
                      hLR: float=0,
                      q: float=0)->float | np.ndarray:
        """Return X-sec avg dynamic head for constant inputs.
        
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
        
        T = aq.T2L
                
        exp = np.exp(-(time - t0) / (aq.mu * T))
        h =hLR +  (h0 - hLR) * exp +  (R + q) * T * (1 - exp)   
        return h if h.size > 1 else h.item()
    
    
    def asymptote(self, R:float, q: float, hLR: float)->float:
        """Return the transient solution for t at infinity (sready X-section average)
        
        havg = hLR + (R + q) c G
        
        """
        aq = self.aq
        return hLR + (R + q) * aq.T2L
    
class Brug(ABC):
    def __init__(self, aq):
        self.aq = aq

class Brug13302(Brug):
    """One-layer transient cross section bounded by surface water at x=p/m b.
    
    Case 133.02 is the case wat which surface water suddeny rise
    with a fixed amount, with no entry resistance.
    """
    def steady(self, x: float | np.ndarray, dh: float) -> float | np.ndarray:
        x = np.atleast_1d(x)
        h = dh * np.ones_like(x)
        return h if h.size > 1 else h.item()
    
    def steady_avg(self, dh:float) -> float:
        return dh
            
    def transient(self, dh: float, time: float | np.ndarray, x: float | np.ndarray, eps=1e-12) -> float | np.ndarray:
        
        assert (   np.isscalar(x) and not np.isscalar(time)
                or np.isscalar(time) and not np.isscalar(x)
                or np.isscalar(x) and np.isscalar(time)
        ), "x and time must be both scalars or (scalar and array) not both arrays"
        
        aq = self.aq
        
        T = 4 * aq.b**2 * aq.mu / (np.pi**2 * aq.kD)
        tau = time / T
        px2b = np.pi * x / (2 * aq.b)
        
        N = int(np.ceil(np.sqrt(T/(time[1] + eps))))
        N = 30
        
        s0 = self.steady(x, dh=dh)        
        F = 4 * dh / np.pi
        
        for n in range(N + 1):
            n2p1 = 2*n + 1
            ds = (
            (-1)**n / n2p1
            * np.cos(n2p1 * px2b)
            * np.exp(-n2p1**2 * tau)
            )
            if n == 0:
                s = ds
            else:
                s += ds
        return s0 - F * s
    
    def transient_avg(self, dh: float, time: float | np.ndarray, eps=1e-12) -> float | np.ndarray:
        """Return cross-section average head."""
                
        aq = self.aq
        
        T = 4 * aq.b**2 * aq.mu / (np.pi**2 * aq.kD)
        tau = time / T        
        
        N = int(np.ceil(np.sqrt(T/(time[1] + eps))))
        N = 30
                
        F = 8 * dh / np.pi**2
        
        for n in range(N + 1):
            n2p1 = 2*n + 1
            ds = (
            1 / n2p1**2            
            * np.exp(-n2p1**2 * tau)
            )
            if n == 0:
                s = ds
            else:
                s += ds
        return dh - F * s
    
    def block_response(self, time: np.ndarray, x:float) -> np.ndarray:
        """Return the block response."""
        if x is None:
            dh = self.transient_avg(dh=1, time=time)
        else:
            dh = self.transient(dh=1, time=time, x=x)
        dh[1:] -= dh[:-1]
        return dh[dh >= 0.001]
    
    def sim_lfilter(self, rch=None, h_summer=0, h_winter=0):
        """Simulate the head using Bruggeman's 133.02.
        
        These solutions are without entry resistance
        
        """        
        tdata = pd.DataFrame(rch, columns=['rch'])
        summer = np.logical_and(tdata.index.month > 3, tdata.index.month < 10)
        hLR = np.zeros(len(tdata)) + h_winter
        hLR[summer] = h_summer
        tdata['hLR'] = hLR
        
        time = (tdata.index - tdata.index[0]) / np.timedelta64(1, 'D')
                
        Br1 = Brug13302(self.aq)        
        
        b1 = Br1.block_response(time=time)        
        tdata['h'] = (  lfilter(b1, 1, tdata['hLR'] - h_winter)             
             + h_winter
             )
        return tdata

    
class Brug13316(Brug):
    """One-layer transient cross section bounded by surface water at x=p/m b.
    
    Case 133.16 is the case where recharge starts at t=0,
    without entry resistance.
    """
    
    def steady(self, x: float | np.ndarray, R: float) -> float | np.ndarray:
        x = np.atleast_1d(x)
        aq=self.aq
        h = R / (2 * aq.kD) * (aq.b**2 - x**2)
        return h if h.size > 1 else h.item()

    def steady_avg(self, R: float) -> float | np.ndarray:        
        aq=self.aq
        h = R * aq.b**2 / (3 * aq.kD)
        return h
    
    def transient(self, R: float, time: float | np.ndarray,
                  x: float | np.ndarray,
                  eps=1e-12) -> float | np.ndarray:
        
        assert (   np.isscalar(x) and not np.isscalar(time)
                or np.isscalar(time) and not np.isscalar(x)
                or np.isscalar(x) and np.isscalar(time)
        ), "x and time must be both scalars or (scalar and array) not both arrays"
        
        aq = self.aq
        
        T = 4 * aq.b**2 * aq.mu / (np.pi**2 * aq.kD)
        tau = time / T
        px2b = np.pi * x / (2 * aq.b)
        
        N = int(np.ceil(np.sqrt(T/(time[1] + eps))))
        N = 30
        
        s0 = self.steady(x, R)
        
        # --- Factor in boek 16, moet 8 zijn       
        F = 16 * R * aq.b**2 / (np.pi**3 * aq.kD)
        
        for n in range(N + 1):
            n2p1 = 2*n + 1
            ds = (
                (-1)**n / n2p1**3
                * np.cos(n2p1 * px2b)
                * np.exp(-n2p1**2 * tau)
            )
            if n == 0:
                s = ds
            else:
                s += ds
        return s0 - F * s
    
    def transient_avg(self, R: float, time: float | np.ndarray,
                  eps=1e-12) -> float | np.ndarray:
        """Return X-section average head."""
        
        aq = self.aq
        
        T = 4 * aq.b**2 * aq.mu / (np.pi**2 * aq.kD)
        tau = time / T        
        
        N = int(np.ceil(np.sqrt(T/(time[1] + eps))))
        N = 30
        
        s0 = self.steady_avg(R=R)
        
        # --- Factor in boek 16, moet 8 zijn       
        F =  32 * R * aq.b**2 / (np.pi**4 * aq.kD)
        
        for n in range(N + 1):
            n2p1 = 2*n + 1
            ds = (
                1 / n2p1**4                
                * np.exp(-n2p1**2 * tau)
            )
            if n == 0:
                s = ds
            else:
                s += ds
        return s0 - F * s
    
    def block_response(self, time: np.ndarray, x:float) -> np.ndarray:
        """Return the block response."""
        if x is None:
             h = self.transient_avg(R=1, time=time)
        else:
            h = self.transient(R=1, time=time, x=x)
        h[1:] -= h[:-1]
        return h[h >= 0.001]

    def sim_lfilter(self, rch=None):
        """Simulate the head using Bruggeman's 133.16.
        
        These solutions are without entry resistance
        
        """        
        tdata = pd.DataFrame(rch, columns=['rch'])        
        tdata['hLR'] = np.zeros(len(tdata))        
                
        time = (tdata.index - tdata.index[0]) / np.timedelta64(1, 'D')
        t0 = tdata.index[0]
                
        Br1 = Brug13316(self.aq)        
        
        b1 = Br1.block_response(time=time)        
        tdata['h'] = lfilter(b1, 1, tdata['rch']) + tdata.loc[t0, 'hLR']
        return tdata
    

class Brug13702(Brug):
    """One-layer transient cross section bounded by surface water at x=p/m b.
    
    Case 137.02 is the case in which the surface water suddenly rises
    with a fixed amount, with entry resistance.
    """    
    def steady(self, x: float | np.ndarray, dh: float) -> float | np.ndarray:
        x = np.atleast_1d(x)
        h = dh * np.zeros_like(x)
        return h if h.size > 1 else h.item()
            
    def transient(self, dh: float, time: float | np.ndarray, x: float | np.ndarray) -> float | np.ndarray:
        
        assert (   np.isscalar(x) and not np.isscalar(time)
                or np.isscalar(time) and not np.isscalar(x)
                or np.isscalar(x) and np.isscalar(time)
        ), "x and time must be both scalars or (scalar and array) not both arrays"
        
        aq = self.aq
             
        N = 30
        
        if np.isclose(aq.w, 0):
            raise ValueError("aq.w must be > 0 to use this function")
        
        eps   = aq.b / (aq.k * aq.w)
        alfas = [root_x_tan_x(eps, k) for k in range(N)]
        
        for n, alpha in zip(range(N + 1), alfas):
            T = aq.b**2 * aq.mu / (aq.kD * alpha**2)             
            ds = (
                (np.sin(alpha) / alpha) / (1 + eps / (alpha**2 + eps**2))
                * np.cos(alpha * x /aq.b)
                * np.exp(-time/ T)
            )
            if n == 0:
                s = ds
            else:
                s += ds
        return dh - 2 * dh * s
    
    def transient_avg(self, dh: float, time: float | np.ndarray) -> float | np.ndarray:
        """Returd X-section averaged head."""
        
        aq = self.aq
             
        N = 30
        
        if np.isclose(aq.w, 0):
            raise ValueError("aq.w must be > 0 to use this function")
        
        eps   = aq.b / (aq.k * aq.w)
        alfas = [root_x_tan_x(eps, k) for k in range(N)]
        
        for n, alpha in zip(range(N + 1), alfas):
            T = aq.b**2 * aq.mu / (aq.kD * alpha**2)             
            ds = (
                np.sin(alpha)**2 / (1 + eps / (alpha**2 + eps**2))                
                * np.exp(-time/ T)
            )
            if n == 0:
                s = ds
            else:
                s += ds
        return dh - dh * s
    
    def block_response(self, time: np.ndarray, x:float) -> np.ndarray:
        """Return the block response."""
        if x is None:
            dh = self.transient_avg(dh=1, time=time)
        else:
            dh = self.transient(dh=1, time=time, x=x)
        dh[1:] -= dh[:-1]
        return dh[dh >= 0.001]
    
    def sim_lfilter(self, rch=None, h_summer=0, h_winter=0):
        """Simulate the head using Bruggeman's 133.02.
        
        This solution is with entry resistance
        
        """        
        tdata = pd.DataFrame(rch, columns=['rch'])
        summer = np.logical_and(tdata.index.month > 3, tdata.index.month < 10)
        hLR = np.zeros(len(tdata)) + h_winter
        hLR[summer] = h_summer
        tdata['hLR'] = hLR
        
        time = (tdata.index - tdata.index[0]) / np.timedelta64(1, 'D')
                
        Br1 = Brug13702(self.aq)        
        
        b1 = Br1.block_response(time=time)        
        tdata['h'] = (  lfilter(b1, 1, tdata['hLR'] - h_winter)             
             + h_winter
             )
        return tdata


class Brug13709(Brug):
    """One-layer transient cross section bounded by surface water at x=p/m b.
    
    Case 137.09 is the case in which constant recharge starts at t=0,
    with entry resistance.
    """

    """Sudden rise h of the surface water level."""    
    
    def steady(self, x: float | np.ndarray, R: float) -> float | np.ndarray:
        x = np.atleast_1d(x)
        aq = self.aq
        h = R / (2 * aq.kD) * (aq.b**2 - x**2) + R * aq.b * aq.w / aq.D
        return h if h.size > 1 else h.item()
    
    def steady_avg(self, R: float) -> float | np.ndarray:        
        aq = self.aq
        h = R * aq.b**2 / (3 * aq.kD) + R * aq.b * aq.w / aq.D
        return h

            
    def transient(self, R: float, time: float | np.ndarray, x: float | np.ndarray) -> float | np.ndarray:
        
        assert (   np.isscalar(x) and not np.isscalar(time)
                or np.isscalar(time) and not np.isscalar(x)
                or np.isscalar(x) and np.isscalar(time)
        ), "x and time must be both scalars or (scalar and array) not both arrays"
        
        aq = self.aq
        
        if np.isclose(aq.w, 0):
            raise ValueError("aq.w must be > 0 to use this function")

        
        s0 = self.steady(x=x, R=R)
        
        F = 2 * R * aq.b**2 / aq.kD        
        N = 30
        
        eps   = aq.b / (aq.k * aq.w)
        alfas = [root_x_tan_x(eps, k) for k in range(N)]
        
        for n, alpha in zip(range(N + 1), alfas):
            T = aq.b**2 * aq.mu / (aq.kD * alpha**2)             
            ds = (
                (np.sin(alpha) / alpha**3) / (1 + eps / (alpha**2 + eps**2))
                * np.cos(alpha * x /aq.b)
                * np.exp(-time/ T)
            )
            if n == 0:
                s = ds
            else:
                s += ds
        return s0 - F * s

    def transient_avg(self, R: float, time: float | np.ndarray) -> float | np.ndarray:
        """Return X-section average head."""
               
        aq = self.aq
        
        if np.isclose(aq.w, 0):
            raise ValueError("aq.w must be > 0 to use this function")

        
        s0 = self.steady_avg(R=R)
        
        F = R * aq.b**2 / aq.kD        
        N = 30
        
        eps   = aq.b / (aq.k * aq.w)
        alfas = [root_x_tan_x(eps, k) for k in range(N)]
        
        for n, alpha in zip(range(N + 1), alfas):
            T = aq.b**2 * aq.mu / (aq.kD * alpha**2)             
            ds = (
                (np.sin(alpha)**2 / alpha**2) / (1 + eps / (alpha**2 + eps**2))                
                * np.exp(-time/ T)
            )
            if n == 0:
                s = ds
            else:
                s += ds
        return s0 - F * s

    def block_response(self, time: np.ndarray, x:float) -> np.ndarray:
        """Return the block response."""
        if x is None:
            h = self.transient_avg(R=1, time=time)
        else:
            h = self.transient(R=1, time=time, x=x)
        h[1:] -= h[:-1]
        return h[h >= 0.001]
    
    def sim_lfilter(self, rch=None):
        """Simulate the head using Bruggeman's 137.09.
        
        This solutions are with entry resistance
        
        """        
        tdata = pd.DataFrame(rch, columns=['rch'])        
        tdata['hLR'] = np.zeros(len(tdata))        
                
        time = (tdata.index - tdata.index[0]) / np.timedelta64(1, 'D')
        t0 = tdata.index[0]
                
        Br1 = Brug13709(self.aq)        
        
        b1 = Br1.block_response(time=time)        
        tdata['h'] = lfilter(b1, 1, tdata['rch']) + tdata.loc[t0, 'hLR']
        return tdata


class Brug133(Brug):
    """Class to simulate one-layer exactly without entry resistance.
    
    Using Bruggeman(1999 solution 133.02 and 133.16)    
    """
    
    def sim_lfilter(self, rch=None, h_summer=0, h_winter=0):
        """Simulate the head using Bruggeman's 133.02 and 133.16.
        
        These solutions are without entry resistance
        
        """
        
        tdata = pd.DataFrame(rch, columns=['rch'])
        summer = np.logical_and(tdata.index.month > 3, tdata.index.month < 10)
        hLR = np.zeros(len(tdata)) + h_winter
        hLR[summer] = h_summer
        tdata['hLR'] = hLR
        
        time = (tdata.index - tdata.index[0]) / np.timedelta64(1, 'D')
                
        Br1 = Brug13302(self.aq)
        Br2 = Brug13316(self.aq)
        
        b1 = Br1.block_response(time=time)
        b2 = Br2.block_response(time=time)
        tdata['h'] = (  lfilter(b1, 1, tdata['hLR'] - h_winter)
             + lfilter(b2, 1, tdata['rch'])
             + h_winter
             )
        return tdata

class Brug137(Brug):
    """Class to simulate one-layer exactly without entry resistance.
    
    Using Bruggeman(1999 solution 133.02 and 133.16)    
    """
    
    def sim_lfilter(self, rch=None, h_summer=0, h_winter=0):
        """Simulate the head using Bruggeman's 133.02 and 133.16.
        
        These solutions are without entry resistance
        
        """        
        tdata = pd.DataFrame(rch, columns=['rch'])
        summer = np.logical_and(tdata.index.month > 3, tdata.index.month < 10)
        hLR = np.zeros(len(tdata)) + h_winter
        hLR[summer] = h_summer
        tdata['hLR'] = hLR
        
        time = (tdata.index - tdata.index[0]) / np.timedelta64(1, 'D')
                
        Br1 = Brug13702(self.aq)
        Br2 = Brug13709(self.aq)
        
        b1 = Br1.block_response(time=time)
        b2 = Br2.block_response(time=time)
        tdata['h'] = (  lfilter(b1, 1, tdata['hLR'] - h_winter)
             + lfilter(b2, 1, tdata['rch'])
             + h_winter
             )
        return tdata
    

def ex_dupuit_transient(b=50, R=0.001, h0=0, hLR=0, w=0):
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
    """
    # --- Input data
    aq = Aquifer(k=10, D=10, c=200, w=w, mu=0.2, b=b)
        
    title1 = str(aq).replace(", c = 200", "").replace(", b = 50", "")
    title2 = f"R={R} m/d, h0={h0} m, hLR={hLR} m"
    
    time = np.logspace(0, 4, 81)
    
    fig, ax = plt.subplots(figsize=(10, 7.5))
    fig.suptitle("Dupuit case")
    ax.set_title("Head development voor stationary inputs" + "\n" + 
                 title1 + "\n" + title2)
    ax.set(xlabel='time [d]', ylabel='head [m]', xscale='log')
    
    clrs = cycle('brgkmcy')
    for b in [30, 45, 60, 85, 120, 170, 240]:
        clr = next(clrs)
        aq.b = b
        mdl = Dupuit(aq=aq)
        
        T = aq.T1L
        
        h = mdl.transient_avg(time=time, R=R, h0=h0, hLR=hLR)
        hinf = mdl.asymptote(R=R, hLR=hLR)

        ax.plot(time[1:], h[1:], color=clr, label=f'b={aq.b:7.0f} d, muT={aq.mu * T:8.3g} d')
        ax.plot(time[-1], hinf, 'o', mfc=clr)
    ax.grid(True)
    ax.legend(loc="upper left")    


def ex_dupuit_transient_pd(rch=None, b=50, h0=0, h_summer=-0.9, h_winter=-1.1):
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
    """
    # --- Input data
    aq = Aquifer(k=10, D=10, c=200, w=0., mu=0.2, b=b)
        
    title1 = str(aq).replace(", c = 200", "").replace(", b = 50", "")
    title2 = f"h0={h0} m, h_summer={h_summer} m, h_winter={h_winter} m"
    
    # --- Get meteo data
    if rch is None:
        meteo = ggor_meteo.Meteo()
        rch = meteo.recharge
        
    fig, ax = plt.subplots(figsize=(10, 7.8))
    fig.suptitle("Dupuit case")
    ax.set_title("Head driven by meteo data" + "\n" + 
                 title1 + "\n" + title2)
    ax.set(xlabel='time', ylabel='head [m]', xscale='log')
    
    clrs = cycle('brgkmcy')
    for b in [50, 500]:
        clr = next(clrs)
        aq.b = b
        mdl = Dupuit(aq=aq)
        
        T = aq.T1L
        
        tdata = mdl.transient_pd(rch, h_summer=h_summer, h_winter=h_winter)
        
        ax.plot(tdata.index, tdata['h'], color=clr,
                label=f'muT={aq.mu * T:8.3g} d')
    ax.grid(True)
    ax.legend(loc="lower right")
  

def ex_base_case_steady():
    """Show the base case steady and compare with single layer."""
    
    from itertools import cycle
    
    # --- Input data
    aq = Aquifer(k=10, D=10, c=200, w=0, mu=0.2, b=50)
    phi, hLR, R = 0., 0., 0.001

    # --- Points along X-section
    x = np.linspace(-aq.b, aq.b, 101)
    
    # --- Setup figure
    fig, ax = plt.subplots(figsize=(10, 7))
    fig.suptitle("Base case, steady-state")
    
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
    
def ex_base_case_transient0(b=50, R=0.001, h0=0, hLR=0, q=0):
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
    """
    # --- Input data
    aq = Aquifer(k=10, D=10, c=200, w=0, mu=0.2, b=b)
        
    title1 = str(aq).replace(", c = 200", "")
    title2 = f"R={R} m/d, h0={h0} m, hLR={hLR} m, q={q} m/d"
    
    time = np.logspace(0, 4, 81)
    
    fig, ax = plt.subplots(figsize=(10, 7))
    fig.suptitle("Base case, transient 0")
    ax.set_title("Head development voor stationary inputs" + "\n" + 
                 title1 + "\n" + title2)
    ax.set(xlabel='time [d]', ylabel='head [m]', xscale='log')
    
    clrs = cycle('brgkmcy')
    for c in [10, 30, 100, 300, 1000, 3000]:
        clr = next(clrs)
        aq.c = c
        mdl = Base_case(aq=aq)
        
        T = aq.c * aq.G
        
        h = mdl.transient(time=time, R=R, h0=h0, hLR=hLR, q=q)
        hinf = mdl.asymptote(R=R, hLR=hLR, q=q)

        ax.plot(time[1:], h[1:], color=clr, label=f'c={aq.c:7.0f} d, muT={aq.mu * T:8.3g} d')
        ax.plot(time[-1], hinf, 'o', mfc=clr)
    ax.grid(True)
    ax.legend(loc="lower right")    

def ex_base_transient(rch=None, b=50, h0=0, h_summer=-0.9, h_winter=-1.1, q=0):
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
    aq = Aquifer(k=10, D=10, c=200, w=0, mu=0.2, b=b)
        
    title1 = str(aq).replace(", c = 200", "")
    title2 = f"h0={h0} m, h_summer={h_summer} m, h_winter={h_winter} m, q={q} m/d"
    
    # --- Get meteo data
    if rch is None:
        meteo = ggor_meteo.Meteo()
        rch = meteo.recharge
        
    fig, ax = plt.subplots(figsize=(10, 7))
    fig.suptitle("Base case, transient")
    ax.set_title("Head driven by meteo data" + "\n" + 
                 title1 + "\n" + title2)
    ax.set(xlabel='time', ylabel='head [m]', xscale='log')
    
    clrs = cycle('brgkmcy')
    for c in [10, 1000]: # 30, 100, 300, 1000, 3000]:
        clr = next(clrs)
        aq.c = c
        mdl = Base_case(aq=aq)
        
        T = aq.c * aq.G
        
        tdata = mdl.transient_pd(rch, h_summer=h_summer, h_winter=h_winter, q=q)
        
        ax.plot(tdata.index, tdata['h'], color=clr,
                label=f'c={aq.c:7.0f} d, muT={aq.mu * T:8.3g} d')
    ax.grid(True)
    ax.legend(loc="lower right")
    
def compare_limits():
    """Graphically compare the limits (b/lam coth(b/lam) + 1) with (1/3 (b/lam)^3)
    
    The single layer and layer with leakage to an underlying aquifer should become
    the same when
    
    (b/lam coth(b/lam) + 1) --> (1/3 (b/lam)^3) for lam --> infinity or
    b/lam --> 0.
    
    That this is true is shown graphically here.
    """
        
    fig, ax = plt.subplots()

    ax.set_title(r"$\frac{b}{\lambda}\coth\frac{b}{\lambda}-1\to\frac{b^{2}}{3\lambda^{2}}$")
    ax.set_xlabel(r"$b/\lambda$")
    ax.set_ylabel("f(x)")
    
    x = np.logspace(-1, np.log10(4), 21)
    
    ax.plot(x, x / np.tanh(x)-1, label=r"$\frac{b}{\lambda} \coth \left(\frac{b}{\lambda}\right) - 1$")
    ax.plot(x, (x**2) /3, label=r'$ \frac{1}{3} \left(\frac{b}{\lambda}\right)^2$')
    
    ax.grid(True)
    ax.legend()
    
    parts = Path(os.getcwd()).parts
    pth = os.path.join(*parts[:parts.index('GGOR') + 1], 'doc', 'images')
    
    fig.savefig(os.path.join(pth, "limits_1_and_2_layers.png"))
    
    plt.show()

def ex_brug13302():
    aq = Aquifer(k=10, D=10, c=200, w=0, mu=0.15, b=50)
    xs = np.array([0, 0.25, 0.5, 0.58, 0.75, 0.95]) * aq.b
    dh = 0.1
    time = np.linspace(0, 20, 101)
    
    brug = Brug13302(aq)
    dup  = Dupuit(aq)

    fig, ax = plt.subplots(figsize=(10, 7.5))
    
    fig.suptitle("Bruggeman(1999, solution 133.02)")
    ax.set_title(f"Sudden rise {dh} mof the surface water level" +
                 "\n" + str(aq).replace(", w=0","").replace(", c=200",""))
    ax.set(xlabel='t d[]', ylabel='h [m]')
    
    clrs = cycle('brgkmcy')
    for x in xs:
        clr = next(clrs)
        ht = brug.transient(dh=dh, time=time, x=x)
        ha = brug.transient_avg(dh=dh, time=time)
        hd = dup.transient(R=0, time=time, h0=0, hLR=dh)
        ax.plot(time, ht, '-', color=clr, label=f'Brug133.02, x={x:.1f} m')
        ax.plot(time, ha, 'o', mec=clr, mfc='none', label='Brug133.02 transient_avg')
        ax.plot(time, hd, 'x', mec=clr, mfc='none', label='Dupuit transient_avg')
        
    ax.grid(True)
    ax.legend(loc='center')

def ex_brug13702():
    aq = Aquifer(k=10, D=10, c=200, w=1, mu=0.15, b=50)
    xs = np.array([0, 0.25, 0.5, 0.58, 0.75, 0.95, 0.99]) * aq.b
    dh = 0.1
    time = np.linspace(0, 20, 101)
    
    brug1 = Brug13302(aq)
    brug2 = Brug13702(aq)
    dup   = Dupuit(aq)

    fig, ax = plt.subplots(figsize=(10, 7.5))
    
    fig.suptitle("Bruggeman (1999, solution 137.02)")
    ax.set_title(f"Sudden rise {dh} m of surface water. With entry resistance" +
                 "\n" + str(aq).replace(", c=200",""))
    ax.set(xlabel='t d[]', ylabel='h [m]')
    
    clrs = cycle('brgkmcy')
    for x in xs:
        clr = next(clrs)
        h1 = brug1.transient(dh=dh, time=time, x=x)
        h2 = brug2.transient(dh=dh, time=time, x=x)
        hd = dup.transient(time=time, R=0., h0=0, hLR=dh)
        hb = brug2.transient_avg(dh=dh, time=time)
        ax.plot(time, h1, '-', color=clr,
                label=f'h1, x={x:.1f} m, no   ditch resistance')
        ax.plot(time, h2, '.', color=clr,
                label=f'h2, x={x:.1f} m, with ditch resistance')
        ax.plot(time, hd, 'o', mec=clr, mfc='none',
                label="Dupuit with ditch resistance")
        ax.plot(time, hb, 'x', mec=clr, mfc='none',
                label="brug avg with ditch resistance")
    ax.grid(True)
    ax.legend(loc='center')


def ex_brug13316():
    aq = Aquifer(k=10, D=10, c=200, w=0, mu=0.15, b=50)
    xs = np.array([0, 0.25, 0.5, 0.58, 0.75, 0.95, 0.99]) * aq.b
    R = 0.001
    time = np.linspace(0, 20, 101)
    
    brug = Brug13316(aq)
    dup  = Dupuit(aq)

    fig, ax = plt.subplots(figsize=(10, 7.5))
    
    fig.suptitle("Bruggeman (1999, solution 133.16)")
    ax.set_title(f"Constant precipitation of {R} m/d" +
                 "\n" + str(aq).replace(", w=0", "").replace(", c=200",""))
    ax.set(xlabel='t d[]', ylabel='h [m]')
    
    clrs = cycle("brgkmcy")
    for x in xs:
        clr = next(clrs)
        ht = brug.transient(R=R, time=time, x=x)
        hb = brug.transient_avg(R=R, time=time)
        hd = dup.transient(time=time, R=R)
        hdup = dup.steady(x=x, hLR=0, R=R)
        
        ax.plot(time, ht, color=clr, label=f'Brug133.16, x={x:.1f} m')
        ax.plot(time, hb, 's', mec=clr, mfc='none', label='Brug133.16 avg')
        ax.plot(time, hd, 'o', mec=clr, mfc='none', label='Dupuit transient')      
        ax.plot(time[-1], hdup, 'x', mec=clr, mfc='none', label='Dupuit')
        
    ax.grid(True)
    ax.legend(loc='center')


def ex_brug13316a():
    aq = Aquifer(k=10, D=10, c=200, w=0.001, mu=0.15, b=50)
    xs = np.array([0, 0.25, 0.5, 0.58, 0.75, 0.95, 0.99]) * aq.b
    R = 0.001
    time = np.linspace(0, 20, 101)
    
    xs = np.array([0, 0.25, 0.5, 0.75, 0.9, 0.95]) * aq.b
    brug1 = Brug13316(aq)
    brug2 = Brug13709(aq)
    dup = Dupuit(aq)

    fig, ax = plt.subplots(figsize=(10, 7.5))
    
    fig.suptitle("Bruggeman (1999, solution 133.16 and 137.09)")
    ax.set_title(f"Constant precipitation of {R} m/d" +
                 "\n" + str(aq).replace(", c=200",""))
    ax.set(xlabel='t d[]', ylabel='h [m]')
    
    clrs = cycle("brgkmcy")
    for x in xs:
        clr = next(clrs)
        ht1 = brug1.transient(R=R, time=time, x=x)
        ht2 = brug2.transient(R=R, time=time, x=x)
        hdup = dup.steady(x=x, hLR=0, R=R)
        
        ax.plot(time, ht1, color=clr, label=f'Brug133.16, x={x:.1f} m')
        ax.plot(time, ht2, '.', color=clr, label=f'Brug137.09, x={x:.1f} m')
        ax.plot(time[-1], hdup, 'x', mec=clr, mfc='none', label='Dupuit')
        
    ax.grid(True)
    ax.legend(loc='center')


def ex_brug13709():
    aq = Aquifer(k=10, D=10, c=200, w=1, mu=0.15, b=50)
    xs = np.array([0, 0.25, 0.5, 0.58, 0.75, 0.95, 0.99]) * aq.b
    R = 0.001
    time = np.linspace(0, 20, 101)
    
    brug1 = Brug13316(aq)
    brug2 = Brug13709(aq)
    dup  = Dupuit(aq)

    fig, ax = plt.subplots(figsize=(10, 7.5))
    
    fig.suptitle("Bruggeman (137.09)")
    ax.set_title(f"Constant precipitation {R} m/d. With entry resistance." +
                 "\n" + str(aq).replace(", c=200",""))
    ax.set(xlabel='t d[]', ylabel='h [m]')
    
    clrs = cycle("brgkmcy")
    for x in xs:
        clr = next(clrs)
        h1t = brug1.transient(R=R, time=time, x=x)
        h2t = brug2.transient(R=R, time=time, x=x)
        h2a = brug2.transient_avg(R=R, time=time)
        hdt = dup.transient(R=R, time=time)
        hd  = dup.steady(x=x, hLR=0, R=R)
        
        ax.plot(time, h1t, '-', color=clr,
                label=f'Brug133.16, x={x} m, no resis.')
        ax.plot(time, h2t, '.', color=clr,
                label=f'Brug137.09, x={x} m, + resis.')
        ax.plot(time, h2a, 's', mec=clr, mfc='none',
                label= 'Brug137.09 avg, + resis.')
        ax.plot(time, hdt, 'o', mec=clr, mfc='none',
                label='Dupuit + resitance')
        ax.plot(time[-1], hd, 'x', color=clr)
        
    ax.grid(True)
    ax.legend(loc='center')

    
if __name__ == "__main__":
    if False:
        ex_base_case_steady()
    if False:
        ex_dupuit_transient(b=50, R=0.001, h0=0, hLR=0)
        ex_dupuit_transient_pd(rch=None, b=50, h0=0, h_summer=-0.9, h_winter=-1.1)
    if False:
        ex_base_case_transient(b=50, R=0.001, h0=0, hLR=0, q=0)
        ex_base_case_transient(b=500, R=0.001, h0=0, hLR=0, q=0)
        ex_base_case_transient(b=5000, R=0.001, h0=0, hLR=0, q=0)    
        ex_base_case_transient(b=50, R=0.0, h0=3, hLR=0, q=0)
        ex_base_case_transient(b=50, R=0.0, h0=0, hLR=3, q=0)
        ex_base_case_transient(b=50, R=0.0, h0=0, hLR=0, q=0.001)
    if False:
        rch = ggor_meteo.Meteo().recharge
        h_summer, h_winter = -0.9, -1.1
        h0, q = 0, 0
        ex_base_transient(rch, b=50, h0=3, h_summer=0, h_winter=0, q=q)
    if False:
        compare_limits()
    if False:
        ex_brug13316a()
    if True:
        ex_brug13302()
        ex_brug13316()
        ex_brug13702()
        ex_brug13709()
        
        
    plt.show()
    

