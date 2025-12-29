# %% Docstring
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
# %% Dataclass Aquifer
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
    
    @property
    def wbcD(self):
        """Return factor w b /(c D)."""
        return (self.w * self.b) / (self.c * self.D)
    
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


# %% Analytical solution (Dupuit and Base_class)
class AnalyticalSolution(ABC):
    
    def __init__(self, aq: Aquifer) -> None:
        """Instantiate analytic simulator."""
        self.aq = aq
        
    @abstractmethod
    def steady(self, Nx: int, phi: float, hLR: float, R: float, x:float | None=None)->float | np.ndarray:
        """Return steady-state solution (h, x)"""
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
                       q: float=0,
                       eps: float=0.0001) -> np.ndarray:
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
        try:
            h = self.transient(time=time, R=R, hLR=dh, q=q)
        except Exception:
            h = self.transient(time=time, R=R, hLR=dh)            
        h[1:] -= h[:-1]
        h = h[:np.max(np.where(h[1:] >= eps)[0], initial=-1) + 2]
        return h
     
    def sim_by_lfilter(self, rch, h_summer=0, h_winter=0, q=0):
        """Return head simulated by convolution."""
        tdata = pd.DataFrame(rch, columns=['rch'])
        time = np.asarray((tdata.index - tdata.index[0]) / np.timedelta64(1, 'D'))
        
        summer = np.logical_and(rch.index.month >3, rch.index.month < 10)
        hLR = np.zeros(len(rch)) + h_winter
        hLR[summer] = h_summer
        tdata['hLR'] = hLR
        tdata['q'] = q
        
        b_R  = self.block_response(time, R=1)
        b_dh = self.block_response(time, dh=1)
        
        # b_q always works even if q does not exist (then just zeros are returned)
        b_q  = self.block_response(time, q=1)

        
        h = (  lfilter(b_R,  1, tdata['rch'])
             + lfilter(b_dh, 1, tdata['hLR'] - h_winter)
             + lfilter(b_q, 1, tdata['q'])
             + h_winter
        )
        
        tdata['h'] = h
        return tdata

class Dupuit(AnalyticalSolution):
    """
    The Dupuit case is a single layer with constant kD and
    no regional aquifr below it, but entry resistance
    is included.
    
    The transient case is only for X-section-average head.
    """
    
    def steady(self, R:float, hLR:float, x: float | np.ndarray | None=None)->tuple:
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
        if x is None:
            return hLR + R * aq.T1L
        else:
            x = np.atleast_1d(x)    
            h = hLR + R / (2 * aq.kD) * (aq.b**2 - x**2) + R * aq.w * aq.b / aq.D
            return h if len(h) > 1 else h.item()
    

    def transient_pd(self, rch: pd.Series,
            h_summer: float=0, h_winter: float=0, q: float=0) -> pd.DataFrame:
        """Return result of dynamic simulation.

        Parameters
        ----------
        rch: pd.Series [m/d]
            recharge with index = pd.Timestamps
        h0: float
            head at t=0
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
        h0 = tdata['h'].values[0]
                        
        # --- end of the previous day
        t0 = tdata.index[0] - np.timedelta64(1, 'D')
                
        aq = self.aq
            
        T = aq.T1L # inlcludes entry resistance
    
        hcol = list(tdata.columns).index('h')
        
        for it, (t, R, hLR) in enumerate(tdata[['rch', 'hLR']].itertuples()):

            # --- allow timesteps to vary
            dt = (t - t0) / np.timedelta64(1, 'D')
                        
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
        h = np.zeros_like(time)
        
        t0 = time[0]
        
        aq = self.aq
        
        T = aq.T1L # includes entry resistance
                
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
    def steady_q(self, x: float | np.ndarray, R:float, q:float, hLR:float) -> float | np.ndarray:
        aq = self.aq
        bw_cD = aq.b * aq.w / (aq.c * aq.D)
        Rqc = (R + q) * aq.c
        chx = np.cosh(x / aq.lam)
        chb = np.cosh(aq.b / aq.lam)
        h = hLR + Rqc *(aq.G + bw_cD) + Rqc * (1 - bw_cD) * (1 + chx / chb)
        return h if h.size > 1 else h.item()
     
    def steady(self, R: float,
               phi: float | None=None,
               q: float | None=None,
               hLR: float=0,
               x: float | np.ndarray | None=None)->tuple:
        """Return steady state solution of cross section x for given q.
        
        If phi is not None -> compute h based on phi not q
        if phi is None -> compute h base on q not phi
        if x is not None --> return h at position x
        if x is None -> return X-section average head.

        Paramters        
        ---------
        R: float [L/T]
            recharge rate
        phi: float [L] | None
            water level in underlying regioinal aquifer.
        q: float [L/T] | None
            Upward seepage from regional aquifer
        hLR: float [L], default=0
            water level in the ditch
        x: float | np.ndarray [m] | None
            x-coordinates  -b <= x <= b            
        """
        assert (phi is None) != (q is None), (
            """phi and q can not both be None."""
        )        
        assert np.isscalar(phi) != np.isscalar(q), (
            """Either phi or q must be a scalar."""
        )
        aq = self.aq        
                
        if phi is not None:        
            havg = hLR + (R * aq.c + (phi -hLR)) * (aq.G + aq.wbcD) / (aq.G + aq.wbcD + 1)
            q = (phi - havg) / aq.c
        else:
            havg = hLR + (R + q) * aq.c * (aq.G + aq.wbcD)
            phi = havg + q * aq.c
        
        if x is None:
            return havg
        else:        
            x = np.atleast_1d(x)
            chx = np.cosh(x / aq.lam)
            chb = np.cosh(aq.b / aq.lam)
            
            hb = hLR + (R + q) * aq.c * aq.wbcD
            
            h = hb * chx / chb + (havg + (R + q) * aq.c) * (1 - chx / chb)
            
            return h if len(h) > 1 else h.item()
    

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
            
            exp = np.exp(-dt / (aq.mu * T))
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
        h = np.zeros_like(time)
        
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
    
# %% Bruggeman's solutions 133 and 137
class Brug(ABC):
    def __init__(self, aq):
        self.aq = aq
        
    def block_response(self, time: np.ndarray, x:float | None=None, eps=0.0001) -> np.ndarray:
        """Return the block response.
        
        Parameters
        ----------
        time: np.ndarray of float
            time
        x: float | None
            Position in the cross section  -b <=x<=b
            if None, then X-section average value is returned.
        """
        h = self.transient(time=time, R=1, dh=1, x=x)
        h[1:] -= h[:-1]                
        h = h[:np.max(np.where(h[1:] >= eps)[0], initial=-1) + 2]
        return h

    def sim_lfilter(self, rch:pd.Series=None,
                    h_summer:float=0, h_winter:float=0, x:float=None) -> np.ndarray:
        """Simulate the head using Bruggeman's solution 133.02.16, 137.02/09.
        
        Solutions 133.02/16 are without and 137.02/09 with entry resistance.
        
        Parameters
        ----------
        rch: pd.Series wit daily values an pd.timestamp index
            The recharge
        h_summer: float
            Ditch water level during summer (April-September).
        h_winter: float
            Ditch water level during winter.
        x: float | None
            Position in X-section -b<=x<=b
            If None, then X-section average is returned.
        """        
        tdata = pd.DataFrame(rch, columns=['rch'])
        summer = np.logical_and(tdata.index.month > 3, tdata.index.month < 10)
        hLR = np.zeros(len(tdata)) + h_winter
        hLR[summer] = h_summer
        tdata['hLR'] = hLR
        
        time = np.array((tdata.index - tdata.index[0]) / np.timedelta64(1, 'D'))
        
        b1 = self.block_response(time=time, x=x)        
        tdata['h'] = (  lfilter(b1, 1, tdata['hLR'] - h_winter)             
             + h_winter
             )
        return tdata

class Brug13302(Brug):
    """One-layer transient cross section bounded by surface water at x=p/m b.
    
    Case 133.02 is the case wat which surface water suddeny rise
    with a fixed amount, with no entry resistance.
    """
    def steady(self, dh:float, R:float | None=None, x: float | np.ndarray | None=None) -> float | np.ndarray:
        if x is None:
            # --- X-section average value
            return dh
        else:
            # --- Heads at x.
            x = np.atleast_1d(x)
            h = dh * np.ones_like(x)
            return h if h.size > 1 else h.item()
    
            
    def transient(self, time: float | np.ndarray,
                  dh: float,
                  R:float | None=None,
                  x: float | np.ndarray|None=None,
                  N=30,
                  eps=1e-12) -> float | np.ndarray:
        
        assert not (isinstance(x, np.ndarray) and isinstance(time, np.ndarray)), (
                "x and time must be both scalars or (scalar and array) not both arrays"
        )

        aq = self.aq
        T = 4 * aq.b**2 * aq.mu / (np.pi**2 * aq.kD)
        tau = time / T
        
        s0 = self.steady(dh=dh, x=x)        
        
        if x is None:
            # --- Return X-section average head.                 
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
            return s0 - F * s
        else:
            # --- Return head at x        
            px2b = np.pi * x / (2 * aq.b)
            
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

    
class Brug13316(Brug):
    """One-layer transient cross section bounded by surface water at x=p/m b.
    
    Case 133.16 is the case where recharge starts at t=0,
    without entry resistance.
    """
    
    def steady(self, R: float,
               dh:float | None=None,
               x: float | np.ndarray | None=None) -> float | np.ndarray:
        aq=self.aq
        if x is None:
            # --- Return X-section average.            
            h = R * aq.b**2 / (3 * aq.kD)
            return h
        else:
            # --- Return values at x.
            x = np.atleast_1d(x)            
            h = R / (2 * aq.kD) * (aq.b**2 - x**2)
            return h if h.size > 1 else h.item()
    
    def transient(self,
                  time: float | np.ndarray,
                  R:float,
                  dh: float | None=None,
                  x: float | np.ndarray | None=None,
                  N=30,
                  eps=1e-12) -> float | np.ndarray:
        
        assert not (isinstance(x, np.ndarray) and isinstance(time, np.ndarray)), (
                "x and time must be both scalars or (scalar and array) not both arrays"
        )
        
        aq = self.aq
        T = 4 * aq.b**2 * aq.mu / (np.pi**2 * aq.kD)
        tau = time / T
        
        s0 = self.steady(R=R, x=x)
    
        if x is None:
            # --- Return X-section average head.
            
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
        else:
            # --- Return head at position(s) x.        
            px2b = np.pi * x / (2 * aq.b)
                        
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

class Brug13702(Brug):
    """One-layer transient cross section bounded by surface water at x=p/m b.
    
    Case 137.02 is the case in which the surface water suddenly rises
    with a fixed amount, with entry resistance.
    """    
    def steady(self, dh:float, R:float | None=None, x: float | np.ndarray | None=None) -> float | np.ndarray:
        if x is None:
            return dh
        else:
            x = np.atleast_1d(x)
            h = dh + np.zeros_like(x)
        return h if h.size > 1 else h.item()
            
    def transient(self, 
                  time: float | np.ndarray,
                  dh: float,
                  R: float | None=None,
                  x: float | np.ndarray | None=None, N=30) -> float | np.ndarray:
        
        assert not (isinstance(x, np.ndarray) and isinstance(time, np.ndarray)), (
                "x and time must be both scalars or (scalar and array) not both arrays"
        )        
        aq = self.aq
         
        if np.isclose(aq.w, 0):
            # raise ValueError("aq.w must be > 0 to use this function")
            aq.w = 0.0001
    
        brug_eps   = aq.b / (aq.k * aq.w)
        alfas = [root_x_tan_x(brug_eps, k) for k in range(N)]
        
        s0 = self.steady(dh=dh, x=x)
        
        F = 2 * s0
        
        if x is None:
            for n, alpha in zip(range(N + 1), alfas):
                T = aq.b**2 * aq.mu / (aq.kD * alpha**2)             
                ds = (
                    (np.sin(alpha)/alpha)**2 / (1 + brug_eps / (alpha**2 + brug_eps**2))                
                    * np.exp(-time/ T)
                )
                if n == 0:
                    s = ds
                else:
                    s += ds
            return s0 - F * s
        else:
            for n, alpha in zip(range(N + 1), alfas):
                T = aq.b**2 * aq.mu / (aq.kD * alpha**2)             
                ds = (
                    (np.sin(alpha) / alpha) / (1 + brug_eps / (alpha**2 + brug_eps**2))
                    * np.cos(alpha * x / aq.b)
                    * np.exp(-time/ T)
                )
                if n == 0:
                    s = ds
                else:
                    s += ds
        return s0 - F * s


class Brug13709(Brug):
    """One-layer transient cross section bounded by surface water at x=p/m b.
    
    Case 137.09 is the case in which constant recharge starts at t=0,
    with entry resistance.
    """
    def steady(self, R:float, dh:float | None=None, x: float | np.ndarray | None=None) -> float | np.ndarray:
        if x is None:
            # --- Return average for cross section
            aq = self.aq
            h = R * aq.b**2 / (3 * aq.kD) + R * aq.b * aq.w / aq.D
            return h
        else:
            # --- Return values at x
            x = np.atleast_1d(x)
            aq = self.aq
            h = R / (2 * aq.kD) * (aq.b**2 - x**2) + R * aq.b * aq.w / aq.D
            return h if h.size > 1 else h.item()
                
    def transient(self, time: float | np.ndarray,
                  R: float,
                  dh:float | None=None,
                  x: float | np.ndarray | None=None,
                  N: float=30) -> float | np.ndarray:
        
        assert not (isinstance(x, np.ndarray) and isinstance(time, np.ndarray)), (
                "x and time must both be arrays")
                
        aq = self.aq
        
        if np.isclose(aq.w, 0):
            # raise ValueError("aq.w must be > 0 to use this function")
            aq.w = 1e-4
        
        brug_eps   = aq.b / (aq.k * aq.w)
        alfas = [root_x_tan_x(brug_eps, k) for k in range(N)]
        
        s0 = self.steady(R=R, x=x)
        
        F = 2 * R * aq.b**2 / aq.kD
        
        if x is None:
            # --- Return X-section average heads
            for n, alpha in zip(range(N + 1), alfas):
                T = aq.b**2 * aq.mu / (aq.kD * alpha**2)             
                ds = (
                    (np.sin(alpha)**2 / alpha**4) / (1 + brug_eps / (alpha**2 + brug_eps**2))                
                    * np.exp(-time/ T)
                )
                if n == 0:
                    s = ds
                else:
                    s += ds
            return s0 - F * s          
        else:
            # --- Return head values at x
            for n, alpha in zip(range(N + 1), alfas):
                T = aq.b**2 * aq.mu / (aq.kD * alpha**2)             
                ds = (
                    (np.sin(alpha) / alpha**3) / (1 + brug_eps / (alpha**2 + brug_eps**2))
                    * np.cos(alpha * x /aq.b)
                    * np.exp(-time/ T)
                )
                if n == 0:
                    s = ds
                else:
                    s += ds
            return s0 - F * s

class Brug133(Brug):
    """Class to simulate one-layer exactly without entry resistance.
    
    Using Bruggeman(1999 solution 133.02 and 133.16)    
    """
    
    def sim_lfilter(self, rch=None, h_summer=0, h_winter=0, x=None):
        """Simulate the head using Bruggeman's 133.02 and 133.16.
        
        These solutions are without entry resistance
        
        Parameters
        ----------
        rch: pd.Series
            recharge [m/d], daily values with pd.timestamp index
        h_summer: float
            Ditch level during summer (April-September)
        h_winter: float
            Ditch water level during winter
        x: float | None
            Position in X-section  -b <= x <= b
            if None, the X-section average value will be returned.
            (For steady state X-sectio average is at x=b/sqrt(3)=0.58b)
        """
        tdata = pd.DataFrame(rch, columns=['rch'])
        summer = np.logical_and(tdata.index.month > 3, tdata.index.month < 10)
        hLR = np.zeros(len(tdata)) + h_winter
        hLR[summer] = h_summer
        tdata['hLR'] = hLR
        
        time = np.asarray((tdata.index - tdata.index[0]) / np.timedelta64(1, 'D'))
                
        Br1 = Brug13302(self.aq)
        Br2 = Brug13316(self.aq)
        
        b1 = Br1.block_response(time=time, x=None)
        b2 = Br2.block_response(time=time)
        tdata['h'] = (
               lfilter(b1, 1, tdata['hLR'] - h_winter)
             + lfilter(b2, 1, tdata['rch'])
             + h_winter
             )
        return tdata

class Brug137(Brug):
    """Class to simulate one-layer exactly without entry resistance.
    
    Using Bruggeman(1999 solution 137.02 and 137.09)    
    """
    def sim_lfilter(self, rch=None, h_summer=0, h_winter=0, x=None):
        """Simulate the head using Bruggeman's 137.02 and 137.09.
        
        These solutions are with entry resistance
        
        Parameters
        ----------
        rch: pd.Series
            recharge [m/d], daily values with pd.timestamp index
        h_summer: float
            Ditch level during summer (April-September)
        h_winter: float
            Ditch water level during winter
        x: float | None
            Position in X-section  -b <= x <= b
            if None, the X-section average value will be returned.
            (For steady state X-sectio average is at x=b/sqrt(3)=0.58b)
        """
        tdata = pd.DataFrame(rch, columns=['rch'])
        summer = np.logical_and(tdata.index.month > 3, tdata.index.month < 10)
        hLR = np.zeros(len(tdata)) + h_winter
        hLR[summer] = h_summer
        tdata['hLR'] = hLR
        
        time = np.asarray((tdata.index - tdata.index[0]) / np.timedelta64(1, 'D'))
                
        Br1 = Brug13702(self.aq)
        Br2 = Brug13709(self.aq)
        
        b1 = Br1.block_response(time=time, x=x)
        b2 = Br2.block_response(time=time, x=x)
        tdata['h'] = (
               lfilter(b1, 1, tdata['hLR'] - h_winter)
             + lfilter(b2, 1, tdata['rch'])
             + h_winter
             )
        return tdata
    
# %% Examples 
def ex_dupuit_transient(b=50, R=0.001, h0=0, hLR=0, w=0):
    """Show head development for steady inputs together with asymptote.

    This is done for w=w and w=0, for different b.
    
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
    # --- Input data for w=w and w=0
    aqww = Aquifer(k=10, D=10, c=200, w=w, mu=0.2, b=b)
    aqw0 = Aquifer(k=10, D=10, c=200, w=0, mu=0.2, b=b)
        
    title1 = str(aqww).replace(", c=200", "").replace(", b=50", "")
    title2 = f"R={R} m/d, h0={h0} m, hLR={hLR} m"
    
    time = np.logspace(0, 4, 81)
    
    fig, ax = plt.subplots(figsize=(10, 7.5))
    fig.suptitle("Dupuit case")
    ax.set_title("Head development voor stationary inputs" + "\n" + 
                 title1 + "\n" + title2)
    ax.set(xlabel='time [d]', ylabel='head [m]', xscale='log')
    
    # --- Compute and show for different b
    clrs = cycle('brgkmcy')
    for b in [30, 45, 60, 85, 120, 170, 240]:
        clr = next(clrs)
        aqww.b = b
        aqw0.b = b
        mdlww = Dupuit(aq=aqww)
        mdlw0 = Dupuit(aq=aqw0)
        
        Tww = aqww.T1L
        Tw0 = aqw0.T1L
        
        hww = mdlww.transient(time=time, R=R, h0=h0, hLR=hLR)
        hinfww = mdlww.asymptote(R=R, hLR=hLR)
        hw0 = mdlw0.transient(time=time, R=R, h0=h0, hLR=hLR)
        hinfw0 = mdlw0.asymptote(R=R, hLR=hLR)

        ax.plot(time, hww, '-',  color=clr, label=f'b={aqww.b:5.0f} d, muT={aqww.mu * Tww:6.3g} d, w={aqww.w} d')
        ax.plot(time[-1], hinfww, 'o', mfc=clr)
        ax.plot(time, hw0, '--', color=clr, label=f'b={aqw0.b:5.0f} d, muT={aqw0.mu * Tw0:6.3g} d, w={aqw0.w} d')
        ax.plot(time[-1], hinfw0, 'o', mfc=clr)

    ax.grid(True)
    ax.legend(loc="upper left")    

def ex_dupuit_transient_pd(rch=None, bs=None, h0=0, h_summer=-0.9, h_winter=-1.1):
    """Show transient head development driven by meteo
    
    The computations are done for two different b.
    
    Parameters
    ----------
    tdata: pd DataFrame with fields RH (precip) and EV24 (evapotranspiraton)
        The meteo data
    bs: np.ndarray of float [L]
        half-widths of the X-section to test
    h0: float
        Initial head at t=0
    h_winter: float
        Ditch water level during winter.
    h_summer: float
        Ditch water level during summer.
    """
    # --- Input data
    aq = Aquifer(k=10, D=10, c=200, w=0., mu=0.2, b=50)
    title1 = str(aq).replace(", c=200", "").replace(", b=50", "")
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
    # --- For given values of b
    for b in bs:
        clr = next(clrs)
        
        # --- Accomodate b
        aq.b = b
        mdl = Dupuit(aq=aq)
        
        T = aq.T1L # includes entry resistance
        
        tdata = mdl.transient_pd(rch, h_summer=h_summer, h_winter=h_winter)
        
        ax.plot(tdata.index, tdata['h'], color=clr, lw=0.5,
                label=f'b={b} m, muT={aq.mu * T:8.3g} d')
    ax.grid(True)
    ax.legend(loc="lower right")
  
def ex_sim_lfilter_dupuit(rch=None, bs=None, h0=0, h_summer=0.9, h_winter=-1.1):
    """Show transient head development driven by meteo
    
    The computation is done in two ways:
    
    1. By convolution
    2. By simulation
    
    The results are the same. One can  use rch *= 0 (see below), to
    simulate with zero recharge bu with varying hLR.
    
    Parameters
    ----------
    tdata: pd DataFrame with fields RH (precip) and EV24 (evapotranspiraton)
        The meteo data
    bs: np.ndarray of float [L]
        half-widths of the X-section
    h0: float
        Initial head at t=0
    h_winter: float
        Ditch water level during winter.
    h_summer: float
        Ditch water level during summer.
    """
    # --- Input data
    aq = Aquifer(k=10, D=10, c=200, w=0., mu=0.2, b=50)
        
    title1 = str(aq).replace(", c=200", "").replace(", b=50", "")
    title2 = f"h0={h0} m, h_summer={h_summer} m, h_winter={h_winter} m"
    
    # --- Get meteo data
    if rch is None:
        meteo = ggor_meteo.Meteo()
        rch = meteo.recharge
        
    # test
    # rch *= 0.
        
    fig, ax = plt.subplots(figsize=(10, 7.8))
    fig.suptitle("Dupuit case by simulation and by convolution")
    ax.set_title("Head driven by meteo data, using lfilter" + "\n" + 
                 title1 + "\n" + title2)
    ax.set(xlabel='time', ylabel='head [m]', xscale='linear')
    
    clrs = cycle('brgkmcy')
    for b in bs:
        clr = next(clrs)
        
        aq.b = b
        mdl = Dupuit(aq=aq)
        T = mdl.aq.T1L

        # --- Compuation by ordinary simulation
        tdata = mdl.transient_pd(rch, h_summer=h_summer, h_winter=h_winter)
        
        ax.plot(tdata.index, tdata['h'], '.', color=clr, lw=0.5, ms=2,
                label=f'By simulation,  b={b} m, muT={aq.mu * T:8.3g} d')

        # --- Computation by convolution (using scipy.signal.lfilter)
        tdata = mdl.sim_by_lfilter(rch, h_summer=h_summer, h_winter=h_winter)
        
        ax.plot(tdata.index, tdata['h'], color=clr, lw=0.5,
                label=f'By convolution, b={b} m, muT={aq.mu * T:8.3g} d')

    ax.grid(True)
    ax.legend(loc="lower right")

def ex_sim_lfilter_base_case(rch=None, bs=None, h0=0, h_summer=0.9, h_winter=-1.1, q=0):
    """Show transient head development driven by meteo.
    
    The computation is done bye
    
    1. Direct simulation.
    2. Convolution using scipy.signal.lfilter
    
    You may use the rch *=0 below to see the effect of hLR alone.
    
    Parameters
    ----------
    tdata: pd DataFrame with fields RH (precip) and EV24 (evapotranspiraton)
        The meteo data
    bs: np.ndarray of float [L]
        half-widths of the X-section
    h0: float
        Initial head at t=0
    h_winter: float
        Ditch water level during winter.
    h_summer: float
        Ditch water level during summer.
    """
    # --- Input data
    aq = Aquifer(k=10, D=10, c=200, w=0., mu=0.2, b=50)
        
    title1 = str(aq).replace(", c=200", "").replace(", b=50", "")
    title2 = f"h0={h0} m, h_summer={h_summer} m, h_winter={h_winter} m"
    
    # --- Get meteo data
    if rch is None:
        meteo = ggor_meteo.Meteo()
        rch = meteo.recharge
        
    # test    
    # rch *= 0
        
    fig, ax = plt.subplots(figsize=(10, 7.5))
    fig.suptitle("Base-case computed by direct simulation and by convolution")
    ax.set_title("Head driven by meteo data, using lfilter" + "\n" + 
                 title1 + "\n" + title2)
    ax.set(xlabel='time', ylabel='head [m]')
    
    clrs = cycle('brgkmcy')
    for b in np.atleast_1d(bs):
        clr = next(clrs)
        
        aq.b = b
        mdl = Base_case(aq=aq)
        T = mdl.aq.T2L

        # --- Direct simulation.
        tdata = mdl.transient_pd(rch, h_summer=h_summer, h_winter=h_winter, q=q)
        ax.plot(tdata.index, tdata['h'], '.', color=clr, lw=0.5, ms=2,
                label=f'By simulation,  b={b} m, muT={aq.mu * T:8.3g} d')

        # --- Convolution using lfilter.
        tdata = mdl.sim_by_lfilter(rch, h_summer=h_summer, h_winter=h_winter, q=q)
        ax.plot(tdata.index, tdata['h'], color=clr, lw=0.5,
                label=f'By convolution, b={b} m, muT={aq.mu * T:8.3g} d')

    ax.grid(True)
    ax.legend(loc="lower right")
  
def ex_2cases(rch=None, test=False, c=None, bs=None, w=None, h0=0, h_summer=0.9, h_winter=-1.1, q=0):
    """Show transient head development driven by meteo in two cases.
    
    Compare two cases:
    1. Dupuit
    2. Base-case
    
    It is shown that, for the same recharge, the head computed by simulation and
    by convolution (lfilter) are exactly the same (for large c and same w).
    
    Parameters
    ----------
    tdata: pd DataFrame with fields RH (precip) and EV24 (evapotranspiraton)
        The meteo data
    bs: np.ndarray of float [L]
        half-widths of the X-section
    h0: float
        Initial head at t=0
    h_winter: float
        Ditch water level during winter.
    h_summer: float
        Ditch water level during summer.
    """
    # --- Input data
    aq = Aquifer(k=10, D=10, c=c, w=w, mu=0.2, b=50)
        
    title1 = str(aq).replace(", b=50", "")
    title2 = f"h0={h0} m, h_summer={h_summer} m, h_winter={h_winter} m"
    
    # --- Get meteo data
    if rch is None:
        meteo = ggor_meteo.Meteo()
        rch = meteo.recharge
       
    rch1 = rch.copy() 
    rch2 = rch.copy()
    
    # --- Test, varies recharge constant during seasons.
    # --- For large c, both cases yield the same results
    
    summer = np.logical_and(rch.index.month > 3, rch.index.month < 10)
    if test==1:
        rsummer = 0.001
        h_summer, h_winter = 0, 0   
        r = np.zeros(len(rch))
        r[summer] = 0.001
        rch1[:] = r.copy()
        rch2[:] = r.copy()
        title3 = f"Recharge varies by season, hLR constant between 0 and {rsummer} m/d."
    elif test==2:
        rch1 *= 0.
        rch2 *= 0.
        title3 = "No recharge only hLR varies by season."
    else:
        title3 = "Recharge and hLR vary continuously"
        
    fig, ax = plt.subplots(figsize=(10, 9.5))
    fig.suptitle("Base case")
    ax.set_title("Head driven by meteo data, using lfilter" + "\n" + 
                 title1 + "\n" + title2 + "\n" + title3)
    ax.set(xlabel='time', ylabel='head [m]')
    
    for b in np.atleast_1d(bs):
        
        aq.b = b
        mdl1 = Dupuit(aq=aq)
        mdl2 = Base_case(aq=aq)
        
        T1 = mdl1.aq.T1L
        T2 = mdl2.aq.T2L

        if True:
            # --- Dupuit directly.           
            tdata1 = mdl1.transient_pd(rch1, h_summer=h_summer, h_winter=h_winter, q=q)
            ax.plot(tdata1.index, tdata1['h'], '.', color='b', lw=0.5, ms=2,
                    label=f'Dupuit    direct  muT={aq.mu * T1:8.3g} d')
            
            # --- Base-case directly
            tdata2 = mdl2.transient_pd(rch2, h_summer=h_summer, h_winter=h_winter, q=q)
            ax.plot(tdata2.index, tdata2['h'], '.', color='r', lw=0.5, ms=2,
                    label=f'Base_case direct  muT={aq.mu * T2:8.3g} d')
        
        if True:
            # --- Dupuit by convolution
            tdata1 = mdl1.sim_by_lfilter(rch1, h_summer=h_summer, h_winter=h_winter, q=q)
            ax.plot(tdata1.index, tdata1['h'], color='b', lw=1.5,
                    label=f'Dupuit    lfilter  muT={aq.mu * T1:8.3g} d')
            
            # --- Base-case by convolution
            tdata2 = mdl2.sim_by_lfilter(rch2, h_summer=h_summer, h_winter=h_winter, q=q)
            ax.plot(tdata2.index, tdata2['h'], color='r', lw=0.5,
                    label=f'Base-case lfilter  muT={aq.mu * T2:8.3g} d')

    ax.grid(True)
    ax.legend(loc="lower right")
  

def ex_base_case_steady(phi=0, hLR=0, R=0.001, w=0):
    """Show the base case steady and compare with single layer."""
    
    from itertools import cycle
    
    # --- Input data
    aq = Aquifer(k=10, D=10, c=200, w=w, mu=0.2, b=50)
    
    # --- Points along X-section
    x = np.linspace(-aq.b, aq.b, 101)
    
    # --- Setup figure
    fig, ax = plt.subplots(figsize=(10, 7))
    fig.suptitle("Base case, steady-state")
    
    ax.set_title('Analytical base-case, steady state' + "\n"
                 + str(aq).replace(', c=200', '')
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
        dup = Dupuit(aq=aq)
    
        # --- head along X-section and in center x=0
        hx = mdl.steady(phi=phi, hLR=hLR, R=R, x=x)
        h0 = mdl.steady(phi=phi, hLR=hLR, R=R, x=0)
        hm = mdl.steady(phi=phi, hLR=hLR, R=R)
        
        ax.plot(x, hx, color=clr, label=f'c={c} d')
        ax.plot(0, h0, 'o', ms=8, mec=clr, mfc='none')
        ax.plot(x, np.zeros_like(x) + hm, color=clr, label='hm')
    
    # --- Steady analytical solution single layer, with no leakage.    
    ax.plot(x, dup.steady(R=R, hLR=hLR, x=x), '.', color='k', lw=0.5, label='steady one-layer')
    ax.plot(x, np.zeros_like(x) + dup.steady(R=R, hLR=hLR), 'o', ms=10, mec='k', mfc='none', lw=0.5, label='steady one-layer, mean')    
    ax.grid(True)
    ax.legend(loc='upper right')
    plt.show()
    
def ex_base_case_steady_1(phi=0., hLR=0., R=0.01, c=200, w=1):
    """Show the base case equivalence if specified by phi or q + entry resistance.
    
    The X-section average head and the head along the section are computed and shown,
    with arbitrary entry resistance w >= 0.
    
    The section is computed for given phi and for given q.
    The computation is done using the extended analytical fomulas
    derived from which the x-section average head was eliminated
    (implimented right in this exmaple).
    And it is done for the implementation in the class "Base-case"
    where the compuation is split in two parts, with the
    X-section average head first computed and the head along the
    X-section thereafter using the result.

    The result shows that the same result is achieved.
    For the case of c -> infty, just compare with the
    Dupuit solution.
    """
    
    from itertools import cycle
    
    # --- Input data
    aq = Aquifer(k=10, D=10, c=c, w=w, mu=0.2, b=50)
    
    # --- Points along X-section
    x = np.linspace(-aq.b, aq.b, 26)
    
    # --- Setup figure
    fig, ax = plt.subplots(figsize=(10, 7))
    fig.suptitle("Base case, steady-state")
    
    ax.set_title('Analytical base-case, steady state' + "\n"
                 + str(aq)
                 )
    ax.set(xlabel='x [m]', ylabel='h [m]')
    
    # --- Simulate and show for different c-values.
    clrs = cycle('brgkmcy')
    
    # --- Manage graph color.
    clr = next(clrs)

    # --- Replace c, this gives new aquifer and new model.
    aq.c = c
    mdl = Base_case(aq=aq)    
    dup = Dupuit(aq=aq)    
    
    # --- Derived analytical stuff
    hfm = hLR + (phi - hLR + R * aq.c) * (aq.G + aq.wbcD) / (aq.G + aq.wbcD + 1)
    
    q = (phi - hfm) / aq.c
    
    hqm = hLR + (R + q) * aq.c * (aq.G + aq.wbcD)

    chx = np.cosh(x / aq.lam) / np.cosh(aq.b / aq.lam)
    hfx = hLR + (phi - hLR + R * aq.c) * (
        1 + (aq.wbcD / (aq.G + aq.wbcD + 1) - 1) * chx
    )
    hqx = hLR + (R + q) * aq.c * ((aq.G +aq.wbcD + 1) -(aq.G + 1) * chx)

    ax.plot(x, np.zeros_like(x) + hfm, '+', ms=10, mfc='none', label='hfm', )
    ax.plot(x, np.zeros_like(x) + hqm, 'x', ms=10, mfc='none', label='hqm')
    ax.plot(x, hfx, 's', ms=10, mfc='none', label='hfx')
    ax.plot(x, hqx, '*', ms=10, mfc='none', label='hqx')
    
    hdup = dup.steady(R=R, hLR=hLR, x=x)
    hdum = dup.steady(R=R, hLR=hLR)
    ax.plot(x, hdup, '.-', label='Dupuit x (the same for c=infty)')
    ax.plot(x, np.zeros_like(x) + hdum, '.-', label='Dupuit mean')
    

    # --- head along X-section and in center x=0
    hx = mdl.steady(R=R, phi=phi, hLR=hLR, x=x)
    hm = mdl.steady(R=R, phi=phi, hLR=hLR)
    q = (phi - hm) / aq.c
    hq = mdl.steady(R=R, hLR=hLR, q=q, x=x)        
    haq = mdl.steady(R=R, q=q, hLR=hLR)
    
    ax.plot(x, hx, color=clr, label=f'h_phi, c={c} d')
    ax.plot(x, hm + np.zeros_like(x), '-', color=clr, label=f'hm_phi, c={c} d')
    
    clr = next(clrs)
    ax.plot(x, hq, '.', color=clr, label=f'h_q, c={c} d')
    ax.plot(x, haq + np.zeros_like(x), '.', color=clr, label=f'hm_q, c={c} d')
    
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
        
    title1 = str(aq).replace(", c=200", "")
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

def ex_base_case_transient_pd(rch=None, cs=None, b=50, h0=0, h_summer=-0.9, h_winter=-1.1, q=0):
    """Show transient head development driven by meteo
    
    Parameters
    ----------
    rch: pd.Series of daily recharge)
        The meteo data
    cs: np.ndarray of floats
        c values to test
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
        
    title1 = str(aq).replace(", c=200", "")
    title2 = f"h0={h0} m, h_summer={h_summer} m, h_winter={h_winter} m, q={q} m/d"
    
    # --- Get meteo data
    if rch is None:
        meteo = ggor_meteo.Meteo()
        rch = meteo.recharge
        
    fig, ax = plt.subplots(figsize=(10, 7.5))
    fig.suptitle("Base case, transient")
    ax.set_title("Head driven by meteo data" + "\n" + 
                 title1 + "\n" + title2)
    ax.set(xlabel='time', ylabel='head [m]', xscale='log')
    
    clrs = cycle('brgkmcy')
    for c in cs: # 30, 100, 300, 1000, 3000]:
        clr = next(clrs)
        aq.c = c
        mdl = Base_case(aq=aq)
        
        T = aq.T2L # includes entry resistance
        
        tdata = mdl.transient_pd(rch, h_summer=h_summer, h_winter=h_winter, q=q)
        
        ax.plot(tdata.index, tdata['h'], color=clr, lw=0.5,
                label=f'c={aq.c:7.0f} d, muT={aq.mu * T:.1f} d')
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

def ex_brug13302_13702_Dupuit(xs_b=None, w=None, dh=0.1):
    """Compare Bruggeman's solution 133.102 with 137.02 and Dupuit.
    
    Solutions 133.02 (wihtout entry resistance) and 137.02 (with
    entry resistan ce) are for a sudden change of head at the
    boundaries, i.e. at x=+/- b. De solution of Dupuit can be
    made to simulate that by setting recharge R=0.
    
    Parameters
    ----------
    xs_b: np.ndarray of floats
        relative coordinates, xs_b = x/b
    w: float
        entry resistance [d] to be applied to 137.02 and Dupuit.
    dh: float
        Sudden head change at x=0.    
    """
    aq = Aquifer(k=10, D=10, c=200, w=w, mu=0.15, b=50)
    xs = np.atleast_1d(xs_b) * aq.b
    
    time = np.linspace(0, 20, 101)
    
    brug1 = Brug13302(aq)
    brug2 = Brug13702(aq)
    dup   = Dupuit(aq)

    fig, ax = plt.subplots(figsize=(10, 7.5))
    
    fig.suptitle(f"Bruggeman 133.02 (w=0), 137.02 (w={w}) and Dupuit (w={w} d)")
    ax.set_title(f"Sudden rise {dh} m of surface water. With/without entry resistance" +
                 "\n" + str(aq).replace(", c=200",""))
    ax.set(xlabel='t d[]', ylabel='h [m]')
    
    for x in xs:        
        h1t = brug1.transient(time=time, dh=dh, x=x)
        h1a = brug1.transient(time=time, dh=dh)
        
        h2t = brug2.transient(time=time, dh=dh, x=x)
        h2a = brug2.transient(time=time, dh=dh)
        
        hd1 = dup.transient(time=time, R=0., h0=0, hLR=dh)
        hds = dup.steady(R=0, hLR=dh)
        
        ax.plot(time, h1t, '-', color='blue',
                label=f'Brug133.02, w={0} d, x={x:.1f} m')
        ax.plot(time, h1a, '-', color='blue',
                label=f'Bug133.02, w={0} d, x=None')

        ax.plot(time, h2t, '.', color='red',
                label=f'Brug137.02, w={w} d, x={x:.1f} m')
        ax.plot(time, h2a, 'o', mec='red', mfc='none',
                label=f'Brug137.02, w={w} d, x=None')
        
        ax.plot(time, hd1, 'o', mec='green', mfc='none',        
                label=f"Dupuit tr, w={w} d, muT={dup.aq.mu * dup.aq.T1L:.1f} d")
        ax.plot(time[[0, -1]], [hds, hds], '--', color='green',
                label=f"Dupuit st, w={w} d")
    ax.grid(True)
    ax.legend(loc='center')


def ex_brug133_137_vs_dupuit(rch=None, test=None, b=100, h_summer=-0.9, h_winter=-1.1, w=0.001):
    """Compare results 133.02/16, 137.02/09 with Dupuit
    
    This combines varyring recharge with varying hLR
    
    Compare the X-section average values
    
    Compute using lfilter
    """
    aq = Aquifer(k=10, D=10, c=200, w=w, mu=0.15, b=b)
    
    title2 = str(aq).replace(", c=200", "")
    
    # --- Get meteo data
    if rch is None:
        meteo = ggor_meteo.Meteo()
        rch = meteo.recharge

    summer = np.logical_and(rch.index.month > 3, rch.index.month < 10)
    if test==1:
        R = 0.001
        r = np.zeros(len(rch))
        r[summer] = R
        rch[:] = r
        h_summer, h_winter = 0, 0
        title3=f"Recharge varies by season between 0 and {R} m/d"
    if test==2:
        rch *= 0
        title3 = "No recharge, hLR varies by season."
    else:
        title3 = "Recharge varies continuously and hLR by season."
    
    brug133 = Brug133(aq)
    brug137 = Brug137(aq)
    dup  = Dupuit(aq)

    fig, ax = plt.subplots(figsize=(10, 9.5))
    
    fig.suptitle("Compare Bruggeman 133.02/16, 137.02/09 with Dupuit (all lfilter))")
    ax.set_title(title2 + "\n" + title3)
    ax.set(xlabel='t d[]', ylabel='h [m]')
    
    tdata3 = brug133.sim_lfilter(rch=rch, h_summer=h_summer, h_winter=h_winter)
    tdata7 = brug137.sim_lfilter(rch=rch, h_summer=h_summer, h_winter=h_winter)        
    hd = dup.transient_pd(rch, h_summer, h_winter)
    
    ax.plot(rch.index, tdata3['h'], lw=0.5, label=f"Brug133   (w={0} d)")
    ax.plot(rch.index, tdata7['h'], lw=0.5, label=f"Brug137   (w={w} d")
    ax.plot(rch.index, hd['h'],    '.', lw=0.5, ms=2, label=f"Dupuit tr (w={w} d)")
        
    ax.grid(True)
    ax.legend(loc='center')
   
    
def ex_brug13716_13309_avg_dupuit(xs_b=None, w=None):
    """Compare Bruggeman's 133.16 with 133.09 and Dupuit.
    
    Does this for the x where the steady head equals the
    X-section average head without entry resistance.
    
    The average head occurs wehere x=b/sqrt(3)=0.577 b
     For convenience, just specify the xs and the w
     
     Bruggeman 133.16 is without entry resistance.
     Bruggeman 137.09 is without entry resistance.
     
     Parameters
     ----------
     xs_b: np.ndarray of float
        array of x/b values.
    w: float
        entry resistance for case Bruggeman 133.09 and Dupuit.     
    """
    
    aq = Aquifer(k=10, D=10, c=200, w=w, mu=0.15, b=50)
    xs = np.atleast_1d(xs_b) * aq.b
    R = 0.001
    time = np.linspace(0, 20, 101)
    
    brug1 = Brug13316(aq)
    brug2 = Brug13709(aq)
    dup  = Dupuit(aq)

    fig, ax = plt.subplots(figsize=(10, 7.5))
    
    fig.suptitle(f"Bruggeman 133.16 (w=0), 137.09 (w={w}) and Dupuit (w={w} d)")
    ax.set_title(f"Constant precipitation {R} m/d. With/without entry resistance." +
                 "\n" + str(aq).replace(", c=200",""))
    ax.set(xlabel='t d[]', ylabel='h [m]')
    
    for x in xs:        
        h1t = brug1.transient(time=time, R=R, x=x)
        h1a = brug1.transient(time=time, R=R)
    
        h2t = brug2.transient(time=time, R=R, x=x)
        h2a = brug2.transient(time=time, R=R)
    
        hdt = dup.transient(time=time, R=R)
        hd  = dup.steady(R=R, hLR=0, x=x)
        
        ax.plot(time, h1t, '-', color='blue',
                label=f'Brug133.16, w={0} d, x={x:.1f} m')
        ax.plot(time, h1a, 's', mec='blue', mfc='none',
                label= f'Brug137.16, w={0} d, X-sec. avg')

        ax.plot(time, h2t, '.', color='red',
                label=f'Brug137.09, w={w} d, x={x:.1f} m')
        ax.plot(time, h2a, '*', mec='red', mfc='none',
                label= f'Brug137.09, w={w} d, X-sec. avg')
        
        ax.plot(time[[0, -1]], [hd, hd], '--', color='green',
                label=f"Dupuit st, w={w} d, x={x:.1f} m")
        ax.plot(time, hdt, 'o', mec='green', mfc='none',
                label=f'Dupuit tr, w={w} d, X-sec. avg')
        
    ax.grid(True)
    ax.legend(loc='center')

 # %% __main__   
if __name__ == "__main__":
    if False:
        ex_base_case_steady(phi=0, hLR=0, R=0.001, w=1)
        ex_base_case_steady_1(phi=0.0, hLR=0., R=0.01, c=200, w=1)
    if False:
        ex_dupuit_transient(b=50, R=0.001, h0=0, hLR=0, w=1)
        ex_dupuit_transient_pd(rch=None, bs=[50, 250], h0=0, h_summer=-0.9, h_winter=-1.1)
    if False:
        ex_base_case_transient_pd(cs=[50, 250], b=50,   h0=0, h_summer=0, h_winter=0, q=0)
        ex_base_case_transient_pd(cs=[50, 250], b=500,  h0=0, h_summer=0, h_winter=0, q=0)
        ex_base_case_transient_pd(cs=[50, 250], b=5000, h0=0, h_summer=0, h_winter=0, q=0)    
        ex_base_case_transient_pd(cs=[50, 250], b=50,   h0=3, h_summer=0, h_winter=0, q=0)
        ex_base_case_transient_pd(cs=[50, 250], b=50,   h0=0, h_summer=0, h_winter=0, q=0)
        ex_base_case_transient_pd(cs=[50, 250], b=50,   h0=0, h_summer=0, h_winter=0, q=0.001)
    if False:
        compare_limits()
    if False:               
        ex_brug13302_13702_Dupuit(xs_b=[np.sqrt(1/3)], w=1.0, dh=0.1)
        ex_brug13716_13309_avg_dupuit(xs_b=[np.sqrt(1/3)], w=0.5)    
    if False:
        ex_sim_lfilter_dupuit(rch=None, bs=[50, 150, 300], h0=0, h_summer=-0.9, h_winter=-1.1)
        ex_sim_lfilter_base_case(bs=[50, 150, 250],h0=0, h_summer=-0.9, h_winter=-1.1, q=0)
        ex_2cases(test=1, bs=[50, 150, 250], h0=0, c=50000, w=0.5, h_summer=-0.9, h_winter=-1.1, q=0)
        ex_2cases(test=2, bs=[50, 150, 250], h0=0, c=50000, w=0.5, h_summer=-0.9, h_winter=-1.1, q=0)
        ex_2cases(test=0, bs=[50, 150, 250], h0=0, c=50000, w=0.5, h_summer=-0.9, h_winter=-1.1, q=0)
    if True:  
        ex_brug133_137_vs_dupuit(rch=None, test=2, b=250, h_summer=-0.9, h_winter=-1.1, w=5)
        ex_brug133_137_vs_dupuit(rch=None, test=1, b=250, h_summer=-0.9, h_winter=-1.1, w=5)
        ex_brug133_137_vs_dupuit(rch=None, test=0, b=250, h_summer=-0.9, h_winter=-1.1, w=5)
    try:
        plt.show()
    except Exception:
        pass
    print('Done')
    
# %%
    

