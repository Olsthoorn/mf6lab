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
            
        T = aq.T1L
    
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
        h = np.zeros_like(time)
        
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
        wb_cD = aq.w * aq.b / (aq.c * aq.D)
                
        if phi is not None:        
            havg = hLR + (R * aq.c + (phi -hLR)) * (aq.G + wb_cD) / (aq.G + wb_cD + 1)
            q = (phi - havg) / aq.c
        else:
            havg = hLR + (R + q) * aq.c * (aq.G + wb_cD)
            phi = havg + q * aq.c
        
        if x is None:
            return havg
        else:        
            x = np.atleast_1d(x)
            chx = np.cosh(x / aq.lam)
            chb = np.cosh(aq.b / aq.lam)
            
            hb = hLR + (R + q) * aq.c * wb_cD
            
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
            h = dh * np.zeros_like(x)
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
            raise ValueError("aq.w must be > 0 to use this function")
    
        brug_eps   = aq.b / (aq.k * aq.w)
        alfas = [root_x_tan_x(brug_eps, k) for k in range(N)]
        
        s0 = self.steady(dh=dh, x=x)
        
        if x is None:
            for n, alpha in zip(range(N + 1), alfas):
                T = aq.b**2 * aq.mu / (aq.kD * alpha**2)             
                ds = (
                    np.sin(alpha)**2 / (1 + brug_eps / (alpha**2 + brug_eps**2))                
                    * np.exp(-time/ T)
                )
                if n == 0:
                    s = ds
                else:
                    s += ds
            return s0 - s0 * s
        else:
            for n, alpha in zip(range(N + 1), alfas):
                T = aq.b**2 * aq.mu / (aq.kD * alpha**2)             
                ds = (
                    (np.sin(alpha) / alpha) / (1 + brug_eps / (alpha**2 + brug_eps**2))
                    * np.cos(alpha * x /aq.b)
                    * np.exp(-time/ T)
                )
                if n == 0:
                    s = ds
                else:
                    s += ds
        return s0 - 2 * s0 * s


class Brug13709(Brug):
    """One-layer transient cross section bounded by surface water at x=p/m b.
    
    Case 137.09 is the case in which constant recharge starts at t=0,
    with entry resistance.
    """

    """Sudden rise h of the surface water level."""    
    
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
            raise ValueError("aq.w must be > 0 to use this function")
        
        brug_eps   = aq.b / (aq.k * aq.w)
        alfas = [root_x_tan_x(brug_eps, k) for k in range(N)]
        
        s0 = self.steady_avg(R=R, x=x)
        
        if x is None:
            # --- Return X-section average heads
            F = R * aq.b**2 / aq.kD        

            for n, alpha in zip(range(N + 1), alfas):
                T = aq.b**2 * aq.mu / (aq.kD * alpha**2)             
                ds = (
                    (np.sin(alpha)**2 / alpha**2) / (1 + brug_eps / (alpha**2 + brug_eps**2))                
                    * np.exp(-time/ T)
                )
                if n == 0:
                    s = ds
                else:
                    s += ds
            return s0 - F * s          
        else:
            # --- Return head values at x        
            F = 2 * R * aq.b**2 / aq.kD        

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
  
def ex_sim_lfilter_dupuit(rch=None, b=50, h0=0, h_summer=0.9, h_winter=-1.1):
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
        
    # test
    rch *= 0.
        
    fig, ax = plt.subplots(figsize=(10, 7.8))
    fig.suptitle("Dupuit case")
    ax.set_title("Head driven by meteo data, using lfilter" + "\n" + 
                 title1 + "\n" + title2)
    ax.set(xlabel='time', ylabel='head [m]', xscale='log')
    
    clrs = cycle('brgkmcy')
    for b in [250]:
        clr = next(clrs)
        aq.b = b
        mdl = Dupuit(aq=aq)
        T = mdl.aq.T1L

        tdata = mdl.sim_by_lfilter(rch, h_summer=h_summer, h_winter=h_winter)
        ax.plot(tdata.index, tdata['h'], color=clr, lw=0.5,
                label=f'lfilter  muT={aq.mu * T:8.3g} d')

        tdata = mdl.transient_pd(rch, h_summer=h_summer, h_winter=h_winter)
        ax.plot(tdata.index, tdata['h'], '.', color=clr, lw=0.5,
                label=f'direct  muT={aq.mu * T:8.3g} d')
    ax.grid(True)
    ax.legend(loc="lower right")

def ex_sim_lfilter_base_case(rch=None, b=50, h0=0, h_summer=0.9, h_winter=-1.1, q=0):
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
        
    # test    
    rch *= 0
        
    fig, ax = plt.subplots(figsize=(10, 7.8))
    fig.suptitle("Base case")
    ax.set_title("Head driven by meteo data, using lfilter" + "\n" + 
                 title1 + "\n" + title2)
    ax.set(xlabel='time', ylabel='head [m]', xscale='log')
    
    clrs = cycle('brgkmcy')
    for b in [250]:
        clr = next(clrs)
        aq.b = b
        mdl = Base_case(aq=aq)
        T = mdl.aq.T2L
        
        tdata = mdl.sim_by_lfilter(rch, h_summer=h_summer, h_winter=h_winter, q=q)
        ax.plot(tdata.index, tdata['h'], color=clr, lw=0.5,
                label=f'lfilter  muT={aq.mu * T:8.3g} d')

        tdata = mdl.transient_pd(rch, h_summer=h_summer, h_winter=h_winter, q=q)
        ax.plot(tdata.index, tdata['h'], '.', color=clr, lw=0.5,
                label=f'direct  muT={aq.mu * T:8.3g} d')
    ax.grid(True)
    ax.legend(loc="lower right")
  
def ex_2cases(rch=None, c=50, b=50, h0=0, h_summer=0.9, h_winter=-1.1, q=0):
    """Show transient head development driven by meteo in two cases.
    
    First case is Dupuit, the second is base-case. It is shown that
    for the same recharge, the head computed by simulation and
    by convolution (lfilter) are exactly the same.
    
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
    aq = Aquifer(k=10, D=10, c=c, w=0., mu=0.2, b=b)
        
    title1 = str(aq).replace(", b = 50", "")
    title2 = f"h0={h0} m, h_summer={h_summer} m, h_winter={h_winter} m"
    
    # --- Get meteo data
    if rch is None:
        meteo = ggor_meteo.Meteo()
        rch = meteo.recharge
       
    rch1 = rch.copy() 
    rch2 = rch.copy()
    # test
    summer = np.logical_and(rch.index.month > 3, rch.index.month < 10)
    r = np.zeros(len(rch))
    r[summer] = 0.001
    rch1[:] = r.copy()
    rch2[:] = r.copy()
    
    fig, ax = plt.subplots(figsize=(10, 7.8))
    fig.suptitle("Base case")
    ax.set_title("Head driven by meteo data, using lfilter" + "\n" + 
                 title1 + "\n" + title2)
    ax.set(xlabel='time', ylabel='head [m]', xscale='log')
    
    for b in [250]:        
        aq.b = b
        mdl1 = Dupuit(aq=aq)
        mdl2 = Base_case(aq=aq)
        
        T1 = mdl1.aq.T1L
        T2 = mdl2.aq.T2L
        
        if True:
            tdata1 = mdl1.sim_by_lfilter(rch1, h_summer=h_summer, h_winter=h_winter, q=q)
            ax.plot(tdata1.index, tdata1['h'], color='b', lw=1.5,
                    label=f'Dupuit lfilter  muT={aq.mu * T1:8.3g} d')
            tdata2 = mdl2.sim_by_lfilter(rch2, h_summer=h_summer, h_winter=h_winter, q=q)
            ax.plot(tdata2.index, tdata2['h'], color='r', lw=0.5,
                    label=f'Base-case lfilter  muT={aq.mu * T2:8.3g} d')
        if True:
            tdata1 = mdl1.transient_pd(rch1, h_summer=h_summer, h_winter=h_winter, q=q)
            ax.plot(tdata1.index, tdata1['h'], '.', color='b', lw=0.5,
                    label=f'Dupuit direct  muT={aq.mu * T1:8.3g} d')
            tdata2 = mdl2.transient_pd(rch2, h_summer=h_summer, h_winter=h_winter, q=q)
            ax.plot(tdata2.index, tdata2['h'], '.', color='r', lw=0.5,
                    label=f'Base_case direct  muT={aq.mu * T2:8.3g} d')

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
    
    ax.set_title('Analytical base-case, steady state' + "\n"
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
    
def ex_base_case_steady_1():
    """Show the base case equivalence if specified by phi or q + entry resistance."""
    
    from itertools import cycle
    
    # --- Input data
    aq = Aquifer(k=10, D=10, c=200, w=0.5, mu=0.2, b=50)
    phi, hLR, R = 0., 0., 0.1

    # --- Points along X-section
    x = np.linspace(-aq.b, aq.b, 26)
    
    # --- Setup figure
    fig, ax = plt.subplots(figsize=(10, 7))
    fig.suptitle("Base case, steady-state")
    
    ax.set_title('Analytical base-case, steady state' + "\n"
                 + str(aq).replace(' c=200,', '')
                 )
    ax.set(xlabel='x [m]', ylabel='h [m]')
    
    # --- Simulate and show for different c-values.
    clrs = cycle('brgkmcy')
    for c in [100]:
        
        # --- Manage graph color.
        clr = next(clrs)

        # --- Replace c, this gives new aquifer and new model.
        aq.c = c
        mdl = Base_case(aq=aq)
    
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
        ht = brug.transient(time=time, dh=dh, x=x)
        ha = brug.transient_avg(time=time, dh=dh)
        hd = dup.transient(time=time, R=0, h0=0, hLR=dh)
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
        h1 = brug1.transient(time=time, dh=dh, x=x)
        h2 = brug2.transient(time=time, dh=dh, x=x)
        hd = dup.transient(time=time, R=0., h0=0, hLR=dh)
        hb = brug2.transient(time=time, dh=dh)
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
        ht = brug.transient(time=time, R=R, x=x)
        hb = brug.transient(time=time, R=R)
        hd = dup.transient(time=time, R=R)
        hdup = dup.steady(x=x, hLR=0, R=R)
        
        ax.plot(time, ht, color=clr, label=f'Brug133.16, x={x:.1f} m')
        ax.plot(time, hb, 's', mec=clr, mfc='none', label='Brug133.16 avg')
        ax.plot(time, hd, 'o', mec=clr, mfc='none', label='Dupuit transient')      
        ax.plot(time[-1], hdup, 'x', mec=clr, mfc='none', label='Dupuit')
        
    ax.grid(True)
    ax.legend(loc='center')

def ex_brug133_137_dupuit(rch=None, h_summer=-0.9, h_winter=-1.1, w=0.001):
    """Compare results 133.02/16, 137.02/09 with Dupuit
    
    Compare the X-section average values
    
    Compute using lfilter
    """
    aq = Aquifer(k=10, D=10, c=200, w=w, mu=0.15, b=50)
    
        # --- Get meteo data
    if rch is None:
        meteo = ggor_meteo.Meteo()
        rch = meteo.recharge

    summer = np.logical_and(rch.index.month > 3, rch.index.month < 10)
    r = np.zeros(len(rch))
    r[summer] = 0.001
    rch[:] = r
    
    brug133 = Brug133(aq)
    brug137 = Brug13302(aq)
    dup  = Dupuit(aq)

    fig, ax = plt.subplots(figsize=(10, 7.5))
    
    fig.suptitle("Compare Bruggeman 133.02/16, 137.02/09 with Dupuit (all lfilter))")
    ax.set_title("Varying precipitation")
    ax.set(xlabel='t d[]', ylabel='h [m]')
    
    tdata3 = brug133.sim_lfilter(rch=rch, h_summer=h_summer, h_winter=h_winter)
    tdata7 = brug137.sim_lfilter(rch=rch, h_summer=h_summer, h_winter=h_winter)        
    hd = dup.transient_pd(rch, h_summer, h_winter)
    
    ax.plot(rch.index, tdata3['h'], label="Brug133 (- entry resistance)")
    ax.plot(rch.index, tdata7['h'], label="Brug137 (+ entry resistance)")
    ax.plot(rch.index, hd['h'],     label="Dupuit transient")
        
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
        ht1 = brug1.transient(time=time, R=R, x=x)
        ht2 = brug2.transient(time=time, R=R, x=x)
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
        h1t = brug1.transient(time=time, R=R, x=x)
        h2t = brug2.transient(time=time, R=R, x=x)
        h2a = brug2.transient(R=R, time=time)
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

 # %% __main__   
if __name__ == "__main__":
    if True:
        # ex_base_case_steady()
        ex_base_case_steady_1()
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
    if False:
        ex_brug13302()
        ex_brug13316()
        ex_brug13702()
        ex_brug13709()
    if False:
        ex_sim_lfilter_dupuit(h0=0, h_summer=1.0, h_winter=0)
        ex_sim_lfilter_base_case(h0=0, h_summer=1.0, h_winter=0, q=0)
        ex_2cases(h0=0, c=50, h_summer=0, h_winter=0, q=0)
        ex_brug133_137_dupuit(rch=None, h_summer=-0.9, h_winter=-1.1, w=0.001)
    try:
        plt.show()
    except Exception:
        pass
    print('Done')
    
# %%
    

