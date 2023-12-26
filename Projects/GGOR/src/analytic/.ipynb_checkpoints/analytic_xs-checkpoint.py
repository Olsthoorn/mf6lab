#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analytic computation of head in multi aquiers of cross section between x=0
and x=b, where dphi/dx = 0 at x=0 and phi is fixed beyond the ditch resistance
at x=b. Teh multi-aquifer system has nLay aquifers and nLay aquitards, one
between each pair of aquifers, one at the top of the upper aquier and one
below the bottom aquifer.

Each aquifer has a hydraulic resistance at x=b. Outside x=b, the head is fixed.
The head is also fixed above the top aquitard and below the bottom aquitard.
The specify heads are thus nLay + 2 values per time, where the first head is
the head maintained above the top aquifer, the last one is that maintained
below the bottom aquifer and the in between ones are the heads maintained
beyong the ditch resistance at x=b in each aquifer.

The heads can be simulated as a function of x, 0<=x<=b for the steady-state
situation using method phix. The steady head averaged over the cross section
can be computed with the method phim and the dynamic head averaged over the
cross section can be computed with the mehtod phit.

Aquifer class defines the multi-aquifer system and computes properties from
it's input.

Analytic_nlay class defines the analytic solution and contains the method to
compute the heads and plot them.

Created on Thu Jul  6 00:25:13 2017

@author: Theo
"""
import numpy as np
from matplotlib import pyplot as plt
import scipy.linalg as lg
import pandas as pd
#import pdb

def dots(*args):
    A = args[0]
    for arg in args[1:]:
        A = A.dot(arg)

class Multi_aquifer(object):
    """Aquifer holds the aquifer properties and can compute other ones
    from basic properties kD and c
    parameters:
    -----------
    kD = [m2/d] list of kD values of the aquifers (nLay)
    c  = [d] list of resistances of confining layers (nLay+1)
    w  = [d] list of nLay resistances between aquifer and fully penetrating ditches
    k  = [m/s] list of nLay conductivities of the aquifers
    b  = [m] half-width of cross section, scalar
    """
    def __init__(self,
                 kD=None, c=None, D=None, d=None, w=None, S=None, b=None, **kwargs):
        '''Return a multi_aquifer object.

        parameters
        ----------
        kD : np.ndarray, dim (naquif) [L2/T]
          Tranmissivities of the aquifers.
        c : np.ndarray, dim (naquif) [T]
          Vertical resistance [T] of the aquifers. First c on top of system.
        D : np.ndarray, dim (len(naquif)) [L]
          Thickness of the aquifers.
        d : np.ndarray, dim (len(c)) [L]
          Thickness of the aquitards.
        w : np.ndarray, dim (naquif) [T]
          Entry resistance of layer, implemented as inflow resistance
          in the form of a vertical wall if resitance w [d] in de aquifer.
        S : np.ndarray, dim (naquif) [-]
          Storage coefficient or Specific yield.
        '''
        try:
            self.nLay = len(D)
            self.kD =np.array(kD)
            self.c = np.array(c) # vertical resistance between aquifers + top and bottom
            self.D = np.array(D) # Thickness of aquifers
            self.d = np.array(d) # Thickness of aquitards
            self.w = np.array(w) # Entry resistance of ditches
            self.S = np.array(S) # Storage coefficients of the aquifers
            self.b = float(b)    # half-width of the cross section.

        except:
            raise Exception("kD, c, D, d, w, S and b must all be specified.")

        # There is always an aquitard at the top and bottom of the system
        assert len(self.c) == self.nLay + 1,\
            "len(c)={} != len(kD) + 1= {}".format(len(c), len(kD) + 1)
        assert len(self.d) == self.nLay + 1, "len(d)={} != len(c)  ={}".format(len(d), len(c))

        # The length of the following must equal the number of layers
        assert len(self.w) == self.nLay, "len(w)={} != len(kD) = {}".format(len(w), len(kD))
        assert len(self.D) == self.nLay, "len(D)={} != len(kD) = {}".format(len(D), len(kD))
        assert len(self.S) == self.nLay, "len(S)={} != len(kD) = {}".format(len(S), len(kD))

        assert np.isscalar(b), "cross section half-width b must be a float."

        # Below are precomputed properties for fast access
        self.kD    = np.matrix(np.diag(self.kD))
        self.S     = np.matrix(np.diag(self.S))
        self.S_m1  = self.S.I
        self.k     = np.diag(self.kD)/self.D
        self.I     = np.matrix(np.eye(self.nLay))
        self.kD_m1 =  self.kD.I
        self.Kw_m1 = np.matrix(np.diag(1./(self.k * self.w)))
        self.Kw    = self.Kw_m1.I
        self.A_m2  = self.get_sysmat_A_m2() # compute system matrix
        self.A_m1  = np.matrix(lg.sqrtm(self.A_m2))
        self.A     = self.A_m1.I
        self.A2    = self.A_m2.I
        self.F     = self.Kw * self.A_m1 * lg.sinhm(self.b * self.A_m1)\
                                + lg.coshm(self.b * self.A_m1)
        self.F_m1  = self.F.I


        # Properties neede for transient simulation of the average
        # head in each layer of the cross section.
        #self.B     = self.kD * self.A_m2 * (self.I - self.A / self.b *
        #                np.matrix(lg.sinhm(self.b * self.A_m1)) * self.F.I).I
        self.B     = self.I - self.A / self.b * lg.sinhm(self.b * self.A_m1) * self.F_m1
        self.B_m1  = self.B.I
        self.G     = self.S_m1 * self.kD * self.A_m2 * self.B_m1
        self.D, self.V  = lg.eig(self.G)
        self.D = np.matrix(np.diag(self.D))
        self.V = np.matrix(self.V)


    def get_sysmat_A_m2(self):
        '''Return the system matrix of the multi-aquifer system.

        This is the big Lambda ^ -2 in the derivation.
        In the code A instead of big Lambda is used.
        The _m2 stands for "to the power -2"
        '''
        c  = self.c
        kD = np.diag(self.kD)
        v0 = -1/(c[1:-1] * kD[1:])
        v1 = +1/(c[ :-1] * kD) + 1/(c[1:] * kD)
        v2 = -1/(c[1:-1] * kD[:-1])
        return np.matrix(np.diag(v0, -1) + np.diag(v1, 0) +np.diag(v2, 1))



class Analytic_nlay(object):
    """ Analytic multilayer groundwater flow solution for a cross section for
    0 <= x <= b.

    The groundwater system is a multilayer aquifer system consisting of nLay
    aquifers and nLay+1 aquitards, one between each pair of adjacent aquifers,
    one boven the top aquifer and one below the bottom one.

    Heads are fixed above the top and below the bottom aquitard and in all
    aquifers at x=b, outside a resistance w, which is at x=b.

    Steady-state and transient heads can be computed and plotted. The steady-
    state can be plotted as a function of x; and the transient heads as a
    function of time, but averaged over the width of the cross section.

    This multilayer groundwater-flow solution makes intense use of linear
    algebra and all involved arrays are cast to np.matrix instead of
    np.ndarray to allow using linear algebra in a natural way.

    @TO 170712
    """
    def __init__(self, aquifer):
        '''Return Analytic_nLay object

        parameters
        ----------
            aquifer : Multi_aquifer instance (see this module)
        '''
        self.aq = aquifer
        self.nLay = self.aq.nLay
        return


    def hqx(self, N=None, hLR=None, ht=None, hb=None, n=51):
        '''Return heads in the cross section.

        This function was fully tested. 20181123

        parameters
        ----------
        b: float
            half-width of cross section
        N: vector (nLay).
            Uniform infiltration rate in the layers
        hLR : vector (nLay,)
            presecrived ditch level (at x=0)
        ht: float
            head above top aquifer
        hb: float
            head below bottom aquifer
        n: int
            number of points from 0 to b (half-width of the cross section)
        returns
        -------
            h: head array (nlay, nx)
            q: discharge array (nlay, nx)
            x: array of x-avalues between 0 adn half-width b
        '''

        x = np.linspace(0, self.aq.b, n)

        g = np.zeros(self.nLay)
        if ht is not None:
            g[0] = ht / self.aq.c[0]
        if hb is not None:
            g[-1]= hb / self.aq.c[-1]

        h   = np.zeros((self.nLay, len(x)))
        q   = np.zeros((self.nLay, len(x)))

        hLR = np.array(hLR)[:, np.newaxis]
        Npg = (np.array(N) + np.array(g))[:, np.newaxis]

        for i, x_ in enumerate(x):
            h[:, i:i+1]= hLR \
                + (self.aq.I - lg.coshm(x_ * self.aq.A_m1) * self.aq.F_m1) \
                * self.aq.A2 * self.aq.kD_m1 * Npg

            q[:, i:i+1] = \
                self.aq.kD * self.aq.A_m1 * lg.sinhm(x_ * self.aq.A_m1) \
                * self.aq.F_m1 * self.aq.A2 * self.aq.kD_m1 * Npg
        return (h, q, x)


    def hqm(self, N=None, hLR=None, ht=None, hb=None):
        '''Return mean head and outflow for all layers in the cross section.

        This function was fully tested. 20181123

        parameters
        ----------
        hLR : vector (nLay,)
            presecrived ditch level (at x=0)
        N: vector (nLay).
            Uniform infiltration rate in the layers
        htb: vector of len 2
            head above top aquifer and below bottom aquifer
        returns
        -------
            h: head array (nlay, nx)
            q: discharge array (nlay, nx)
        '''

        # right-hand vector representing head at top of topaquifer and
        # below bottom of bottom aquifer

        g    = np.zeros(self.nLay)
        if ht is not None:
            g[0] = ht / self.aq.c[0]
        if hb is not None:
            g[-1]= hb / self.aq.c[-1]

        h   = np.zeros((self.nLay, len(x)))
        q   = np.zeros((self.nLay, len(x)))

        hLR = np.array(hLR)[:, np.newaxis]
        Npg = (np.array(N) + np.array(g))[:, np.newaxis]

        h = hLR + (self.aq.I - self.aq.A / self.aq.b
               * lg.sinhm(self.aq.b * self.aq.A_m1) * self.aq.F_m1) \
            * self.aq.A2 * self.aq.kD_m1 * Npg

        q = self.aq.kD * self.aq.A_m1 * lg.sinhm(self.aq.b * self.aq.A_m1) \
            * self.aq.F_m1 * self.aq.A2 * self.aq.kD_m1 * Npg

        return np.array(h)[:,0], np.array(q)[:,0]


    def hmt(self, t, N=None, hLR=None, ht=None, hb=None, **kwargs):
        """ Returns transient head, aveaged over the cross section
        and discharge to the ditches.

        parameters
        ----------
        data: pd.DataFrame of time data
            index is timestamps
        t: [ m ] np.ndarray of shape (nt,)
            a vector of times to compute the average head in x-section
        hLR: [ m ] np.ndarray of shape (nLay, nt)
            specified heads at x = b
        N: [m/d], np.ndarray of shape (nLay, nt)
            specified injection in the aquifers.
        ht: [ m ] np.ndaray of shape (nt,)
            maintained head above top aquifer and below bottom aquifer.
        hb: [ m ] np.ndarray of shape (nt,)
            maintained head below the bottom aquifer
        """

        assert t.ndim == 1, "t must be a vector with times"
        nt = len(t)

        if ht is not None:
            if np.isscalar(ht):
                ht = np.zeros(nt) * ht
            else:
                assert np.all(hb.shape == t.shape), 'hb must have same shape as t'
        if hb is not None:
            if np.isscalar(hb):
                hb = np.zeros(nt) * hb
            else:
                assert np.all(hb.shape == t.shape), 'ht must have same shape as t'

        nLay   = self.aq.nLay

        assert np.all(hLR.shape == (nLay, nt)), 'hLR must have shape (nLay, nt)'
        assert np.all(N.shape == (nLay, nt)), 'N must have shape (nLay, nt)'

        Npg = self.get_Npg(N, ht, hb)

        hLR = np.matrix(hLR)
        Npg = np.matrix(Npg)

        D, D_m1, V, V_m1, S_m1 = self.aq.D, self.aq.D.I,\
                    self.aq.V, self.aq.V.I, self.aq.S_m1       # Eigen values and left eigen vectors

        # Coefficient for discharge to multilayer diches
        #TAm2 = -(self.aq.kD * self.aq.A_m2 * ((self.aq.A / self.aq.b)
        #        * lg.sinhm(self.aq.b * self.aq.A_m1) * self.aq.F_m1).I)
        shmf = lg.sinhm(self.aq.b * self.aq.A_m1) * self.aq.F_m1
        K = self.aq.kD * self.aq.A_m1 * shmf * \
            (self.aq.I - self.aq.A / self.aq.b * shmf).I / self.aq.b

        #Arrays to be filled in
        hmt = np.matrix(np.zeros((nLay, nt)))
        hmt[:, 0] = hLR[:, 0] # start at ditch level
        qditch = np.matrix(np.zeros_like(hmt))

        for i, dt in enumerate(np.diff(t)):

            hLR_ = hLR[:, i] # (nLay, 1) array/matrix
            hmt_ = hmt[:, i] # same
            e = lg.expm(-D * dt)

            '''
            hmt[:, i+1] = hLR_ + V * e * V_m1 * (hmt_ - hLR_) + \
            V * (self.aq.I - e) * D_m1 * V_m1 * S_m1 * Npg[:, i]
            '''
            hmt[:, i+1] = hLR_ + V * D_m1 * e * D * V_m1 * (hmt_ - hLR_) + \
            V * D_m1 * (self.aq.I - e) * V_m1 * S_m1 * Npg[:, i]

            hmt05 = 0.5*(hmt[:, i+1] + hmt[:, i])
            qditch[:,i+1] = K * (hmt05 - hLR_)

        return np.array(hmt), np.array(qditch) # nLay, nt



    def get_Npg(self, N=None, ht=None, hb=None):
        """
        Returns array Npg = recharge + leakage effect top and bottom aquitards.

        parameters
        ----------
        N: ndarray (Nlay, Nt)
            uniform recharge of all layers and all times
        ht: ndarray(nt)
            head boven the top aquifer. None of not present.
        hb: ndarray(nt)
            head below the bottom aquifer. None if absent.
        """

        assert N.shape[0] == self.aq.nLay

        if ht is not None:
            N[0, :]  += ht / self.aq.c[ 0]
        if hb is not None:
            N[-1, :] += hb / self.aq.c[-1]

        Npg = N.copy() # Don't overwrite original
        return Npg


    def plotx(self, x, h, q, **kwargs):
        """Plots the steady-state heads in the cross section

        parameters
        ----------
        x: [ m ] np.ndarray of len nx
            x coordinates along cross section between 0 and b, e.g. np.linspace(0, b, n)
        h: [ m ] np.ndarray of shape (nt, nLay + 2)
            Specified heads,  (only row 1 is used)
        q: [m/d] np.ndarray, of hape (nt, nLay)
            Specified injection (only row 1 is used)
        """
        plt.plot(x.T, self.phix(x, h, q, **kwargs))


    def plott(self, t, h, q, **kwargs):
        """Plots the transient head, aveated over the cross section.

        parameters
        ----------
        t: [ m ], np.ndarray, of shape (nt,)
            a vector of times to compute the average head in x-section.
        h: [ m ], np.ndarray of shape (nt, nLay + 2)
            specified heads,
        q: [m/d], np.ndarray of shape (nt, nLay)
            specified injection
        """
        plt.plot(t, self.phit(t, h, q, **kwargs))


    def plotm(self, h, q, **kwargs):
        """Plots steady-state head, averaged over the cross section

        parameters
        ----------
        h: [ m ], np.ndarray of shape (nt, nLay + 1)
            specified heads, (only row 1 is used)
        q: [m/d], np.ndarray of shape (t, nLay)
        specified injection (only row 1 is used)
        """
        b = self.aquifer.b
        for iL, fi in enumerate(self.phim(h, q)):
            plt.plot([0, b], [fi, fi], label="layer{}".format(iL))


if __name__ == '__main__':

    """Runs a multilayer example"""

    props = {'kD': [150., 150, 150], #Transmissivities of aquifers
             'c': [25000, 25, 25, 250000], # Resistance of aquitards
             'w': [0.01, 0.01, 0.01], # Entry resistances
             'D': [25., 25., 25.], # Thickness of quifers
             'S':[0.15, 0.1, 0.1], # Storage coefficients
             'd': [2., 2., 2., 2.],
             'b': 100, # half width of xsection
             } # Thickness of aquitards (not used)

    # instantiate Aquifer object
    aq = Multi_aquifer(kD=props['kD'], c=props['c'], w=props['w'],
                          D=props['D'], d=props['d'], S=props['S'],
                          b=props['b'])

    # file containing daily values of precipitation and evapotranspiration
    tne_file = 'PE-00-08.txt'

    # read into pandas
    tne = pd.read_table(tne_file, header=None, parse_dates=True, dayfirst=True,
                        delim_whitespace=True, index_col=0)

    tne.columns = ['N', 'E']
    tne.index.name = 'date'

    # add column with time in days, not dates
    tne['t'] = np.asarray(np.asarray(tne.index - tne.index[0],
                           dtype='timedelta64[D]') + 1, dtype=float)

    # compute q in m/d (injection in first layer only = recharge) (nt, nLay)
    b = 100.
    q = np.zeros((len(tne), aq.nLay)) # shape (nt, nLay)
    q[:, 0] = (tne.N - tne.E) / 1000.

    # given heads [nt, nLay + 2]
    h = np.zeros((len(tne), aq.nLay + 2)) # shape (nt, nLay + 2)

    # instantiate solution object
    multiaq = Analytic_nlay(aquifer=aq)

    #%% Steady, head in cross section
    hLR = np.array([0.0, 0.0, -0.0])
    N0  = np.array([0.002, 0.0, 0.0])
    ht  = 0.
    hb  = -1.

    # head and discharge as function of x in cross section
    h, q, x = multiaq.hqx(hLR=hLR, N=N0, ht=ht, hb=hb, n=51)

    # average head in cross section and discharge at x=b
    hm, qd = multiaq.hqm(hLR=hLR, N=N0, ht=ht, hb=hb)

    fig, ax = plt.subplots(2, 1, figsize=(10, 5))
    ax[0].set_title('Multilayer aquifer cross section')
    ax[1].set_xlabel('m from left to center')
    ax[0].set_ylabel('head')
    ax[1].set_ylabel('flow m2/d')
    ax[0].grid()
    ax[1].grid()
    for i, (h_, q_) in enumerate(zip(h, q)):
        ax[0].plot(x, h_, label='layer{}'.format(i))
        ax[1].plot(x, q_, label='layer{}'.format(i))

        ax[0].plot(x[[0,-1]], np.array([hm[i], hm[i]]), '--', label='layer{}'.format(i))
        ax[1].plot(x[-1], qd[i], '.', label='layer{}'.format(i))

    h_ = N0[0] / 2 * (np.max(x)**2 - x**2) / props['kD'][0]
    ax[0].plot(x, h_, '.', label='simple')

    ax[0].legend()
    ax[1].legend()


    #%% Trransient simulation average head in cross section

    nLay = aq.nLay
    t = (tne.index - tne.index[0]).astype('timedelta64[D]').values
    N = (tne['N'] - tne['E']).values

    # Make transient simulate a step response
    nt = 25
    N = N0[:, np.newaxis] * np.ones((1, nt))
    t = t[:nt]
    hLR = np.zeros_like(N)
    ht = np.zeros_like(N[0])
    hb = np.zeros_like(N[0])
    hb += -1

    ht, qd = multiaq.hmt(t, N=N, hLR=hLR, ht=ht, hb=hb)

    fig, ax3 = plt.subplots(2, 1, figsize=(10, 5))
    ax3[0].set_title('Multiaquifer transient, average head in cross section')
    ax3[0].set_ylabel('Head')
    ax3[1].set_ylabel('Q_ditch / b [m/d]')
    ax3[1].set_xlabel('tijd in d')
    ax3[0].set_ylim((0, 0.06))
    ax3[0].grid()
    ax3[1].grid()
    for i, (ht_, qd_) in enumerate(zip(ht, qd)):
        ax3[0].plot(t, ht_, label='layer{}'.format(i))
        ax3[1].plot(t, qd_, label='layer{}'.format(i))
    ax3[0].legend()
    ax3[1].legend()

    #%%

    #phi = solution.phix(x, h, q)
    #print("phi:\n", phi)
'''
    # plot head in cross section steady
    fig, ax= plt.subplots()
    ax.set_title('Head in cross section')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('head [m]')
    solution.plotx(x, np.mean(h, axis=0), np.mean(q, axis=0), label="test")
    solution.plotm(np.mean(h, axis=0), np.mean(q, axis=0), q)  # average head in cross secton (steady)
    plt.legend()
    plt.show()

    # plot transient heads (averaged over the cross section)
    fig, ax = plt.subplots()
    ax.set_title("Head as function of time")
    ax.set_xlabel('t [d]')
    ax.set_ylabel('phi [m]')
    ax.grid()
    solution.plott(tne.t, h, q)
    plt.show()
'''
