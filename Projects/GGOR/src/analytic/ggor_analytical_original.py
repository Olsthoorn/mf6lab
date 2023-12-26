#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 22 10:58:55 2018

Analyical transient simulations of cross sections between parallel ditches.

The problem considered is a cross section of a parcel between side
ditches with the same properties, subjected to varying recharge and water
levels in the ditches.

Computed are the cross-section averaged head, the discharge to the ditches.
Input are time-varying net recharge, ditch level and upward seepage.

Note that all input data are considered average values from the previous
time till the current one and the head is computed at the current time.
The times are the index of the pd.DataFrame with the precipitation and the
evapotranspiration data.

Alos, note that snow is not considered. Neither is interceptioin and storage
in the unsaturatied zone.
Interception and storage in the unsaturated zone may be implemented by a
prefilter on the meteo data. This is currently not implemented.

The meteo data are obtained from KNMI weather stations. The data can be
downloaded on the fly or be read from an already downloaded file. IN both
cases are the original precipitation and evapotranspiratond ata converted from
0.1 mm/d to m/d.

A number of analytic solutions is possible. A class is defined to encapsulate
each one of them. To run, the class instance delegates the actual computation
to a dedicated function.

Names of the different analytial solutions:
    L + [1|2] + [q]f] [+ W]
    This tells that the solution has 1 or 2 computed layers and the boundary
    condition in the underlying regional aquifer is eigher seepge (q) or head
    (f) and that it has or does not have entry/outflow resistance at the ditch.

    So we can name solutions as follows
    L1q, L2q, L1f L2f, L1qw, L1fw, L2f2, L2qw

@author: Theo
"""
import numpy as np
import pandas as pd
import scipy.linalg as la
import matplotlib.pyplot as plt
from KNMI import knmi
import GGOR.src.numeric.gg_mflow_1parcel as gn

def gen_testdata(data, **kwargs):
    """Return copy of time_data with altered input columns suitable for testing.

    The tuples consist of a set floats to be used repeatly for interval time:
        interval, val1, val2, val3, ...
        Changes occur after each interval

    Parameters
    ----------
    data: pd.DataFrame witih datetime index
        input for modelling
    RH: tuple of floats for net pecipitation variation
        interval, val1, val2, ...
    EV24: tuple of floats for net evapotranspiration variation
        interval, val1, val2, ...
    hLR: tuple of floats for ditch level variation
        interval, val1, val2, ...
    q: tuple of floats for vertical seepage variation
        interval, val1, val2, ...
    h1: tuple of floats for head in regional aquifer variation
        interval, vavl1, val2, val3, ...
        Note that the top layer is h0

    Example
    -------
    gen_testdata(data, RH=(200, 0.02, 0.01, 0.03, 0.01)),
                       EV24 = (365, 0., -0.001)
                       q=(150, -0.001, 0.001, 0, -0.003)
                       )
        This fills the data colums 'RH', EV24' and 'q' with successive
        values repeatly at each interval of resp. 200, 365 and 150 days.
    """
    data = data.copy() # leave time_data intact
    for W in kwargs:
        # index array telling which of the tuple values to pick
        I = np.asarray((data.index - data.index[0]) / np.timedelta64(1, 'D')
                          // kwargs[W][0] % (len(kwargs[W]) - 1) + 1, dtype=int)
        data[W] = np.array(kwargs[W])[I] # immediately pick the right values

    return data

def newfig(title='title?', xlabel='xlabel?', ylabel='ylabel?', xscale='linear', yscale='linear',
           xlim=None, ylim=None, size_inches=(12, 11), **kwargs):
    """Return ax for new plot."""
    fig, ax = plt.subplots()
    fig.set_size_inches(size_inches)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.grid()
    return ax

def newfig2(titles=['title1', 'title2'], xlabel='time',
            ylabels=['heads [m]', 'flows [m2/d]'],
            xscale='linear',
            yscale=['linear', 'linear'],
            sharex=True,
            sharey=False,
            xlims=None,
            ylims=None,
            size_inches=(12, 6),
            **kwargs):
    """Return ax[0], ax[1] for new plot."""
    fig, ax = plt.subplots(2, 1, sharex=sharex, sharey=sharey)
    fig.set_size_inches(size_inches)
    for a, title, ylabel in zip(ax, titles, ylabels):
        a.set_title(title)
        a.set_xlabel(xlabel)
        a.set_ylabel(ylabel)
        a.grid()
    ax[0].set_xscale(xscale)
    ax[1].set_xscale(xscale)
    ax[0].set_yscale(yscale[0])
    ax[1].set_yscale(yscale[1])
    if xlims is not None:
        ax[0].set_xlim(xlims[0])
        ax[1].set_xlim(xlims[1])
    if ylims is not None:
        ax[0].set_ylim(ylims[0])
        ax[1].set_ylim(ylims[1])
    return ax

#%% The properties dictionary

# Properties is a dict that holding the properties required by
# the analytical solution in name: value pairs.
# It is the user's responsibility to make sure that all required
# properties are contained in the dict.

props = {'b': 250, # distance between the parallel ditches
         'dx': 1.0, # only use for numerical model
         'hdr': 0.1, #m drain elevation
         'cdr': 1000,  #d drainage resistance (for drainage only)
         'c': (500, 50., 30.), # vertical resistance between layers (aquifers)
         'bd': 1.0, # width of each ditch
         'w': np.array([[1.0, 1.0], # always > 0, first col > second col
                        [3., 1.],
                        [np.inf, np.inf],
                        [np.inf, np.inf]]), # entree and outflow resistance [d]
         'z': np.array([[  0,  -8.],
                        [ -10, -30],
                        [-35,  -50],
                        [-52,  -60]]), # layer boundary elevations
         'S': np.array([[0.1, 0.001],
                        [0.1, 0.001],
                        [0.1, 0.001],
                        [0.1, 0.001]]), # specific yield and storage coeff. of layers
         'k': np.array([[ 2, 0.2],
                        [10,  5],
                        [50, 10],
                        [25,  5]]), # horizontal and vertical conductivity of layers
         'hLR': (0.2, -0.2, -0.2, -0.3), # steady state ditch level for each layer
         'q' : (0.002, 0., 0., 0.), # injection in each aquifer m/d steady state
         'zp': -0.2, # summer level in ditches, to used in dynamic simulations
         'wp': -0.2, # winter level in ditches, to used in dynamic simulations
           }

def set_hLR(data=None, props=None):
    """Add column hLR to data inline.

    Parameters
    ----------
    data: pd.DataFrame
        input data frame with datetime index.
    props: dict
        properties of aquifer, aquiclude, ditches.
    """
    # Add or update column "summer" (bool) in data
    data['summer'] = [m in {4, 5, 6, 7, 8, 9} for m in [d.month for d in data.index]]
    data['hLR'] = props['wp']
    data['hLR'].loc[data['summer']] = props['zp']


def check_cols(data=None, cols=None):
    """Check presence of input columns in pd.DataFrame.

    Parameters
    ----------
    data: pd.DataFame
        DataFrame with input colums to check
    cols: list of str
        names of columns that must be present
    """
    r = set(cols).difference(data.columns)
    if len(r) > 0:
        r = ["'" + k + "'" for k in r]
        raise ValueError("Missing the columns {{{}}} in DataFrame".format(', '.join(r)))

def getDT(data):
    """Return stepsize array from data.

    Dt is the diff of the index of data prepended wit one value
    equal to the value of the first Dt, so that Dt has the same
    lenght as the index.

    Parameters
    ----------
    data: pd.DataFrame
        the time dependent data
    """
    Dt = np.diff((data.index - data.index[0]) / np.timedelta64(1, 'D'))
    Dt = np.hstack((Dt[0], Dt))
    return Dt

# %% Stand-alone simulator

# This simulator is meant to work indpendently of the class solution
# It must therefore contain properties, solution name and the different
# methods that are designed for specific solutions.
def single_Layer_transient(solution_name, props=None, data=None):
    """Return results tim_data with simulation results included in its columns.

    Parameters
    ----------
    solution_name: str, solution name, one of
        'l1f' : Single layer with given head in regional aquifer
        'l1q' : Single layer with upward seepage no ditch resistance
        'l2q  : Two aquifer system with given seepage from lower aquifer'
        'l2qw': Same but with ditch resistance
    data: pd.DataFrame
        required fields: hLR, RH, EV24, q or h1
            hLR: left-right ditch-water level
            RH: precipitation
            EV24: evapotranspiration
            q: net seepage=, upward positive
            h1: head in the regional aquifer
        RH and EV24 must be in m/d!

    Returns
    -------
    data: pd.DataFrame
        this is the input DataFrame with extra columns, which contain
        the results of the simulation, namely:
            h0 for groundwater
            qd: the discharge to the two ditches
            sto: the storage during each timestep
    The seepage is not added; it is obtained from RH - EV24 - qd - sto.
    In case the solution includes a resistant top layer, the
    seepage through this toplayer is also included as qt.
    Notice that all the flows are in m/d, that is, averaged
    over the entire cross section and the repective timestep.

    The resistance between cover layer and regional aquifer is concentrated
    in the given cover-layer resistance 'c' at the bottom of the cover
    and is constant given in the properties.
    """
    data = data.copy() # leave original intact

    check_cols(data, ['RH', 'EV24', 'hLR'])

    # Extract data for single-layer problem
    b, c, w = props['b'], props['c'][0], props['w'][0]
    mu, S, k =  props['S'][0, 0], props['S'][1, 1], props['k'][0, 0]
    cdr, hdr = props['cdr'], props['hdr']

    D = np.abs(np.diff(props['z'], axis=1))[0][0]

    # Initialize head in top layer and regional aquifer
    h0 = np.zeros(len(data) + 1) # initialize h0, first aquifer
    h1 = np.zeros(len(data) + 1) # initialize h1, second aquifer (phi)
    hm = np.zeros(len(data))
    h0[0] = data['hLR'].iloc[0]
    h1[0] = data['hLR'].iloc[0]

    Dt = getDT(data)
    qs0 = np.zeros(len(data))
    qv0  = np.zeros(len(data))
    qdr  = np.zeros(len(data))
    qb0 = np.zeros(len(data))
    tol = 1e-3
    if solution_name == 'L1f': # for given Phi (f)
        check_cols(data, ['h1'])
        data['q'] = np.NaN
        for i, (dt, t1, hlr, RH, EV24, phi, hdr) in enumerate(zip(
                                Dt,                  # time step
                                np.cumsum(Dt),       # t1 (end of time step)
                                data['hLR'].values,  # ditch water level
                                data[  'RH'].values, # precipirtation
                                data['EV24'].values, # evapotranpiration
                                data[ 'h1'].values,  # phi, head in regional aqufier
                                data['hdr'].values)): # (tile) drainage elevation
            N = RH - EV24
            t = t1 - dt   # initialize t0, begining of time step
            hh = np.array([h0[i], 0])
            loop_counter = 0
            while t < t1:
                loop_counter += 1
                if loop_counter == 100:
                    print("Warning iteration {:i} no convergence!".format(i))
                    break

                w_    = w[0] if hh[0] <= hlr else w[1] # in or outflow resistance depedning on h - hLE
                if   hh[0] > hdr + tol:   cdr_ = np.inf
                elif hh[0] < hdr - tol:   cdr_ = cdr
                else:
                    lam    = np.sqrt(k * D * c)
                    Lamb   = 1 / ((b / lam) / np.tanh(b / lam) + (w_ / D) * (c / b))
                    rising = (phi - hh[0]) + (N * c - Lamb * (N * c - (hlr - phi))) > 0
                    if rising:      cdr_ = np.inf
                    else:           cdr_ = cdr
                ch     = c  / (1 + c / cdr_)    # c_hat (see theory)
                Th     = mu * ch
                phih   = (phi + (c / cdr_) * hdr) / (1 + c / cdr_)
                lamh   = np.sqrt(k * D * ch)
                Lambh  = 1 / ((b / lamh) / np.tanh(b / lamh) + (w_ / D) * (ch / b))
                B      = N * ch - Lambh * (N * ch - (hlr - phih))
                r      = (hh[0] - phih - B) / (hdr - phih - B)
                if r > 1:
                    dtau = Th * np.log(r) # estimate time step until crossing hdr
                    if t + dtau > dt:
                        dtau = t1 - t # limit to within current time step
                else:
                    dtau = t1 - t
                e = np.exp(-dtau / Th)
                f = (1 - e) / (dtau / Th)
                hh[1] = phih + (hh[0] - phih) * e + (
                    N * ch - Lambh * (N * ch - (hlr - phih))) * (1 - e)
                hm     = phih + (hh[0] - phih) * f + (
                    N * ch - Lambh * (N * ch - (hlr - phih))) * (1 - f)
                qs0[i] += mu * np.diff(hh) / dtau * (dtau / dt)
                qv0[i] += (phi - hm) / c   * f * dtau / dt
                qdr[i] += (hm  - hdr) / cdr * f * dtau / dt
                hh[0] = hh[1]
                t += dtau
            qb0[i] += N + qv0[i] - qdr[i] - qs0[i]
            h0[i+1] = hh[1]

        # Gather results and store in data columns
        data['h0']    = h0[1:]
        data['qs0']   = qs0
        data['qv0']   = qv0
        data['qdr']   = qdr
        qN  = (data['RH'] - data['EV24']).values
        data['qb0']   = qN + data['qv0'] - data['qdr'] - data['qs0']
        data['qsum0'] = qN + data['qv0'] - data['qdr'] - data['qb0'] -data['qs0']
        data['qv1']   = -data['qv0']
        data['qb1']   = 0.       # no ditches
        data['qs1']   = S * np.diff(np.hstack((data['h1'].values, data['h1'][-1]))) /Dt
        data['q1']    = data['qs1'] + data['qv0'] + data['qb1']    # undefined because phi is given
        data['qsum1'] = data['q1'] - data['qs1'] - data['qv0'] - data['qb1']
        print('sum(abs(qsum0)) = {}, sum(abs(qsum1)) = {} m/d'.
              format(data['qsum0'].abs().sum(), data['qsum1'].abs().sum()))

    elif solution_name == 'L1q': # for given q (f)
        check_cols(data, ['q'])
        for i, (dt, t1, hlr, RH, EV24, q, hdr) in enumerate(zip(
                                Dt,                  # time step
                                np.cumsum(Dt),       # t1 (end of time step)
                                data['hLR'].values,  # ditch water level
                                data[  'RH'].values, # precipirtation
                                data['EV24'].values, # evapotranpiration
                                data[ 'q'].values,   # injection in lower aquifer
                                data['hdr'].values)): # (tile) drainage elevation
            N = RH - EV24
            hh = np.array([[h0[i], 0],
                           [h0[i], 0]])
            if i == 400:
                print('i =' , i)
            for iter in range(2): # iterate over phi because q is given
                loop_counter = 0
                phi = hh[iter][0] + q * c
                t = t1 - dt
                while t < t1:
                    loop_counter += 1
                    if loop_counter == 100:
                        print("Warning iteration {:i} no convergence!".format(i))
                        break

                    w_    = w[0] if hh[iter][0] <= hlr else w[1] # in or outflow resistance depedning on h - hLE
                    if   hh[iter][0] > hdr + tol:   cdr_ = np.inf
                    elif hh[iter][0] < hdr - tol:   cdr_ = cdr
                    else:
                        lam    = np.sqrt(k * D * c)
                        Lamb   = 1 / ((b / lam) / np.tanh(b / lam) + (w_ / D) * (c / b))
                        rising = (phi - hh[iter][0]) + (N * c - Lamb * (N * c - (hlr - phi))) > 0
                        if rising:      cdr_ = np.inf
                        else:           cdr_ = cdr
                    ch     = c  / (1 + c / cdr_)    # c_hat (see theory)
                    Th     = mu * ch
                    phih   = (phi + (c / cdr_) * hdr) / (1 + c / cdr_)
                    lamh   = np.sqrt(k * D * ch)
                    Lambh  = 1 / ((b / lamh) / np.tanh(b / lamh) + (w_ / D) * (ch / b))
                    B      = N * ch - Lambh * (N * ch - (hlr - phih))
                    r      = (hh[iter][0] - phih - B) / (hdr - phih - B)
                    if r > 1:
                        dtau = Th * np.log(r) # estimate time step until crossing hdr
                        if t + dtau > dt:
                            dtau = t1 - t # limit to within current time step
                    else:
                        dtau = t1 - t
                    e = np.exp(-dtau / Th)
                    f = (1 - e) / (dtau / Th)
                    hh[iter][1] = phih + (hh[iter][0] - phih) * e + (
                        N * ch - Lambh * (N * ch - (hlr - phih))) * (1 - e)
                    hm     = phih + (hh[iter][0] - phih) * f + (
                        N * ch - Lambh * (N * ch - (hlr - phih))) * (1 - f)
                    qs0[i] += mu * np.diff(hh[iter]) / dtau * (dtau / dt)
                    qv0[i] += (phi - hm) / c   * f * dtau / dt
                    qdr[i] += (hm  - hdr) / cdr * f * dtau / dt
                    hh[iter][0] = hh[iter][1]
                    t += dtau
                qb0[i] += N + qv0[i] - qdr[i] - qs0[i]
                if iter == 0:
                    hh[iter + 1][0] = hh[iter][1]
            # results below divided by 2 due to iteration iter
            h0[i+1] = hh[:, 1].mean() # mean of the two iterations at t1
            h1[i+1] = (h0[i] + h0[i+1]) / 2 + q * c
            qs0[i] /= 2 # due to iter, 2 loops
            qv0[i] /= 2
            qdr[i] /= 2
            qb0[i] /= 2

        # Gather results and store in data columns
        data['h0']    = h0[1:]
        data['h1']    = h1[1:]
        data['hdr']   = hdr
        data['qs0']   = qs0
        data['qv0']   = qv0
        data['qdr']   = qdr
        qN  = (data['RH'] - data['EV24']).values
        data['qb0']   = qN + data['qv0'] - data['qdr'] - data['qs0']
        data['qsum0'] = qN + data['qv0'] - data['qdr'] - data['qb0'] -data['qs0']
        data['qs1']   = S * np.diff(h1) / Dt
        data['qv1']   = -data['qv0']
        data['qb1']   = 0.       # no ditches
        data['q1']    = data['q']    # undefined because phi is given
        data['qsum1'] = data['q1'] + data['qv1'] - data['qs1'] - data['qb1']
        print('sum(abs(qsum0)) = {} m/d, sum(abs(qsum1))'.
              format(data['qsum0'].abs().sum(), data['qsum1'].abs().sum()))
    elif solution_name == 'xxL1q': # For given q
        # Direct computation with given q
        lamb = np.sqrt(k * D * c)
        for i, (dt, hlr, RH, EV24, q) in enumerate(zip(Dt,
                                data['hLR'].values,
                                data['RH'].values,
                                data['EV24'].values,
                                data['q'].values)):
            n = RH - EV24
            w_ = w[0] if h0[i] < hlr else w[1]
            Lamb[i] = 1 / ((b / lamb) / np.tanh(b / lamb) + (w_ / D) * (c  / b))
            T = mu * c / Lamb[i]
            e = np.exp(-dt / T)
            h0[i + 1] =  hlr + (h0[i] - hlr) * e + (c * (n + q) * (1 - Lamb[i]) / Lamb[i]) * (1 - e)
    else:
        raise ValueError("Don't recognize solution name <{}>".
                         format(solution_name))

    return data

def single_layer_steady(solution_name, props=None, data=None, dx=1.0):
    """Return simulated results for cross section.

    Parameters
    ----------
    solutions_name: str, the name of the case/solution
        'L1q' : Single layer given seepage.
        'L1h' : Single layer with given head in regional aquifer.
    props: dict
        required properties of the system
    data: dict or pd.DataFrame
        Data for w hich to simulate the steady situation.
        If dict, it must contain 'RH', 'EV24', 'q' or ''h1' depending on the case.
        if pd.DataFrame, it must contain columns 'RH', 'EV24', 'q' or 'h1', depending on
        the specific case.
        If dict, the average situation is computed as a steady-state cross section.
    dx: float
        Desired x-axis step size in cross section. The width is known from
        the properties.
    """
    k, c, b, w = props['k'][0, 0], props['c'][0], props['b'], props['w'][0].mean()
    D = np.diff(props['x'], axis=1)[0]
    lam = np.sqrt(k * D * c)

    # X coordinages
    x = np.hstack((np.arange(0, b, dx), b)).unique() # make sure b is included
    x = np.hstack((x[::-1]), x[1:]) # generate for total section -b < x < b

    if isinstance(data, pd.DataFrame):
        check_cols(data, ['RH', 'EV24', 'hLR', 'h1'])
        if solution_name == 'L1q':
            check_cols(data, ['h0', 'q'])

        data = dict(data.mean()) # turn DataFrame into a Series
        N = data['RH'] - data['EV24']

    Lamb = 1 / (b  / lam / np.tanh(b  / lam) + (b /c) / (D / w))
    phi = data['h1'] if not (solution_name in ['L1q']) else data['h0'] - data['q'] * c

    hx = phi + N * c - (N * c - (data['hLR'] - phi)
                    ) * Lamb * (b / lam) * np.cosh(x / lam) / np.sinh(b / lam)

    return hx

def multi_layer_steady(b=None, dx=None, kh=None, z=None, c=None, w=None,
                       hdr=None, cdr=None, q=None, hLR=None, plot=True, **kwargs):
    """Return steady multilayer solution for symmetrial cross section.

    With all data in a dict kwargs, call this function using
    multi_layer_steady(**kwargs).
    No data pd.DatFrame input is needed, because not transient.
    If plot==True, then a graph is made showing the heads in
    the corss section for all layers and printing the water
    balance.

    Parameters
    ----------
    b: float, half-width of cross section [L].
    dx: step size for plot [L].
    kh: vector of conductivities for respective aquifers [L/T].
    z: vector of successive layer elevations (first at top) [L].
    c: Vector of aquitard resistances below the layers. Default np.inf [T].
    w: Vector of intry resistances [T].
    hdr: Elevation of drainage [L].
    cdr: Drainage resistance [T].
    q: Recharge as a vector with one value per aquifer [L/T]. q[0] must be
        RH - EV24, the net recharge to the top layer.
    hLR: Prescribed heads at x=b, vector with one value for each aquifer [L]

    @TO 20200615
    """
    def sysmat(c, kD):
        """Return system matrix."""
        dL  = 1 / (c[:-1] * kD)
        dR  = 1 / (c[1: ] * kD)
        Am2 = -np.diag(dL[1:], k=-1) - np.diag(dR[:-1], k=1) + np.diag(dL + dR, k=0)
        Am1 = la.sqrtm(Am2)
        A1  = la.inv(Am1)
        A2  = A1 @ A1
        return Am2, Am1, A1, A2

    vect = lambda x: x[:, np.newaxis]
    x = np.hstack((np.arange(0, b, dx), b))
    D = np.abs(np.diff(z))
    kD = kh * D
    c = np.hstack((cdr, c, np.ones(len(kD) -len(c)) * np.inf))
    hLR = vect(np.array(hLR))
    q = vect(np.array(q))
    g = np.zeros_like(q)
    g[0]  = hdr / (cdr * kD[0])

    Am2, Am1, A1, A2 = sysmat(c, kD)

    T   = np.diag(kD)
    Tm1 = la.inv(T)
    w = np.array(w)
    kh = np.array(kh)
    Kw  = np.diag(kh * w)

    coshmb = la.coshm(Am1 * b)
    sinhmb = la.sinhm(Am1 * b)
    F = la.inv(coshmb + Kw @ Am1 @ sinhmb)
    I = np.eye(len(kD))

    # Compute Phi for all x-values
    Phi = np.zeros((len(kD), len(x)))
    hq = A2 @ Tm1 @ (q + g)

    for i, xx in enumerate(x):
        coshmx = la.coshm(Am1 * xx)
        Phi[:, i] = (coshmx @F @ hLR + (I - coshmx @ F) @ hq).ravel()

    Phim = hq - A1 / b @ sinhmb @ F @ (hq - hLR) # average head in the layers
    qb   = T @ Am1 / b @ sinhmb @ F @ (hq - hLR) # discharge to ditches
    ql   = T @ Am2 @ Phim - g

    if plot:
        ax = newfig(f'Phi in all layers hdr={hdr:6.2f}, cdr={cdr:6.2f}', 'x [m]', 'Phi [m]',  size_inches=(14,8))
        clrs = []
        for iL, p in enumerate(Phi):
            h=ax.plot(x, p, label=
            f"layer={iL}, kD={kD[iL]:6.1f}, c_={c[iL + 1]:6.1f}, hLR={hLR[iL, 0]:6.2f}, " +
            f"w={w[iL]:.4g}, qb={qb[iL,0]:.4g} ql={ql[iL,0]:.4g}, q={q[iL,0]:.4g}, " +
            f"g={g[iL,0]:4g}, q - qb - ql = {q[iL,0] - qb[iL,0] - ql[iL,0]:.4g} m/d")
            clrs.append(h[-1].get_color())
        ax.hlines(Phim, x[0], x[-1], colors=clrs, ls='dashed')

        # Water budget
        WB = np.hstack((q, -qb, -(T @ Am2 @ Phim - g), q -qb -(T @ Am2 @ Phim - g)))
        print("\n\nWater balance for each layer (to layer is positive, from layer negative):")
        print(' '.join([f"{x:10s}" for x in ['supplied', '-ditch', '-leakage', 'sum [m/d]' ]]))
        print(' '.join([f"{x:10s}" for x in ['q', '-qb', '-ql', 'sum' ]]))
        for wb in WB:
            print(' '.join([f"{x:10.7f}" for x in wb]))
        print("\n")
    return Phi

def multi_layer_transient(data=None, b=None, k=None, z=None, c=None, w=None,
            cdr=None, hdr=None, S=None, check=True, **kwargs):
    """Return steady 2layer dynamic solution for symmetrial cross section.

    The code can deal with multiple layers. But for sake of limiting the input
    the imlementation is limited to two aquifer and 3 aquitards, the
    first of which is the drainage resistance and the last one is infinite by
    default.

    It checkes wheather water infitlrates or exfiltrates from ditch and adjusts
    the ditch resistance accordingly.

    It does not do dynamic drainage by switching the drainage resistance
    in accordance with the head in the top layer. This would lead to insta0
    bilities without explicitly precomputing the moments on which the head
    in the top layer crosses hdr. It may be added in the future.

    The resistance of the top aquitard (cdr) or of the aquitard below the bottom
    layer must not be infinite, or the system matrix will be singular,
    and no results can be computed.

    It can do ditches in all layers. Set w = np.inf for layers that do not have
    penetrating ditches.

    If plot == True, then a plot of the heads is shown and the water budget
    is printed for the top layer and the second layer. This printing is then
    done by selecting the proper columns.

    Parameters
    ----------
    data: a pd.DataFrame holding the dynamic data
        RH  : Preciptaion injected in the first aquifer
        EV24: Evapotranspiration from the first aquifer
        qv : seepage injected in the second aquifer
        hLR: Prescribed heads at x=b, vector with one value for each aquifer [L]
            these are the ditch levels in both the first and second aquifer
            as these are assumed the same [L].
    b: Float, half-width of cross section [L].
    k: Array [nLay, 2] of conductivities. First col kh, second col kv [L/T].
    z: Array [nLay, 2] first col top of layer, second bottom of layers [L].
    c: Vector of aquitard resistances below the layers [nLay].
       Set the last equal to inf. Or it's done by default. No leakage from below.
       The top c is that on top of the top aquifer, i.e. normally cdr.
    w: Ditch resistances [nlay, 2] (first col inflow, second col outflow) [T].
    S: storage coeff (nLay, 2) first col Sy,  second col S [-].
    cdr: Float. Resitance of top layer
    hdr: Float. Fixed head above top resistance layer.

    Returns
    -------
    data: pd.DataFrame
        Results, data copied from input and augmented with simulation results.

    @TO 20200701
    """
    def check_shape(nm, var, nLay=2, ncol=2):
        """Verify parameter shape.

        Parameters
        ----------
        nm: str
            parameter name
        p: np.array
            parameter
        nLay: int
            number of rows required
        """
        var = np.asarray(var)
        assert var.shape[0] >= nLay, AssertionError(f"{nm} must have at least {nLay} rows.")
        assert var.shape[1] >= ncol, AssertionError(f"{nm} must have at least {ncol} columns.")
        return var[:nLay, :ncol] if ncol > 1 else var[:nLay, ncol]


    def sysmat(c, kD):    # System matrix
        """Return the system matrix and related matrices.

        Parameters
        ----------
        c: vector of 2 floats
            resistance of the layer on top of each aquifer (top = drainage)
        kD: vector of 2 floats
            transmissivity of the layers
        """
        if len(c) != len(kD) + 1:
            raise ValueError(f"len c={len(c)} != len(kD) + 1={len(kD)+1}!")
        dL  = 1 / (c[:-1] * kD) # left  diagonal
        dR  = 1 / (c[1: ] * kD) # right diagional
        Am2 = -np.diag(dL[1:], k=-1) - np.diag(dR[:-1], k=1) + np.diag(dL + dR, k=0)
        Am1 = la.sqrtm(Am2)
        A1  = la.inv(Am1)
        A2  = A1 @ A1
        return Am2, Am1, A1, A2

    nLay = 2
    data = data.copy() # keep original data intact
    idx = data.index

    c  = np.hstack((cdr, np.array(c)[:nLay -1], np.inf))
    w  = check_shape('w', w)
    # ditch resistance, first col inffiltration, second exfiltration
    assert np.all(w[:, 0] >= w[:, 1]), AssertionError("All w[:,0] must be >= w[:, 1]")

    S  = check_shape('S', S); S=np.hstack((S[0, 0], S[1:, 1]))
    kh = check_shape('k', k, ncol=1)
    z  = check_shape('z', z)
    D  = np.abs(np.diff(z, axis=1))[:, 0]
    kD = kh * D

    Am2, Am1, A1, A2 = sysmat(c, kD)

    T   = np.diag(kD) # transmissivity
    Tm1 = la.inv(T)

    # unifrom injection into each aquifer
    q   = np.zeros((nLay, len(data)))
    q[0, :] = (data['RH'] - data['EV24']).values[:]
    q[1, :] =  data['q'].values[:]

    # Leakage from top aquifer i.e. hdr / (cdr * kD0)
    g = np.zeros((nLay, len(data)))
    g[0, :]  = hdr / cdr # this may be made time variable
    #g[-1, :] = hLast / c[-1]  # could be used, nut we take c[-1] = inf

    coshmb = la.coshm(Am1 * b)
    sinhmb = la.sinhm(Am1 * b)
    I = np.eye(len(kD))

    # Compute Phi for all x-values
    G    = A2 @ Tm1 @ np.diag(S)
    E, V = la.eig(G)
    Em1  = la.inv(np.diag(E))
    Vm1 = la.inv(V)

    # Initialize time steps and heads
    Dt  = np.hstack(((idx[1] - idx[0]) / np.timedelta64(1, 'D'), np.diff((idx - idx[0]) / np.timedelta64(1, 'D'))))
    Phi = np.zeros((2, len(data) + 1))
    Phi[:, 0] = data['hLR'].iloc[0]
    qb = np.zeros_like(q) # from layers to boundary at x=b

    # Loop over time steps
    for it, (dt, hlr) in enumerate(zip(Dt, data['hLR'])):
        hLR = hlr * np.ones((nLay, 1))

        # infiltration or exfiltration ?
        w_ = (w[:, 0] * (Phi[:, it] < hLR[:, 0]) + w[:,1] * (Phi[:, it] >= hLR[:, 0]))[:, np.newaxis]
        Kw = np.diag(kh * w_)
        F = la.inv(coshmb + Kw @ Am1 @ sinhmb)

        e = la.expm(-Em1 * dt)
        hq = A2 @ Tm1 @ (q[:, it:it+1] + g[:, it:it+1])
        qb[:, it:it+1] = T @ Am1 / b @ sinhmb @ F @ (hq - hLR)

        # steady state ultimate solution for t->inf
        hss = A2 @ Tm1 @ (q[:,it:it+1] + g[:, it:it+1] - qb[:,it:it+1])

        # Compute head
        Phi[:, it + 1 : it + 2] = V @ e @ Vm1 @  Phi[:, it : it + 1] + V @ (I - e) @ Vm1 @ hss

    Phim = (Phi[:,:-1] + Phi[:, 1:]) / 2
    qs = S[:, np.newaxis] * np.diff(Phi, axis=1) / Dt  # from top layer into storage

    # leakage through aquitards fram all layers
    ql = np.zeros_like(q)
    for it in range(len(Dt)):
        ql[:,it:it+1] = T @ Am2 @ Phim[:, it:it+1] - g[:, it:it+1]

    Phi = Phi[:, 1:] # cut off first day (before first data in index)

    # Store results
    data['h0'] = Phi[0]  # head top layer
    data['h1'] = Phi[1]  # head bot layer
    data['hdr'] = hdr
    data['qs0'] = qs[0]  # into storage top layer
    data['qs1'] = qs[1]  # into stroage in bot layer
    data['q0']  = q[0]   # injection into top layer (RH - EV24)
    data['q1']  = q[1]   # injectino into bot layer
    data['ql0'] = ql[0]  # leakage from top layer
    data['ql1'] = ql[1]  # leakage from bog layer
    #data['qb00'] = qb[0]  # to ditch from top layer
    data['qb0'] = data['q0'] - data['qs0'] - data['ql0']
    #data['qb10'] = qb[1]  # to dicht from second layer
    data['qb1'] = data['q1'] - data['qs1'] - data['ql1']
    data['qdr'] = (Phim[0] - hdr) / cdr      # to drain from top layer
    data['qv0'] = (Phim[1] - Phim[0]) / c[1] # to top layer from bot layer
    data['qv1'] = -data['qv0']               # to bot layer from top layer
    data['sumq0'] = data['q0'] - data['qs0'] - data['qb0'] - data['ql0'] # water balance top layer
    data['sumq1'] = data['q1'] - data['qs1'] - data['qb1'] - data['ql1'] # water balance bot layer
    data['sumq01'] = data['q0'] - data['qs0'] - data['qb0'] - data['qdr'] + data['qv0']
    data['sumq11'] = data['q1'] - data['qs1'] - data['qb1'] - data['qv0']

    # note that it must also be valid that q0 = qv0 - qde

    if check:
        """Check the water budget."""

        print('\n\nFirst water balance:')
        # show the water balance
        ttl_ = ['supply', 'storage', 'toDitch', 'leakage', 'sumq' ]
        ttl0 = ['q0', 'qs0', 'qb0', 'ql0', 'sumq0']
        ttl1 = ['q1', 'qs1', 'qb1', 'ql1', 'sumq1']

        mxcols = pd.options.display.max_columns
        pd.options.display.max_columns = len(ttl_) + 1

        print("\nWater balance first layer")
        print('                  ', ' '.join([f"{k:10s}" for k in ttl_]))
        print(data[ttl0])

        print("\nWater balance second layer")
        print('                  ', ' '.join([f"{k:10s}" for k in ttl_]))
        print(data[ttl1])
        print()
        pd.options.display.max_columns = mxcols

        print('\n\nSecond water balance:')
        # show the water balance
        ttl0_ = ['supply', 'storage', 'toDitch', 'drn', 'leak', 'sumq' ]
        ttl0  = ['q0',     'qs0',     'qb0',     'qdr', 'qv0', 'sumq01']
        ttl1_ = ['supply', 'storage', 'toDitch', 'drn', 'leak', 'sumq' ]
        ttl1  = ['q1',     'qs1',     'qb1',     'qv1',        'sumq11']

        print("\nWater balance first layer")
        print('                  ', ' '.join([f"{k:10s}" for k in ttl0_]))
        print(data[ttl0])

        print("\nWater balance second layer")
        print('                  ', ' '.join([f"{k:10s}" for k in ttl1_]))
        print(data[ttl1])
        print()
        pd.options.display.max_columns = mxcols


    return data # return updated DataFrame


def getGXG(data=None, startyr=None, nyr=8):
    """Return GXG coputed over 8 hydrological years starting at startyr.

    The GXG is computed from the 14th and 28th of the month groundwater
    head values in the period startyr/04/01 through endyr/03/31, so called
    hydrological years, in total nyr. endYr - startYr = nyr -1

    We add boolean columns 'GHG', 'GLG', 'GVG' to the DataFrame data indicating
    which dates within the nyr hydrological years contribute to the GXG.
    To get the GLG, just data['h0'].loc[data['GLG']].

    Parameters
    ----------
    data: pd.DataFrame with datetime index and column 'h0'
        input data to compute the GXG
    startyr: int
        year of time series
    nyr: int
        number of years, such that startyr + nYr = endYr

    Returns
    -------
    GXg object
    """
    AND = np.logical_and
    OR  = np.logical_or

    data['GLG'] = False # boolean column showing which dates contribute to the GLG
    data['GHG'] = False # same for GHG
    data['GVG'] = False # same for GVG
    h0 = data['h0'] # simplifies some lines and the return line below
    dt = data.index[0] - np.datetime64(data.index[0].date()) # if index is offset by some time within the day
    for iyr in range(nyr): # run over the hydrological years and set the boolean columns of data
        y1 = startyr + iyr
        t1 = np.datetime64(f'{y1    }-04-01') + dt # start hydrological year
        t2 = np.datetime64(f'{y1 + 1}-03-28') + dt # end hydrologial year

        # dGXG = boolean indicating measurement dates 14th and 28th of each month
        dGXG = AND(AND(data.index >= t1, data.index <= t2), data.index.day % 14 == 0)
        # dGVG boolean, indicating dates that contribute to spring level
        dGVG = AND(dGXG, OR(
                        AND(data.index.month == 3, data.index.day % 14 == 0),
                        AND(data.index.month == 4, data.index.day == 14)
                        ))

        data.loc[h0[dGXG].nlargest( 3).index, 'GHG'] = True # set GHG
        data.loc[h0[dGXG].nsmallest(3).index, 'GLG'] = True # set GLG
        data.loc[dGVG, 'GVG'] = True # set GVG

    # Return the actual GLG, GHG and GVG as a tuple, actual points are retrieved by the boolean columns
    return h0.loc[data['GLG']].mean(), h0.loc[data['GHG']].mean(), h0.loc[data['GVG']].mean()


# def plot_2layer_heads_and_flows(data=None, props=None):
#     """Plot heads and flows.

#     Parameters
#     ----------
#     data: pd.DataFrame
#         time varying data and results
#     props: dict
#         aquifer and aquitard properties
#     """
#     ax = newfig2([f"Phi in all layers hdr={props['hdr']:6.2f}, cdr={props['cdr']:6.2f}", 'Flows [L/T]'], 'time',
#                   ['Phi [m]', 'q[m/d]'],  size_inches=(14,12))
#     ax[0].plot(data.index, data['h0'], label='top layer')
#     ax[0].plot(data.index, data['h1'], label='bottom layer')
#     ax[0].legend(loc='best')

#     ax[1].plot(data.index, data['RH'], label='P')
#     ax[1].plot(data.index, -data['EV24'], label='E')
#     ax[1].plot(data.index, data['qdr'], label='qdrn')
#     ax[1].plot(data.index, data['qb0'], label='qditch')
#     ax[1].plot(data.index, data['qv0'], label='qv0')
#     ax[1].plot(data.index, data['qs0'], label='qs0')
#     ax[1].legend(loc='best')
#     return ax

def plot_watbal(ax=None, data=None, titles=None, xlabel=None, ylabels=['m/d', 'm/d'],
                size_inches=(14, 8), sharex=True, sharey=True,
                single_layer=False, **kwargs):
    """Plot the running water balance.

    Parameters
    ----------
    ax: list
        ax[0], ax[1] are the axes for plotting the top and bottom layer.
    titles: list of 2 strings
    xlabel: str
    ylabels: list of 2 strings
    size_inches: tuple of two
        Figure size. Only applied when a new figure is generated (ax is None).
    kwargs: dict with extra parameters passed to newfig2 if present

    Returns
    -------
    ax: list of two axes
        plotting axes
    """
    LBL = { 'q1' : {'leg': 'q-in' ,   'clr': 'green',    'sign': +1},
            'RH' : {'leg': 'RCH',   'clr': 'green',    'sign': +1},
            'EV24':{'leg': 'EVT',   'clr': 'gold',     'sign': -1},
            'DRN': {'leg': 'DRN',   'clr': 'lavender', 'sign': +1},
            'RIV': {'leg': 'DITCH', 'clr': 'magenta',  'sign': +1},
            'qb0': {'leg': 'DITCH', 'clr': 'indigo',   'sign': -1},
            'qb1': {'leg': 'DITCH', 'clr': 'indigo',   'sign': -1},
            'qv0': {'leg': 'LEAK',  'clr': 'gray',     'sign': +1},
            'qv1': {'leg': 'LEAK',  'clr': 'gray',     'sign': +1},
            'qs0': {'leg': 'STO',   'clr': 'cyan',     'sign': -1},
            'qs1': {'leg': 'STO',   'clr': 'cyan',     'sign': -1},
            'qdr': {'leg': 'DRN',   'clr': 'blue',     'sign': -1},
            }

    if ax is None:
        ax = newfig2(titles, xlabel, ylabels, sharey=True, size_inches=size_inches, **kwargs)
        ax[0].set_title('titles[0]')
        ax[1].set_title('titles[1]')
    elif isinstance(ax, plt.Axes):
        ax = [ax]
    for a, title, ylabel in zip(ax, titles, ylabels):
        a.set_title(title)
        a.set_xlabel(xlabel)
        a.set_ylabel(ylabel)
        a.grid(True)

    check_cols(data, ['RH', 'EV24', 'qs0', 'qb0', 'qv0', 'qdr'])
    W0 = data[['RH', 'EV24', 'qs0', 'qb0', 'qv0', 'qdr']].copy()
    for k in W0.columns:
        W0[k] *= LBL[k]['sign']
    C0 = [LBL[k]['clr'] for k in W0.columns]
    Lbl0 = [LBL[k]['leg'] for k in W0.columns]

    index = W0.index
    W0 = np.asarray(W0)
    ax[0].stackplot(index, (W0 * (W0>0)).T, colors=C0, labels=Lbl0)
    ax[0].stackplot(index, (W0 * (W0<0)).T, colors=C0) # no labels

    ax[0].legend(loc='best', fontsize='xx-small')

    # Check water balance layer 1
    #np.sum(W0 * (W0 > 0) + W0 * (W0 < 0), axis=1)

    if not single_layer:

        check_cols(data, ['q1', 'qs1', 'qb1', 'qv1'])
        W1 = data[['q1', 'qs1', 'qb1', 'qv1']].copy()
        for k in W1.columns:
            W1[k] *= LBL[k]['sign']

        C1 = [LBL[k]['clr'] for k in W1.columns]
        Lbl1 = [LBL[k]['leg'] for k in W1.columns]

        index = W1.index
        W1 = np.asarray(W1)

        ax[1].stackplot(index, (W1 * (W1>0)).T, colors=C1, labels=Lbl1)
        ax[1].stackplot(index, (W1 * (W1<0)).T, colors=C1) # no labels

        ax[1].legend(loc='best', fontsize='xx-small')

        ax[0].set_xlim(index[[0, -1]])

        # Check water balance layer 1
        #W1 * (W1 > 0) - W1 * (W1 < 0)

    # make sure y-axes are shared
    ax[0].get_shared_y_axes().join(ax[0], ax[1])

    return ax

def plot_heads(ax=None, data=None, title=None, xlabel='time', ylabel=['m'],
           size_inches=(14, 8), loc='best', **kwargs):
    """Plot the running heads in both layers.

    Parameters
    ----------
    ax: plt.Axies
        Axes to plot on.
    title: str
        The title of the chart.
    xlabel: str
        The xlabel
    ylabel: str
        The ylabel of th chart.
    size_inches: tuple of two
        Width and height om image in inches if image is generated and ax is None.
    kwargs: Dict
        Extra parameters passed to newfig or newfig2 if present.

    Returns
    -------
    The one or two plt.Axes`ax
    """
    if ax is None:
        ax = [newfig(title, xlabel, ylabel, size_inches=size_inches, **kwargs)]
    else:
        ax.grid(True)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

    check_cols(data, ['h0', 'hLR', 'hdr', 'h1'])
    ax.plot(data.index, data[ 'h0'], 'k', lw=1, label='head top layer')
    ax.plot(data.index, data['hLR'], 'r', lw=1, label='hLR (ditches)')
    ax.plot(data.index, data['hdr'], 'g', lw=1, label='zDrainage')
    ax.plot(data.index, data[ 'h1'], 'b', lw=1, label='h1 (phi)')
    ax.legend(loc=loc)
    return ax


def plotGXG(ax=None, data=None, startyr=None, nyr=8):
    """Plot GXG on an existing axes.

    Parameters
    ----------
    ax: plt.Axes
        An existing axes to plot on.
    data: pd.Series
        The head data.
    startyr: int
        The first hydrological year to use.
    nyr: int
        The number of hydrological years to includein GXG computation.
    """

    if not startyr:
        startyr = data.index[0].year + 1

    glg, ghg, gvg = getGXG(data=data, startyr=startyr, nyr=nyr)

    p = data

    # plot the values of h0 pertaining to glg, ghg or gvg
    ax.plot(p.index[p['GLG']], p['h0'].loc[p['GLG']], 'ro', label='GLG data', mfc='none')
    ax.plot(p.index[p['GHG']], p['h0'].loc[p['GHG']], 'go', label='GHG data', mfc='none')
    ax.plot(p.index[p['GVG']], p['h0'].loc[p['GVG']], 'bo', label='GVG data', mfc='none')

    # plot the glg, ghg and gvg values as dotted line over the evaluation period
    t0 = np.datetime64(f'{startyr}-04-01')
    t1 = np.datetime64(f'{startyr + nyr}-03-31')
    ax.plot([t0, t1], [glg, glg], 'r--', label='GLG')
    ax.plot([t0, t1], [ghg, ghg], 'g--', label='GHG')
    ax.plot([t0, t1], [gvg, gvg], 'b--', label='GVG')

    ax.legend(loc='lower left')
    return


def plot_hydrological_year_boundaries(ax=None, startyr=None, nyr=8):
    """Plot hydrological year boundaries on a given axis.

    Parameters
    ----------
    ax: plt.Axes
        an existing axes with a datatime x-axis
    staryr: int
        year of first hydrological year
    nyr: int
        numbef of years to include
    """
    if isinstance(ax, plt.Axes): ax = [ax]

    for a in ax:
        for iy in np.arange(nyr + 1):
            t = np.datetime64(f'{startyr + iy}-04-01')
            a.axvline(t, color='gray', ls=':')
            a.axvline(t, color='gray', ls=':')

#%% Solution class

class Solution:
    """Analytic solution base class.

    Specific analytic solution classes should be derived from this (see below).

    @TO 2020-07-01
    """

    def __init__(self, props=None):
        """Return an instance of an analytical solution only storing name and properties.

        Parameters
        ----------
        props: dict
            a dict containing the properrties. The necessary properties
            are given in the example tamplate in this class. Not all
            properties are used by all solutions. Unsued properties
            may be omitted from the actual properties.
        """
        self.name = str(self.__class__).split('.')[-1].split("'")[0]
        self.props = props
        return


    def check_props(self, props=None):
        """Return verified properties.

        Parameters
        ----------
        props: dict
            properties of the section to be simulated
        """
        missing_params = set(self.required_params) - set(props.keys())

        if missing_params:
            raise ValueError('Missing required properties: ({})'.
                             format(missing_params))
        return props


    def check_cols(self, data=None):
        """Verify input data for presence of required columns.

        Parameters
        ----------
        data: pd.DataFrame
            Input data with required columns which are
            'RH', 'EV24', 'hLR', 'q' or 'h1' dep. on self.name.

        Returns
        -------
        None

        """
        missing_cols = set(['RH', 'EV24', 'hLR']).difference(data.columns)

        if missing_cols:
            raise KeyError("{" + " ,".join([f"'{k}'" for k in missing_cols]) + "} are missing.")


    def sim(self, data=None):
        """Compute and store head and flows in added columns of the input data.

        Parameters
        ----------
        data: pd.DataFrame with all required time series in columns
        required columns: 'hLR', 'RH','EV24','q'|'h1'
            meaning:
            hLR: [m] ditch water level,
            RH: [m/d] precipitation
            EV24: [m/d] evap
            q: [m/d] upward seepage,
            h1: [m]head in regional aquifer.
            h0: [m] head above shallow aquifer

        Returns
        -------
        data with extra or overwritten columns:
            'h0','qd', 'qs', 'q2'[, q0]
        h0: [m] simulated head in cover layer
        qd: [m/d] discharge via ditches
        qsto: [m/d] stored
        if applicable for the particular solution:
            q2: [m/d] computed seepage from regional aquifer
            q0: [m/d] computed seepage from overlying layer with constant head

        All point variable (heads) are valid at the timestamps; the flow
        values are average values during each time step.

        The first head is hLR
        The length of the first time step is assumed equal to that
        of the second time step.

        The head at the beginning of the first time step is assumed
        equal to that of ditches during the first time step.
        """
        self.data = single_Layer_transient(solution_name=self.name,
                                          props=self.props,
                                          data=data)
        return


    def plot(self, titles=['heads', 'flows layer 0', 'flows layer 1'],
             xlabel='time', ylabels=['m', 'm/d', 'm/d'],
             size_inches=(14, 8), **kwargs):
        """Plot results of 2 -layer analytical simulation.

        Parameters
        ----------
        titles: 2-list of 2-list of titles
            titles for the two head graphs, titles for the two flow graphs.
        xlabel: str
            xlabel
        ylabels: 2-list of 2-list of titiles
            y-axis titles for the head graphs and for the flow graphs.
        size_inches: 2 tuple (w, h)
            Size of each of the two figures.
        kwargs: dict
            Additional paramters to pass.

        Returns
        -------
        None, however, the axes of the two head and two flow plots
        are stored in self.ax as a [2, 2] arrray of axes. Note that
        self.ax[:, 0] are the head axes and self.ax[:, 1] are the flow axes.
        """
        #self.ax = plot_2layer_heads_and_flows(data=self.data, props=self.props)

        titles = [f'({self.name}) ' + title for title in titles]

        fig, self.ax = plt.subplots(3, 1, sharex=True)
        fig.set_size_inches(size_inches)

        plot_heads(ax=self.ax[0], data=self.data, title=titles[0],
                    xlabel='time', ylabels=ylabels[0])

        if self.name == 'modflow':
            gn.plot_watbal(ax=self.ax[1:], data=self.data, titles=titles[1:], xlabel=xlabel,
                        ylabels=ylabels[1:], **kwargs)

        else:
            plot_watbal(ax=self.ax[1:], data=self.data, titles=titles[1:], xlabel=xlabel,
                        ylabels=ylabels[1:], single_layer=False, **kwargs)

        startyr = self.data.index[0].year + 1

        plotGXG(ax=self.ax[0], data=self.data, startyr=startyr, nyr=8)

        plot_hydrological_year_boundaries(ax=self.ax, startyr=startyr, nyr=8)


    def plotGXG(self, ax=None, startyr=None, nyr=8):
        """Plot the points contributing to GLG, GHG and GVG respectively.

        Parameters
        ----------
        startyr: int
            year of first hydrological year to include in GXG computation
        nyr: int
            number of years to include in computation of GXG. (Default = 8)
        """
        plotGXG(ax=self.ax[0], data=self.data['h0'], startyr=startyr, nyr=nyr)


    def plot_hydrological_year_boundaries(self, startyr=None, nyr=8):
        """Plot the boundaries of the hydrological years on self.ax.

        Parameters
        ----------
        startyr: int
            The first hydraulical year (April 1 - March 31).
        nyr: int
            The number of years to include.
        """
        for ax in self.ax.ravel():
            plot_hydrological_year_boundaries(ax=ax, startyr=startyr, nyr=nyr)


    def getGXG(self, startyr=None, nyr=8):
        """Add boolean coluns GLG, GHG and GVG columns to self.data.

        These columns indicate which lines/dates pertain to the GXG.
        Get the values using self.data['h0'].iloc['GHG'] etc.

        Parameters
        ----------
        startyr: int
            year of first hydrological year to include in GXG computation
        nyr: int
            number of years to include in computation of GXG. (Default = 8)

        Returns
        -------
        gxg: 3-tuple of floats
            (glg, ghg, gvg)
        """
        glg, ghg, gvg = getGXG(data=self.data, startyr=startyr, nyr=nyr)
        return glg, ghg, gvg

# Specific analytical solution as classes derived from base class "Solution".
class L1f(Solution):
    """Return analytical solution with given phi in regional aquifer.

    The solution has drainage, ditches in top aquifer only with different
    entry and exit ditch resistance. Regional head is given.
    """

    def __init__(self, props=None):
        super().__init__(props=props)
        self.name = 'L1f'

class L1q(Solution):
    """Return analytical solution wiith given seepage from regional aquifer."""

    def __init__(self, props=None):
        super().__init__(props=props)
        self.name = 'L1q'

class L1(Solution):
    """One layer aka Kraaijenhoff vd Leur (Carslaw & Jaeger (1959, p87))."""

    def __init__(self, props=None):
        super().__init__(props=props)


class L2(Solution):
    """Return analytic two-layer solution."""

    def __init__(self, props=None):
        super().__init__(props=props)
        self.name = 'L2'

    def sim(self, data=None):
        """Simulate 2-layer system using multilayer analytical solution."""
        self.data = multi_layer_transient(data=data, **self.props)

class Lnum(Solution):
    """Return numeric solution using MODFLOW."""

    def __init__(self, props=None):
        super().__init__(props=props)
        self.name = 'modflow'

    def sim(self, data=None):
        """Simulate 2-layer system using multilayer analytical solution."""
        self.data = gn.modflow(props=props, data=data)


if __name__ == '__main__':
    # home = '/Users/Theo/GRWMODELS/python/GGOR/'

    data = knmi.get_weather(stn=240, start='20100101', end='20191231')
    set_hLR(data=data, props=props)

    q = 0.005 # md/d
    ayear = 365 # days
    hyear = 182 # days
    dh1 = props['c'][0] * q # phi change as equivalent to q

    data = gen_testdata(data=data,
                          RH  =(2 * ayear, 0.0, 0.002 * 0., 0.002 * 0.),
                          EV24=(2 * ayear, 0.0, 0.0, 0.0),
                          #hLR =(1 * hyear, 0.0, 0.0,  -0.0, 0., 0., ),
                          q   =(2 * ayear, 0.0, q * 0., 0. -q * 0.),
                          h1  =(2 * ayear, 0.0, -dh1 * 0., -dh1 * 0., 0., dh1 * 0.),
                          hdr =(5 * hyear, props['hdr'] * 0., 0.)
                          )

    ttiles = ['heads', 'flows layer 0', 'flows layer 1']
    ylabels = ['m', 'm/d', 'm/d']
    size_inches = (14, 10)

    if False: # analytic with given head in regional aquifer
        l1f = L1f(props=props)
        l1f.sim(data)
        l1f.plot(titles=titles, ylabels=ylabels, size_inches=size_inches)
    elif False: # analytic with given seepage from regional aquifer
        l1q = L1q(props=props)
        l1q.sim(data=data)
        l1q.plot(titles=titles, ylabels=ylabels, size_inches=size_inches)
    elif False: # Analytic two layers, with ditches in both aquifers
        l2 = L2(props=props)
        l2.sim(data)
        l2.plot(titles=titles, ylabels=ylabels, size_inches=size_inches)
    elif True: # numerical with dichtes in both aquifers
        mf = Lnum(props=props)
        mf.sim(data)
        mf.plot(titles=titles, ylabels=ylabels, size_inches=size_inches)
    else:
        print("?? Nothing to do!!")
