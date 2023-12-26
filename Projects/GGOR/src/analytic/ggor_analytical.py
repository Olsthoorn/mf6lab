#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%% Info
"""
Created on Tue May 22 10:58:55 2018.

Analyical transient simulations of cross sections between parallel ditches.

The problem considered is a cross section of a parcel between side
ditches with the same properties, subjected to varying recharge and water
levels in the ditches.

Computed are the cross-section averaged head, the discharge to the ditches.
Input are time-varying net recharge, ditch level and upward seepage.

Note that all input tdata are considered average values from the previous
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

@Theo Olsthoorn 20200901

Version 4. Contains true circumference of ditches and possibly penetration of ditch into regional aquifer.

@ 20200925
"""
# Modules
import numpy as np
import pandas as pd
import scipy.linalg as la
import matplotlib.pyplot as plt
import pdb

import os
PYTHON = os.path.expanduser('~/GRWMODELS/python/')
TOOLS  = os.path.expanduser('~/GRWMODELS/python/tools')

import sys
sys.path.insert(0, TOOLS)
sys.path.insert(0, PYTHON)

from KNMI import knmi
import GGOR.src.numeric.gg_mflow as gt
#import GGOR.src.numeric.gg_mflow_1parcel as gn

NOT = np.logical_not
AND = np.logical_and
OR  = np.logical_or

ayear = 365 # days

# The required properties dictionary
props_required = {'b': 'Section half width [m]',
         'O_parcel': 'Parcel circumference [m]. Not used analytically.',
         'A_parcel': 'Parcel circumference [m2], Not used analytically.',
         'AHN': 'Ground sruface elevation [m above national datum]',
         'd_drain': 'Drain depth below AHN [m]',
         'c_drain': 'Aeal drainage resistance [d]',
         'c_CB': 'Conf. bed resistance between the two layers (aquifers} [d]',
         'b_ditch': 'Ditch half-width [m]',
         'd_ditch': 'Ditch depth from AHN [m]',
         'wo_ditch': 'Outflow resistance of ditch [d]',
         'wi_ditch': 'Inflow resistance of ditch [d]',
         'n_trench': 'Numb. of trenches in section [-]. Not used analytically.',
         'd_trench': 'Trench depth relative to AHN [m]. Not used analytically.',
         'D1': 'Thickness of top aquifer (coverlayer) [m]',
         'D_CB': 'Thickness of confing bed [m]',
         'D2': 'Thickness of bot aquifer (regional) [m]',
         'sy': 'Specific yield of top aquifer [-]',
         'S2': 'Storage ceofficient of regional aquifer [-]',
         'kh': 'Hydraulic cond. of top aquifer [m/d]',
         'kv': 'Vertical hydr. cond. of top aquifer [m/d]',
         'kh2': 'Hydraulic cond. of regional aquifer [m/d]',
         'kv2': 'Vertical hydr. cond. of regional aquifer [m/d].',
         'h_summer': 'Summer ditch level [m like AHN]',
         'h_winter': 'Winter ditch level [m like AHN]',
         'q_up' : 'Upeward positive seepage from regional aquif into top aquif [m/d]',
          }

# Functions
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
            yscales=['linear', 'linear'],
            sharex=True,
            sharey=False,
            xlims=None,
            ylims=None,
            size_inches=(12, 6),
            **kwargs):
    """Return ax[0], ax[1] for new plot."""
    fig, ax = plt.subplots(2, 1, sharex=sharex, sharey=sharey)
    fig.set_size_inches(size_inches)
    for a, title, ylabel, yscale in zip(ax, titles, ylabels, yscales):
        a.set_title(title)
        a.set_xlabel(xlabel)
        a.set_ylabel(ylabel)
        a.set_xscale(xscale)
        a.set_yscale(yscale)
        a.grid()
    if xlims is not None:
        ax[0].set_xlim(xlims[0])
        ax[1].set_xlim(xlims[1])
    if ylims is not None:
        ax[0].set_ylim(ylims[0])
        ax[1].set_ylim(ylims[1])
    return ax


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

def getDT(tdata):
    """Return stepsize array from tdata.

    Dt is the diff of the index of tdata prepended wit one value
    equal to the value of the first Dt, so that Dt has the same
    lenght as the index.

    Parameters
    ----------
    tdata: pd.DataFrame
        the time dependent data
    """
    Dt = np.diff((tdata.index - tdata.index[0]) / np.timedelta64(1, 'D'))
    Dt = np.hstack((Dt[0], Dt))
    return Dt

def add_var_to_tseries(tdata=None, var='q_up', props=None):
    """Return time series variable var added.

    Add column var to series getting the value or monthly value from dict
    or series props. The var must be in props like 'var' or 'var01', var02,
    ... 'var12' in case the variable is given by monthly values.

    Used for specified heads ('phi')and seepage ('q_up')

    Parameters
    ----------
    tdata: time-series data for recharge
        the date and time of the series
    props: dict (pd.Series with props as keys)
        properties of the parcel.
    """
    if var in props:
        tdata.loc[:, var] = props[var]
    else:
        tdata[var] = np.nan
        # replace by specified monthly values
        month = np.array([t.month for t in tdata.index], dtype=int)
        for m in range(1, 13):
            try:
                var = '{:}{:02d}'.format(var, m)
                tdata.loc[m == month, var] = props[var]
            except:
                raise ValueError("Missing vairable '{}' and or monthly var '{} in parcel_data"
                                 .format(var[:-2], var))
    return tdata

def prop2tvalues(tindex=None, var=None, props=None):
    """Return array of time values with var specified in props.

    variable 'var' must be in props or correspondign monthly values must
    be in props specified as var01, var02 ... var12.
    A LookupError is raised if neither is present.

    Parameters
    ----------
    tindex: pd.DateTime index
        times corresponding to values
    var: str
        variable name in props
    props: dict or series
        properties from which var to pick
    """
    x = np.zeros(len(tindex))
    if var in props:
        x[:] = props[var]
    else:
        month = np.array([t.month for t in tdata.index], dtype=int)
        for m in range(1, 13):
            x[month == m] = props[f'{var:}{m:02d}']
    return x



#% Stand-alone simulator
# This simulator is meant to work indpendently of the class solution
# It must therefore contain properties, solution name and the different
# methods that are designed for specific solutions.
def single_Layer_transient(solution_name=None, parcel_data=None, tdata=None,
                           use_w_not_c=False):
    """Return results tim_data with simulation results included in its columns.

    Parameters
    ----------
    solution_name: str, solution name, one of
        'l1f' : Single layer with given head in regional aquifer
        'l1q' : Single layer with upward seepage no ditch resistance
        'l2q  : Two aquifer system with given seepage from lower aquifer'
        'l2qw': Same but with ditch resistance
    props: dict
        required spacial propeties.
    tdata: pd.DataFrame
        required fields: hLR, RH, EV24, q_up or h1
            hLR: left-right ditch-water level
            RH: precipitation
            EV24: evapotranspiration
            q_up: net seepage=, upward positive
            h1: head in the regional aquifer
        RH and EV24 must be in m/d!
    use_w_not_c: bool
        triggers use of specified analytical ditch resistance w instead of
        computing the analytical ditch resistance from the real c, the
        ditch circumference and the resistance from pertial penetration of
        the ditch.

    Returns
    -------
    tdata: pd.DataFrame
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
    def get_omega(AHN=None, hlr=None, D1=None, D_CB=None, d_dtich=None, b_ditch=None):
        """Return half ditch circumference from circumstances.

        Parameters
        ----------
        AHN: float
            ground surface elevation [m]
        hlr: float
            elevation of water level in ditch [m]
        D1: float
            thickness of top aquifer (cover layer) [m]
        D_CB: float
            thickness of intermediate aquitard [m]
        d_ditch: depth of ditch bottom below ground surface
        b_dtich: half-width of ditch [m]

        Returns
        -------
        omega1, omgega2: (float,float)
            half circumference of ditch in top and bottom aquifer
        """
        z_cov_bot = AHN - D1        # elev of bottom of first aquifer
        z_reg_top = AHN - D1 - D_CB # elevatin of top of second aquifer

        z_ditch_bot1   = max(z_cov_bot, AHN - d_ditch) # elev ditch bottom in aquifer1
        z_ditch_bot2   =                AHN - d_ditch  # elev ditch bottom in aquifer2

        omega1 = 0. if hlr < z_ditch_bot1 else \
            (hlr - z_ditch_bot1) + min(b_ditch, z_ditch_bot1 - z_cov_bot)
        omega2 = 0. if z_reg_top < z_ditch_bot2 else \
            b_ditch + z_reg_top - z_ditch_bot2
        return omega1, omega2


    def w_from_c(cd=None, D=None, omega=None, kh=None, kv=None):
        """Return analytical ditch resistance from real circumstances.

        Parameters
        ----------
        cd: float
            ditch bottom resistrance [d]
        D:  float
            aquifer thickness [m]
        omega: float
            half circumference of ditch [m]
        kh: float
            horizontal aquifer conductivity [m/d]
        kv: float
            vertical aquifer conductivity
        """
        if omega == 0.:
            return np.inf
        else:
            w = cd * (D / omega) + 2 * D / (
                np.pi * np.sqrt(kh * kv)) * np.log(D * np.sqrt(kh / kv) / omega)
            return w

    def w_ghb_riv_from_wi_wo(wi=None, wo=None):
        """Return modflow equivalent ghb and riv resistance from wo wi.

        Parameters
        ----------
        wi: float
            analytical resistance for flow from ditch to aquifer (i=in) [d]
        wo: float
            analytical resistance for flow from aquifer to ditch (o=out):

        Returns
        -------
        w_ghb, w_riv: (float, float)
            modflow equivalent resistances [d]
        """
        wtol = 1e-3 # d
        w_ghb = wi
        if np.isinf(wo) or wi < wo + wtol:
            w_riv = np.inf
        else:
            w_riv = wo / (1 - wo / wi)
        return w_ghb, w_riv


    def update_phi(q_up=None, ht=None, hlr=None, c=None, w=None, b=None, D=None):
        """Return update of phi in case seepage flow from below is given.

        Parameters
        ----------
        q_up: float
            prescribed seeage (upward postiive)
        ht: float
            current head in top aquifer [m]
        hlr: float
            water level in dtich [m]
        c: float
            resistance of intermediate aquitard [m]
        w: float
            analytic resistance of ditch [m/d]
        b: float
            half-width of cross section [m]
        D: float
            thickness of second aquifer
        """
        cb = w * b / D
        return (ht + c  / cb * hlr + c * q_up) / (1 + c / cb)


    def get_B(N=None, qdr=None, phi=None, hlr=None, c=None, Lam=None):
        """Compute the B variable."""
        B      = (N - qdr) * c - ((N - qdr) * c - (hlr - phi)) * Lam
        return B


    def move_on(t=None, dtau=None, ht=None, T=None, phi=None, B=None):
        e = np.exp(-dtau / T)     # At and of current time step.
        f = (1 - e) / (dtau / T)  # Exact aveage over current time step.
        hend    = phi + (ht - phi) * e + B * (1 - e)
        dhdt    = (phi - ht + B) * e / T
        h0mean  = phi + (ht - phi) * f + B * (1 - f)
        return t + dtau, hend, dhdt, h0mean


    def get_runoff(tau=None, N=None, hlr=None, h0=None, phi=None, hdr=None, Lam=None, c=None, T=None):
        """Return drainage required to lower the water table to hdr by end of time step."""
        e = 0. if tau == np.inf else np.exp(-tau / T)
        qdr = N + ((hlr - phi) /c * Lam + (
            (h0 - phi) / c * e  - (hdr - phi) / c) / (1 - e)) / (1 - Lam)
        return qdr


    # Only with given seepagre, do we iterate (once) to update phi to a beter average during the time step.
    # However, even 1 (no) iterations would be fine if the time step is small relative to the system time T = mu c
    # Hence one may set NITER to any positive int. This should push phi during the time step to its
    # the real average value yielding the desired seepage q_up on avearge during the time step. Advice, stick
    # to NITER=2, it's sufficient and takes least coputation time.
    NITER = 1 if solution_name == 'L1f' else 2

    tdata = tdata.copy() # leave original intact

    check_cols(tdata, ['RH', 'EV24', 'summer'])

    # Geenrate heads object and a cbc object.
    HDS = HDS_obj(tdata=tdata, parcel_data=parcel_data, solution_name=solution_name)
    CBC = CBC_obj(tdata=tdata, parcel_data=parcel_data, solution_name=solution_name)

    cbc = CBC.data

    Dt   = getDT(tdata)

    htol = 1e-3 # d, head tolerance

    nper = len(tdata)

    # We need heads one longer than number of time steps to get a start.
    h0 = np.zeros(nper + 1) # h0[i] end of prev. time step, start of current one
    h1 = np.zeros(nper + 1) # h1[i] end of prev. time step, start of current one

    print('Sumulating {}, NITER={}'.format(solution_name, NITER))

    for ip in range(len(parcel_data)):

        print('Parcel {}'.format(ip), end='')

        # Get parcel properties
        props = parcel_data.iloc[ip]

        AHN = props['AHN'] # Ground surface elevation
        b, c   = props['b'], props['c_CB']

        if use_w_not_c:
            # both layers the same wo and wi
            wo1, wo2 = props['wo_ditch'], props['wo_ditch2']
            wi1, wi2 = props['wi_ditch'], props['wi_ditch2']
        else:
            co, ci = props['co_ditch'], props['ci_ditch']
            d_ditch, b_ditch = props['d_ditch'], props['b_ditch']

        mu, S  = props['sy'], props['S2']
        kh1, kv1, kh2, kv2 = props['kh'], props['kv'], props['kh2'], props['kv2']
        hdr    = AHN - props['d_drain']
        phi, D1, D_CB, D2  = props['phi'], props['D1'], props['D_CB'], props['D2']

        lam = kh1 * D1 * c


        tdata['hLR'] = props['h_winter']
        tdata.loc[tdata['summer'], 'hLR'] = props['h_summer']

        # add column q_up to tdata, using data from props
        tdata = add_var_to_tseries(tdata, var='q_up', props=props)
        tdata = add_var_to_tseries(tdata, var='phi' , props=props)

        # Initialize
        h0[0] = tdata['hLR'].iloc[0]
        h1[0] = phi

        was_summer = not tdata['summer'].iloc[0]
        # Iterate this parcel over all time steps
        for it, (dt, t1, hlr, prec, ev24, q_up, phi, is_summer) in enumerate(zip(
                                Dt,                   # time step
                                np.cumsum(Dt),        # t1 (end of time step)
                                tdata[ 'hLR'].values, # ditch water level
                                tdata[  'RH'].values, # precipirtation
                                tdata['EV24'].values, # evaoptranspiration
                                tdata['q_up'].values, # upward seepage series (used if solution is 'L1q')
                                tdata[ 'phi'].values, # head in lower aquifer (used if solution is 'L1f')
                                tdata['summer'].values # season indicator
                                )):
            # Initialize before each time step
            N = prec - ev24 # Current net recharge

            change_season = is_summer != was_summer;   was_summer = is_summer

            if change_season:
                if use_w_not_c:
                    pass
                else:
                    omega1, omega2 = get_omega(AHN=AHN, hlr=hlr,
                            D1=D1, D_CB=D_CB, d_dtich=d_ditch, b_ditch=b_ditch)
                    wo1 = w_from_c(cd=co, D=D1, omega=omega1, kh=kh1, kv=kv1)
                    wi1 = w_from_c(cd=ci, D=D1, omega=omega1, kh=kh1, kv=kv1)
                    wo2 = w_from_c(cd=co, D=D2, omega=omega2, kh=kh2, kv=kv2)
                    wi2 = w_from_c(cd=ci, D=D2, omega=omega2, kh=kh2, kv=kv2)

                # Ditch restance to resistance used in MODFLOW's GHB and RIV packages
                w_ghb1, w_riv1 = w_ghb_riv_from_wi_wo(wi=wi1, wo=wo1)
                w_ghb2, w_riv2 = w_ghb_riv_from_wi_wo(wi=wi2, wo=wo2)

            w1 = wi1 if h0[it] > hlr else wo1
            Lam  = 1 / ((b / lam) / np.tanh(b / lam) + (w1 / D1) * (c / b))

            w2 = wi2 if hlr > h1[it] else wo2
            if solution_name == 'L1q':
                phi = update_phi(q_up=q_up, ht=h0[it], hlr=hlr, c=c, w=w2, b=b, D=D2)

            T = mu / (1/c + 1/w1 * D1 / b)

            qroff = 0.
            B = get_B(N=N, qdr=qroff, phi=phi, hlr=hlr, c=c, Lam=Lam)
            rising = (phi - h0[it] + B) > 0

            for iter in range(NITER): # Use update phi in case q_up is specified (L1q instead of L1f)
                ht = h0[it]
                if iter == NITER - 1:
                    qdr = 0.
                    qv0 = 0.
                    qb0 = 0.
                    qb1 = 0.
                    qghb1, qghb2 = 0., 0.
                    qriv1, qriv2 = 0., 0.

                # Because the head change during a time step is monotonous, due to constant
                # boundary conditions during the time step, each time step may be
                # split into at most two sub steps. The first being
                # upto the moment that a rising water table intersects the drainage level
                # within the time step, and the rest of the time step when it will
                # be at the drainage level. Each time step in which there's no hit
                # will be convered in a single step.
                # Water tables above the drainage level cannot occur.
                t = t1 - dt
                while t < t1:
                    if ht > hdr - htol:
                        if ht > hdr + htol: # ht > drainage level --> let drain during system time to drainge level
                            qroff = get_runoff(tau=T, N=N, hlr=hlr, h0=ht, phi=phi, hdr=hdr, Lam=Lam, c=c, T=T)
                        elif rising: # het at drainage level, if rising keep it there
                            qroff = get_runoff(tau=np.inf, N=N, hlr=hlr, h0=ht, phi=phi, hdr=hdr, Lam=Lam, c=c, T=T)
                        B = get_B(N=N, qdr=qroff, phi=phi, hlr=hlr, c=c, Lam=Lam)
                        dtau = t1 - t
                    else: # ht below drainage level, normal situations, check for hit of drainage level.
                        qroff = 0.
                        B = get_B(N=N, qdr=qroff, phi=phi, hlr=hlr, c=c, Lam=Lam)
                        r = (phi - ht + B) / (phi - hdr + B)
                        future_hit  = r > 1
                        if future_hit and rising: # we hit hdr in the future
                            thit = t + T * np.log(r)
                            if thit <= t1:  # water table will reach hdr now!
                                dtau = thit - t
                            else:
                                dtau = t1 - t
                        else:
                            dtau = t1 - t

                    t, ht, dhdti, h0mean = move_on(
                                t=t, dtau=dtau, ht=ht, T=T, phi=phi, B=B)

                    if iter == NITER - 1:
                        qv0 += (phi - h0mean) / c * dtau / dt              # exact for fixed phi:
                        qdr += qroff  * dtau / dt                          # exact for drainage
                        qb0 += (N - qdr - (hlr - phi) / c) * Lam * dtau/dt #exact for fixed phi
                        qb1 += (phi - hlr) * D2 / w2 / b * dtau / dt       # ziede theorie

                        # Split ditch flow over ghb en riv (riv only works when outflow to ditch)
                        qghb1 += -qb0 if qb0 < 0. else  -qb0 / (1. + w_ghb1/ w_riv1)
                        qriv1 +=   0. if qb0 < 0. else  -qb0 / (w_riv1 / w_ghb1 + 1.)
                        qghb2 += -qb1 if qb1 < 0. else  -qb1 / (np.inf if np.isinf(w_ghb2) else 1. + w_ghb2 / w_riv2)
                        qriv2 +=   0. if qb1 < 0. else  -qb1 / (np.inf if np.isinf(w_riv2) else 1. +  w_riv2 / w_ghb2)

                # Update phi for the next iterataion
                if solution_name == 'L1q':
                    phi = update_phi(q_up=q_up, ht=0.5 * (h0[it] + ht), hlr=hlr, c=c, w=w2, b=b, D=D2)

            #if ht > 0.18:
            #    print('stop hier')

            # results below divided by 2 due to iteration iter
            h0[it + 1] = ht
            h1[it + 1] = phi
            # Getting the water budget terms in the way used by MODFLOW.
            # Modflow's water budget parts summed always yiedsl zero for each cell.
            cbc['RCH'][0, ip, it] =   prec
            cbc['EVT'][0, ip, it] =  -ev24
            cbc['STO'][0, ip, it] =  -mu * (h0[it + 1] - h0[it]) / dt
            cbc['DRN'][0, ip ,it] =  -qdr
            cbc['GHB'][0, ip, it] = qghb1
            cbc['RIV'][0, ip, it] = qriv1
            cbc['FLF'][0, ip, it] =  +qv0

            # To the dtiches:
            cbc['STO'][1, ip, it] =  -S  * (h1[it + 1] - h1[it]) / dt
            cbc['GHB'][1, ip, it] = qghb2
            cbc['RIV'][1, ip, it] = qriv2
            cbc['FLF'][1, ip, it] =  -qv0
            cbc['WEL'][1, ip, it] = qv0 - qghb2 - qriv2 - cbc['STO'][1, ip, it]
            if it % 100 == 0: print('.', end='')
            # END it
        # Gather results and store in tdata columns
        HDS.data[0, ip, :] = h0[1:]
        HDS.data[1, ip, :] = h1[1:]
        print('{} time steps'.format(it))

    HDS.GXG = GXG_obj(HDS)
    print('Done.')

    return tdata, HDS, CBC


def single_layer_steady(solution_name, props=None, tdata=None, dx=1.0, verbose=False, use_w_not_c=None):
    """Return simulated results for cross section the mean situation in tdata.

    Parameters
    ----------
    solutions_name: str, the name of the case/solution
        'L1q' : Single layer given seepage.
        'L1h' : Single layer with given head in regional aquifer.
    props: dict
        required properties of the system. Must contain the static data required by
        the model including 'h_summer' and 'h_winter'', 'phi' and 'q_up'.'
    tdata: dict or pd.DataFrame or pd.Series
        Data for w hich to simulate the steady situation.
        Must contain 'RH', 'EV24' and the index must be a DateTime. Generally
        tdata is one record of a larger time series, for which a steady
        state result is desired.
    dx: float
        Desired x-axis step size in cross section. The width is known from
        the properties.
    use_w_not_c: bool
       Telling whether to use the analytical ditch resistance w or the real one c
       with the effect of partially penentrating ditches.
    """
    htol = 1e-5 # m
    MAXITER = 100
    if isinstance(props, pd.DataFrame):
        props = props.iloc[0]

    kh, kv, D, c, b = props['kh'], props['kv'], props['D1'], props['c_CB'], props['b']

    if use_w_not_c:
        wo, wi =  props['wo_ditch'], props['wi_ditch']
    else:
        co, ci, omega = props['co_ditch'], props['ci_ditch'], props['ditch_omega1']
        wpp1 =   2 / (np.pi * np.sqrt(kh * kv)) * np.log(D * np.sqrt(kh / kv)/ props['ditch_omega1'])
        wo = co * D / omega + wpp1
        wi = ci * D / omega + wpp1

    lam = np.sqrt(kh * D * c)

    assert isinstance(tdata, pd.DataFrame), "tdata must be a pd.DataFrame not a {}".format(type(tdata))

    # X coordinates
    x = np.unique(np.hstack((np.arange(0, b, dx), b))) # make sure b is included
    x = np.unique(np.hstack((-x[::-1], x[1:]))) # generate for total section -b < x < b

    tdata['q_up'] = props['q_up']
    tdata['phi' ] = props['phi']
    tdata['summer'] = [t.month >= 4 and t.month < 9 for t in tdata.index]
    tdata['hLR'] = props['h_winter']
    tdata.loc[tdata['summer'], 'hLR'] = props['h_summer']
    t0 = tdata.index[0]

    tdata = dict(tdata.loc[t0])

    N = tdata['RH'] - tdata['EV24']

    hx_prev = np.ones_like(x) *  tdata['hLR']

    success = False
    if solution_name == 'L1':
        for iter in range(MAXITER):
            w = wo if hx_prev.mean() > tdata['hLR'] else wi
            hx = tdata['hLR'] + N / (2 * kh * D) * (b ** 2 - x ** 2) + N * b * w  / D
            if np.mean(np.abs(hx - hx_prev) < htol):
                success = True
                break
            else:
                hx_prev = hx
    else:
        for iter in range(MAXITER):
            w = wo if hx_prev.mean() > tdata['hLR'] else wi
            Lamb = 1 / (b  / lam / np.tanh(b  / lam) + (b /c) / (D / w))
            phi = tdata['phi'] if solution_name == 'L1f' else np.mean(hx_prev) + tdata['q_up'] * c

            hx = phi + N * c - (N * c - (tdata['hLR'] - phi)
                            ) * Lamb * (b / lam) * np.cosh(x / lam) / np.sinh(b / lam)
            if verbose:
                print('iter={}, std(hx - hk_prev) = {}'.format(iter, np.std(hx - hx_prev)))
            if np.std(hx - hx_prev) < htol:
                hx = np.vstack((hx, phi * np.ones_like(x)))
                success = True
                break
            else:
                hx_prev = hx

    if not success:
        print("No sucess for solution name = {}".format(solution_name))
    return hx, x, success


def sysmat(c, kD):    # System matrix
    """Return the system matrix and related matrices for mutlti-layer solutions.

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
    Am2 = -np.diag(dL[1:].ravel(), k=-1) - np.diag(dR[:-1].ravel(), k=1)\
                                            + np.diag((dL + dR).ravel(), k=0)
    Am1 = la.sqrtm(Am2)
    A1  = la.inv(Am1)
    A2  = A1 @ A1
    return Am2, Am1, A1, A2


def multi_layer_steady(props=None, tdata=None, dx=1., plot=True, use_w_not_c=None, **kwargs):
    """Return steady multilayer solution for symmetrial cross section.

    With all tdata in a dict kwargs, call this function using
    multi_layer_steady(**kwargs).
    No tdata pd.DatFrame input is needed, because not transient.
    If plot==True, then a graph is made showing the heads in
    the corss section for all layers and printing the water
    balance.

    Parameters
    ----------
    props: dict
        properties of the parcel
    tdata: dict
        time data for the moment to compute the steady state solution
    use_w_not_c: bool
       Telling whether to use the analytical ditch resistance w or the real one c
       with the effect of partially penentrating ditches.

    @TO 20200615, 20200908
    """
    WINF = 1e9 # in case w is inf (e.g. no ditch present in reg. aquifer)
    nlay = 2

    b  = props['b']
    kh = np.array([props['kh'], props['kh2']])[:, np.newaxis]
    kv = np.array([props['kv'], props['kv2']])[:, np.newaxis]
    D  = np.array([props['D1'], props[ 'D2']])[:, np.newaxis]
    hdr, cdr, c = props['AHN'] - props['d_drain'], props['c_drain'], props['c_CB']
    kD = kh * D
    #pdb.set_trace()
    if use_w_not_c:
        assert props['wi_ditch'] >= props['wo_ditch'],\
            AssertionError("wi_ditch must >= wo_ditch")
        assert props['wi_ditch2'] >= props['wo_ditch2'],\
            AssertionError("wi_ditch2 must >= wo_ditch2")
        wo = np.array([props['wo_ditch'], props['wo_ditch2']])[:, np.newaxis]
        wi = np.array([props['wi_ditch'], props['wi_ditch2']])[:, np.newaxis]
    else:
        assert props['ci_ditch'] >= props['co_ditch'],\
                AssertionError("ci_ditch must >= co_ditch]")
        omega = np.array([props['ditch_omega1'], props['ditch_omega2']])[:, np.newaxis]
        wpp = 2 / np.pi / np.sqrt(kh * kv) * np.log(D * np.sqrt(kh / kv) / omega)
        wo = props['co_ditch'] * D / omega + wpp # wpp = due to partial penetration
        wi = props['ci_ditch'] * D / omega * wpp

    tdata['q_up'] = props['q_up']
    tdata['summer'] = [t.month >= 4 and t.month <= 9 for t in tdata.index]
    tdata['hLR'] = props['h_winter']
    tdata.loc[tdata['summer'], 'hLR'] = props['h_summer']
    tdata = tdata.iloc[0]

    x = np.unique(np.hstack((np.arange(0, b, dx), b)))
    x = np.unique(np.hstack((-x[::-1], x)))

    c = np.array([cdr, c, np.inf])[:, np.newaxis]
    hLR = tdata['hLR'] * np.ones((nlay, 1))
    q = np.array([tdata['RH'] - tdata['EV24'], tdata['q_up']])[:, np.newaxis]
    g = np.zeros_like(q)
    g[0]  = hdr / (cdr * kD[0])

    w = (wo if q[0] + q[1] > 0 else wi) * np.ones((nlay, 1))
    w[np.isinf(w)] = WINF

    Am2, Am1, A1, A2 = sysmat(c, kD)

    Kw  = np.eye(nlay) * (kh * w)
    T   = np.eye(nlay) * kD
    Tm1 = la.inv(T)

    coshmb = la.coshm(Am1 * b)
    sinhmb = la.sinhm(Am1 * b)
    F = la.inv(coshmb + Kw @ Am1 @ sinhmb)
    I = np.eye(nlay)

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
    return Phi, x


def multi_layer_transient(solution_name=None, parcel_data=None, tdata=None,  check=True, **kwargs):
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
    tdata: a pd.DataFrame holding the dynamic tdata
        RH  : Preciptaion injected in the first aquifer
        EV24: Evapotranspiration from the first aquifer
        qv : seepage injected in the second aquifer
        hLR: Prescribed heads at x=b, vector with one value for each aquifer [L]
            these are the ditch levels in both the first and second aquifer
            as these are assumed the same [L].
    b: Float, half-width of cross section [L].
    k: Array [nlay, 2] of conductivities. First col kh, second col kv [L/T].
    z: Array [nlay, 2] first col top of layer, second bottom of layers [L].
    c: Vector of aquitard resistances below the layers [nlay].
       Set the last equal to inf. Or it's done by default. No leakage from below.
       The top c is that on top of the top aquifer, i.e. normally cdr.
    w: Ditch resistances [nlay, 2] (first col inflow, second col outflow) [T].
    S: storage coeff (nlay, 2) first col Sy,  second col S [-].
    cdr: Float. Resitance of top layer
    hdr: Float. Fixed head above top resistance layer.

    Returns
    -------
    tdata: pd.DataFrame
        Results, tdata copied from input and augmented with simulation results.

    @TO 20200701
    """
    wtol = 1e-3 # [d] minimum difference between wi and wo
    nlay = 2
    tdata = tdata.copy() # keep original tdata intact
    c_bot = 1e4 # large resisance at bottom of system to make sysem matrix invertible

    Dt  = np.hstack(((tdata.index[1] - tdata.index[0]) / np.timedelta64(1, 'D'),
                     np.diff((tdata.index - tdata.index[0]) / np.timedelta64(1, 'D'))))

    HDS = HDS_obj(tdata=tdata, parcel_data=parcel_data, solution_name=solution_name)
    CBC = CBC_obj(tdata=tdata, parcel_data=parcel_data, solution_name=solution_name)

    print('Simularing {}'.format('multi-layer transient.'))

    for ip in range(len(parcel_data)):
        print('Parcel {}'.format(ip), end='')
        props = parcel_data.iloc[ip]

        b = props['b']
        c  = np.array([props['c_drain'], props['c_CB'], c_bot])
        hdr = props['AHN'] - props['d_drain']

        S  = np.array([props['sy'], props[ 'S2']])
        kh = np.array([props['kh'], props['kh2']])
        kv = np.array([props['kv'], props['kv2']])
        D  = np.array([props['D1'], props[ 'D2']])
        kD = kh * D

        if use_w_not_c:
            wo  = np.array([props['wo_ditch'], props['wo_ditch2']])[:, np.newaxis]
            wi  = np.array([props['wi_ditch'], props['wi_ditch2']])[:, np.newaxis]
            assert np.all(wi >= wo), AssertionError("wi must >= wo, [parcel {}]".format(ip))
        else:
            assert props['ci_ditch'] >= props['co_ditch'],\
                    AssertionError("ci must >= co, [parcel {}]".format(ip))
            omega = np.array([props['ditch_omega1'], props['ditch_omega2']])[:, np.newaxis]
            wpp = 2 / np.pi / np.sqrt(kh * kv)  * np.log(D * np.sqrt(kh / kv) / omega)
            wo  = props['co_ditch'] * D / omega + wpp
            wi  = props['ci_ditch'] * D / omega + wpp

        wghb = wi
        wriv = np.inf if wi < wo + wtol else wi * wo / (wi - wo)

        Am2, Am1, A1, A2 = sysmat(c, kD)

        T   = np.diag(kD) # transmissivity
        Tm1 = la.inv(T)

        # unifrom injection into each aquifer
        q       = np.zeros((nlay, len(tdata)))
        q[0, :] = (tdata['RH'] - tdata['EV24']).values[:]
        q[1, :] = prop2tvalues(tdata.index, var='q_up', props=props)

        # Leakage from top and bottom aquifer to outside the model. i.e. hdr / (cdr * kD0)
        g = np.zeros((nlay, len(tdata)))
        g[0,  :]  = hdr / props['c_drain'] # this may be made time variable
        #g[-1, :]  = 0.  / props['c_bot']   # in fact a dummy (happens to be zero

        coshmb = la.coshm(Am1 * b)
        sinhmb = la.sinhm(Am1 * b)
        I = np.eye(len(kD))

        # Compute Phi for all x-values
        G    = A2 @ Tm1 @ np.diag(S)
        E, V = la.eig(G)
        Em1  = la.inv(np.diag(E))
        Vm1  = la.inv(V)

        hLR = np.ones(len(tdata)) * props['h_winter']
        hLR[tdata['summer']] = props['h_summer']

        # Initialize heads
        Phi = np.zeros((2, len(tdata) + 1))
        Phi[0, 0] = hLR[0]

        qb = np.zeros_like(q) # from layers to boundary at x=b

        # Loop over time steps
        for it, (dt, hlr) in enumerate(zip(Dt, hLR)):
            # infiltration or exfiltration ?
            w1 = (wi * (Phi[:, it] < hlr) + wo * (Phi[:, it] >= hlr))[:, np.newaxis]
            Kw = np.diag(kh * w1)
            F = la.inv(coshmb + Kw @ Am1 @ sinhmb)

            e = la.expm(-Em1 * dt)
            hq = A2 @ Tm1 @ (q[:, it:it+1] + g[:, it:it+1])
            qb[:, it:it+1] = T @ Am1 / b @ sinhmb @ F @ (hq - np.ones((nlay, 1)) * hlr)

            # steady state ultimate solution for t->inf
            hss = A2 @ Tm1 @ (q[:,it:it+1] + g[:, it:it+1] - qb[:,it:it+1])

            # Compute head
            Phi[:, it + 1 : it + 2] = V @ e @ Vm1 @  Phi[:, it : it + 1] + V @ (I - e) @ Vm1 @ hss

            if it % 100 == 0: print('.', end='')

        Phim = (Phi[:,:-1] + Phi[:, 1:]) / 2
        qs = S[:, np.newaxis] * np.diff(Phi, axis=1) / Dt  # from top layer into storage

        # leakage through aquitards fram all layers
        ql = np.zeros_like(q)
        for it in range(len(Dt)):
            ql[:,it:it+1] = T @ Am2 @ Phim[:, it:it+1] - g[:, it:it+1]

        Phi = Phi[:, 1:] # cut off first day (before first tdata in index)

        HDS.data[:, ip, :] = Phi

        # Store results
        qb = q - qs - ql      # Flow to ditches in both layers
        # Divide over GHB (in) and RIV according to MODFLOW
        L_in0  = Phi[0] <= hLR
        L_in1  = Phi[1] <= hLR
        L_out0 = Phi[0] >  hLR
        L_out1 = Phi[1] >  hLR
        CBC.data['GHB'][0, ip, L_in0 ] = wriv / (wriv + wghb) * qb[0][L_in0 ]
        CBC.data['GHB'][0, ip, L_out0] = wriv / (wriv + wghb) * qb[0][L_out0]
        CBC.data['GHB'][1, ip, L_in1 ] = wriv / (wriv + wghb) * qb[1][L_in1 ]
        CBC.data['GHB'][1, ip, L_out1] = wriv / (wriv + wghb) * qb[1][L_out1]
        CBC.data['RIV'][0, ip, L_in0 ] = wghb / (wriv + wghb) * qb[0][L_in0 ]
        CBC.data['RIV'][0, ip, L_out0] = wghb / (wriv + wghb) * qb[0][L_out0]
        CBC.data['RIV'][1, ip, L_in1 ] = wghb / (wriv + wghb) * qb[1][L_in1 ]
        CBC.data['RIV'][1, ip, L_out1] = wghb / (wriv + wghb) * qb[1][L_out1]
        CBC.data['STO'][:, ip, :] = qs
        CBC.data['RCH'][0, ip, :] =  tdata['RH'].values
        CBC.data['EVT'][0, ip, :] = -tdata['EV24'].values
        CBC.data['DRN'][0, ip, :] =  ql[0] - ql[1]
        CBC.data['FLF'][0, ip, :] = -ql[1]
        CBC.data['FLF'][1, ip, :] =  ql[1]     #tdata['ql1'] = ql[1]  # leakage from bo1 layer
        CBC.data['WEL'][1, ip, :] =  q[1]
        print(it)

    HDS.GXG = GXG_obj(HDS)

    print('Done!')

    return tdata, HDS, CBC # return updated DataFrame


def plot_hydrological_year_boundaries(ax=None, tindex=None, **kwargs):
    """Plot hydrological year boundaries on a given axis.

    Parameters
    ----------
    ax: plt.Axes
        an existing axes with a datatime x-axis
    tindex: DateTime index
        tindex to use for this graph.
    kwargs: dict
        passed on to plt.axvline
    """
    years = np.unique(np.array([t.year for t in tindex]))

    if not 'color' in kwargs:
        kwargs['color'] = 'gray'

    if isinstance(ax, plt.Axes): ax = [ax]
    for a in ax:
        for yr in years:
            t = np.datetime64(f'{yr}-03-14')
            if t > tindex[0] and t < tindex[-1]:
                a.axvline(t, ls=':', **kwargs)

class HDS_obj:
    """Object to store and retrieve the heads like the on of Modflow."""

    def __init__(self, tdata=None, parcel_data=None, solution_name=None):
        """Return HDS_obj.

        Parameters
        ----------
        hds: ndarray
            heads in shape (nlay, 1, nper)
        tindex: index like pd.Datrame.index
            datetimes for hds
        """
        self.solution_name = solution_name
        self.nparcel = len(parcel_data)
        self.nlay = 2 # for GGOR
        self.nper = len(tdata)
        self.tindex = tdata.index
        self.shape = (self.nlay, self.nparcel, self.nper)
        self.data = np.zeros((self.nlay, self.nparcel, self.nper), dtype=float)

        #self.GXG = GXG_obj(self)

    def plot_gxg(self, ax=None, selection=[0, 1, 3, 4, 5], **kwargs):
        """Plot GXG as markers.

        Parameters
        ----------
        ax: plt.Axes
            axes to plot on, may contain graphs.
        selection: None, int, sequence
            for which parcels to plot the GXG.
        **kwargs
            further arguments passed on to plt.plot()
        """
        selection = gt.selection_check(selection)

        self.GXG.plot(ax=None, selection=[0, 1, 3, 4, 5], **kwargs)


    def plot(self, titles=['top', 'bottom'], xlabel='time',
             ylabels=['head [m]', 'head [m]'], sharex=True, sharey=True,
             selection=None, size_inches=(14, 8), GXG=False):
        """Plot the heads.

        Parameters
        ----------
        selection: int, sequence or slice
            which of the parcels to plot, none means all (first 5)

        """
        selection = gt.selection_check(selection)

        ax = newfig2(titles=titles, xlabel=xlabel, ylabels=ylabels,
                     sharex=sharex, sharey=sharey,
                     size_inches=size_inches)

        clrs = 'brgkmcy'
        for sel in range(self.nparcel)[selection]:
            clr = clrs[sel % len(clrs)]
            ax[0].plot(self.tindex, self.data[0][sel], clr, label='parcel {}'.format(sel))
            ax[1].plot(self.tindex, self.data[1][sel], clr, label='parcel {}'.format(sel))

        ax[0].legend(loc='best', fontsize='xx-small')
        ax[1].legend(loc='best', fontsize='xx-small')

        plot_hydrological_year_boundaries(ax[0], tindex=self.tindex, color='darkgray')
        plot_hydrological_year_boundaries(ax[1], tindex=self.tindex, color='darkgray')

        # plot GXG on graphs of top aquifer only.
        if GXG:
            self.GXG.plot(ax=ax[0], selection=selection, colors=clrs)

        return ax


class GXG_obj:
    """Generate GXG object.

    This object hold the GXG (GLG, GVG, GHG)  i.e. the lowest, hightes and spring
    groundwater head information and their long-term averaged values based on
    the number of hydrologi al years implied in the given time_data.
    (A hydrological year runs form March14 through March 13 the next year, but
     the GXG are based on values of the 14th and 28th of each month only.)

    self.gxg is a recarray with all the individual records. (nyear * 9, nparcel)
    self.GXG is a recarray with the long-time averaged values (nparcel).

    @TO 2020-08-31
    """

    def __init__(self, HDS_obj=None):
        """Initialize GXG object.

        Parameters
        ----------
        time_data: pd.DataFrame
            time_data, we only need its index
        avgHds: np.nd_array shape = (nt, nz, nParcel)
            The xsection-averaged heads for all parcels and all times.
            Heads aveaged along the x-axis taking into account cel width
            and ignoring inactive cells.
        """
        # Hand recordings to compute GXG from (14th and 28th of each month)
        hand   = np.array([t.day % 14 == 0 for t in HDS_obj.tindex], dtype=bool)
        # parcel-average heads (time=row, parcel=col)
        self.ahds    = HDS_obj.data[:, :, hand][0].T   # only top layer
         # Time index of head recordings (14th and 28th of each month
        self.tindex  = HDS_obj.tindex[hand]
        self.nparcel = self.ahds.shape[1]
        self.tshift = self.tindex[0] - np.datetime64(self.tindex[0].date())

        hyear_index = np.array([t.year -1 if t.month < 4 else t.year
                            for t in self.tindex])
        hyears = np.unique(hyear_index)[2:-1] # to be processed

        # Format to store the gxg data in a recarray
        gxg_dtype = [('t', pd.Timestamp), ('hd', float), ('hyear', int),
                 ('l', bool), ('h', bool), ('v', bool)]

        # The gxg recarray has 9 records per hyear (3 glg, 3 ghg, 3gvg and nparcel layers)
        # These are tine individual values contribution to the GXG
        self.gxg = np.zeros((len(hyears) * 9, self.nparcel), dtype=gxg_dtype)

        T = (True, True, True)
        F = (False, False, False)

        for iyr, hyear in enumerate(hyears):
            ah = self.ahds[  hyear_index == hyear]
            td = self.tindex[hyear_index == hyear]
            Ias = np.argsort(ah, axis=0)  # Indices of argsort along time axis

            # Make sure hydrological years start at March 14!!
            #assert td.index[0].month ==3 and td.index[0].day == 14, "hyears must start at 14th of March"

            hyr = (hyear, hyear, hyear)

            for ip in range(self.nparcel):
                Iglg = Ias[ :3, ip]
                Ighg = Ias[-3:, ip]
                #Igvg = resfers to original self.tindex
                Igvg = OR(OR(self.tindex == np.datetime64('{}-{:02d}-{:02d}'.
                            format(hyear, 3, 14)) + self.tshift,
                         self.tindex == np.datetime64('{}-{:02d}-{:02d}'.
                            format(hyear, 3, 28)) + self.tshift),
                         self.tindex == np.datetime64('{}-{:02d}-{:02d}'.
                            format(hyear , 4, 14)) + self.tshift)
                # The three lowest values
                self.gxg[iyr * 9 + 0:iyr * 9 + 3, ip] = np.array(
                    [(t, hd, yr, l, h, v) for t, hd, yr, l, h, v in zip(
                        td[Iglg], ah[Iglg, ip], hyr, T, F, F)], dtype=gxg_dtype)
                # The three highest values
                self.gxg[iyr * 9 + 3:iyr * 9 + 6, ip] = np.array(
                    [(t, hd, yr, l, h, v) for t, hd, yr, l, h, v in zip(
                        td[Ighg], ah[Ighg, ip], hyr, F, T, F)], dtype=gxg_dtype)
                # The three spring values
                self.gxg[iyr * 9 + 6:iyr * 9 + 9, ip] = np.array(
                    [(t, hd, yr, l, h, v) for t, hd, yr, l, h, v in zip(
                        self.tindex[Igvg], self.ahds[Igvg, ip], hyr, F, F, T)], dtype=gxg_dtype)

        # Comptue and store the long-term averaged values, the actual GXG
        dtype = [('id', int), ('GLG', float), ('GHG', float), ('GVG', float)]
        self.GXG = np.ones(self.nparcel, dtype=dtype)
        for ip in range(self.nparcel):
            self.GXG[ip] = (
                ip,
                self.gxg[self.gxg[:, ip]['v'], ip]['hd'].mean(),
                self.gxg[self.gxg[:, ip]['l'], ip]['hd'].mean(),
                self.gxg[self.gxg[:, ip]['h'], ip]['hd'].mean())


    def plot(self, ax=None, selection=[0, 1, 3, 4, 5], **kwargs):
        """Plot GXG.

        Parameters
        ----------
        selection : list
            list if indices to select the parcels for plotting
        nmax: int
            maximum number of graphs to plot
        """
        selection = gt.selection_check(selection)

        if 'colors' in kwargs:
            clrs = kwargs.pop('colors', 'brgkmcy')

        for iclr, ip in enumerate(range(self.nparcel)[selection]):
            g = self.gxg.T[ip]
            clr = clrs[iclr % len(clrs)]
            ax.plot(g['t'][g['v']], g['hd'][g['v']], clr, marker='o',
                    mfc='none', ls='none', label='vg [{}]'.format(ip))
            ax.plot(g['t'][g['l']], g['hd'][g['l']], clr, marker='v',
                    mfc='none', ls='none', label='lg [{}]'.format(ip))
            ax.plot(g['t'][g['h']], g['hd'][g['h']], clr, marker='v',
                    mfc='none', ls='none', label='hg [{}]'.format(ip))

        hyears = np.unique(self.gxg.T[0]['hyear'])
        t = (pd.Timestamp('{}-{:02d}-{:02d}'.format(hyears[ 0], 3, 14)),
             pd.Timestamp('{}-{:02d}-{:02d}'.format(hyears[-1], 2, 28)))

        lw = 0.5
        for iclr, ip in enumerate(range(self.nparcel)[selection]):
            clr = clrs[iclr % len(clrs)]
            ax.hlines(self.GXG['GVG'][self.GXG['id']==ip], *t, clr,
                      ls='solid'  , lw=lw, label='GVG parcel {}'.format(ip))
            ax.hlines(self.GXG['GLG'][self.GXG['id']==ip], *t, clr,
                      ls='dashed' , lw=lw, label='GLG parcel {}'.format(ip))
            ax.hlines(self.GXG['GHG'][self.GXG['id']==ip], *t, clr,
                      ls='dashdot', lw=lw, label='GHG parcel {}'.format(ip))

        ax.legend(loc='best', fontsize='xx-small')
        return ax


class CBC_obj:
    """CBC_obj to store  row by row  flows.

    Stores nlay * nparcel * nper flows for each cbc_label.
    """

    def __init__(self, tdata=None, parcel_data=None, solution_name=None):
        """Gerate a CBC instantiation for the cross sections.

        Parameters
        ----------
        tdata: pd.DataFrame
            time data, (RH, EV24, 'summer')
        parcel_data: pd.DataFrame
            parcel properties
        """
        self.solution_name = solution_name
        self.tindex = tdata.index
        self.nper = len(tdata)
        self.nparcel = len(parcel_data)
        self.nlay = 2
        self.ncol = 1
        self.shape = (self.nlay, self.nparcel, self.nper)
        self.labels = [k for k in gt.watbal_label]
        self.Arel = parcel_data.A_parcel / np.sum(parcel_data.A_parcel)

        self.data = dict()
        for lbl in gt.watbal_label:
            self.data[lbl] = np.zeros(self.shape)

        #dtype = [(lbl, float, (self.nlay, self.nparcel, self.nper))
        #                                         for lbl in self.labels]
        #self.data = np.zeros(1, dtype=dtype)

    def plot(self, titles=['top', 'bottom'], xlabel='time',
             ylabels=['flux [m/s]', 'flux [m/d]'], sharex=True, sharey=False,
             size_inches=(14, 8), selection=None, ax=None, **kwargs):
        """Plot the running water balance.

        Parameters
        ----------
        selection: set or sequence
            the set of parcels to include in the water budget
            ax[0], ax[1] are the axes for plotting the top and bottom layer.
        titles: list of 2 strings
        xlabel: str
        ylabels: list of 2 strings
        size_inches: tuple of two
            Figure size. Only applied when a new figure is generated (ax is None).
        ax: list optional
        kwargs: dict with extra parameters passed to newfig2 if present

        Returns
        -------
        ax: list of two axes
            plotting axes
        """
        selection = gt.selection_check(selection)

        if ax is None:
            ax = newfig2(titles, xlabel, ylabels, sharex=sharex, sharey=sharey,
                         size_inches=size_inches, **kwargs)
            ax[0].set_title(titles[0])
            ax[1].set_title(titles[1])
        else:
            try:
                for a, title, ylabel in zip(ax, titles, ylabels):
                    a.set_title(title)
                    a.set_xlabel(xlabel)
                    a.set_ylabel(ylabel)
                    a.grid(True)
            except:
                raise ValueError('ax must be a list of two plt.Axis objectes.')

        C0   = [gt.watbal_label[k]['clr'] for k in self.data]
        Lbl0 = [gt.watbal_label[k]['leg'] for k in self.data]

        V0 = np.zeros((len(gt.watbal_label), len(self.tindex)))
        V1 = np.zeros((len(gt.watbal_label), len(self.tindex)))

        Arel = self.Arel.values[selection, np.newaxis]
        for i, lbl in enumerate(gt.watbal_label):
            V0[i] = np.sum(self.data[lbl][0, selection, :] * Arel, axis=-2)
            V1[i] = np.sum(self.data[lbl][1, selection, :] * Arel, axis=-2)

        ax[0].stackplot(self.tindex, (V0 * (V0>0)), colors=C0, labels=Lbl0)
        ax[0].stackplot(self.tindex, (V0 * (V0<0)), colors=C0) # no labels
        ax[1].stackplot(self.tindex, (V1 * (V1>0)), colors=C0, labels=Lbl0)
        ax[1].stackplot(self.tindex, (V1 * (V1<0)), colors=C0) # no labels

        ax[0].legend(loc='best', fontsize='xx-small')
        ax[1].legend(loc='best', fontsize='xx-small')

        plot_hydrological_year_boundaries(ax[0], tindex=self.tindex)
        plot_hydrological_year_boundaries(ax[1], tindex=self.tindex)

        return ax

    def get_budget(self, selection=None):
        """Return water budget a pd.DataFrame."""
        if selection is None:
            selection = slice(0, self.nparcel, 1)
        elif isinstance(selection, int):
            selection = slice(0, selection, 1)
        elif isinstance(selection, {tuple, list, np.ndarray, range, slice}):
            pass
        else:
            raise ValueError("selection must be None (for all), an int or a sequence")

        selection = gt.selection_check(selection)

        Arel = self.Arel.values[selection, np.newaxis]

        budget = pd.DataFrame(
            {lbl + f'{i:}':
             np.sum(self.data[lbl][i, selection, :] * Arel, axis=-2)
             for i in (0, 1)
             for lbl in gt.watbal_label.keys()})
        budget.index.name = self.solution_name
        return budget


# Solution class
class Solution:
    """Analytic solution base class.

    Specific analytic solution classes should be derived from this (see below).

    @TO 2020-07-01
    """

    def __init__(self, dirs=None, parcel_data=None):
        """Return an instance of an analytical solution only storing name and properties.

        Parameters
        ----------
        dirs: DirStruct obj
            knows the directory structure of the GGOR project and its case
        props: dict
            a dict containing the properrties. The necessary properties
            are given in the example tamplate in this class. Not all
            properties are used by all solutions. Unsued properties
            may be omitted from the actual properties.
        """
        self.name = str(self.__class__).split('.')[-1].split("'")[0]
        self.case = os.path.basename(dirs.case)
        self.parcel_data = parcel_data
        return


    def check_props(self):
        """Return verified properties.

        Parameters
        ----------
        props: dict
            properties of the section to be simulated
        """
        missing_params = set(self.required_params) - set(self.parcel_data.columns)

        if missing_params:
            raise ValueError('Missing required properties: ({})'.
                             format(missing_params))
        return

    @property
    def cbc(self):
        """Return self.CBC.data."""
        return self.CBC.data

    @property
    def hds(self):
        """Return self.HDS.data."""
        return self.HDS.data

    def check_tcols(self, tdata=None):
        """Verify input data for presence of required columns.

        Parameters
        ----------
        tdata: pd.DataFrame
            Input data with required columns which are 'RH', 'EV24'
        """
        missing_cols = set(['RH', 'EV24', 'summer', 'hand']).difference(tdata.columns)

        if missing_cols:
            raise KeyError("{" + " ,".join([f"'{k}'" for k in missing_cols]) + "} are missing.")


    def sim(self, tdata=None, use_w_not_c=False):
        """Compute and store head and flows in added columns of the input tdata.

        Parameters
        ----------
        tdata: pd.DataFrame with all required time series in columns
        required columns: 'RH','EV24'
            RH: [m/d] precipitation
            EV24: [m/d] evap
        use_w_not_c: bool
            triggers use of prescribed analytical ditch resistance `w` instead
            of computing it from prescribed ditch resistance c, the ditch
            circumference and the extra resistance due to partial ditch
            penetration.

        Returns
        -------
        tdata with extra or overwritten columns:
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
        pdb.set_trace()
        self.tdata, self.HDS, self.CBC = single_Layer_transient(
                solution_name=self.name,
                parcel_data = self.parcel_data,
                tdata=tdata,
                use_w_not_c=use_w_not_c)
        return


    def plot_heads(self, case=None, xlabel='time', ylabels=['m', 'm'],
                   selection=None, size_inches=(14, 8), GXG=False, **kwargs):
        """Plot results of 2 -layer analytical simulation.

        Parameters
        ----------
        case: str
            naem of the current case (used in chart titles)
        xlabel: str
            xlabel
        ylabels: 2-list of 2-list of titiles
            y-axis titles for the head graphs and for the flow graphs.
        selection: None, int, slie or sequence
            of which parcels to plot the heads (max 5)
        size_inches: 2 tuple (w, h)
            Size of each of the two figures.
        GXG: boolean
            tells whether or not the GXG data are to be plotted.
        kwargs: dict
            Additional paramters to pass.

        Returns
        -------
        None, however, the axes of the two head and two flow plots
        are stored in self.ax as a [2, 2] arrray of axes. Note that
        self.ax[:, 0] are the head axes and self.ax[:, 1] are the flow axes.
        """
        #self.ax = plot_2layer_heads_and_flows(tdata=self.tdata, props=self.props)
        case = self.case if case is None else case

        titles = ['Heads top aquifer', 'Heads bottom aquifer']
        titles = [f'{title}, case {case}, (Solution {self.name})' for title in titles]

        titles = [f'({self.name}) ' + title for title in titles]

        ax = self.HDS.plot(titles, xlabel, ylabels, selection=selection,
                      size_inches=size_inches, GXG=GXG, **kwargs)

        return ax


    def plot_cbc(self, case=None, xlabel='time',
             ylabels=['flux [m/d]', 'flux [m/d]'], sharex=True, sharey=False,
             size_inches=(14, 8), selection=None, ax=None, **kwargs):
        """Plot results of 2 -layer analytical simulation.

        Parameters
        ----------
        case: str
            naem of the current case (used in chart titles)
        xlabel: str
            xlabel
        ylabels: 2-list of 2-list of titiles
            y-axis titles for the head graphs and for the flow graphs.
        size_inches: 2 tuple (w, h)
            Size of each of the two figures.
        selection: None, int or sequence, slice, range
            The selection of parcels over whicht to compute the water balance.
        kwargs: dict
            Additional paramters to pass.

        Returns
        -------
        None, however, the axes of the two head and two flow plots
        are stored in self.ax as a [2, 2] arrray of axes. Note that
        self.ax[:, 0] are the head axes and self.ax[:, 1] are the flow axes.
        """
        #self.ax = plot_2layer_heads_and_flows(tdata=self.tdata, props=self.props)

        case = self.case if case is None else case

        titles = ['Fluxes top aquifer.', 'Fluxes bottom aquifer.']
        titles = [f'{title}, case {case}, (Solution {self.name})' for title in titles]

        self.CBC.plot(titles=titles, xlabel=xlabel,
             ylabels=ylabels, sharex=sharex, sharey=sharey,
             size_inches=size_inches, selection=selection, ax=ax, **kwargs)

    def get_budget(self):
        """Return wat budget at given time or index.

        Parameters
        ----------
        timestamp: pd.Timestamp, np.datetime64 or int
            datrtime or index at which water budget is requested.
        """
        return self.CBC.get_budget()

# Specific analytical solution as classes derived from base class "Solution".
class L1(Solution):
    """Class for simulating one layer cross sections without given underlying regional head.

    The solution has drainage, ditches, entry and exit ditch resistance.
    """

    def __init__(self, dirs=None, parcel_data=None):
        """Return a Solution object to simulate 1D transient flow without.

        Parameters
        ----------
        dirs: DirStruct obj
            knows the directory structure of the GGOR project and its case
        parcel_data: pd.DataFrame
            the properties of the parcels
        """
        super().__init__(dirs=dirs, parcel_data=parcel_data)
        self.name = 'L1'

class L1f(Solution):
    """Class for simulating one layer cross sections with given underlying regional head.

    The solution has drainage, ditches in top aquifer only with different
    entry and exit ditch resistance. Regional head is given.
    """

    def __init__(self, dirs=None, parcel_data=None):
        """Return a Solution object to simulate 1D transient flow with given underlying regional head.

        Parameters
        ----------
        dirs: DirStruct obj
            knows the directory structure of the GGOR project and its case
        parcel_data: pd.DataFrame
            the properties of the parcels
        """
        super().__init__(dirs=dirs, parcel_data=parcel_data)
        self.name = 'L1f'

class L1q(Solution):
    """Class for simulating one layer cross sections with given upward regional seepage flux."""

    def __init__(self, dirs=None, parcel_data=None):
        """Return Solution object for simulating 1 layer cross sections with given regional seepage flux.

        Parameters
        ----------
        dirs: DirStruct obj
            knows the directory structure of the GGOR project and its case
        parcel_data: pd.DataFrame
            the properties of the parcels.
        """
        super().__init__(dirs=dirs, parcel_data=parcel_data)
        self.name = 'L1q'

class L1(Solution):
    """One layer aka Kraaijenhoff vd Leur (Carslaw & Jaeger (1959, p87))."""

    def __init__(self, dirs=None, parcel_data=None):
        """Return a Solution object to analytically simulate transient 1-layer flow.

        Parameters
        ----------
        dirs: DirStruct obj
            knows the directory structure of the GGOR project and its case
        parcel_data: pd.DataFrame
            the properties of the parcels
        """
        super().__init__(dirs=dirs, parcel_data=parcel_data)


class L2(Solution):
    """Return analytic two-layer solution."""

    def __init__(self, dirs=None, parcel_data=None):
        """Return a solution object to analytically simulate transient 2-layer flow.

        Parameters
        ----------
        dirs: DirStruct obj
            knows the directory structure of the GGOR project and its case
        parcel_data: pd.DataFrame
            the properties of the parcels
        """
        super().__init__(dirs=dirs, parcel_data=parcel_data)
        self.name = 'L2'

    def sim(self, tdata=None, use_w_not_c=False):
        """Simulate 2-layer system using multilayer analytical solution.

        tdata: pd.DataFrame
            time data, with columns RH (precip), EV24 (evapotranspiration)
            the index are are datetimes.
        use_w_not_c: bool
            triggers use of presscribed analytical ditch resistance `w` instead
            of computing it from the real ditch resistance, the ditch circum-
            ference and the effect of conctraction of streamlines due to
            partially penetration of the ditch(es)
        """
        self.tdata, self.HDS, self.CBC = multi_layer_transient(
                            tdata=tdata, parcel_data=self.parcel_data,
                            use_w_not_c=use_w_not_c)

class Lnum(Solution):
    """Return numeric solution using MODFLOW."""

    def __init__(self, dirs=None, parcel_data=None):
        """Generate a Lnum object (numerical simulation).

        Parameters
        ----------
        dirs: DirStruct obj
            knows the directory structure of the GGOR project and its case
        parcel_data: pd.DataFrame
            parcel properties
        """
        super().__init__(dirs=dirs, parcel_data=parcel_data)
        self.name = 'modflow'
        self.dirs = dirs

    def sim(self, tdata=None, use_w_not_c=False):
        """Numerically simulate 2-layer system using MODFLOW.

        tdata: pd.DataFrame
            time data, with columns RH (precip), EV24 (evapotranspiration)
            the index are are datetimes.
        use_w_not_c: bool
            triggers use of presscribed analytical ditch resistance `w` instead
            of computing it from the real ditch resistance, the ditch circum-
            ference and the effect of conctraction of streamlines due to
            partially penetration of the ditch(es)
        """
        par, spd, bdd, gr = gt.run_modflow(dirs=self.dirs, parcel_data=self.parcel_data,
                     tdata=tdata, laycbd=(1, 0), dx=1., use_w_not_c=False)
        return par, spd, bdd, gr

#%% __main__

if __name__ == '__main__':
    # home = '/Users/Theo/GRWMODELS/python/GGOR/'

    test=True
    size_inches = (14, 10)
    use_w_not_c = False

    # Parameters to generate the model. Well use this as **kwargs
    GGOR_home = os.path.expanduser('~/GRWMODELS/python/GGOR') # home directory
    case = 'AAN_GZK'

    #GGOR directory structure
    dirs = gt.Dir_struct(GGOR_home, case=case)

    #Get the meteo tdata from an existing file or directly from KNMI
    meteo_data = knmi.get_weather(stn=240, start='20100101', end='20191231',
                                  folder=dirs.meteo)

    # Add columns "summer' and "hyear" to it"
    tdata = gt.handle_meteo_data(meteo_data, summer_start=4, summer_end=10)


    if test:
        # Change the tdata such that we have a smooth input
        tdata = gt.gen_testdata(tdata=tdata,
                              RH  =(270, 0.0, 0.002, 0.000),
                              EV24=(360, 0.0, 0.002, 0.000),
                              )

        parcel_data = gt.get_test_parcels(os.path.join(
                                dirs.case, 'pdata_test.xlsx'), 'parcel_tests1').iloc[:4]
    else:
        # Bofek data, coverting from code to soil properties (kh, kv, sy)
        # The BOFEK column represents a dutch standardized soil type. It is used.
        # Teh corresponding values for 'kh', 'kv' and 'Sy' are currently read from
        # and excel worksheet into a pd.DataFrame (table)
        bofek = pd.read_excel(os.path.join(dirs.bofek, "BOFEK eenheden.xlsx"),
                              sheet_name = 'bofek', index_col=0, engine="openpyxl")

        # Create a GGOR_modflow object and get the upgraded parcel_data from it
        parcel_data = gt.GGOR_data(defaults=gt.defaults, bofek=bofek, BMINMAX=(5, 250),
                                   GGOR_home=GGOR_home, case=case).data

   #%%
    parcel_data = parcel_data.iloc[:4]
    obj = None
    GXG = False
    scenario = 4
    if scenario == 1: # analytic with given head in regional aquifer
        l1f = L1f(dirs=dirs, parcel_data=parcel_data)
        l1f.sim(tdata=tdata, use_w_not_c=use_w_not_c)
        l1f.plot_heads(GXG=GXG)
        l1f.plot_cbc()
        budget = l1f.CBC.get_budget()
        obj = L1f
    if scenario == 2: # analytic with given seepage from regional aquifer
        l1q = L1q(dirs=dirs, parcel_data=parcel_data)
        l1q.sim(tdata=tdata, use_w_not_c=use_w_not_c)
        l1q.plot_heads(GXG=GXG)
        l1q.plot_cbc()
        budget = l1q.CBC.get_budget()
        obj = L1q
    if scenario == 3: # Analytic two layers, with ditches in both aquifers
        l2 = L2(diss=dirs, parcel_data=parcel_data)
        l2.sim(tdata=tdata, use_w_not_c=use_w_not_c)
        l2.plot_heads(GXG=GXG)
        l2.plot_cbc()
        budget = l2.CBC.get_budget()
        obj = L2
    if scenario == 4: # numerical with dichtes in both aquifers
        mf = Lnum(dirs=dirs, parcel_data=parcel_data)
        par, spd, bdd, gr = mf.sim(tdata=tdata, use_w_not_c=use_w_not_c)
        heads = gt.Heads_obj(dirs, tdata=tdata, IBOUND=par['IBOUND'], gr=gr)
        watbal = gt.Watbal_obj(dirs,
                               tdata=tdata,
                               parcel_data=parcel_data,
                               IBOUND=par['IBOUND'],
                               gr=gr)
        selection = [0, 1, 2, 3]
        titles=['Parcel averaged heads', 'Parcel averaged heads']
        ax = heads.plot(tdata=tdata,
               parcel_data=parcel_data,
               selection=selection,
               titles=titles,
               size_inches=(14, 8), GXG=GXG)
        ax = watbal.plot(parcel_data=parcel_data,
                         tdata=tdata,
                         sharey=True)

    print('---- All done ! ----')

    if obj is not None:
        d   = obj.CBC.data
        hds = obj.HDS.data
        ax = newfig2(['test HDS', 'test CBC'], 'tijd', ['m', 'm/d'])
        ax[0].plot(obj.tdata.index, hds[0, 0, :], label='h[L0, ip0]')
        ax[0].plot(obj.tdata.index, hds[1, 0, :], label='h[L1, ip0]')

        ax[0].legend()
        ax[1].plot(obj.tdata.index, d['STO'][0, 0, :], color='cyan', label='STO[L0, ip0]')
        ax[1].plot(obj.tdata.index, d['RCH'][0, 0, :], color='green',  label='RCH[L0, ip0]')
        ax[1].plot(obj.tdata.index, d['EVT'][0, 0, :], color='yellow',  label='EVT[L0, ip0]')
        ax[1].plot(obj.tdata.index, d['FLF'][0, 0, :], color='darkgray', label='FLF[L0, ip0]')
        ax[1].plot(obj.tdata.index, d['DRN'][0, 0, :], color='lightgray', label='DRN[L0, ip0]')
        ax[1].plot(obj.tdata.index, d['GHB'][0, 0, :], color='purple', label='GHB[L0, ip0]')
        ax[1].plot(obj.tdata.index, d['RIV'][0, 0, :], color='magenta', label='RIV[L0, ip0]')
        ax[1].legend()
        plt.show()