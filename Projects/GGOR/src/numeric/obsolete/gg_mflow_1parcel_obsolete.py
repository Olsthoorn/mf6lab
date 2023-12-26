#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Implement modflow simulation of one GGOR cross section.

This file generates a GGOR model voor a single parcel using the same data
and input as the analytical model in ggor_analytical.py. This is to compare
the results of the analytical and the numerical model.

Hence, the data setting the properties of the analytical model and its meteo is
used to set up the numerical model

TODO: as per 170825
    verify with analytical solution
    implement greppels
    implement ditch shape
    document
@TO 2020-06-03
"""

# IMPORTS
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from KNMI import knmi
import flopy
import flopy.modflow as fm
import pandas as pd
from fdm.mfgrid import Grid
#import pdb

import GGOR.src.numeric.gg_mflow as gg

print('sys version {}'.format(sys.version))
print('numpy version {}'.format(np.__version__))
print('flopy version {}'.format(flopy.__version__))

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


# functions
def set3D(layvals, dims=None):
    """Return 3D array filled with data stored in self.data according to what.

    Parameters
    ----------
        layvas: sequence
            The layer values to be used
        dims: sequence
            the dimensions of the 3D array to be generated
    """
    return np.asarray(layvals)[:, np.newaxis, np.newaxis] * np.ones(tuple(dims[1:]))


def setBoundary(what=None, gr=None, props=None, tdata=None):
    """Return dictionary for boundary of given type.

    Parameters
    ----------
    what: str, one of
         'WEL': WEL package used to simulate vertical seepage.
         'DRN': DRN package used to simulate tile drainage (or surface runoff.
         'GHB': GHB package used to simulate ditch in- and outflow.
         'RIV': RIV package used to simulate ditch outflow (together with GHB).
    gr : fdm_tools.mfgrid.Grid
        grid object
    props: dict with stationary tdata
        physical properties of this parcel
    tdata : pd.DataFrame with time tdata in columns 'HR', 'EV24', 'hLR', 'summer'
    """
    boundary_dict = {}
    if what=='WEL': # in all of layer 2 (bottom aquifer.
        # Get the indices of all cells within these layers
        I = gr.NOD[-1, :, :].ravel()
        A = gr.DX * gr.DY # their cell area
        lrc = np.array(gr.I2LRC(I))  # their lrc index tuples
        Q_up = A.ravel()[I] * props['q_up']

        # Gnerate basic recarray prototype to store the large amount of tdata efficiently
        dtype = flopy.modflow.ModflowWel.get_default_dtype()
        spd = np.recarray(len(I), dtype=dtype)

        spd['k'] = lrc[:, 0]
        spd['i'] = lrc[:, 1]
        spd['j'] = lrc[:, 2]

        # fill the trainsient tdata
        for isp in [0]: # only first stress period, rest is the same.
            spd['flux'] = Q_up
            boundary_dict[isp] = spd.copy()

    elif what == 'DRN':
        I = gr.NOD[0, :, 1:-1].ravel() # drains in first layer
        lrc  = np.asarray(gr.I2LRC(I), dtype=int)
        cond = gr.DX.ravel()[I] * gr.DY.ravel()[I] / props['c_drain']

        dtype = flopy.modflow.ModflowDrn.get_default_dtype()
        spd = np.recarray(len(I), dtype=dtype)

        spd['k'] = lrc[:, 0]
        spd['i'] = lrc[:, 1]
        spd['j'] = lrc[:, 2]
        spd['elev'] = props['AHN'] - props['d_drain']
        spd['cond'] = cond

        isp = 0
        boundary_dict  = {isp: spd} # ionly first isp, rest is the same.

    elif what == 'GHB': # entry resistance (the highest resistance)
        w = props['wi_ditch']
        I = gr.NOD[:, :, 0].ravel()
        lrc  = np.array(gr.I2LRC(I), dtype=int)

        cond = (gr.DZ.ravel()[I] * gr.DY.ravel()[I]).ravel() / w

        dtype = flopy.modflow.ModflowGhb.get_default_dtype()
        spd = np.recarray(len(I), dtype=dtype)

        spd['k'] = lrc[:, 0]
        spd['i'] = lrc[:, 1]
        spd['j'] = lrc[:, 2]
        spd['cond'] = cond
        for isp, hlr in enumerate(tdata['hLR']):
            spd['bhead'] = hlr
            boundary_dict[isp] = spd.copy()
    elif what == 'RIV': # aprallel exit resistance.
        # set ditch bottom elevation equal to ditch level at each time step.
        I = gr.NOD[:, :, 0].ravel()
        lrc  = np.asarray(gr.I2LRC(I), dtype=int)
        wi, wo = props['wi_ditch'], props['wo_ditch']

        assert wi >= wo, "Ditch entry resist. must be >= ditch entry resistance!"

        w = wo * wi / (wi - wo) if wi > wo else np.inf

        cond = gr.DZ.ravel()[I] * gr.DY.ravel()[I] / w

        dtype = flopy.modflow.ModflowRiv.get_default_dtype()
        spd = np.recarray(len(I), dtype=dtype)

        spd['k'] = lrc[:, 0]
        spd['i'] = lrc[:, 1]
        spd['j'] = lrc[:, 2]
        spd['cond'] = cond
        for isp, hlr in enumerate(tdata['hLR']):
            spd['stage'] = hlr
            spd['rbot']  = hlr
            boundary_dict[isp] = spd.copy()
    return boundary_dict


# THIS MODEL
def modflow(props=None, dx=1.0, tdata=None):
    """Set up and run a modflow model simulating a cross section over time.

    Modflow is setup and run and the results are added as section averaged
    heads and water budget components in columns of the output DataFrame
    with the same index as the input and the results as added columns.
    The input DataFrame remains intact.

    Parameters
    ----------
    props: dict
        the properties that defined the situation of the model
    dx: float
        cell width to be used in the model grid
    tdata: pd.DataFrame
        the time-varying tdata.

    Returns
    -------
    pd.DataFrame copy if input with added columns.

    @TO 2020-07-21
    """
    tdata = tdata.copy() # Don't overwrite input tdata
    tdata['hLR'] = props['h_winter']
    tdata['hLR'].loc[tdata['summer']] = props['h_winter']
    tdata['q_up'] = props['q_up']

    modelname  = 'GGOR_1parcel'

    print("FLOPY MODFLOW model: <<{}>>".format(modelname))

    home = '/Users/Theo/GRWMODELS/python/GGOR/'

    if os.name == 'posix':
        executable = os.path.join(home, 'bin/mfusg.mac')
    else:
        executable = os.path.join(home, 'bin/mfusg_64.exe')

    exe_name   = flopy.mbase.which(executable)

    # MODEL DOMAIN AND GRID DEFINITION, x runs from -b to 0 (left size of cross section)
    xGr = np.hstack((-props['b'], np.arange(-props['b'], 0, dx), 0.))
    yGr = [-0.5, 0.5]
    zGr = props['AHN'] - np.cumsum(
                    np.array([0., props['D1'], props['D_CB'], props['D2']]))

    #LAYCBD = np.ones(len(zGr) // 2 - 1)
    LAYCBD = np.array([1, 0])
    gr = Grid(xGr, yGr, zGr, LAYCBD=LAYCBD)

    # MODEL DATA and PARAMETER VALUES
    # Varying ditch level
    kh = np.array([props['kh'], props['kh2']])
    kv = np.array([props['kv'], props['kv2']])
    sy = np.array([props['sy'], props[ 'S2']])
    S  = np.array([props['sy'], props[ 'S2']])

    hstrt = props['h_summer'] if tdata['summer'].iloc[0] else props['h_winter']

    # A string looks up tdata in the database
    # a number uses this number for all parcels
    HK     = gr.const(kh, lay=True)
    VKA    = gr.const(kv, lay=True)
    SY     = gr.const(sy, lay=True)
    SS     = gr.const( S, lay=True)

    # Vertical hydraulic conductivity of aquitards
    VKCB   = gr.const( (gr.dz[gr.Icbd] / np.array(props['c_CB'])), cbd=True)

    # Include layer number in IBOUND (actually layer number + 1)
    IBOUND = gr.const(1, lay=True)

    STRTHD = gr.const(hstrt * np.ones((gr.nlay, 1)), lay=True)

    # All layers to confined, with constant D and kD
    LAYTYP = np.zeros(gr.nlay, dtype=int)

    # If so, then make sure SS complies with SY
    if LAYTYP[0] == 0: SS[0] = SY[0] / gr.DZ[0]

    # STRESS PERIOD DATA
    NPER   = len(tdata)
    PERLEN = np.diff(tdata.index) / np.timedelta64(1, 'D') # time step lengths
    PERLEN = np.hstack((PERLEN[0], PERLEN))
    NSTP   = np.ones(NPER, dtype=int)
    STEADY = np.ones(NPER, dtype=int) == 0 # all transient

    RECH = {isp: tdata['RH'].values[isp] for isp in range(NPER)}   # Recharge
    EVTR = {isp: tdata['EV24'].values[isp] for isp in range(NPER)} # Evapotranspiration

    # Boudnaries, no fixed heads.
    GHB  = setBoundary(what='GHB', gr=gr, props=props, tdata=tdata)
    RIV  = setBoundary(what='RIV', gr=gr, props=props, tdata=tdata)
    DRN  = setBoundary(what='DRN', gr=gr, props=props, tdata=tdata)
    KWEL = setBoundary(what='WEL', gr=gr, props=props, tdata=tdata)

    # What to save on the route?
    OC   = {(isp, istp-1): ['save head', 'save budget', 'print budget'] for isp, istp in zip(range(NPER), NSTP)}

    # MODEL AND packages ADDED TO IT
    mf  = fm.Modflow(modelname, exe_name=exe_name)

    dis = fm.ModflowDis(mf, gr.nlay, gr.ny, gr.nx,
                        delr=gr.dx, delc=gr.dy, top=gr.Ztop[0], botm=gr.Zbot,
                        laycbd=list(gr.LAYCBD),
                        nper=NPER, perlen=PERLEN, nstp=NSTP, steady=STEADY)
    bas = fm.ModflowBas(mf, ibound=IBOUND, strt=STRTHD)
    lpf = fm.ModflowLpf(mf, hk=HK, vka=VKA, chani=np.ones(gr.nlay) * 1e-20, sy=SY, ss=SS,
                        laytyp=LAYTYP, vkcb=VKCB, ipakcb=53)
    ghb = fm.ModflowGhb(mf, stress_period_data=GHB, ipakcb=53)
    riv = fm.ModflowRiv(mf, stress_period_data=RIV, ipakcb=53)
    drn = fm.ModflowDrn(mf, stress_period_data=DRN, ipakcb=53)
    wel = fm.ModflowWel(mf, stress_period_data=KWEL, ipakcb=53)
    rch = fm.ModflowRch(mf, nrchop=3, rech=RECH, ipakcb=53)
    evt = fm.ModflowEvt(mf, nevtop=3, evtr=EVTR, ipakcb=53)
    oc  = fm.ModflowOc( mf, stress_period_data=OC, compact=True)
    #pcg = fm.ModflowPcg(mf, mxiter=200, iter1=200, hclose=0.001, rclose=0.001)
    sms = fm.ModflowSms(mf) #, mxiter=200, iter1=200, hclose=0.001, rclose=0.001)

    # This is to prevent irritating pyflake warning for not using these items
    for chk in [dis, bas, lpf, ghb, riv, drn, wel, rch, evt, oc, sms]:
        if chk is None:
            print('chk reporting one of the modflow items is None')

    # Write the model input files and run MODFLOW
    mf.write_input()
    success, mfoutput = mf.run_model(silent=False, pause=False)

    print('Running success = {}'.format(success))
    if not success:
        raise Exception('MODFLOW did not terminate normally.')

#%% __main__
if __name__ == '__main__':

    test=False

    # Parameters to generate the model. Well use this as **kwargs
    GGOR_home = os.path.expanduser('~/GRWMODELS/python/GGOR') # home directory
    case = 'AAN_GZK'

    #GGOR directory structure
    dirs = gg.Dir_struct(GGOR_home, case=case)

    #Get the meteo tdata from an existing file or directly from KNMI
    meteo_data = knmi.get_weather(stn=240, start='20100101', end='20191231')

    # Add columns "summer' and "hyear" to it"
    tdata = gg.handle_meteo_data(meteo_data, summer_start=4, summer_end=10)

    if test:
        parcel_data = gg.get_test_parcels(os.path.join(
                                dirs.case, 'pdata_test.xlsx'), 'parcel_tests')
    else:
        # Bofek data, coverting from code to soil properties (kh, kv, sy)
        # The BOFEK column represents a dutch standardized soil type. It is used.
        # Teh corresponding values for 'kh', 'kv' and 'Sy' are currently read from
        # and excel worksheet into a pd.DataFrame (table)
        bofek = pd.read_excel(os.path.join(dirs.bofek, "BOFEK eenheden.xlsx"),
                              sheet_name = 'bofek', index_col=0, engine="openpyxl")

        # Create a GGOR_modflow object and get the upgraded parcel_data from it
        parcel_data = gg.GGOR_data(defaults=gg.defaults, bofek=bofek, BMINMAX=(5, 250),
                                   GGOR_home=GGOR_home, case=case).data


    # MODFLOW input arrays are int the returned dicts
    par, spd, bdd, gr =  gg.run_modflow(dirs=dirs, parcel_data=parcel_data, time_data=meteo_data)

    #%% Get the modflow-computed heads and cell by cell flows

    heads = gg.Heads_obj(dirs, IBOUND=par['IBOUND'], gr=gr)

    watbal = gg.Watbal_obj(dirs,
                           IBOUND=par['IBOUND'],
                           parcel_data=parcel_data,
                           time_data=meteo_data,
                           gr=gr)

    #%% Open HDS file and plot heads (with GXG)
    if test:
        for tst in set(parcel_data['Test'].values):
            selection = list(parcel_data.index[parcel_data['Test'] == tst])
            test_vals_str = '{}'.format(', '.join(
                [str(tv) for tv in parcel_data[tst].iloc[selection]]))

            titles=['Parcel averaged heads, testvar {} in [{}]'.format(tst, test_vals_str),
                    'Parcel averaged heads, testvar {} in [{}]'.format(tst, test_vals_str)]
            ax = heads.plot(time_data=meteo_data,
                           parcel_data=parcel_data,
                           selection=selection,
                           titles=titles,
                           size_inches=(14, 8))
            if False:
                ax = watbal.plot(parcel_data=parcel_data,
                                 time_data=meteo_data,
                                 selection=selection[0],
                                 sharey=True)
    else:
        selection = [0]
        ax = heads.plot(time_data=meteo_data,
               parcel_data=parcel_data,
               selection=selection,
               titles=titles,
               size_inches=(14, 8))

        ax = watbal.plot(parcel_data=parcel_data,
                         time_data=meteo_data,
                         selection=None,   # over all parcels
                         sharey=True)

#%%
    print('---- All done ! ----')
#%%
    ax = watbal.plot(parcel_data=parcel_data,
                         time_data=meteo_data,
                         selection=None,   # over all parcels
                         sharey=True)

