# -*- coding: utf-8 -*-

# Note that the PYTHONPATH is set in mflab/.env
# Use VENV flopy by setting it


"""Simultate GGOR for nparcels using MODFLOW.

* The parcels (area) data are in a shape file.
* The attributes are in the dbf
* The .dbd file is read into pd.DataFrame.
* The meteo is read from existing file or obtained from KNMI site.
* Modflow is used to simultaneously simulate all the parcels dynamically.
* The results are shown for selected parcels (hds, GXG)
* The running water budget is shown for all the parcels combined.

Different scenarios can be dealt with as cases. Scenarios are used to simulate
a regular case or for testing the behavior of the model in given test-circumstances.

A scenario for a test has for instance max 5 parcels. The the number of head
time series to be plotted for verification is then also limited to 5.
Hence a DataFrame needs to be generated from the data that define the properties of these parcels.
This DataFrame can be read from an Excel workbook that bears the name of the case, which allows full control over the parcel properties of the test case.

The next data for a case is the meteo. These data too can be read from
an excel sheet from the same workbook. The required fields are in the
example workbook named "test_basic.xlsx".

@ TO 2020-09-06

The GGOR has been converted to mf6 but using the structured grid, which makes best sense
because we simulate each parcel as a single row of cells in the grid.

Below, the data for each modflow module to be adapted from the defaults, i.e. from those
in the Excel workbook are specified below. Before saving the Modflow files and running
Modflow 6, they well be used to update the default values for each module.

@ TO 2025-07-02
"""
# %% --- imports
import os
import sys
import numpy as np
import time
from pathlib import Path
import logging
from mf6_bootstrap import activate # noqa: F401

activate(verbose=True) # noqa: RUF100

# --- Imports are not missing after activate(__file__)
# ruff: noqa: E402

from settings import props      # pyright: ignore[reportMissingImports]      
from timing import log_timed    # pyright: ignore[reportMissingImports]
from mf6tools import Dirs       # pyright: ignore[reportMissingImports]
import ggor_tools as ggt        # pyright: ignore[reportMissingImports]
from logging_setup import configure_logging # pyright: ignore[reportMissingImports]

# --- setting up the logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
configure_logging()

# --- Rest of the script
case_name = Path(__file__).parent.parent.parts[-1]

start_script = time.perf_counter()

logger.info("Running module as a script")

dirs = Dirs()
dirs.meteo = os.path.join(Path(dirs.proj).parent, 'data', 'meteo')
dirs.bofek = os.path.join(Path(dirs.proj).parent, 'data', 'bofek')


# %% --- tdis ======  time discretization
with log_timed(logger, 'tdata generated and pickled'):
    tdata = ggt.get_tdata(dirs=dirs,  stn=240,
                  start='20100101', end='20191231',
                  folder=dirs.meteo)
    tdata_file = os.path.join(dirs.data, 'tdata.pkl')
    tdata.to_pickle(tdata_file)
    
with log_timed(logger, "Parcel_data generated and pickled"):
    parcel_data = ggt.get_parcel_data(dirs=dirs, defaults=ggt.defaults, BMINMAX=(5, 250))
    pdata_file = os.path.join(dirs.data, 'parcel_data.pkl')
    parcel_data.to_pickle(pdata_file)

start_date_time = str(tdata.index[0])

nper, nstep, tsmult = len(tdata), 1, 1.0

# --- Because we change boundary conditions, nper > 0
dt = np.diff(tdata.index - tdata.index[0]) / np.timedelta64(1, 'D')
dt = np.hstack((dt[0], dt)) # Assume dt[0] ==  t[1] - t[0])

period_data = [[sp_time, nstep, tsmult] for sp_time in dt]

Simtdis = {'perioddata': period_data,
           'nper': len(period_data),
           'start_date_time': start_date_time,
           'time_units': props['time_units'],
           }

with log_timed(logger, "grid generated"):
    gr = ggt.grid_from_parcel_data(parcel_data=parcel_data,
                                   dx=props['dx'])


# %% --- Gwfdis ======
IDOMAIN = gr.const(1, dtype=int)
for iy, b in zip(range(gr.ny), parcel_data['b']):
    IDOMAIN[:, iy, gr.Xm[iy] > b] = 0

IDOMAIN[gr.DZ < props['minDz']] = -1 # Don't need this in GGOR

Gwfdis = {'gr': gr,
     'idomain': IDOMAIN,
     'length_units': props['length_units']}

# %% --- Gwfsto ======= Storage coefficients

S = ggt.set3D(parcel_data[['sy',  'S2']], gr.shape)
Sy = ggt.set3D(parcel_data[['sy', 'sy']], gr.shape)

Gwfsto = {    
    'ss':  S,
    'sy': Sy,
    'iconvert': props['icelltype'],
    'storagecoefficient': True, # Interpret as S instead of Ss where
    }
    
# %% --- Gwfnpf ===== Cell properties
parcel_data['kc'] = parcel_data['D_CB'] / parcel_data['c_CB']

Gwfnpf = {
    'k': ggt.set3D(parcel_data[['kh', 'kh2']], gr.shape),
    'k22': 1.e-20,                  # no flow along y=axis
    'k33': ggt.set3D(parcel_data[['kv', 'kv2']], gr.shape),
    'icelltype': props['icelltype'],
}

# %% --- Gwfic ===== Initial heads
Gwfic = {
    'strt': ggt.set3D(parcel_data['h_winter'], gr.shape)
}

# %% --- Gwfrcha ===== Recharge   
rch_spd = {isp: tdata['RH'].iloc[isp] for isp in range(len(tdata))}

Gwfrcha = {  
    'recharge': rch_spd,
    'readasarrays': True,
    'print_input': True,    
}

# %% --- Gwfevta ================================
with log_timed(logger, "Gwfevta gnerated"):
    Gwfevta = {
        'readasarrays': True,
        'ievt': None,
        'surface': {0: ggt.set3D(parcel_data['AHN'] - parcel_data['ET_surfd'], shape=gr.shape)[0]},
        'depth': {0: ggt.set3D(parcel_data['ET_exdp'], shape=gr.shape)[0]},
        'rate': {isp: tdata['EV24'].iloc[isp] for isp in range(len(tdata))}
    }

# %% --- Prepare boundary conditions ==========
mon = np.array([dt.month for dt in tdata.index])
day = np.array([dt.day for dt in tdata.index])

Isp_start_summer = np.where(np.logical_and(mon ==  4, day == 1))[0]
Isp_start_winter = np.where(np.logical_and(mon == 10, day == 1))[0]

# %% --- wel ===== used to model given seepage from regional aquifer
active = IDOMAIN > 0
Iwel = gr.NOD[-1][active[-1]]

Q = (parcel_data['q_up'].values[:, np.newaxis] * gr.Area) [np.newaxis, : ,:] * gr.const(1)
Q[:-1, : ,:] = 0.0

spd = ([(lrc, Qw) for lrc, Qw in 
            zip(gr.lrc_from_iglob(Iwel,  astuples=True), Q.ravel()[Iwel])])

# Only the first SP needs data, as long as seepage is constant.
# TODO: this likely changes in the future, to monthly seepage values.
# TODO: Specify how monthly seepage values for all parcels are inported.
# TODO: This import may require an extended data file originating from a regional model.
# TODO: its dimensions should be (nparcel x months in time series)
stress_period_data = {0: spd}

Gwfwel = {
    'stress_period_data': stress_period_data,
    'maxbound': len(Iwel)
}    

# %% --- Gwfghb ===== Is used for flow from and toward ditches
    
# --- First cell of top and bottom layer (always). Conduction of bottom layer may be zero.
Ighb = gr.NOD[[0, -1], :, 0].flatten()

LRC_ghb = gr.lrc_from_iglob(Ighb, astuples=True)

# --- Get conductance for the GHB (connection soil Ditch in top and bottom layer)
# --- i.e. one value per parcel in the top layer and one value in the bottom layer.

# --- Use analytic ditch resistance with layer thickness and no partial penetration
condGHB = ggt.get_cond_GHB(pdata=parcel_data, gr=gr)

# --- Get the GHB heads for summer and winter for these ditches
hw = np.vstack((parcel_data['h_winter'], parcel_data['h_winter']))
hs = np.vstack((parcel_data['h_summer'], parcel_data['h_summer']))

ghb_winter = [(lrc, head, cond) for lrc, head, cond in
                        zip(LRC_ghb, hw.ravel(), condGHB.ravel())]
ghb_summer = [(lrc, head, cond) for lrc, head, cond in
                        zip(LRC_ghb, hs.ravel(), condGHB.ravel())]

# --- Input for the first stress period
stress_period_data = {0: ghb_summer if tdata.iloc[0]['summer'] is True else ghb_winter}

# --- Only generate input when summer changes to winter or vice versa
for isp in Isp_start_summer:
    stress_period_data[isp] = ghb_summer
for isp in Isp_start_winter:
    stress_period_data[isp] = ghb_winter

Gwfghb = {
    'stress_period_data': stress_period_data,
    'maxbound': len(LRC_ghb),
}

# %% --- Gwfriv ==============Is used for extra flow toward ditch (lower flow to than from)
# Riv cells are the same as GHB cells
LRC_riv = LRC_ghb

condRIV = ggt.get_RIV_Cond(pdata=parcel_data, gr=gr)

riv_winter = [(lic, stage, cond, rbot) for lic, stage, cond, rbot in
                            zip(LRC_riv, hw.ravel(), condRIV.ravel(), hw.ravel())]
riv_summer = [(lic, stage, cond, rbot) for lic, stage, cond, rbot in
                        zip(LRC_riv, hs.ravel(), condRIV.ravel(), hs.ravel())]

# --- Input for the first stress period
stress_period_data = {0: riv_summer if tdata.iloc[0]['summer'] is True else riv_winter}

# --- Only generates input when summer changes to winter and vice versa.
for isp in Isp_start_summer:
    stress_period_data[isp] = riv_summer
for isp in Isp_start_winter:
    stress_period_data[isp] = riv_winter

Gwfriv = {
    'stress_period_data': stress_period_data,
    'maxbound': len(LRC_riv),
}

# %% --- Gwfdrn ============== Drn is used for drains, surface runoff and trenches
Idrn = gr.NOD[0, :, 1:][active[0, :, 1:]]
LRC_drn = gr.lrc_from_iglob(Idrn, astuples=True)

# --- Get drain elevation and drain conductance as Ny * Nx array (top layer)
elevation = ggt.get_drain_elev_with_trenches(pdata=parcel_data, gr=gr, d_drn=props['drain_depth']).ravel()[Idrn]
condDRN   = ggt.get_cond_DRN(pdata=parcel_data, gr=gr).ravel()[Idrn]

# --- Stress period data
#dtype = flopy.modflow.ModflowDrn.get_default_dtype()
#spd   = np.recarray(gr.nod, dtype=dtype)

spd = [(lrc, elev, cond) for lrc, elev, cond in zip(LRC_drn, elevation, condDRN)]

# --- Input is required only for the first SP, because drain data are constant.
stress_period_data = {0: spd}

Gwfdrn = {
    'stress_period_data': stress_period_data,
    'maxbound': len(LRC_drn),
}

# %% --- Gwfoc ==== Output control for flow model
Gwfoc = {'head_filerecord':   os.path.join(dirs.SIM, "{}Gwf.hds".format(case_name)),
         'budget_filerecord': os.path.join(dirs.SIM, "{}Gwf.cbc".format(case_name)),
         'saverecord': [("HEAD", "FREQUENCY", props['oc_frequency']),
                        ("BUDGET", "FREQUENCY", props['oc_frequency'])],
}

logger.info(f"mf_adapt finished in {time.perf_counter() - start_script:.2f} seconds")

# --- Pickling the parcel_data geopandas.GeoDataFrame
pdata_pkl = os.path.join(dirs.data, 'parcel_data.pkl')
parcel_data.to_pickle(pdata_pkl)
logger.info(f"parcel_data pickled to {pdata_pkl}")

# --- Pickling the tdata pd.DataFrame
tdata_pkl = os.path.join(dirs.data, 'tdata.pkl')
logger.info(f"tdata pickled to {tdata_pkl}")
tdata.to_pickle(tdata_pkl)

if __name__ == '__main__':
    print('---- All done mf_adapt ! ----')
