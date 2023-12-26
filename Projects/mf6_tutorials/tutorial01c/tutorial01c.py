#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 15 16:51:50 2021

First Flopy mdf6 tutorial form the flopy documentation

This is the same probleme as tutorial01, with but the chd with the model
interior is replaced by a well with a time series as extraction input and also
the recharge has been defined as a timearrayseries.Further the adding of packages
to the model has been separated from the specificiation of the input data.
This implies that the adding of packages now is essentially independent of
the specificiation of the input data. More packages may have to be added
in a similar manner.

See for the examples on how to use the ts (time series) and tas (time array series)

https://github.com/modflowpy/flopy/blob/develop/examples/Notebooks/flopy3_mf6_obs_ts_tas.ipynb

@author: Theo 2021-22-16
"""

import numpy as np
import matplotlib.pyplot as plt
import flopy
#import inspect
import os
import sys

tools = os.path.abspath('~/GRWMODELS/python3/tools')
if not tools in sys.path:
    sys.path.insert(0, tools)

from fdm.mfgrid import Grid

sys.path.append("../common")

attribs = lambda obj: [o for o in dir(obj) if not o.startswith('_')]


import etc

#%%
tutorials = '/Users/Theo/GRWMODELS/python/mf6_tutorials/'
exe_dir   = '/Users/Theo/GRWMODELS/mfLab/trunk/bin/'

#%% Simulation
name      = 'tutorial01b'
sim_name  = name
sim_ws    = os.path.join(tutorials, name)
exe_path  = os.path.join(exe_dir, 'mf6.mac')
version   = 'mf6'

if not os.path.isdir(sim_ws):
    os.mkdir(sim_ws)

os.chdir(sim_ws)

#%% Tdis, time dicretization
time_units='Days'
start_date_time = '2021-05-22 18:00'
perioddata = [(20., 100., 1.)]
nper = len(perioddata)

transient=[False]

#%%$ Modflow grid
N, Nlay, H = 100, 1, 50

x = np.linspace(0, 1000, N + 1)
y = np.linspace(1000, 0, N + 1)
z = np.linspace(0, -H, Nlay + 1)

gr = Grid(x, y, z)

#%% Storage package
ss = 1.0e-5
sy = 0.15

#%% Layer properties
k = 600.
icelltype = 1

#%% Initial conditions
h1 = 0.
strt = gr.const(h1)

#%% Constant head at outer ring of cells in firsta nd last layers
chd_rec = []
for layer in [0]:
    for col in range(N):
        chd_rec.append(((layer,    0, col), h1))
        chd_rec.append(((layer, N -1, col), h1))
    for row in range(1, N - 1):
        chd_rec.append(((layer, row,    0), h1))
        chd_rec.append(((layer, row, N -1), h1))


#%% Well
wel_rec = []

# Using the time series "flow" in the current position
wel_rec.append(((0, int(N / 4), int(N / 4)), 'Qwell'))

# Time series
well_flow = [(0, -600), (5, -200), (7, -400), (9, -300), (10, 0), (100, 0)]

# Define the time series in a dict
ts_dict = {'timeseries':well_flow,
           'time_series_namerecord': 'Qwell',
           'interpolation_methodrecord':'stepwise',
           'filename':'Q_well.ts',
           'sfacrecord': 1.0}


#%% Recharge as a tie series array
tas = {0: 0.002, 1: -0.002, 2: 0.040, 3: 0.001, 4: 0., 5: 0.002,
       6: 0.0, 7: -0.002, 8: -0.004, 20: 0.001}

tas_name             = 'KNMI_Rotterdam'
tas_fname            = tas_name + '.tas'
interpolation_method = 'stepwise'


#%% Oc output control
headfile    = "{}.hds".format(name)
budgetfile  = "{}.cbb".format(name)
saverecord  = [("HEAD", "ALL"), ("BUDGET", "ALL")]
printrecord =("HEAD", "LAST"),


#%% Adding the packages this should normally stay untouched for any model
# Except for which packagres to add. All the data are specified above.
sim=flopy.mf6.MFSimulation(sim_name=name, exe_name=exe_path, version=version)

tdis = flopy.mf6.ModflowTdis(sim, time_units=time_units,
            start_date_time=start_date_time, nper=nper, perioddata=perioddata)

# Solver
ims = flopy.mf6.ModflowIms(sim, pname='ims', complexity="SIMPLE")

# Add the groundwater flow model
gwf = flopy.mf6.ModflowGwf(sim, modelname=name, model_nam_file="{}.nam".format(name),
    version=version, exe_name=exe_path, newtonoptions=None,)

# Dis package
dis = flopy.mf6.ModflowGwfdis(gwf,
            nlay=gr.nlay, nrow=gr.nrow, ncol=gr.ncol, delr=gr.dx, delc=gr.dy,
                               top=gr.z[0], botm=gr.z[1:])

# Storage
sto = flopy.mf6.ModflowGwfsto(gwf, storagecoefficient=True, iconvert=1,
                              ss=ss, sy=sy, transient=transient,
                              save_flows=True
                              )
# Initial conditions
ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

# Layer rproperties
npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=icelltype, k=k, save_flows=True)

# Ghb
chd = flopy.mf6.ModflowGwfchd(gwf,
            maxbound=len(chd_rec), stress_period_data=chd_rec, save_flows=True)

# Well
wel = flopy.mf6.ModflowGwfwel(gwf,
            maxbound=len(wel_rec), timeseries=ts_dict,
            stress_period_data=wel_rec, save_flows=True)

# recharge
tas.update(filename=tas_fname,
           time_series_namerecord=tas_name,
           interpolation_methodrecord=interpolation_method)

rcha = flopy.mf6.ModflowGwfrcha(gwf, save_flows=True, readasarrays=True,
                                timearrayseries=tas,
                                recharge='timearrayseries ' + tas_name)

# Output control ('OC') package
oc = flopy.mf6.ModflowGwfoc(gwf,
        saverecord        = saverecord,
        head_filerecord   =[headfile],
        budget_filerecord =[budgetfile],
        printrecord       = printrecord,
    )

# Write simulation files and run
# This should not have to be touched.
sim.write_simulation()

# Run
success, buff = sim.run_simulation()
if not success:
    raise Exception("MODFLOW 6 did not terminate normally.")


#%% Post-process head results
hds = flopy.utils.binaryfile.HeadFile(headfile)

kstpkper=(10, 0)

h = hds.get_data(kstpkper=kstpkper)
ax = etc.newfig("Head contours at kstpkper = {}".format(str(kstpkper)),
                "x [m]", "y [m]")
c = ax.contour(gr.xm, gr.ym, h[0], np.arange(-1, 0, 0.01), linestyles='-', colors="black")
ax.clabel(c, fmt="%.3f")

#%%
times = hds.get_times()
hall = hds.get_alldata()[:, 0, 25, 25]

ax = etc.newfig("Head vs time at cell {}".format(str(wel_rec[0][0])),
                "t [days]", "head [m]")
ax.plot(times, hall, label="tseries at {}".format(str(wel_rec[0][0])))
ax.legend()
plt.show()

