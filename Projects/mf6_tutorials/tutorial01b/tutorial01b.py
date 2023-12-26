#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 15 16:51:50 2021

First Flopy mdf6 tutorial form the flopy documentation

This is the same probleme as tutorial01, with but the chd with the model
interior is replaced by a well with a time series as extraction input.

We further introduce a time-series array for the reharge

The variable well extraction is computed with one stress period and the
extraction of the well as a time series.

Note that the time series must cover the entire simulation period in
order to be interpolated.

See for the examples on how to use the ts (time series) and tas (time array series)

https://github.com/modflowpy/flopy/blob/develop/examples/Notebooks/flopy3_mf6_obs_ts_tas.ipynb

@author: Theo 2021-05-16
"""

import numpy as np
import matplotlib.pyplot as plt
import flopy

import os
import sys

tools = os.path.abspath('~/GRWMODELS/python3/tools')
if not tools in sys.path:
    sys.path.insert(0, tools)

from fdm.mfgrid import Grid

sys.path.append("../common")

attribs = lambda obj: [o for o in dir(obj) if not o.startswith('_')]


import etc

#%% Charactersiti s of hte model

name      = 'tutorial01b'
tutorials = '/Users/Theo/GRWMODELS/python/mf6_tutorials/'
ws        = os.path.join(tutorials, name)
exe_dir   = '/Users/Theo/GRWMODELS/mfLab/trunk/bin/'
exe_name  = os.path.join(exe_dir, 'mf6.mac')

if not os.path.isdir(ws):
    os.mkdir(ws)

os.chdir(ws)

model_nam_file = "{}.nam".format(name)

#%%
N, Nlay, H = 100, 1, 50

x = np.linspace(0, 1000, N + 1)
y = np.linspace(1000, 0, N + 1)
z = np.linspace(0, -H, Nlay + 1)

gr = Grid(x, y, z)


h1, Qw, k = 0., -600., 25.0

well_flow = [(0, -600), (5, -200), (7, -400), (9, -300), (10, 0), (100, 0)]
start_date_time = '2021-05-16 12:30'
perioddata = [(20., 100., 1.)]

#%% Add simulation object
sim=flopy.mf6.MFSimulation(sim_name=name, exe_name=exe_name, version='mf6')

#%% packages added to the sim object

tdis = flopy.mf6.ModflowTdis(sim, pname='tdis', time_units='Days',
            start_date_time=start_date_time, nper=1, perioddata=perioddata)

ims = flopy.mf6.ModflowIms(sim,
                           pname='ims', complexity="SIMPLE")

gwf = flopy.mf6.ModflowGwf(sim,
                           modelname=name, model_nam_file=model_nam_file)

#%% fwf model object
dis = flopy.mf6.ModflowGwfdis(gwf,
            nlay=gr.nlay, nrow=gr.nrow, ncol=gr.ncol, delr=gr.dx, delc=gr.dy,
                               top=gr.z[0], botm=gr.z[1:])

sto = flopy.mf6.ModflowGwfsto(gwf, storagecoefficient=True, iconvert=1,
                              ss=1.0e-5, sy=0.15, transient=[False],
                              save_flows=True
                              )
#%% Initial conditions
ic = flopy.mf6.ModflowGwfic(gwf,
                            pname="ic", strt=gr.const(h1))

#%% Laye rproperties
npf = flopy.mf6.ModflowGwfnpf(gwf,
                              icelltype=1, k=k, save_flows=True)

#%% Stress packages

# Constant head at outer ring of cells in firsta nd last layers
chd_rec = []
for layer in [0]:
    for col in range(N):
        chd_rec.append(((layer,    0, col), h1))
        chd_rec.append(((layer, N -1, col), h1))
    for row in range(1, N - 1):
        chd_rec.append(((layer, row,    0), h1))
        chd_rec.append(((layer, row, N -1), h1))

chd = flopy.mf6.ModflowGwfchd(gwf,
            maxbound=len(chd_rec), stress_period_data=chd_rec, save_flows=True)


#%% Well
wel_rec = []

# The name flow indicates that a tiem series will be used for its input.
# The time series specifies the name in its time_series_namerecord (see next)
wel_rec.append(((0, int(N / 4), int(N / 4)), 'flow'))

wel = flopy.mf6.ModflowGwfwel(gwf,
            maxbound=len(wel_rec), timeseries=None,
            stress_period_data=wel_rec, save_flows=True)

# The time series can be introuced by specifying the attributes for the
# package ".ts" object.
wel.ts.initialize(filename='well_flow.ts',
                  timeseries=well_flow,
                  time_series_namerecord='flow',
                  interpolation_methodrecord='stepwise',
                  sfacrecord=1.0)

#%% recharge

tsa = {0: 0.002, 1: -0.002, 2: 0.040, 3: 0.001, 4: 0., 5: 0.002,
       6: 0.0, 7: -0.002, 8: -0.004, 20: 0.001}

# Three ways to use the timeseriesarray
method = 3

if method == 3:
    #method1 pass the tas array to package constructor, no tas file is generated
    rcha = flopy.mf6.ModflowGwfrcha(gwf, save_flows=True, readasarrays=True,
              timearrayseries=tsa, recharge='timearrayseries KNMI_Rotterdam')

    rcha.tas.time_series_namerecord='KNMI_Rotterdam'
    rcha.tas.interpolation_methodrecord='stepwise'
elif method == 2:
    # Initialize time array series through rcha.tas.initialize, a file is geneated
    rcha = flopy.mf6.ModflowGwfrcha(gwf, save_flows=True, readasarrays=True,
                                   recharge='timearrayseries KNMI_Rotterdam')

    rcha.tas.initialize(filename='KNMI_Rotterdam.tas',
                       tas_array=tsa,
                       time_series_namerecord ='KNMI_Rotterdam',
                       interpolation_methodrecord = 'stepwise',
                       )
elif method == 3:
    tas = tsa.copy()
    tas.update(filename='KNMI_Rotterdam.tas',
               time_series_namerecord ='KNMI_Rotterdam',
               interpolation_methodrecord = 'stepwise',
              )

    # Pass timearrayseries a dictionary of anything that could be passed to rcha.tas.initialize
    rcha = flopy.mf6.ModflowGwfrcha(gwf, save_flows=True, readasarrays=True,
                                    timearrayseries=tas,
                                    recharge='timearrayseries KNMI_Rotterdam')

#%% Output control ('OC') package
headfile   = "{}.hds".format(name)
budgetfile = "{}.cbb".format(name)

oc = flopy.mf6.ModflowGwfoc(gwf,
        saverecord        =[("HEAD", "ALL"), ("BUDGET", "ALL")],
        head_filerecord   =[headfile],
        budget_filerecord =[budgetfile],
        printrecord       =[("HEAD", "LAST")],
    )

#%%
sim.write_simulation()

#%%
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

