# %% Analyzing output of Pennik Series 2 tests
# TO 090223
# %%
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import flopy
import mf_adapt
from mf_adapt import gr
from mf6lab.mf6contourtools import get_contour_levels
from fdm.mf6_face_flows import get_struct_flows # (flowjas, grb_file=None, verbose=False)
import pandas as pd

dirs = mf_adapt.dirs
grb_file = os.path.join(dirs.GWF, mf_adapt.sim_name + 'Gwf.dis.grb')

def addPatch(pgcoords, alpha=0.2, fc='none', ec='none', ax=None):
    codes = np.zeros(len(pgcoords), dtype=int) + path.Path.LINETO
    codes[0] = path.Path.MOVETO
    codes[-1] = path.Path.CLOSEPOLY
    pth = path.Path(pgcoords, codes)
    ptch = patches.PathPatch(pth, alpha=alpha, fc=fc, ec=ec)
    ax.add_patch(ptch)
    return
               

# %% load the unformatted files with the heads, the concentration and the budget terms
sim = flopy.mf6.MFSimulation.load(
    sim_name=mf_adapt.sim_name, version='mf6', sim_ws=dirs.SIM)

gwf = sim.get_model('{}Gwf'.format(sim.name).lower()) # list(sim.model_names)[0])

# Meaurements
hdata = dict()
for tseries in mf_adapt.tseries:
    name = tseries['time_series_namerecord']
    hdata[name] = pd.read_csv(os.path.join(dirs.Data , name + '.csv'), header=None, index_col=0)
    hdata[name].index = np.datetime64(mf_adapt.start_date_time) + hdata[name].index * (np.datetime64('2020-04-01') - np.datetime64('2018-01-01'))
    hdata[name][1] *= mf_adapt.FEET

headsObj = gwf.output.head()
budObj   = gwf.output.budget()

Qchd = [np.abs(q['q']).sum() / 2 for q in budObj.get_data(text='chd')]

hmin, hmax, hInactive = np.unique(headsObj.get_alldata())[[0, -2, -1]]

#  %% Steam function:
flowjas = budObj.get_data(text='FLOW-JA')

fflows = get_struct_flows([flowjas], grb_file=grb_file, verbose=False)

datetimes = np.datetime64(mf_adapt.start_date_time) + np.array(headsObj.times) * np.timedelta64(1, 'D')

fig, ax = plt.subplots()
fig.set_size_inches(12, 5)
ttl = ax.set_title("{}".format(sim.name))
ax.grid(True)
ax.set_xlabel = 'r [m]'
ax.set_ylabel = 'h [m]'

rs = [500, 1000, 2000, 4000]
Ix = gr.Ix(rs)

h1 = headsObj.get_alldata()[:, 1, 0, Ix]
h3 = headsObj.get_alldata()[:, 3, 0, Ix]


#for h, r  in zip(h1.T, rs):
#    ax.plot(datetimes, h, label='WVP1, r= {:.0f} m'.format(r))
#for h, r  in zip(h3.T, rs):
#    ax.plot(datetimes, h, label='WVP2, r= {:.0f} m'.format(r))

for name in hdata:
    ax.plot(hdata[name].index, hdata[name].values, '.-', label=name)

ax.legend()

plt.show()

