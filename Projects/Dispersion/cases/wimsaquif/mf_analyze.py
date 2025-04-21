# %% Analyzing output of wimsaquif
# TO 20250410
# %%
import os
import numpy as np
import matplotlib.pyplot as plt
import flopy

import mf_adapt
from src.mf6contourtools import get_contour_levels
from src import mf6tools
from fdm.mf6_face_flows import get_struct_flows # (flowjas, grb_file=None, verbose=False)

use_models, use_packages  = mf6tools.get_models_and_packages_from_excel(
                                                mf_adapt.params_wbk, sheet_name='NAM')

cmap = 'binary_r'

dirs     = mf_adapt.dirs
grb_file = os.path.join(dirs.GWF, mf_adapt.sim_name + 'Gwf.dis.grb')
gr       = mf_adapt.gr
pr       = mf_adapt.pr # properties from settings
case_name = pr['k_field_pars']['name']    
# %% load the unformatted files with the heads, the concentration and the budget terms
sim = flopy.mf6.MFSimulation.load(
    sim_name=mf_adapt.sim_name, version='mf6', sim_ws=dirs.SIM)

time_units = sim.tdis.time_units.get_data()

gwf = sim.get_model('{}Gwf'.format(sim.name).lower()) # list(sim.model_names)[0])

# %% load the simulation results
headsObj = gwf.output.head()
budObj   = gwf.output.budget()

# %% the date times at the end of each time step
kstpkper = np.array(headsObj.get_kstpkper())
sim_time =  np.array(headsObj.get_times())

#  %% Get the flow right face for computing the stream function:
flowjas = budObj.get_data(text='FLOW-JA')

# fflows: dict with face flows (keys: 'fff', 'frf', 'flf') as in MF5 
fflows = get_struct_flows([flowjas], grb_file=grb_file, verbose=False)
pmin, pmax = np.inf, -np.inf
for frf in fflows['frf']:
    pmin = min(pmin, np.min(frf.sum(axis=0)))
    pmax = max(pmax, np.max(frf.sum(axis=0)))

# %% Get suitable levels for contouring the stream function psi [m2/d]
hmin, hmax, hInactive = np.unique(headsObj.get_alldata())[[0, -2, -1]]

# cLevels = get_contour_levels(max(cmin, pr['cFresh']), min(cmax, pr['cSalt']), 20)[1:] # don't show zero to keep it transparant.
pLevels = get_contour_levels(pmin, pmax, 20)
hLevels = get_contour_levels(hmin, hmax, 20)

dpsi = np.diff(pLevels)[0]

startdatetime = np.datetime64(sim.tdis.start_date_time.get_data().replace('t', ' '))
datetimes = startdatetime + np.array(headsObj.times) * np.timedelta64(1, 'm')

# %% Set up initial plot of cross section with contouring
fig, ax = plt.subplots(figsize=(10,6))
ttl =  ax.set_title(mf_adapt.section_name + 'steamlines and headcontours\n' + pr['k_field_pars']['k_field_str'])
ax.set(title=ttl, xlabel='x [m]', ylabel='z [m]')

# logo:
ax.text(0.85, 0.05, str(np.datetime64('today')),
        fontsize=10, fontweight='normal', va='bottom', transform=fig.transFigure) # ax.transAxis

frame = -1

# Data for the the first frame
psi = gr.psi_row(fflows['frf'][frame], row=0)

# Save psi, for later use by Modpath
np.save(os.path.join(dirs.GWF, "psi.npy"), psi)

h = np.squeeze(headsObj.get_data(kstpkper=kstpkper[frame]))
caxH = ax.contour( gr.xc, gr.zc,  np.ma.array(  h, mask=(h==hInactive)), levels=hLevels)
caxP = ax.contour( gr.Xp, gr.Zpx(), psi, levels=pLevels, colors='k', linewidths=0.5) 

Qchd = budObj.get_data(text='chd', kstpkper=kstpkper[frame])[0]['q']
Qthrough = Qchd[Qchd > 0].sum()

# Update title

ttl.set_text("Stream function (case " + case_name + f') Qmdl = {Qthrough:.1f} {pr['Qdim']}\n'
             + pr['k_field_pars']['k_field_str'])

ax.grid()
ax.legend(loc='lower left')
    
plt.savefig(os.path.join(dirs.images, 'stream_func_' + case_name + '.png'))

plt.show()
