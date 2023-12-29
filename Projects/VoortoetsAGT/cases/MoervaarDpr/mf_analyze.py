# %% Analyzing outpu

# %%
import os
import numpy as np
import matplotlib.pyplot as plt
import flopy
import mf_adapt
from etc import newfig
from src.mf6contourtools import get_contour_levels
from fdm.mf6_face_flows import get_struct_flows # (flowjas, grb_file=None, verbose=False)

dirs = mf_adapt.dirs
grb_file = os.path.join(dirs.GWF, mf_adapt.sim_name + 'Gwf.dis.grb')

# %% load the unformatted files with the heads, the concentration and the budget terms
sim = flopy.mf6.MFSimulation.load(
    sim_name=mf_adapt.sim_name, version='mf6', sim_ws=dirs.SIM)

gwf = sim.get_model('{}Gwf'.format(sim.name).lower()) # list(sim.model_names)[0])

headsObj = gwf.output.head()
budObj   = gwf.output.budget()

h = headsObj.get_alldata()

gr = mf_adapt.gr

#  %% Steam function:
flowjas = budObj.get_data(text='FLOW-JA')

fflows = get_struct_flows([flowjas], grb_file=grb_file, verbose=False)

datetimes = np.datetime64(mf_adapt.start_date_time) + np.array(headsObj.times) * np.timedelta64(1, 'D')

iper = -1

# %% Get psi and suitable levels
P =  mf_adapt.rch * gr.dx.sum() / 4 # estimate of psi extremes
psi = gr.psi_row(fflows['frf'][iper], row=0)
levels = get_contour_levels(-P, P, 100)
dpsi = np.diff(levels)[0]
rch = mf_adapt.rch

htop = np.array([h[iper, iz, 0, ix] for ix, iz in enumerate(mf_adapt.Iz[0])])

Zpx = gr.Zpx_limited_to_water_table(htop)
layer_patches = gr.layer_patches_x(row=0) # Get layer patches

# Plot the cross section
title = "{} voor rch = {:.1f} mm/d. Tussen 2 stroomlijnen = {:.2f} m2/d = rch over {:.0f} m".format(
    sim.name, rch * 1000, dpsi, dpsi / rch)
ax = newfig(title, 'x along section [m]', 'elevation [m]', figsize=(15, 8))

# Set patch colors and add to the axis before anything else

# Contour the stream lines
ax.contour(gr.Xp, Zpx, psi, levels=levels, lw=0.5, label='Psi')

# Plot the layers (we dont have to do this)
for iL, zp in enumerate(gr.Z):
    ax.plot(gr.xm, zp[0], 'k-', lw=0.5) # , label=f'layer {iL}')

# Plot the water table
ax.plot(gr.xm, htop, 'b-', lw=1, label='water tafel')

# Get layer names and add to the legend
colors = mf_adapt.lay['Color'].values
codes  = mf_adapt.lay['Code'].values
names  = mf_adapt.lay['Name'].values

for p, clr, code, name in zip(layer_patches, colors, codes, names):
    p.set_fc(clr)
    p.set_ec('k')
    p.set_lw(0.25)
    p.set_alpha(1.0)
    p.set_label(code + ' ' + name)
    ax.add_patch(p)

# TODO:
ax.text(0.85, 0.05, str(np.datetime64('today')),
        fontsize=10, fontweight='normal', va='bottom', transform=plt.gcf().transFigure) # ax.transAxis

# add logo and date.
# More text about Psi.
ax.legend(loc='lower left')

plt.savefig(os.path.join(dirs.images, sim.name + '.png'))

plt.show()
# %%
