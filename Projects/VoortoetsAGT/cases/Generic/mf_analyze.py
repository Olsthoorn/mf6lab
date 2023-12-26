# %% Analyzing outpu

# %%
import os
import numpy as np
import matplotlib.pyplot as plt
import flopy
import mf_adapt
from mf_adapt import gr
from etc import newfig
from mf6lab.mf6contourtools import get_contour_levels
from fdm.mf6_face_flows import get_struct_flows # (flowjas, grb_file=None, verbose=False)
from etc import color_cycler

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

psi = gr.psi_row(fflows['frf'][iper], row=0)
levels = get_contour_levels(-5, 5, 80)

htop = np.array([h[iper, iz, 0, ix] for ix, iz in enumerate(mf_adapt.Iz[0])])

Zpx = gr.Zpx_limited_to_water_table(htop)
layer_patches = gr.layer_patches_x(row=0) # Get layer patches

# Plot the cross section
ax = newfig("{}".format(sim.name), 'x along section [m]', 'elevation [m]', figsize=(12, 5))

colors = ['yellow', 'orange', 'gold', 'lightskyblue',
          'lightsalmon', 'violet', 'chocolate', 'yellowgreen']

# Set patch colors and add to the axis before anything else
for p, clr in zip(layer_patches, color_cycler(colors)):
    p.set_fc(clr)
    p.set_ec('k')
    p.set_lw(0.25)
    p.set_alpha(1.0)
    ax.add_patch(p)

# Contour the stream lines
ax.contour(gr.Xp, Zpx, psi, levels=levels, lw=0.5, label='Psi')

# Plot the layers (we dont have to do this)
for iL, zp in enumerate(gr.Z):
    ax.plot(gr.xm, zp[0], 'k-', lw=0.5, label=f'layer {iL}')

# Plot the water table
ax.plot(gr.xm, htop, 'b-', lw=1, label='water tafel')

ax.legend()

plt.show()

# %%
