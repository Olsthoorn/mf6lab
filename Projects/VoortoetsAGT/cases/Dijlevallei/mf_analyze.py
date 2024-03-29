# === Analyze and visualize simulation output =====
import os
import numpy as np
import matplotlib.pyplot as plt
import flopy

# === Local imports =====
import mf_adapt
from etc import newfig
from src.mf6contourtools import get_contour_levels
from fdm.mf6_face_flows import get_struct_flows # (flowjas, grb_file=None, verbose=False)

# === From mf_adapt =====
dirs = mf_adapt.dirs
sim_name = mf_adapt.sim_name
start_date_time = mf_adapt.start_date_time
gr = mf_adapt.gr
rch = mf_adapt.rch
Iz = mf_adapt.Iz
lay = mf_adapt.lay

# === load the model with simulation results. =====
sim = flopy.mf6.MFSimulation.load(sim_name=sim_name,
                                  version='mf6',
                                  sim_ws=dirs.SIM)

# === load simulation results =====
gwf = sim.get_model('{}Gwf'.format(sim.name).lower()) # list(sim.model_names)[0])
headsObj = gwf.output.head()
budObj   = gwf.output.budget()

# === Use head from last time step (steady) =====
h = headsObj.get_data(kstpkper=headsObj.get_kstpkper()[-1])

datetimes = np.datetime64(start_date_time) + np.array(headsObj.times) * np.timedelta64(1, 'D')

# === Get suitable levels for the stream function
P =  rch * gr.dx.sum() # estimate of psi extremes
levels = get_contour_levels(-P, P, 200)
dpsi = np.diff(levels)[0] # flow between successive stream lines

# === Plot the cross section =====
title = "{} voor rch = {:.1f} mm/d. Tussen 2 stroomlijnen = {:.2f} m2/d = rch over {:.0f} m".format(
    mf_adapt.section_name, rch * 1000, dpsi, dpsi / rch)

ax = newfig(title, 'x along section [m]', 'elevation [m]', figsize=(15, 8))

# ==== Stream function and contouring it. =====
flowjas = budObj.get_data(text='FLOW-JA')
grb_file = os.path.join(dirs.GWF, sim_name + 'Gwf.dis.grb')
fflows = get_struct_flows([flowjas], grb_file=grb_file, verbose=False)
psi = gr.psi_row(fflows['frf'][-1], row=0) # Stream function

# === Water table =====
htop = np.array([h[iz, 0, ix] for ix, iz in enumerate(Iz)])
ax.plot(gr.xm, htop, 'b-', lw=1, label='water tafel')

Zpx = gr.Zpx_limited_to_water_table(htop) # Adapt grid to water table
ax.contour(gr.Xp, Zpx, psi, levels=levels, lw=0.5, label='Psi')

# === Plot the layers ===== (we dont have to do this)
for iL, zp in enumerate(gr.Z):
    ax.plot(gr.xm, zp[0], 'k-', lw=0.5) # , label=f'layer {iL}')

# === Plot the layers as colored patches =====
# ==== layer colors, codes and names =====
colors = lay['Color'].values
codes  = lay['Code'].values
names  = lay['Name'].values

layer_patches = gr.layer_patches_x(row=0) # Get layer patches
for p, clr, code, name in zip(layer_patches, colors, codes, names):
    p.set_fc(clr)
    p.set_ec('k')
    p.set_lw(0.25)
    p.set_alpha(1.0)
    p.set_label(code + ' ' + name)
    ax.add_patch(p)

ax.legend(loc='lower left')

# === Logo ======
ax.text(0.85, 0.05, str(np.datetime64('today')),
        fontsize=10, fontweight='normal', va='bottom', transform=plt.gcf().transFigure) # ax.transAxis

plt.savefig(os.path.join(dirs.images, sim.name + '.png'))

plt.show()
