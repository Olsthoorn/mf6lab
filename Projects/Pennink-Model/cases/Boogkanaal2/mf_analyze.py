# === Analyze and visualize simulation output =====
import os
import numpy as np
import matplotlib.pyplot as plt
import flopy
from PIL import Image


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
Iz = mf_adapt.Iz
pr = mf_adapt.pr
k33 = mf_adapt.k33
IDOMAIN = mf_adapt.IDOMAIN

photo = Image.open(os.path.join(dirs.photos, pr['photo']))

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
hMin, hMax, hInact = np.unique(h)[[0, -2, -1]]
hlevels = get_contour_levels(hMin, hMax, 50)
dh = np.diff(hlevels)[0]

datetimes = np.datetime64(start_date_time) + np.array(headsObj.times) * np.timedelta64(1, 'D')

# === Get suitable levels for the stream function
P =  budObj.get_data(text='chd')[0]['q']
plevels = get_contour_levels(np.sum(P[P<0]), np.sum(P[P>0]), 50)
dpsi = np.diff(plevels)[0] # flow between successive stream lines

# === Plot the cross section =====
title = "{}, dh = {:.3f} m, dpsi = {:.2f} m2/d, Qtotal = {:.2} m2/d".format(
    mf_adapt.section_name, dh, dpsi, P[P > 0].sum())

ax = newfig(title, 'x along section [m]', 'elevation [m]', figsize=(10, 10))

ax.imshow(photo, extent=pr['extent'], alpha=0.75)

# ==== Stream function and contouring it. =====
flowjas = budObj.get_data(text='FLOW-JA')
grb_file = os.path.join(dirs.GWF, sim_name + 'Gwf.dis.grb')
fflows = get_struct_flows([flowjas], grb_file=grb_file, verbose=False)
psi = gr.psi_row(fflows['frf'][-1], row=0) # Stream function

# === Water table =====
htop = np.array([h[iz, 0, ix] for ix, iz in enumerate(Iz)])
ax.plot(gr.xm, htop, 'b-', lw=1, label='water tafel', zorder=5)

Zpx = gr.Zpx_limited_to_water_table(htop) # Adapt grid to water table
ax.contour(gr.Xp, Zpx, psi, levels=plevels, lw=0.5, label='Psi',  colors='orange')

# === Contour the heads ====
caxh = ax.contour(gr.XM_on_zplanes[:, 0, :],
                  gr.Z_limited_to_water_table(htop[np.newaxis, :])[:, 0, :],
                  gr.heads_on_zplanes(h, fflows['flf'][-1], gr.const(pr['k33']))[:, 0, :],
                  levels=hlevels, colors='w')

# === Plot the layers ===== (we dont have to do this)
# for iL, zp in enumerate(gr.Z):
# f   ax.plot(gr.xm, zp[0], 'k-', lw=0.5) # , label=f'layer {iL}')
          
ax.plot(*pr['frame'].T, 'k', lw=0.5, label='frame')

for p in pr['canalPatches']:
    ax.add_patch(p)

ax.legend(loc='lower left')

# === Logo ======
ax.text(0.70, 0.05, str(np.datetime64('today')),
        fontsize=10, fontweight='normal', va='bottom', transform=plt.gcf().transFigure) # ax.transAxis

plt.savefig(os.path.join(dirs.images, sim.name + '.png'))

plt.show()
