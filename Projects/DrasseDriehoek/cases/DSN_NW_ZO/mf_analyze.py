# %% Analyzing outpu

# %%
import os
import numpy as np
import matplotlib.pyplot as plt
import flopy
import mf_adapt
from etc import newfig
from src.mf6contourtools import get_contour_levels
from fdm.mf6_face_flows import get_struct_flows
from settings import props as pr

# === From mf_adapt =====
dirs = mf_adapt.dirs
sim_name = mf_adapt.sim_name
gr = mf_adapt.gr
Iz = mf_adapt.Iz # Layer number of top active cells
start_date_time = mf_adapt.start_date_time

# === Flopy's grid file =====

# === load the unformatted files with the heads, budget and possible concentration =====
sim = flopy.mf6.MFSimulation.load(sim_name=sim_name,
                                  version='mf6',
                                  sim_ws=dirs.SIM)

# === Groundwater FLow Object =====
gwf = sim.get_model('{}Gwf'.format(sim.name).lower()) # list(sim.model_names)[0])

headsObj = gwf.output.head()
budObj   = gwf.output.budget()

h = headsObj.get_data(kstpkper=headsObj.get_kstpkper()[-1])
datetimes = np.datetime64(start_date_time) + np.array(headsObj.times) * np.timedelta64(1, 'D')

# === Get strcutured flows flf, fff and frf first (we only need frf) ====
flowjas = budObj.get_data(text='FLOW-JA')
grb_file = os.path.join(dirs.GWF, sim_name + 'Gwf.dis.grb') # Flopy's grid file
fflows = get_struct_flows([flowjas], grb_file=grb_file, verbose=False)

#  === Stream function ===== Stream function is called psi.
psi = gr.psi_row(fflows['frf'][-1], row=0)

# === Get psi and suitable stream functioin levels =====
P =  pr['rch'] * gr.dx.sum()    # estimate extremes of the stream function
levels = get_contour_levels(-P, P, 200)

# === Flow between successive stream lines =====
dpsi = np.diff(levels)[0]

# === Plot the cross section =====
title = "{} voor rch = {:.1f} mm/d. Tussen 2 stroomlijnen = {:.2f} m2/d = rch over {:.0f} m".format(
    mf_adapt.section_name, pr['rch'] * 1000, dpsi, dpsi / pr['rch'])

ax = newfig(title, 'Lijnafstand [m]', 'mTAW', figsize=(15, 8))

# === plot the water table =====
htop = np.array([h[iz, 0, ix] for ix, iz in enumerate(Iz)]) # Water table
ax.plot(gr.xm, htop, 'b-', lw=1, label='water tafel') # plot water table

# ==== Contour the stream lines =====
Zpx = gr.Zpx_limited_to_water_table(htop) # Adapt psi contour grid to water table
ax.contour(gr.Xp, Zpx, psi, levels=levels, lw=0.5, label='Psi')

# === Plot the layers ===== (we dont have to do this)
for iL, zp in enumerate(gr.Z):
    ax.plot(gr.xm, zp[0], 'k-', lw=0.5) # , label=f'layer {iL}')

# === Plot colored layers (using layer patches( =====
layer_patches = gr.layer_patches_x(row=0) # Get layer patches

# === Colors, codes and names of the layers =====
colors = mf_adapt.lay['Color'].values
codes  = mf_adapt.lay['Code'].values
names  = mf_adapt.lay['Name'].values

# === set patch properties =====
for p, clr, code, name in zip(layer_patches, colors, codes, names):
    p.set_fc(clr)
    p.set_ec('k')
    p.set_lw(0.25)
    p.set_alpha(1.0)
    p.set_label(code + ' ' + name)
    ax.add_patch(p)

# === put legend at a suitable place =====
ax.legend(loc='lower left')

# === Logo at lower richt of figure =====
ax.text(0.85, 0.05, str(np.datetime64('today')),
        fontsize=10,
        fontweight='normal',
        va='bottom',
        transform=plt.gcf().transFigure)

# === Save the picture =====
plt.savefig(os.path.join(dirs.images, sim.name + '.png'))

# === Finish plot and show it =====
plt.show()
# %%
