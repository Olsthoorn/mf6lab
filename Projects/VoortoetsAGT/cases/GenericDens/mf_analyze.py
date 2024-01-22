# Analyzation and visualization of simulation results

# i=== imports =====
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib as mpl
import flopy

# === local imports =====
import mf_adapt
from etc import newfig
from src.mf6contourtools import get_contour_levels
from src import mf6tools
from fdm.mf6_face_flows import get_struct_flows # (flowjas, grb_file=None, verbose=False)

# ==== constants =====
SAVE_ANIMATION = True
cmap = 'RdYlBu_r' # 'seismic', 'viridis_r'

# === From mf_adapt ====
dirs     = mf_adapt.dirs
sim_name = mf_adapt.sim_name
gr       = mf_adapt.gr_new
lay      = mf_adapt.lay
pr       = mf_adapt.pr
IDOMAIN  = mf_adapt.IDOMAIN

Iglob_wt = gr.Iglob_from_lrc(np.vstack((mf_adapt.Iz, np.zeros(gr.nx, dtype=int), gr.NOD[0, 0])).T)

# === load the the model with the simulation results ===== 
sim = flopy.mf6.MFSimulation.load(
    sim_name=mf_adapt.sim_name, version='mf6', sim_ws=dirs.SIM)

time_units = sim.tdis.time_units.get_data()

gwf = sim.get_model('{}Gwf'.format(sim.name).lower()) # list(sim.model_names)[0])
gwt = sim.get_model('{}Gwt'.format(sim.name).lower())

# === load the simulation results =====
headsObj = gwf.output.head()
budObj   = gwf.output.budget()
concObj  = gwt.output.concentration()

# === the date times at the end of each time step =====
datetimes = np.datetime64(mf_adapt.start_date_time) + np.array(headsObj.times) * np.timedelta64(1, 'D')

kstpkper = np.array(headsObj.get_kstpkper())

# %% Get the range needed for contouring
hmin, hmax, hInactive = np.unique(headsObj.get_alldata())[[0, -2, -1]] # the last one may be 1e30 (inactive)
cmin, cmax, cInactive = np.unique(concObj.get_alldata())[[ 0, -2, -1]] # the last one may be 1e30 (inactive)

#  === flow right face for computing the stream function =====
flowjas = budObj.get_data(text='FLOW-JA')
grb_file = os.path.join(dirs.GWF, sim_name + 'Gwf.dis.grb')
fflows = get_struct_flows([flowjas], grb_file=grb_file, verbose=False)

# === Contour levels =====
P   =  pr['rch'] * gr.dx.sum() # estimate of extremes for stream function psi
pLevels = get_contour_levels(-P, P, 200)

cLevels = np.linspace(pr['cFresh'], pr['cSalt'], 50)

# === flow between successive stream lines =====
dpsi = np.diff(pLevels)[0]

startdatetime = np.datetime64(sim.tdis.start_date_time.get_data().replace('t', ' '))
datetimes = startdatetime + np.array(concObj.times) * np.timedelta64(1, 'm')

title = "{} voor rch = {:.1f} mm/d. Tussen 2 stroomlijnen = {:.2f} m2/d = rch over {:.0f} m".format(
    mf_adapt.section_name, pr['rch'] * 1000, dpsi, dpsi / pr['rch'])

# === setup initial plot for animation =====

frame = 0

# === Get initial values for animation =====
psi = gr.psi_row(fflows['frf'][frame], row=0)
c   = concObj.get_data(kstpkper=kstpkper[frame])[:, 0, :]
c[c < pr['cFresh']] = pr['cFresh']
c[c > pr['cSalt']]  = pr['cSalt']

ax  = newfig("", 'Lijnafstand [m]', 'mTAW', figsize=(15, 8))
fig = plt.gcf()
ttl = ax.set_title(title)

# === Add layer patches ===== with colors, codes and layer names
colors = lay['Color'].values
codes  = lay['Code'].values
names  = lay['Name'].values

layer_patches = mf_adapt.gr.layer_patches_x(row=0)
for p, clr, code, name in zip(layer_patches, colors, codes, names):
    p.set_fc(clr)
    p.set_ec('k')
    p.set_lw(0.25)
    p.set_alpha(1.0)
    p.set_label(code + ' ' + name)
    ax.add_patch(p)

# === Plot the layers ==== (we don't have to do this) thin black line around them.
for iL, zp in enumerate(gr.Z):
    ax.plot(gr.xm, zp[0], 'k-', lw=0.5) # , label=f'layer {iL}')

ax.legend(loc='lower left')

# === Add logo =====
ax.text(0.85, 0.05, str(np.datetime64('today')),
        fontsize=10, fontweight='normal', va='bottom', transform=plt.gcf().transFigure) # ax.transAxis

# === Initial water table =====
hwt = headsObj.get_data(kstpkper=kstpkper[frame]).ravel()[Iglob_wt]
wt = ax.plot(gr.xm, hwt, 'b-', lw=1, label='water tafel')

Zpx = gr.Zpx_limited_to_water_table(hwt)
Zmx = gr.Zmx_limited_to_water_table(hwt)

caxP = ax.contour( gr.Xp, Zpx, psi, levels=pLevels, linewidths=0.5) 
caxC = ax.contourf(gr.XM[:, 0, :], Zmx,  np.ma.array( c, mask=np.isnan(c)),
                   levels=cLevels, cmap=mpl.colormaps[cmap])

def update(frame):
    """Animation update function."""
    global kstpkper, ttl, wt, caxP, caxC, ax, cmap

    # Update title
    ttl.set_text(title + ' frame={}, datetime={}'.format(frame, datetimes[frame]))
    
    # Update water table
    hwt = headsObj.get_data(kstpkper=kstpkper[frame]).ravel()[Iglob_wt]
    wt[0].set_ydata(hwt)
    
    # Update psi
    psi = gr.psi_row(fflows['frf'][frame])
    for coll in caxP.collections:
        coll.remove()
    Zpx = gr.Zpx_limited_to_water_table(hwt)
    caxP = ax.contour(gr.Xp, Zpx, psi, levels=pLevels)

    c = concObj.get_data(kstpkper=kstpkper[frame])[ :, 0, :]
    c[c < pr['cFresh']] = pr['cFresh']
    c[c > pr['cSalt']]  = pr['cSalt']
    
    for coll in caxC.collections:
        coll.remove()
    Zmx = gr.Zmx_limited_to_water_table(hwt)
    caxC = ax.contourf(gr.XM[:, 0, :], Zmx,  np.ma.array( c, mask=np.isnan(c)),
                       levels=cLevels, cmap=mpl.colormaps[cmap])
    
    mf6tools.show_animation_progress(frame, nbreak=50, ntot=len(kstpkper))
    
    return ttl, wt, caxP, caxC # Return affected artists (caxH removed)

ani = animation.FuncAnimation(fig=fig, func=update,                              
                              frames=len(kstpkper),
                              interval=30,
                              repeat=False)

if SAVE_ANIMATION:
    videoname = os.path.join(dirs.images, sim.name + '.mp4')
    print("Saving video to: {}".format(videoname))
    ani.save(filename = videoname, writer='ffmpeg')

plt.savefig(os.path.join(dirs.images, sim.name + '.png'))

plt.show()
# %%
