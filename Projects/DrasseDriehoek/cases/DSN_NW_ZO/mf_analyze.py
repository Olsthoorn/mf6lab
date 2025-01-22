# %% Analyzing output of Modflow

# %%
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import flopy
import mf_adapt
import settings
from src.mf6contourtools import get_contour_levels
from fdm.mf6_face_flows import get_struct_flows
from src import mf6tools

# ==== constants =====
SAVE_ANIMATION = True

# === From mf_adapt =====
dirs = mf_adapt.dirs
sim_name = mf_adapt.sim_name
gr = mf_adapt.gr
Iz = mf_adapt.Iz # Layer number of top active cells
start_date_time = mf_adapt.start_date_time

# === Flopy's grid file =====
Iglob_wt = gr.Iglob_from_lrc(np.vstack((mf_adapt.Iz, np.zeros(gr.nx, dtype=int), gr.NOD[0, 0])).T)

# === load the model with the unformatted simlation results =====
sim = flopy.mf6.MFSimulation.load(sim_name=sim_name,
                                  version='mf6',
                                  sim_ws=dirs.SIM)

time_units = sim.tdis.time_units.get_data()

gwf = sim.get_model(f'{sim.name.lower()}Gwf') # [sim.model_names][0])

# === load the simulation results =====
headsObj = gwf.output.head()
budObj   = gwf.output.budget()

# === the date times at the end of each time step =====
datetimes = np.datetime64(mf_adapt.start_date_time) + np.array(headsObj.times) * np.timedelta64(1, 'D')

kstpkper = np.array(headsObj.get_kstpkper())

# === Get the range needed for contouring =====
hmin, hmax, hInactive = np.unique(headsObj.get_alldata())[[0, -2, -1]] # the last one may be 1e30 (inactive)

#  === flow right face for computing the stream function =====
flowjas = budObj.get_data(text='FLOW-JA')
grb_file = os.path.join(dirs.GWF, sim_name + 'Gwf.dis.grb')
fflows = get_struct_flows([flowjas], grb_file=grb_file, verbose=False)

# === rch ====
modelArea = gr.Area.sum()
rch = np.array([r['q'].sum() / modelArea for r in budObj.get_data(text='rch')])

# === Get psi and suitable stream-function levels =====
P =  settings.props['rch'] * gr.dx.sum()    # estimate extremes of the stream function
pLevels = get_contour_levels(-P, P, 200)

# === Flow between successive stream lines =====
dpsi = np.diff(pLevels)[0]

startdatetime = np.datetime64(sim.tdis.start_date_time.get_data().replace('t', ' '))
datetimes = startdatetime + np.array(headsObj.times) * np.timedelta64(1, 'D')

# === setup initial plot for animation =====

frame = 0

title = "{} voor rch = {:.1f} mm/d. Tussen 2 stroomlijnen = {:.2f} m2/d = rch over {:.0f} m".format(
    mf_adapt.section_name, settings.props['rch'] * 1000, dpsi, dpsi / settings.props['rch'])

settings.fig_kwargs['title'] = title

# === Get initial values for animation =====
psi = gr.psi_row(fflows['frf'][frame], row=0)

# === Plot colored layers =====
ax = settings.plot_cross_section(gr=gr, props=settings.lay, IzL=settings.IzL, IzR=settings.IzR, fault_ix=settings.fault_ix, **settings.fig_kwargs)

fig = plt.gcf()
title = mf_adapt.section_name
ttl = ax.set_title(title)

# === Logo at lower richt of figure =====
ax.text(0.85, 0.05, str(np.datetime64('today')),
        fontsize=10,
        fontweight='normal',
        va='bottom',
        transform=plt.gcf().transFigure)

# === Initial water table =====
hwt = headsObj.get_data(kstpkper=kstpkper[frame]).ravel()[Iglob_wt]
wt = ax.plot(gr.xm, hwt, 'b-', lw=1, label='water tafel')

Zpx = gr.Zpx_limited_to_water_table(hwt)
Zmx = gr.Zmx_limited_to_water_table(hwt)

caxP = ax.contour( gr.Xp, Zpx, psi, levels=pLevels, linewidths=0.5)

def update(frame):
    """Animation update function."""
    global kstpkper, ttl, wt, caxP, ax, cmap

    # Update title
    ttl.set_text(title + 
                 f'. dpsi={dpsi:.2f} m2/d, rch={rch[frame] * 1e3:.1f} mm/d = inf over {dpsi / rch[frame]:.0f} m, datetime={datetimes[frame]}, frame={frame}')
    
    # Update water table
    hwt = headsObj.get_data(kstpkper=kstpkper[frame]).ravel()[Iglob_wt]
    wt[0].set_ydata(hwt)
    
    # Update psi
    psi = gr.psi_row(fflows['frf'][frame])
    caxP.set_paths([])
    Zpx = gr.Zpx_limited_to_water_table(hwt)
    caxP = ax.contour(gr.Xp, Zpx, psi, levels=pLevels)
        
    mf6tools.show_animation_progress(frame, nbreak=50, ntot=len(kstpkper))
    
    return ttl, wt, caxP # Return affected artists (caxH removed)

ani = animation.FuncAnimation(fig=fig, func=update,                              
                              frames=len(kstpkper),
                              interval=30,
                              repeat=False)

if SAVE_ANIMATION:
    videoname = os.path.join(dirs.images, sim.name + '.mp4')
    print("Saving video to: {}".format(videoname))
    ani.save(filename = videoname, writer='ffmpeg')

# === Save the picture =====
plt.savefig(os.path.join(dirs.images, sim.name + '.png'))

# === Finish plot and show it =====
plt.show()
# %%
