# %% Analyzing outpu

# %%
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib as mpl
import flopy
import mf_adapt
from etc import newfig
from src.mf6contourtools import get_contour_levels
from fdm.mf6_face_flows import get_struct_flows # (flowjas, grb_file=None, verbose=False)
from src import mf6tools

use_models, use_packages  = mf6tools.get_models_and_packages_from_excel(
                                                mf_adapt.params_wbk, sheet_name='NAM')

# === local settings ===
SAVE_ANIMATION = True

cmap = 'RdYlBu_r' # 'seismic'

dirs     = mf_adapt.dirs
grb_file = os.path.join(dirs.GWF, mf_adapt.sim_name + 'Gwf.dis.grb')
gr       = mf_adapt.gr
pr       = mf_adapt.pr # properties from settings

# === load the unformatted files with the heads, the concentration and the budget terms
sim = flopy.mf6.MFSimulation.load(
    sim_name=mf_adapt.sim_name, version='mf6', sim_ws=dirs.SIM)

time_units = sim.tdis.time_units.get_data()

gwf = sim.get_model('{}Gwf'.format(sim.name).lower()) # list(sim.model_names)[0])
gwt = sim.get_model('{}Gwt'.format(sim.name).lower())

# === Load the simulation results =====
# headsObj = gwf.output.head()
budObj   = gwf.output.budget()
concObj  = gwt.output.concentration()
densObj = flopy.utils.binaryfile.HeadFile(os.path.join(dirs.SIM, mf_adapt.sim_name + 'Gwf.rho'), text="density")

# === The date times at the end of each time step
kstpkper = np.array(concObj.get_kstpkper())
sim_time =  np.array(concObj.get_times())

# === Get the flf, fff, frf for computing the stream function =====
flowjas = budObj.get_data(text='FLOW-JA')
fflows = get_struct_flows([flowjas], grb_file=grb_file, verbose=False)

# === Get suitable levels for contouring the stream function psi [m2/d]
P = np.ceil(np.abs(mf_adapt.Qin))
pLevels = np.linspace(-P, P, 51)
dpsi = np.diff(pLevels)[0]

cLevels = np.linspace(pr['cFresh'], pr['cSalt'], 51)
cLevels2 = np.array([0.01, 0.1, 0.5, 0.9, 0.99]) * pr['cSalt']

startdatetime = np.datetime64(pr['start_date_time'])
datetimes = startdatetime + np.array(concObj.times) * np.timedelta64(1, 's')

# === Set up initial plot for animation =====

title = "{} Qin = {:.3f} m2/d. Stream lines: dPsi = {:.3f} m2/d".format(mf_adapt.section_name, mf_adapt.Qin, dpsi)

ax  = newfig("", 'x [m]', 'z [m]', figsize=(15, 8))
fig = plt.gcf()

ttl = ax.set_title(title)

# === Add logo =====
ax.text(0.85, 0.05, str(np.datetime64('today')),
        fontsize=10, fontweight='normal', va='bottom', transform=fig.transFigure)

# ax.legend(loc='lower left')

frame = 0

psi = gr.psi_row(fflows['frf'][frame], row=0)
c   = concObj.get_data(kstpkper=kstpkper[frame])[:, 0, :]
c[c < pr['cFresh']] = pr['cFresh']
c[c > pr['cSalt' ]] = pr['cSalt']

caxP = ax.contour( gr.Xp, gr.Zpx(), psi, levels=pLevels, linewidths=0.5) 
caxC = ax.contourf(gr.XM[:, 0, :], gr.ZM[:, 0, :], c,
                   levels=cLevels, cmap=mpl.colormaps[cmap])

if 'Gwtdsp' in use_packages:
    diff_txt = "With diffusion, Diffc = {} m2/d".format(pr['disp']['diffc'])
else:
    diff_txt = "Without diffusion"

fontdict = {'fontsize': 15,
            'fontweight': 'bold',
            'color': 'yellow',
            'backgroundcolor': 'blue',
            'va': 'bottom'}

ax.text(0.2, 0.7, diff_txt, **fontdict, transform=fig.transFigure)


def update(frame):
    """Animation update function."""
    global kstpkper, ttl, caxP, caxC, ax

    # Update title
    ttl.set_text(title + ' frame={}, sim time={:.3f} d'.format(frame + 1, sim_time[frame]))
        
    # Update psi
    psi = gr.psi_row(fflows['frf'][frame])
    for coll in caxP.collections:
        coll.remove()
    caxP = ax.contour(gr.Xp, gr.Zpx(), psi, levels=pLevels, linewidths=0.5)

    # Update c
    c = concObj.get_data(kstpkper=kstpkper[frame])[ :, 0, :]
    c[c < pr['cFresh']] = pr['cFresh']
    c[c > pr['cSalt']]  = pr['cSalt']
    
    for coll in caxC.collections:
        coll.remove()
    caxC = ax.contourf(gr.XM[:, 0, :], gr.ZM[:, 0, :],  c,
                       levels=cLevels, cmap=mpl.colormaps[cmap])
    
    # Show progress:
    print('.', end="")
    if (frame + 1) % 50 == 0 or (frame + 1) == len(kstpkper):
        print('{} frames\n'.format(frame + 1))
        if frame + 1 == len(kstpkper):
            caxC2 = ax.contour(gr.XM[:, 0, :], gr.ZM[:, 0, :],
                       concObj.get_data(kstpkper=kstpkper[-1])[ :, 0, :],
                       levels=cLevels2, colors='black', linewidths=1.0)
            ax.clabel(caxC2, caxC2.levels, inline=True, fmt=r'%2.2f', fontsize=10) 
    
    return ttl, caxP, caxC   # Return affected artists (caxH removed)
# =================================================================
ani = animation.FuncAnimation(fig=fig, func=update,                              
                              frames=len(kstpkper),
                              interval=30,
                              repeat=False)

# Add Contours

if SAVE_ANIMATION:
    videoname = os.path.join(dirs.images, sim.name + '.mp4')
    print("\nCounting frames to video in the makeing: {}\n".format(videoname))
    print("Saving video to: {}".format(videoname))
    ani.save(filename = videoname, writer='ffmpeg')

plt.savefig(os.path.join(dirs.images, sim.name + '.png'))

plt.show()

