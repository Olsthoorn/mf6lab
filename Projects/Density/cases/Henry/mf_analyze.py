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


SAVE_ANIMATION = True

dirs     = mf_adapt.dirs
grb_file = os.path.join(dirs.GWF, mf_adapt.sim_name + 'Gwf.dis.grb')
gr       = mf_adapt.gr
pr       = mf_adapt.pr # properties from settings

# %% load the unformatted files with the heads, the concentration and the budget terms
sim = flopy.mf6.MFSimulation.load(
    sim_name=mf_adapt.sim_name, version='mf6', sim_ws=dirs.SIM)

time_units = sim.tdis.time_units.get_data()

gwf = sim.get_model('{}Gwf'.format(sim.name).lower()) # list(sim.model_names)[0])
gwt = sim.get_model('{}Gwt'.format(sim.name).lower())

# %% load the simulation results
headsObj = gwf.output.head()
budObj   = gwf.output.budget()
concObj  = gwt.output.concentration()
densObj = flopy.utils.binaryfile.HeadFile(os.path.join(dirs.SIM, mf_adapt.sim_name + 'Gwf.rho'), text="density")

# %% the date times at the end of each time step
kstpkper = np.array(headsObj.get_kstpkper())

#  %% Get the flow right face for computing the stream function:
flowjas = budObj.get_data(text='FLOW-JA')

# fflows is dict with face flows (keys: 'fff', 'frf', 'flf') as in MF5 and before, works only with structured grid)
fflows = get_struct_flows([flowjas], grb_file=grb_file, verbose=False)

# %% Get suitable levels for contouring the stream function psi [m2/d]

P = mf_adapt.pr['qL'] * mf_adapt.gr.dz.sum()

# Contour levels for conc, stream funtion and heads
pLevels = np.unique(get_contour_levels(-P, P, 100))

SALT_STEP = 0.05
cLevels = np.linspace(pr['FRESH'], pr['SALT'], int((pr['SALT'] - pr['FRESH']) / SALT_STEP) + 1)

# Flow between each pair of stream lines
dpsi = np.diff(pLevels)[0]

startdatetime = np.datetime64(sim.tdis.start_date_time.get_data().replace('t', ' '))
datetimes = startdatetime + np.array(concObj.times) * np.timedelta64(1, 's')

title = "{} voor QL = {:.4g} m3/s. Tussen 2 stroomlijnen = {:.4g} m2/s".format(mf_adapt.section_name, P, dpsi)

# %% Set up initial plot of cross section with contouring

frame = 0

# Get initial values to start animation
#h = headsObj.get_data(kstpkper=kstpkper[frame])[:, 0, :]
psi = gr.psi_row(fflows['frf'][frame], row=0)
c   = concObj.get_data(kstpkper=kstpkper[frame])[:, 0, :]
c[c < mf_adapt.FRESH] = mf_adapt.FRESH
c[c > mf_adapt.SALT]  = mf_adapt.SALT

ax  = newfig("", 'x [m]', 'z [m]', figsize=(15, 8))
fig = plt.gcf()
ttl = ax.set_title(title)

# Add logo:
ax.text(0.85, 0.05, str(np.datetime64('today')),
        fontsize=10, fontweight='normal', va='bottom', transform=plt.gcf().transFigure) # ax.transAxis

# ax.legend(loc='lower left')

#caxH = ax.contour( gr.xc, gr.zc,  np.ma.array(  h, mask=np.isnan(h)), levels=hLevels)
caxP = ax.contour( gr.Xp, gr.Zpx(), psi, levels=pLevels, lw=0.5, label='Psi') 
caxC = ax.contourf(gr.XM[:, 0, :], gr.ZM[:, 0, :], c,
                   levels=cLevels, cmap=mpl.colormaps['OrRd'])

# Animation function
def update(frame):
    """Animation update function."""
    global kstpkper, ttl, caxP, caxC, ax

    # Update title
    ttl.set_text(title + ' frame={}, datetime={}'.format(frame, datetimes[frame]))
        
    # Update psi
    psi = gr.psi_row(fflows['frf'][frame])
    for coll in caxP.collections:
        coll.remove()
    caxP = ax.contour(gr.Xp, gr.Zpx(), psi, levels=pLevels)

    c = concObj.get_data(kstpkper=kstpkper[frame])[ :, 0, :]
    
    for coll in caxC.collections:
        coll.remove()
    caxC = ax.contourf(gr.XM[:, 0, :], gr.ZM[:, 0, :],  c,
                       levels=cLevels, cmap=mpl.colormaps['OrRd'])
    print('.', end="")
    if frame % 50 == 0:
        print('\n')
    
    return ttl, caxP, caxC # Return affected artists (caxH removed)

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
