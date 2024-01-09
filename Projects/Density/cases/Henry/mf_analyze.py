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


def h_water_table(h, Iz):
    """Return elevation of the water table.
    
    Parameters
    ----------
    h: heads for a given kstpkperature
        e.g.: h = headsObj(kstpkper=(3, 4))
    Iz: z-index of first non-nan cell
        in mf_adapt.Iz
        
    Example
    -------
    hwt = w_water_table(h, Iz)
    """
    return h.ravel()[Iz]

SAVE_ANIMATION = True

dirs     = mf_adapt.dirs
grb_file = os.path.join(dirs.GWF, mf_adapt.sim_name + 'Gwf.dis.grb')
gr       = mf_adapt.gr_new

# Globel cell numbers for top active cells:
Iz = np.array([gr.NOD[iz, 0, ix] for ix, iz in zip(range(gr.nx), mf_adapt.Iz)], dtype=int)

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

# %% the date times at the end of each time step
datetimes = np.datetime64(mf_adapt.start_date_time) + np.array(headsObj.times) * np.timedelta64(1, 'D')

kstpkper = np.array(headsObj.get_kstpkper())

# %% Get the range needed for contouring
hmin, hmax, hInactive = np.unique(headsObj.get_alldata())[[0, -2, -1]]
cmin, cmax, cInactive = np.unique(concObj.get_alldata())[[ 0, -2, -1]]

#  %% Get the flow right face for computing the stream function:
flowjas = budObj.get_data(text='FLOW-JA')

# fflows is dict with face flows (keys: 'fff', 'frf', 'flf') as in MF5 and before, works only with structured grid)
fflows = get_struct_flows([flowjas], grb_file=grb_file, verbose=False)

# %% Get suitable levels for contouring the stream function psi [m2/d]
rch = mf_adapt.rch
P   =  rch * gr.dx.sum() / 2 # estimate of extremes for stream function psi

# Contour levels for conc, stream funtion and heads
pLevels = get_contour_levels(-P, P, 100)
# hLevels = get_contour_levels(hmin, hmax, 20)
# cLevels = get_contour_levels(cmin, cmax, 20)
# cLevels = cLevels[cLevels >= 0.0]
SALT_STEP = 500.0
cLevels = np.linspace(mf_adapt.FRESH, mf_adapt.SALT, int((mf_adapt.SALT - mf_adapt.FRESH) / SALT_STEP) + 1)

# Flow between each pair of stream lines
dpsi = np.diff(pLevels)[0]

startdatetime = np.datetime64(sim.tdis.start_date_time.get_data().replace('t', ' '))
datetimes = startdatetime + np.array(concObj.times) * np.timedelta64(1, 'm')

title = "{} voor rch = {:.1f} mm/d. Tussen 2 stroomlijnen = {:.2f} m2/d = rch over {:.0f} m".format(
    mf_adapt.section_name, rch * 1000, dpsi, dpsi / rch)

# %% Set up initial plot of cross section with contouring

frame = 0

# Get initial values to start animation
#h = headsObj.get_data(kstpkper=kstpkper[frame])[:, 0, :]
psi = gr.psi_row(fflows['frf'][frame], row=0)
c   = concObj.get_data(kstpkper=kstpkper[frame])[:, 0, :]
c[c < mf_adapt.FRESH] = mf_adapt.FRESH
c[c > mf_adapt.SALT]  = mf_adapt.SALT

ax  = newfig("", 'Lijnafstand [m]', 'mTAW', figsize=(15, 8))
fig = plt.gcf()
ttl = ax.set_title(title)

# Add layer patches with colors, codes and layer names
colors = mf_adapt.lay['Color'].values
codes  = mf_adapt.lay['Code'].values
names  = mf_adapt.lay['Name'].values

layer_patches = mf_adapt.gr.layer_patches_x(row=0)

if True:
    for p, clr, code, name in zip(layer_patches, colors, codes, names):
        p.set_fc(clr)
        p.set_ec('k')
        p.set_lw(0.25)
        p.set_alpha(1.0)
        p.set_label(code + ' ' + name)
        ax.add_patch(p)

# Plot the layers (we don't have to do this) thin black line around them.
for iL, zp in enumerate(gr.Z):
    ax.plot(gr.xm, zp[0], 'k-', lw=0.5) # , label=f'layer {iL}')

# Add logo:
ax.text(0.85, 0.05, str(np.datetime64('today')),
        fontsize=10, fontweight='normal', va='bottom', transform=plt.gcf().transFigure) # ax.transAxis

ax.legend(loc='lower left')

# Initial water table and contouring
hwt = headsObj.get_data(kstpkper=kstpkper[frame]).ravel()[Iz]
wt = ax.plot(gr.xm, hwt, 'b-', lw=1, label='water tafel')

Zpx = gr.Zpx_limited_to_water_table(hwt)
Zmx = gr.Zmx_limited_to_water_table(hwt)

#caxH = ax.contour( gr.xc, gr.zc,  np.ma.array(  h, mask=np.isnan(h)), levels=hLevels)
caxP = ax.contour( gr.Xp, Zpx, psi, levels=pLevels, lw=0.5, label='Psi') 
caxC = ax.contourf(gr.XM[:, 0, :], Zmx,  np.ma.array( c, mask=np.isnan(c)),
                   levels=cLevels, cmap=mpl.colormaps['viridis_r'])

# Animation function
def update(frame):
    """Animation update function."""
    global kstpkper, ttl, wt, caxP, caxC, ax

    # Update title
    ttl.set_text(title + ' frame={}, datetime={}'.format(frame, datetimes[frame]))
    
    # Update water table
    hwt = headsObj.get_data(kstpkper=kstpkper[frame]).ravel()[Iz]
    wt[0].set_ydata(hwt)
    
    # Update psi
    psi = gr.psi_row(fflows['frf'][frame])
    for coll in caxP.collections:
        coll.remove()
    Zpx = gr.Zpx_limited_to_water_table(hwt)
    caxP = ax.contour(gr.Xp, Zpx, psi, levels=pLevels)

    #h = headsObj.get_data(kstpkper=kstpkper[frame])[:, 0, :]
    #for coll in caxH.collections:
    #    coll.remove()
    #caxH = ax.contour( gr.xc, gr.zc,  np.ma.array(  h, mask=np.isnan(h)),
    #                    levels=hLevels)

    c = concObj.get_data(kstpkper=kstpkper[frame])[ :, 0, :]
    c[c < mf_adapt.FRESH] = mf_adapt.FRESH
    c[c > mf_adapt.SALT]  = mf_adapt.SALT
    
    for coll in caxC.collections:
        coll.remove()
    Zmx = gr.Zmx_limited_to_water_table(hwt)
    caxC = ax.contourf(gr.XM[:, 0, :], Zmx,  np.ma.array( c, mask=np.isnan(c)),
                       levels=cLevels, cmap=mpl.colormaps['viridis_r'])
    
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
