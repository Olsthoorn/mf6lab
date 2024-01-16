# %% Analyzing output of Pennik Series 2 tests
# TO 090223
# %%
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
import matplotlib.path as path
from matplotlib.colors import LinearSegmentedColormap
import flopy
import mf_adapt
from src.mf6contourtools import get_contour_levels
from src import mf6tools
from fdm.mf6_face_flows import get_struct_flows # (flowjas, grb_file=None, verbose=False)
from PIL import Image

use_models, use_packages  = mf6tools.get_models_and_packages_from_excel(
                                                mf_adapt.params_wbk, sheet_name='NAM')

SAVE_ANIMATION = True
cmap = 'binary_r'

colors = ['gray', 'black']
colors = ['gray', 'white']
cmap = LinearSegmentedColormap.from_list("GrBl", colors)

dirs     = mf_adapt.dirs
grb_file = os.path.join(dirs.GWF, mf_adapt.sim_name + 'Gwf.dis.grb')
gr       = mf_adapt.gr
pr       = mf_adapt.pr # properties from settings
dirs.photos = (os.path.join(dirs.case, 'photos'))
foto = Image.open(os.path.join(dirs.photos, pr['photo']))


def addPatch(pgcoords, alpha=0.2, fc='none', ec='none', ax=None):
    codes = np.zeros(len(pgcoords), dtype=int) + path.Path.LINETO
    codes[0] = path.Path.MOVETO
    codes[-1] = path.Path.CLOSEPOLY
    pth = path.Path(pgcoords, codes)
    ptch = patches.PathPatch(pth, alpha=alpha, fc=fc, ec=ec)
    ax.add_patch(ptch)
    return
               
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
kstpkper = np.array(headsObj.get_kstpkper())
sim_time =  np.array(headsObj.get_times())

#  %% Get the flow right face for computing the stream function:
flowjas = budObj.get_data(text='FLOW-JA')

# fflows: dict with face flows (keys: 'fff', 'frf', 'flf') as in MF5 
fflows = get_struct_flows([flowjas], grb_file=grb_file, verbose=False)
pmin, pmax = np.inf, -np.inf
for frf in fflows['frf']:
    pmin = min(pmin, np.min(frf.sum(axis=0)))
    pmax = max(pmax, np.max(frf.sum(axis=0)))

# %% Get suitable levels for contouring the stream function psi [m2/d]
hmin, hmax, hInactive = np.unique(headsObj.get_alldata())[[0, -2, -1]]
cmin, cmax, cInactive = np.unique(concObj.get_alldata( ))[[0, -2, -1]]

cLevels = get_contour_levels(max(cmin, pr['cFresh']), min(cmax, pr['cSalt']), 20)[1:] # don't show zero to keep it transparant.
pLevels = get_contour_levels(pmin, pmax, 20)
hLevels = get_contour_levels(hmin, hmax, 20)

# Get top active cell:
Ix, IzWt = gr.NOD[0, 0, :], np.zeros(gr.nx, dtype=int)
while len(Ix) > 0:
    IxNew = []
    for ix in Ix:
        if mf_adapt.IDOMAIN[IzWt[ix], 0, ix] < 0:
            IzWt[ix] += 1
            IxNew.append(ix)
    Ix = IxNew

dpsi = np.diff(pLevels)[0]

startdatetime = np.datetime64(sim.tdis.start_date_time.get_data().replace('t', ' '))
datetimes = startdatetime + np.array(concObj.times) * np.timedelta64(1, 'm')

# %% Set up initial plot of cross section with contouring
fig, ax = plt.subplots()
fig.set_size_inches(12,12)
ax.grid(False)
ax.set_xlabel = 'x [cm]'
ax.set_ylabel = 'y [cm]'
ttl =  ax.set_title(mf_adapt.section_name)

extent = pr['extent']

ax.imshow(foto, extent=extent)

# Box around model (only to check coordinates of extent)
ax.plot((0, pr['L'], pr['L'], 0, 0), (0, 0, pr['H'], pr['H'], 0), 'k')

if False:     
    p1 = patches.PathPatch(mf_adapt.sand, alpha=0.5, fc='brown', ec='none')
    p2 = patches.PathPatch(mf_adapt.canalL, alpha=0.5, fc='Blue', ec='black')
    p3 = patches.PathPatch(mf_adapt.canalR, alpha=0.2, fc='Gray', ec='black')

    for pgcoords, clr in zip([mf_adapt.sand, mf_adapt.canalL, mf_adapt.canalR],
                            ['brown', 'blue', 'green']):
        addPatch(pgcoords, alpha=0.5, fc=clr, ec='black', ax=ax)

# logo:
ax.text(0.85, 0.05, str(np.datetime64('today')),
        fontsize=10, fontweight='normal', va='bottom', transform=fig.transFigure) # ax.transAxis

# ax.legend(loc='lower left')

frame = 0

# Data for the the first frame
psi = gr.psi_row(fflows['frf'][frame], row=0)
c   = concObj.get_data(kstpkper=kstpkper[frame])[:, 0, :]
c[c < pr['cFresh']] = pr['cFresh']
c[np.logical_and(c > pr['cSalt'  ], c <= cmax)] = pr['cSalt']

# caxH = ax.contour( gr.xc, gr.zc,  np.ma.array(  h, mask=(h=hInactive)), levels=hLevels)
caxP = ax.contour( gr.Xp, gr.Zpx(), psi, levels=pLevels, colors='yellow', linewidths=0.5) 
caxC = ax.contourf(gr.XM[:, 0, :], gr.ZM[:, 0, :], np.ma.array( c, mask=(c == cInactive)),
                   levels=cLevels, cmap=cmap, alpha=0.5)

def h_water_table(headsObj, frame, IzWt, gr):
    """Return the water table.
    
    Parameters
    ----------
    headsObj: heads object
        from flopy simulation results
    frame: int
        frame number in movie
    IzWt: ndarray of int
        array of iz indices of active cells
    gr: Grid object
        holds the model network.
    """
    hwt = np.array([headsObj.get_data(kstpkper=kstpkper[frame])[iz, 0, ix] for ix, iz in zip(gr.NOD[0, 0, :], IzWt)])
    hwt[hwt == hInactive] = np.nan
    return hwt

wtbl = ax.plot(gr.xm, h_water_table(headsObj, frame, IzWt, gr), 'b', lw=1, label='water table')[0]

# if 'Gwtdsp' in use_packages:
#     diff_txt = "With diffusion, Diffc = {} m2/d".format(pr['disp']['diffc'])
# else:
#     diff_txt = "Without diffusion"
# 
# fontdict = {'fontsize': 15,
#             'fontweight': 'bold',
#             'color': 'yellow',
#             'backgroundcolor': 'blue',
#             'va': 'bottom'}
#
# ax.text(0.2, 0.7, diff_txt, **fontdict, transform=fig.transFigure) # ax.transAxis)


def update(frame):
    """Animation update function."""
    global ttl, caxH, caxP, caxC, ax, kstpkper
    
    Qchd = budObj.get_data(text='chd', kstpkper=kstpkper[frame])[0]['q']
    Qthrough = Qchd[Qchd > 0].sum()
    
    # Update title
    ttl.set_text(mf_adapt.section_name + ' frame={}, sim time={:.0f} min, Qmdl = {:.1f} {}'.format(
        frame + 1, sim_time[frame], Qthrough, pr['Qdim']))

    # Update Psi    
    psi = gr.psi_row(fflows['frf'][frame])
    for coll in caxP.collections:
        coll.remove()
    caxP = ax.contour(gr.Xp, gr.Zpx(), psi, levels=pLevels, colors='yellow', linewidths=0.5)

    # Update h
    wtbl.set_ydata(h_water_table(headsObj, frame, IzWt, gr))

    # Update c
    c = concObj.get_data(kstpkper=kstpkper[frame])[ :, 0, :]
    c[np.logical_and(c > pr['cSalt'  ], c <= cmax)] = pr['cSalt']
    for coll in caxC.collections:
        coll.remove()
    caxC = ax.contourf(gr.XM[:, 0, :], gr.ZM[:, 0, :],  np.ma.array(c, mask=(c == cInactive)),
                       levels=cLevels, cmap=cmap, alpha=1.0)
        
    mf6tools.show_animation_progress(frame, nbreak=50, ntot=len(kstpkper))
    
    return (ttl, caxP, caxC) # Return affected artists

ani = animation.FuncAnimation(fig=fig, func=update,                              
                              frames=len(concObj.times),
                              interval=120,
                              repeat=False)

if SAVE_ANIMATION:
    videoname = os.path.join(dirs.images, sim.name + '.mp4')
    print("Saving video to: {}".format(videoname))
    ani.save(filename = videoname, writer='ffmpeg')
    
plt.savefig(os.path.join(dirs.images, sim.name + '.png'))

plt.show()
