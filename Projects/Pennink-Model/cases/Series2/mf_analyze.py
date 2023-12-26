# %% Analyzing output of Pennik Series 2 tests
# TO 090223
# %%
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path
import flopy
import mf_adapt
from mf_adapt import gr
from mf6lab.mf6contourtools import get_contour_levels
from fdm.mf6_face_flows import get_struct_flows # (flowjas, grb_file=None, verbose=False)
from PIL import Image
import matplotlib.animation as animation
from matplotlib import colormaps

SAVE_ANIMATION = False

dirs = mf_adapt.dirs
grb_file = os.path.join(dirs.GWF, mf_adapt.sim_name + 'Gwf.dis.grb')

def addPatch(pgcoords, alpha=0.2, fc='none', ec='none', ax=None):
    codes = np.zeros(len(pgcoords), dtype=int) + path.Path.LINETO
    codes[0] = path.Path.MOVETO
    codes[-1] = path.Path.CLOSEPOLY
    pth = path.Path(pgcoords, codes)
    ptch = patches.PathPatch(pth, alpha=alpha, fc=fc, ec=ec)
    ax.add_patch(ptch)
    return

def get_psi_mask(h):
    """Return mask for plotting psi, masking inactive head cells.
    
    Parameters
    ----------
    h: np.ndarray
        heads (containgin 1e30 for inactive cells)

    Returns
    -------
    mask, shape = (h.shape[0] + 1, hshape[1] - 1)
    """
    inactive = 1e30
    h = h.copy()
    h[np.isnan(h)] = inactive
    h = np.minimum(h[:, 1:], h[:, :-1]) # vertical grid lines
    h = np.vstack((h[0], h, h[-1]))
    h = np.minimum(h[:-1, :], h[1:, :]) # horizontal  grid lines
    psi_mask = h == inactive
    return psi_mask
               
# %%
dirs.photos = (os.path.join(dirs.case, 'photos'))
foto = Image.open(os.path.join(dirs.photos, 'Series2_01_p30.jpg'))

# %% load the unformatted files with the heads, the concentration and the budget terms
sim = flopy.mf6.MFSimulation.load(
    sim_name=mf_adapt.sim_name, version='mf6', sim_ws=dirs.SIM)

gwf = sim.get_model('{}Gwf'.format(sim.name).lower()) # list(sim.model_names)[0])
gwt = sim.get_model('{}Gwt'.format(sim.name).lower())

headsObj = gwf.output.head()
budObj   = gwf.output.budget()
concObj  = gwt.output.concentration()

Qchd = [np.abs(q['q']).sum() / 2 for q in budObj.get_data(text='chd')]

kstpkper = np.array(headsObj.kstpkper[40]) - 1

hmin, hmax, hInactive = np.unique(headsObj.get_alldata())[[0, -2, -1]]
cmin, cmax, cInactive = np.unique(concObj.get_alldata())[[ 0, -2, -1]]

#  %% Steam function:
flowjas = budObj.get_data(text='FLOW-JA')

fflows = get_struct_flows([flowjas], grb_file=grb_file, verbose=False)

psi = gr.psi_row(fflows['frf'][0])
psi_min, psi_max = np.unique(psi)[[0, -1]]

crange = get_contour_levels(cmin, cmax, 20); crange = crange[crange > 0.0]
prange = get_contour_levels(*np.unique(psi)[[0, -1]], 20)
hrange = get_contour_levels(hmin, hmax, 20)

startdatetime = np.datetime64(sim.tdis.start_date_time.get_data().replace('t', ' '))
datetimes = startdatetime + np.array(concObj.times) * np.timedelta64(1, 'm')

fig, ax = plt.subplots()
fig.set_size_inches(12,12)
ttl = ax.set_title("Pennink {}".format(sim.name))
ax.grid(False)
ax.set_xlabel = 'x [cm]'
ax.set_ylabel = 'y [cm]'

extent = (-20, 85, -10, 82)
ax.imshow(foto, extent=extent)
ax.plot((0, 65, 65, 0, 0), (0, 0, 65, 65, 0), 'k')

if False:     
    p1 = patches.PathPatch(mf_adapt.sand, alpha=0.5, fc='brown', ec='none')
    p2 = patches.PathPatch(mf_adapt.canalL, alpha=0.5, fc='Blue', ec='black')
    p3 = patches.PathPatch(mf_adapt.canalR, alpha=0.2, fc='Gray', ec='black')

    for pgcoords, clr in zip([mf_adapt.sand, mf_adapt.canalL, mf_adapt.canalR],
                            ['brown', 'blue', 'green']):
        addPatch(pgcoords, alpha=0.5, fc=clr, ec='black', ax=ax)

h = headsObj.get_data(kstpkper=(0,0))[:, 0, :]
c = concObj.get_data(kstpkper=(0,0))[:, 0, :]

cmap = 'binary'

caxH = ax.contour( gr.xc, gr.zc,  np.ma.array(  h, mask=np.isnan(h)), levels=hrange)
caxP = ax.contour( gr.Xp, gr.Zpx, np.ma.array(psi, mask=get_psi_mask(h)), levels=prange) 
caxC = ax.contourf(gr.xc, gr.zc,  np.ma.array( c, mask=np.isnan(c)),
                   levels=crange,
                   cmap=cmap)

time_units = sim.tdis.time_units.get_data()

def update(frame):
    """Animation update function."""
    global caxH, caxP, caxC, Qchd
    ttl.set_text("Pennink {}, Frame={}, t = {} {}, Q = {:.2f} cm3/min"
                 .format(sim.name, frame,
                         datetimes[frame],
                         time_units,
                         Qchd[frame])
                 )
    kstpkper = np.array(concObj.kstpkper[frame]) - 1
    h = headsObj.get_data(kstpkper=kstpkper)[:, 0, :]
    c = concObj.get_data(kstpkper=kstpkper)[ :, 0, :]
    psi = gr.psi_row(fflows['frf'][frame])
    
    for coll in caxP.collections:
        coll.remove()
    caxP = ax.contour(gr.Xp, gr.Zpx, np.ma.array(psi, mask=get_psi_mask(h)), levels=prange)

    for coll in caxH.collections:
        coll.remove()
    caxH = ax.contour( gr.xc, gr.zc,  np.ma.array(  h, mask=np.isnan(h)), levels=hrange)

    for coll in caxC.collections:
        coll.remove()
    caxC = ax.contourf(gr.xc, gr.zc,  np.ma.array( c, mask=np.isnan(c)),
                       levels=crange, cmap=cmap)
    
    return (frame, caxH, caxP, caxC, ttl) # Return affected artists

ani = animation.FuncAnimation(fig=fig, func=update,                              
                              frames=len(concObj.times),
                              interval=30,
                              repeat=False)

if SAVE_ANIMATION:
    videoname = os.path.join(dirs.SIM, sim.name + '.mp4')
    print("Saving video to: {}".format(videoname))
    ani.save(filename = videoname, writer='ffmpeg')
else:
    plt.show()

