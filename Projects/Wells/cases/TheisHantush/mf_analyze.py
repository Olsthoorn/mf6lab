# === Analyzing output =====

# === General imports ======
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import flopy
from functools import partial

# ==== local imports =====
from src import mf6tools
import mf_adapt
from etc import newfig
from src.mf6contourtools import get_contour_levels
from fdm.mf6_face_flows import get_structured_flows_as_dict
from settings import props as pr

SAVE_ANIMATION = False

# === From mf_adapt =====
dirs = mf_adapt.dirs
sim_name = mf_adapt.sim_name
section_name = mf_adapt.section_name
gr = mf_adapt.gr
Iz = mf_adapt.Iz # Layer number of top active cells
start_date_time = mf_adapt.start_date_time
lay = mf_adapt.lay # layer definictions from Excel workbook sheet Lay

# === load simulation =====
sim = flopy.mf6.MFSimulation.load(sim_name=sim_name,
                                  version='mf6',
                                  sim_ws=dirs.SIM)

# === Groundwater FLow Object =====
gwf = sim.get_model('{}Gwf'.format(sim.name).lower()) # list(sim.model_names)[0])
headsObj = gwf.output.head()
budObj   = gwf.output.budget()

kstpkper  = headsObj.get_kstpkper()
sim_times = headsObj.get_times()

#datetimes = np.datetime64(start_date_time) + np.array(headsObj.times) * np.timedelta64(1, 'D')

# === Get strcutured flows flf, fff and frf first (we only need frf) ====
grb_file = os.path.join(dirs.GWF, sim_name + 'Gwf.dis.grb') # Flopy's grid file
fflows = get_structured_flows_as_dict(budObj, grb_file=grb_file)

#  === Stream function ===== Stream function is called psi.
psi = gr.psi_row(fflows[kstpkper[-1]]['frf'], row=0)
P =  abs(pr['Q'])    # estimate extremes of the stream function
p_levels = get_contour_levels(-P, P, 100)
dpsi = np.diff(p_levels)[0]

# === Get h_levels =====
hmin, hmax = np.unique(headsObj.get_alldata())[[0, -1]]
h_levels = get_contour_levels(-5, hmax, n=51)


# === New figure =====
ax = newfig(mf_adapt.section_name, 'distance r [m]', 'elevation [m]', xscale='linear', xlim=(0, 10000), figsize=(15, 8))
fig = plt.gcf()



def init_func(ax=ax, title="", gr=gr, lay=lay):
    
    global ttl, caxp, caxh, wt, hdlines
    
    ttl = ax.set_title(title)
    
    # === plot the water table =====
    wt, = ax.plot(gr.xm, np.zeros_like(gr.xm), 'b-', lw=1, label='water tafel')

    # === head curves ======
    hdlines = []
    for ilay in range(gr.nlay):
        hdlines.append(ax.plot(gr.xm, np.zeros_like(gr.xm), label='head in layer {}'.format(ilay))[0])
      
    # === contours ======  
    caxh = ax.contour([[0, 1], [0, 1]], [[0, 0], [1, 1]], [[0, 1], [2, 3]])
    caxp = ax.contour([[0, 1], [0, 1]], [[0, 0], [1, 1]], [[0, 1], [2, 3]])

    # for ca in caxh.collections:
    #     ca.set_animated = True
    # for ca in caxp.collections:
    #     ca.set_animated = True
        
    # === Plot the layers ===== (we dont have to do this)
    for zp in gr.Z:
        ax.plot(gr.xm, zp[0], 'k-', lw=0.5)

    # === Plot colored layers (using layer patches( =====
    layer_patches = gr.layer_patches_x(row=0) # Get layer patches

    # === Colors, codes and names of the layers =====
    colors = lay['Color'].values
    codes  = lay['Code'].values
    names  = lay['Name'].values

    # === set patch properties =====
    for p, clr, code, name in zip(layer_patches, colors, codes, names):
        p.set_fc(clr)
        p.set_ec('k')
        p.set_lw(0.25)
        p.set_alpha(1.0)
        p.set_label(code + ' ' + name)
        ax.add_patch(p)

    # === put legend at a suitable place =====
    ax.legend(loc='lower right')

    # === Logo at lower richt of figure =====
    ax.text(0.85, 0.05, str(np.datetime64('today')),
            fontsize=10,
            fontweight='normal',
            va='bottom',
            transform=plt.gcf().transFigure)
    
    return ttl, caxp, caxh, wt, hdlines

def update(frame, budObj=None, headsObj=None, kstpkper=None, sim_times=None, p_levels=None, h_levels=None, Iz=None):
                
    global ttl, caxp, caxh, wt, hdlines, section_name, props
    # def update(frame, budObj=None, headsObj=None, kstpkper=None, p_levels=None, h_levels=None, Iz=None,
    #            ttl=None, caxp=None, caxh=None, wt=None, hdlines=None):
    """Animation update function."""    
    # === Update title =====
    Qactual = budObj.get_data(text='wel', kstpkper=kstpkper[frame])[0]['q'].sum()
    ttl.set_text(". {} dpsi = {} m3/d, Q actual = {:.1f} m3/d, frame={}, sim_time = {} d".format(section_name, dpsi, Qactual, frame, sim_times[frame]))

    # === water table =====
    hds = headsObj.get_data(kstpkper=kstpkper[frame])
    hwt = np.array([hds[iz, 0, ix] for ix, iz in enumerate(Iz)]) # Water table
    wt.set(ydata=hwt)
    
    # === layer head lines =====
    for hdline, hLayer in zip(hdlines,  hds[:, 0, :]):
        hdline.set_ydata(hLayer)
        
    # ==== Contour the stream lines =====
    flow_right_face = fflows[kstpkper[frame]]['frf'][0]
    psi = gr.psi_row(flow_right_face)
    for coll in caxp.collections:
        coll.remove()
    Zpx = gr.Zpx_limited_to_water_table(hwt) # Adapt psi contour grid to water table
    caxp = ax.contour(gr.Xp, Zpx, psi, levels=p_levels, linewidths=0.5)

    # === Contour heads =====
    flow_lower_face = fflows[kstpkper[frame]]['flf'][0]
    hZ = gr.heads_on_zplanes(hds, flow_lower_face, gr.const(mf_adapt.lay['k33'].values))
    for coll in caxh.collections:
        coll.remove()
    
    caxh = ax.contour(gr.XM_on_zplanes[:, 0, :], gr.Z_limited_to_water_table(hwt),  hZ,
                       levels=h_levels, colors='k')
    
    #for ca in caxp.collections:
    #    ca.set_animated = True
    #for ca in caxh.collections:
    #    ca.set_animated = True


    # === label heads of last time step =====
    if frame + 1 == len(kstpkper):
        ax.clabel(caxh, h_levels, inline=True, fmt=r'%2.2f', fontsize=10) 

    # === Show progress =====
    mf6tools.show_animation_progress(frame, nbreak=50, ntot=len(kstpkper))
    
    return ttl, caxp, caxh, wt, hdlines


ani = animation.FuncAnimation(fig,
                              partial(update, budObj=budObj, headsObj=headsObj, kstpkper=kstpkper, p_levels=p_levels, h_levels=h_levels, Iz=Iz, sim_times=sim_times),
                              frames=range(len(kstpkper)),
                              init_func=init_func,
                              fargs=None,
                              interval=200,
                              save_count=None,
                              repeat=False,
                              blit=False)

if SAVE_ANIMATION:
    videoname = os.path.join(dirs.images, sim.name + '.mp4')
    print("\nCounting frames to video in the makeing: {}\n".format(videoname))
    print("Saving video to: {}".format(videoname))
    ani.save(filename=videoname, writer='ffmpeg')

plt.savefig(os.path.join(dirs.images, sim.name + '.png'))

# === Finish plot and show it =====
plt.show()
# %%
