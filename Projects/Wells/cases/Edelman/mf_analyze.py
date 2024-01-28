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

SAVE_ANIMATION = True

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
P =  abs(pr['Q'])    # estimate extremes of the stream function
p_levels = get_contour_levels(-P, P, 51)
dpsi = np.diff(p_levels)[0]

# === Get h_levels =====
hmin, hmax = np.unique(headsObj.get_alldata())[[0, -1]]
h_levels = get_contour_levels(hmin, hmax, n=21)


def init_func(ax=None, gr=None, lay=None,
            fflows=None, headsObj=None,
            kstpkper=None, p_levels=None, h_levels=None, Iz=None, sim_times=None
            ):
    
    global ttl, caxp, caxh, wt, hdlines, section_name, props
    
    frame = 0
    Qactual = 0.0
    title = (". {} dpsi = {:.3f} m2/d, Q actual = {:.1f} m2/d, frame={}, sim_time = {:.1f} d"
             .format(section_name, dpsi, Qactual, frame, sim_times[frame]))
    ttl = ax.set_title(title)
    
    # === plot the water table =====
    
    wt, = ax.plot(gr.xm, np.zeros_like(gr.xm), 'b-', lw=1, label='water tafel')

    # === head curves ======
    hdlines = []
    for ilay in range(gr.nlay):
        hdlines.append(ax.plot(gr.xm, np.zeros_like(gr.xm), label='head in layer {}'.format(ilay))[0])
    
    # === contours ======      
    hds = headsObj.get_data(kstpkper=kstpkper[0])
    hwt = np.array([hds[iz, 0, ix] for ix, iz in enumerate(Iz)]) # Water table

    flow_right_face = fflows[kstpkper[0]]['frf'][0]
    psi = gr.psi_row(flow_right_face)

    Zpx = gr.Zpx_limited_to_water_table(hwt) # Adapt psi contour grid to water table
    caxp = ax.contour(gr.Xp, Zpx, psi, levels=p_levels, linewidths=0.5)
        
    flow_lower_face = fflows[kstpkper[0]]['flf'][0]
    hZ = gr.heads_on_zplanes(hds,
                             flow_lower_face,
                             gr.const(mf_adapt.lay['k33'].values),
                             continuous=False)
    caxh = ax.contour(gr.XM_on_zplanes[:, 0, :],
                gr.Z_limited_to_water_table(hwt.reshape(gr.shape[1:]))[:, 0, :],
                hZ[:, 0, :],
                levels=h_levels,
                colors='k')
        
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

def animate(frame, budObj=None, fflows=None, headsObj=None, kstpkper=None,
            sim_times=None, p_levels=None, h_levels=None, Iz=None):
                
    global ttl, caxp, caxh, wt, hdlines, section_name, props
    """Animation update function."""
      
    # === Update title =====
    Qactual = budObj.get_data(text='wel', kstpkper=kstpkper[frame])[0]['q'].sum()
    ttl.set_text(". {} dpsi = {:.3f} m2/d, Q actual = {:.1f} m2/d, frame={}, sim_time = {:.1f} d".format(section_name, dpsi, Qactual, frame, sim_times[frame]))

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
    hZ = gr.heads_on_zplanes(hds, flow_lower_face,
                             gr.const(mf_adapt.lay['k33'].values),
                             continuous=False)
    for coll in caxh.collections:
        coll.remove()
    
    caxh = ax.contour(gr.XM_on_zplanes[:, 0, :],
                gr.Z_limited_to_water_table(hwt.reshape(gr.shape[1:]))[:, 0, :],
                hZ[:, 0, :],
                levels=h_levels,
                colors='k')
    
    # === label heads of last time step =====
    if frame + 1 == len(kstpkper):
        ax.clabel(caxh, h_levels, inline=True, fmt=r'%2.2f', fontsize=10) 

    # === Show progress =====
    mf6tools.show_animation_progress(frame, nbreak=50, ntot=len(kstpkper))
    
    return ttl, caxp, caxh, wt, hdlines


if True:

    # === New figure =====
    ax = newfig(mf_adapt.section_name, 'distance r [m]', 'elevation [m]', xscale='linear',
                xlim=(0, 0.8 * pr['L']), figsize=(15, 8))
    fig = plt.gcf()


    ani = animation.FuncAnimation(fig,
                                partial(animate, budObj=budObj, fflows=fflows, headsObj=headsObj,
                                        kstpkper=kstpkper, p_levels=p_levels, h_levels=h_levels, Iz=Iz, sim_times=sim_times),
                                frames=range(len(kstpkper)),
                                init_func=
                                partial(init_func, ax=ax, gr=gr, lay=lay,
                                        fflows=fflows, headsObj=headsObj,
                                        kstpkper=kstpkper,
                                        p_levels=p_levels,
                                        h_levels=h_levels,
                                        Iz=Iz,
                                        sim_times=sim_times),
                                fargs=None,
                                interval=200,
                                save_count=None,
                                repeat=False,
                                blit=False,
                                )

    if SAVE_ANIMATION:
        videoname = os.path.join(dirs.images, sim.name + '.mp4')
        print("\nCounting frames to video in the makeing: {}\n".format(videoname))
        print("Saving video to: {}".format(videoname))
        ani.save(filename=videoname, writer='ffmpeg')

    plt.savefig(os.path.join(dirs.images, sim.name + '.png'))

    # === Finish plot and show it =====
    plt.show()


if __name__ == '__main__':
    
    def checkFlows(headsObj=None, fflows=None, lay=None, istep=-1):
        
        kstpkper = headsObj.get_kstpkper()
        ksp = kstpkper[istep]
        
        FLF = fflows[ksp]['flf'][0]
        
        assert np.all(FLF[-1].ravel() == 0.), "All flf[-1] must be zero! Sum is {} m3/d".format(FLF[-1].sum())
        
        qvMF6 = -FLF[:-1] / gr.Area[np.newaxis, : ,:] # upward postivie
        
        # Hand compute the flow through the tops and bottoms
        hds = headsObj.get_data(kstpkper=ksp)
        k33 = gr.const(lay['k33'].values)
        cCell = gr.DZ / k33
        c = (cCell[:-1] + cCell[1:]) / 2.
        qvHand = np.diff(hds, axis=0) / c
        
        print("Comparing the average flf computed by Modflow 6 and by 'hand':")
        print("mean qv by MF6  ={:6f} m/d".format(qvMF6.mean()))
        print("mean qv by hand ={:6f} m/d".format(qvHand.mean()))
        print("This verifies the Flow-Loewer-Face implementation.")
        print()
        

    def checkHeadComputation(headsObj=None, fflows=None, lay=None, istep=-1):
        """Check the computation of heads at top and bottom of the model layers.
        
        This head computation is used for head contouring such that the knick occurs as the
        laye boundaries and not in the cell mids, for which the head was computed by MF6.
        
        Contouring with heads at op and bottom of the cells is appropriate because that is
        where the layer properties change.
        
        htop = hcenter + int_zmid^ztop dh/dz dz;  dh/dz = -q / k33
        hbot = hcenter - int_zmid^zbot dh/dz dz:  dh/dz = -q / k33
        
        Compute htop and hbot from two sides.
        First is that qz is uniform above and below zmid but different.
        Second is that qz varies linearly over the thickness of the cell from zbot to ztop.
        """
        kstpkper = headsObj.get_kstpkper()
        ksp = kstpkper[istep]
        
        FLF = fflows[ksp]['flf'][0]
        
        assert np.all(FLF[-1].ravel() == 0.), "All flf[-1] must be zero! Sum is {} m3/d".format(FLF[-1].sum())

        # Hand compute the flow through the tops and bottoms
        hds = headsObj.get_data(kstpkper=ksp) # gr.shape
        k33 = gr.const(lay['k33'].values)     # gr.shape
        cCell = gr.DZ / k33                   # gr.shape
        cFace = (cCell[:-1] + cCell[1:]) / 2

        # Resistance across layer boundaries LB
        # cLB = (cCell[:-1] + cCell[1:]) / 2.  # (nlay - 1, nrow, ncol)
        # qvHand = -np.diff(hds, axis=0) / cLB # (nlay - 1, nrow, ncol)
        
        # qv Lower Face, upward = positive!
        qvLF = -FLF[:-1] / gr.Area[np.newaxis, : ,:] # (nlay -1, nrow, ncol), upward positive

        # Head difference between layer centers
        dh = (-qvLF * cCell[:-1]
              -qvLF * cCell[ 1:]) / 2
        
        dh = -qvLF * cFace
        
        print(
        """"
        Comparing
        heads [1:] computed from hds[:-1] and qv with
        heads [:1] computed from hds[ 1:] and qv
        {}
        {}""".format(
            (hds[ 1:] + dh).mean(axis=1).mean(axis=-1),
            (hds[:-1] - dh).mean(axis=1).mean(axis=-1)
            )
        )
        
        # qv at top and bottom of the layers (vertical = positive)
        qtop = gr.const(0.) # gr.shape
        qbot = gr.const(0.) # gr.shape
        qtop[ 1:] = qvLF   # gr.shape, top zeros
        qbot[:-1] = qvLF   # gr.shape, bot zeros
        
        print(
        """"
        Comparing
        heads at top from hds and qtop with
        heads at bot from hds and qbot with qv constant above and below layer centers
        {}
        {}
        """.format(            
            (hds - qtop * cCell / 2).mean(axis=1).mean(axis=-1),
            (hds + qbot * cCell / 2).mean(axis=1).mean(axis=-1)
            )
        )

        print(
        """"
        Comparing
        heads at top from hds and qtop with
        heads at bot from hds and qbot with qv varying linearly across layers
        {}
        {}
        """.format(            
            (hds - (3 * qtop + 1 * qbot) * cCell / 8).mean(axis=1).mean(axis=-1),
            (hds + (1 * qtop + 3 * qbot) * cCell / 8).mean(axis=1).mean(axis=-1)
            )
        )
        print()

    checkFlows(headsObj=headsObj, fflows=fflows, lay=lay, istep=-1)
    checkHeadComputation(headsObj=headsObj, fflows=fflows, lay=lay, istep=-1)
