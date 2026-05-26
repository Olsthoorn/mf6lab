# %% === Generatting modelled x-sections through the ARK canal 

# %% === Imports

import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt

from itertools import cycle
from typing import Any
from glob import glob
from pprint import pprint
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Patch, PathPatch, Path

from tools.fdm.src.mfgrid import Grid
from tools.fdm.src.fdm3Blom import Fdm3
from tools.etc.etc import logo

from mf6lab.Projects.ARK_RWS.src.ARK_geotop import Dirs


print(sys.executable)

# --- Needed to make figure separate from the notebook and interactive
# matplotlib qt

# --- Notebook name for logo
NOTEBOOK_NAME = "ARK_fdm.ipynb"
# --- Seet the namespace for the relevant directories
dirs = Dirs()

# --- Get the paths and names of  the geotop pdf files in the order they are in dirs.dino
xsec_paths = {i:name for i, name in enumerate(glob(dirs.dino + '*.pdf'))}
xsec_names = {i:os.path.basename(name) for i, name in enumerate(glob(dirs.dino + '*.pdf'))}

# %%

# --- Pickling
def pickleto(var:Any, basename:str, parent:str=dirs.data):
    """Pickle var to os.path.join(dirs.data, basename)"""
    if not basename.endswith('.pkl'):
        basename += ".pkl"

    pkl_file = os.path.join(parent, basename)
    with open(pkl_file, 'wb') as f:
        print(f"Pickled {basename} --> {parent}")        
        pickle.dump(var, f)

# --- Unpickling
def picklefrom(basename:str, parent:str=dirs.data)->Any:
    """Unpickle varname from os.path.join(parent, basename)"""
    if not basename.endswith(".pkl"):
        basename += ".pkl"
          
    pkl_file = os.path.join(parent, basename)
    with open(pkl_file, 'rb') as f:
        print(f"Loaded {basename} <-- {parent}")        
        return pickle.load(f)


def spy(idx_arr):
    """Show where the index labels are in the xsec idx_arr."""
    arr = np.squeeze(idx_arr)
        
    fig, ax = plt.subplots(figsize=(10, 6))
    if np.issubdtype(idx_arr.dtype, np.integer):
        title = "Location of legend indices in xsec.idx_arr"
        classes = np.unique(arr)
        cmap = plt.get_cmap('tab20', len(classes) - 1)
    else:
        title = "imshow of given array"
        cmap = plt.get_cmap('viridis')
    ax.set_title(title)
    mappable = ax.imshow(arr, cmap=cmap, origin='upper')
    fig.colorbar(mappable)
    plt.show()

def extent_patch(extent, **kwargs):
    """Return a rectangle patch object.
    
    Parameters
    ----------
    extent: np.array | tuple | list of 4 floats
        xmin, xmax, ymin, ymax
    kwargs: dict
        more parameter passed on to the Rectangle object.
    
    use it as
    as.add_patch(rect)
    """
    xmin, xmax, ymin, ymax = extent
    return Rectangle((xmin, ymin), xmax - xmin, ymax - ymin, **kwargs)
    
def find_horizonal_intersections(cs, z0):
    """Return horizontal contourline intersections at z0.
    
    Parameters
    ----------
    cs: contour set (returned by plt.contour or ax.contour)
        Contour set with all contour information.
    z0: float
        z-location where contour labels are to be placed along horiontal.
    """
    points = []

    for path in cs.get_paths():
        v = path.vertices
        x, z = v[:, 0], v[:, 1]
        if len(x) < 2:
            continue
        dz0 = z - z0
        mask = dz0[:-1] * dz0[1:] <=0
        idz = np.where(mask)[0]
        if idz.size > 0:
            for i in idz:
                t = (z0 - z[i]) / (z[i+1] - z[i] + 1e-12)
                xi = x[i] + t * (x[i+1] - x[i])
                points.append((xi, z0))
    return points

def find_vertical_intersections(cs, x0):
    """Return vertical contourline intersections at x0.
    
    Parameters
    ----------
    cs: contour set (returned by plt.contour or ax.contour)
        Contour set with all contour information.
    x0: float
        x-location where contour labels are to be place along vertical.
    """
    points = []

    for path in cs.get_paths():
        v = path.vertices
        x, z = v[:, 0], v[:, 1]
        dx0 = x - x0
        mask = dx0[:-1] * dx0[1:] <=0
        idx = np.where(mask)[0]
        if idx.size > 0:
            i = idx[0]
            t = (x0 - x[i]) / (x[i+1] - x[i] + 1e-12)
            zi = z[i] + t * (z[i+1] - z[i])
            if zi < -45. or zi > -5.:
                continue
            points.append((x0, zi))
    return points

# --- Color for empty legend (empty voxel with geo_unit 'none')
WHITE_01 = np.array([1., 1., 1.])

# %% # Thickness and hydraulic resistance of canal bottom

def D_ARK_bot(xsec):
    """Return thickness of resistance layer at canal bottom.
    
    Parameters
    ----------
    xsec: xsec object
        xsec data
    Returns
    -------
    Resistance layer bottom thickness ad cell mids at canal bottom
        
    Uses xsec.c_ARK_bot_dict dictionary with fields
    ------------------------------------------
        D0: float
            Initial resistance layer thickness
        fd: float
            abrased fraction of may be > 1)
    """
    D0 = xsec.c_ARK_bot_dict['D0']
    fd = xsec.c_ARK_bot_dict['fd']

    gr = xsec.gr
    xm = gr.XM[gr.mask_ARK_bot]
    b = xsec.b
    D = D0 + fd * D0 / 2 * (np.cos(2 * np.pi * xm/b) - 1)    
    return D.clip(1e-3, None)

def set_kv_ARK_bot(xsec):
    """Return set_kv_canal_bottom by adjusting kv ARK-bottom
    
    Parameters
    ----------
    xsec: xsec class object with all its data.
        The cross section data
        
    uses xsec.c_ARK_bot_dict dictionary with fields
    ------------------------------------------
        D0: float
            initial thickness of bottom resistance layer
        fd: float
            Thicknes sreduction factor.
        c: float
            Specific resistance of canal bottom [d/m]

    Returns
    -------
    Nothing, resets gr.kv at bottom of ARK.
    """
    c  = xsec.c_ARK_bot_dict['c']
    
    D = D_ARK_bot(xsec)

    gr = xsec.gr
    mask = gr.mask_ARK_bot    

    cbot = c * D
    gr.kv[mask] = 0.5 * gr.DZ[mask] / cbot
    return None

def  get_kv_ARK_bot(xsec):
    """Return kv of ARK bottom."""
    gr = xsec.gr
    return gr.kv[gr.mask_ARK_bot]

def set_c_ARK_bot(xsec):
    """Set the vertical resistance of the ARK bottom.
    Parameters
    ----------
    xsec: xsec class object with all its data.
        The cross section data
    D0: float
        initial thickness of bottom resistance layer
    fd: float
        Thicknes sreduction factor.
    c: float
        Specific resistance of canal bottom [d/m]

    Returns
    -------
    Nothing, resets gr.kv at bottom of ARK.    
    """
    
    set_kv_ARK_bot(xsec)
    return None

def get_c_ARK_bot(xsec):
    """Return resistance of ARK bottom."""
    gr = xsec.gr
    kv = gr.kv
    return gr.DZ[gr.mask_ARK_bot] / kv[gr.mask_ARK_bot]


def c_bot_patch(xsec):
    """Return a patch showing the undulating resistance layer at the ARK-bottom"""
    gr = xsec.gr
    D0 = xsec.c_ARK_bot_dict['D0']
    mask = gr.mask_ARK_bot
    
    D = np.hstack((D0, D_ARK_bot(xsec), D0))
    z = gr.ZM[mask] - 0.5 * gr.DZ[mask]
    zb = np.hstack((z[0], z, z[-1]))
    zt = zb + D
    x = np.hstack((-xsec.b, gr.XM[mask], xsec.b))
    pts = np.vstack((
        np.vstack((x, zb)).T,
        np.vstack((x[::-1], zt[::-1])).T,
        [x[0], zb[0]]
    ))
    codes = np.zeros(len(pts), dtype=int)
    codes[:] = Path.LINETO
    codes[0] = Path.MOVETO
    codes[-1] = Path.CLOSEPOLY
    pth = Path(pts, codes=codes)
    return PathPatch(pth, fc='darkgray', ec='none')    
    
# %% --- Get pickled x-sections

geotop_xsecs = picklefrom("geotop_xsecs.pkl")

# %% Deal with the second cross section only

isec = 1
xsec = geotop_xsecs[xsec_names[isec]]

print(f"Dealing with xsec {isec}:\n{xsec.name}") 

# --- find the ARK it wronly has index 1 in row 6
ix_ARK = np.where(xsec.idx_arr[6] == 1)[0]

spy(xsec.idx_arr)

# %% === Repair the idx_arr for the ARK canal which looks now filled with material 'a'

idx_arr = xsec.idx_arr.copy()

# --- The columns have idx 1 incorrectly
cols = np.where(idx_arr[1] == 1)[0]

# --- Replace by index 0 ('none')
for j in cols:
    rows = idx_arr[:, j] == 1
    idx_arr[rows, j] = 0
    
# --- Check
spy(idx_arr)

# %% === When this works replace the xsec.idx_arr with the corrected version

# --- Replace origional idx_arr by corrected one
xsec.idx_arr = idx_arr

# --- Check
spy(xsec.idx_arr)
print('idx_arr repaired')

# %% Model grid

def set_ARK_xsec(xsec):
    # --- Basic properties
    xsec.ground_elev = -1.3
    xsec.stage  = -0.4  # --- Water level of ARK    
    xsec.d_damw = 0.5   # --- Thickness of sheet piling
    xsec.z_damw = -12.  # --- Bottom elevation of sheet piling
    xsec.hpp    = xsec.ground_elev - 1 # --- Water level in ditches
    xsec.c_drainage = 100. # --- Areal drainage resistance (phi - hpp) = Nc about 0.1 m
    
    # --- Define resistance undulating abrased ARK-bottom
    # D0: initial D, fd: D-fraction erased, c: material resistance d/m
    xsec.c_ARK_bot_dict = dict(D0=1., fd=1., c=100)    

    
    # --- Add legend color and rhow (wet density) to geo_units
    for (_, gu), color in zip(xsec.geo_units.items(), xsec.leg_colors[1:]):
        gu['color'] = color
        gu['rhow'] = gu['n'] * 1000 + (1 - gu['n']) * gu['rho']
            
    # --- where is the ARK in the cross section?
    ixARK = np.where(xsec.idx_arr[10] == 0)[0][0]
    izARK = np.where(xsec.idx_arr[:, ixARK] == 0)[0]
    
    # --- Estimate the canal bottom at dz/2 below the lowest ARK cell center
    xsec.zbot = np.round(xsec.zm[izARK[-1]] - xsec.dz / 2, 1)
    
    # --- ARK water body extent
    ae = np.round(np.array([xsec.x[ixARK], xsec.x[ixARK+1], xsec.zbot, xsec.stage]), 1)
    
    # --- Half width of ARK:
    xsec.b = (ae[1] - ae[0]) / 2

    # --- Centralize xsec.x around xARKmid
    xsec.xARKmid_orig     = 0.5 * (ae[0] + ae[1])
    t = (xsec.xARKmid_orig - xsec.x[0]) / (xsec.x[-1] - xsec.x[0])
    xsec.xy_hart_ARK=np.round(xsec.xyRD[0] + t * np.diff(xsec.xyRD, axis=0))
    
    # --- Keep original xsec.world_extent
    xsec.world_extent_orig = xsec.world_extent
    
    # --- Centralize the x-coordinates around the heart line of the canal (xARKmid_orig)
    xsec.world_extent[:2] -= xsec.xARKmid_orig
    
    # --- Refine away from the edges of the canal  
    slog = np.logspace(0, 2, 10)

    # --- Use b as shorthand of xsec.b
    b = xsec.b
    
    # --- Define cell face x - coordinates between heart line and right side of the canal.
    xc = b - np.hstack((1.0, 1.5, 2.5, np.linspace(5, 45, 9), b))[::-1]

    # --- All x grid coordinates from the heart line of the ARK to the right (x>0)
    x_ = np.hstack((
        xc,                                       # --- inside ARK, right of middle
        b, b + xsec.d_damw,                       # --- sheet piling
        b + slog[slog > xsec.d_damw],             # --- increasing cell wirdth first 100 m 
        np.linspace(0, 2000, 21).clip(200, None)  # --- beyond this to 2000 m 100 m cells
    ))

    # --- For less clutter, round the x-values
    x_ = np.round(x_, 1)

    # --- Mirror around heart line of ARK and remove doubles   
    x = np.unique(np.hstack((-x_[::-1], x_)))
    
    # --- Geneate the model grid, using the new x and the old z below ground surface.
    gr = Grid(x, None, xsec.z[xsec.z <= xsec.ground_elev])
    
    # --- Generate mask arrays to easily get cell id's later.
    # --- mask_ARK water body
    gr.mask_ARK = gr.inblock(xx=(-b, b), yy=None, zz=(xsec.zbot, xsec.stage))
    
    # --- Canal bottom resistance layer (take the lowerst layer of ARK water body
    iARK = np.where(gr.mask_ARK[0, 0, :])[0]           # --- Any ix inside the canal
    izARK= np.where(gr.mask_ARK[:, 0, iARK[0]])[0][-1] #
    
    # --- Get the mask of the bottom layer of the canal water body
    gr.mask_ARK_bot = gr.const(0, dtype=bool)
    gr.mask_ARK_bot[izARK, 0, iARK] = True

    # --- extent of the sheet piling left and right along canal      
    gr.mask_damwL = gr.inblock(xx=(-xsec.b - xsec.d_damw, -xsec.b), zz=(xsec.z_damw, 0))
    gr.mask_damwR = gr.inblock(xx=(+xsec.b, +xsec.b + xsec.d_damw), zz=(xsec.z_damw, 0))
    
    # --- Also get extents
    gr.ARK_extent = np.array([-b, b, xsec.zbot, xsec.stage])
    gr.damwL_extent = np.array([-b, -b - xsec.d_damw, xsec.z_damw, 0])
    gr.damwR_extent = np.array([b, b + xsec.d_damw, xsec.z_damw, 0])

    # --- Top boundary condition outside the ARK canal.
    # --- We use DRN as top boundary condition.
    # --- There is neither evaporation, nor recharge, just seepage from the ARK. So DRN is ok.
    # --- Drains are in the cells with (hpp - 0.25 <= zm <= hpp + 0.25) and outside ARK
    # --- hpp is surface water level outside the canal (Dutch: polder peil)
    gr.mask_DRN = np.logical_and(
        gr.inblock(xx=(gr.x[0], gr.x[-1]), zz=(xsec.hpp - 0.25, xsec.hpp + 0.25)),
        ~gr.mask_ARK
        )

    # --- DRN cell Id's 
    Idrn = gr.NOD[gr.mask_DRN]

    # --- DRN boundary condition
    DRN = np.zeros(len(Idrn), dtype=Fdm3.dtypes['drn'])
    DRN['Ig'] = Idrn
    DRN['C'] = xsec.c_drainage
    DRN['h'] = xsec.hpp
    gr.DRN = DRN

    # --- Other grid Arrays for the model
    IBOUND = gr.const(1, dtype=int)
    IBOUND[gr.mask_damwL] =  0 # --- Inactive
    IBOUND[gr.mask_damwR] =  0 # --- Inactive
    IBOUND[gr.mask_ARK]   = -1 # --- Fixed head
    gr.IBOUND = IBOUND
    
    # --- No fixed Q
    gr.FQ = gr.const(0.)
    
    # --- No ohter fixed heads
    gr.FH = None
    
    # --- Initial heads are all hpp except where ARK cuts into the model grid.
    HI = gr.const(xsec.hpp)
    HI[gr.mask_ARK] = xsec.stage
    gr.HI = HI
    
    # === Filling the hydraulic properties using the legend id of each cell.
    # --- We do this by overlapping the xsec-geotop-grid by the model grid and get the index-id
    # --- for each model cell at once.
    gr.idx_arr = xsec.overlap(gr)
    
    # --- Get the hydraulic properties for the model grid
    # --- The shape of props arrays will be equal to gr.shape
    props  = xsec.get_props(idx_arr=gr.idx_arr)
    
    # --- Make the grid propertie arrays 3D to match the 3D grid of the model, even if ny=1
    for var in props:
        props[var] = props[var][:, np.newaxis, :]
    
    # --- The ARK-water body --> Give it a  very high vertical conductivity
    # --- The horizontal conductivity does not matter here. Leave it zero.
    props['kv'][gr.mask_ARK] = 1000.
    
    # --- Get the wet density of all cells in the waterbody.    
    props['rho_wet'][gr.mask_ARK] = 1000.  

    # --- Add gr property arrays to grid object for later reference
    gr.kh = props['kh']
    gr.kv = props['kv']
    gr.K   = (gr.kh, gr.kh, gr.kv)
    gr.n   = props['n']
    gr.rho = props['rho']
    gr.rhow = props['rho_wet']

    # --- We don't really need S because the problem is steady state
    # --- However, the steady-state solution is appoached transiently for stability reasons.
    gr.S = gr.const(0.001)
    
    # --- Set distance between head lines and stream lines for contouring results:
    gr.dphi = 0.1   # --- m   head steps
    gr.dpsi = 0.1   # --- m2/d steps in the stream function (between pairs of stream lines)

    # --- Add gr object to xsec. So that xsec carries all that is necessary to model it
    xsec.gr = gr
    
    # --- Set the resistance of the bottom of the ARK
    # --- This uses xsec.c_ARK_bot_dict to adapt the lowest row of the canal water body cells.
    set_c_ARK_bot(xsec)
    
    return xsec


def simulate_fdm3(xsec):
    """Setup the fdm3 model, simulate and return results."""
    
    gr = xsec.gr
    
    # --- Instantiate the model
    mdl = Fdm3(gr=xsec.gr, K=gr.K, c=None, S=None, IBOUND=gr.IBOUND, HI=gr.HI, FQ=gr.FQ)

    # --- Run the model
    out = mdl.simulate(DRN=gr.DRN, RIV=None, GHB=None, FDR=None, tm=None, htol=1e-7, maxiter=50, verbose=False)

    # --- Compute the stream function from out['Qx]
    out['psi'] = gr.psi_row(out['Qx'])
    return out

# %% === Plot the cross section

def plot_setup(xsec, xlim=None):
    # --- Setup the figure. Adjust width to get a good size of the xsec.
    fig, ax = plt.subplots(figsize=(13, 6))

    # --- Adjust the right edge of teh axes to allow space for the legend's xbox
    fig.subplots_adjust(right=0.75)

    # --- Make axes title
    gr = xsec.gr
    _cd = xsec.c_ARK_bot_dict
    ttl = (
        f"Stijghoogten, stroomlijnen en drukoverschot [dPhi={gr.dphi} m dPsi={gr.dpsi} m2/d]" +
        f", hartlijn ARK ={xsec.xy_hart_ARK[0]}" + '\n' +    
        fr"Weerstandlaag kanaalbodem: $D_0$={_cd['D0']} m, " +
        f"uitschuurfractie={100*_cd['fd']}%, " +
        f"spec. weerst. ARK-bodem = {_cd['c']} d/m"
    )
    # --- place titles
    fig.suptitle(xsec.name)
    ax.set_title(ttl, fontsize=10)
    ax.set_xlabel('x van hartlijn ARK')
    ax.set_ylabel('z [NAP]')
    ax.grid()
    
    # --- Add logo to reference the picture in the future
    logo(fig, NOTEBOOK_NAME)
    
    ax.set_xlim(xlim)
    return fig, ax
    
    
def get_phi_and_psi_levels(xsec, out=None):
    """Return head and stream function levels."""
    gr = xsec.gr    
    psi = out['psi']
    # --- phi and psi levels for contourin heads and stream lines    
    phi_levels = np.arange(np.floor(xsec.hpp), np.ceil(xsec.stage), gr.dphi)
    psi_levels = np.arange(np.floor(psi.min()), np.ceil(psi.max()), gr.dpsi)
    return phi_levels, psi_levels


def add_head_contours_and_stream_lines(xsec, ax=None, out=None):
    """compute and put the head contours and the stream lines on the axes."""
    gr = xsec.gr
    psi = out['psi']
    # --- Contour heads and streamlines
    phi_levels, psi_levels = get_phi_and_psi_levels(xsec, out=out)
    cs_phi = ax.contour(gr.xm, gr.zm, out['Phi'][:, 0, :], colors='b', linewidths=0.5, levels=phi_levels)
    cs_psi = ax.contour(gr.x[1:-1], gr.z, psi, colors='r', linewidths=0.5, levels=psi_levels)

    h_lbl_pts = find_horizonal_intersections(cs_phi, z0=-35)
    v_lbl_pts = find_vertical_intersections(cs_phi, x0=0)
    lbl_pts = h_lbl_pts + v_lbl_pts
    ax.clabel(cs_phi, fmt='%.2f', manual=lbl_pts, inline=False, fontsize=8)

    # ax.clabel(cs_phi, levels=cs_phi.levels, fontsize=10, fmt='%.2f')
    
    # --- Add head just below the resistance layer
    phi7 = out['Phi'][np.where(gr.zm <= -7)[0][0], 0, :]
    ax.plot(gr.xm, phi7, '--', color='k', lw=1)

    
def add_geologic_background(xsec, ax=None, out=None):
    """Add the geologi background."""
    gr = xsec.gr
    
    # --- Get the color of each cell
    arr_RGB = xsec.leg_color_array(gr.idx_arr) / 255.

    # --- Use pcolormesh to color the cells. gr.idx_arr is here just a dummy for it's shape
    pcm = ax.pcolormesh(gr.X[:,0,:], gr.Z_corners[:,0,:], gr.idx_arr, shading='flat', alpha=1.0)

    # --- Wipe out the data array (gr.idx_arr) inside pcolormesh, as we'll set facecolor manually next.
    pcm.set_array(None)

    # --- Set the facecolor manually
    pcm.set_facecolor(arr_RGB.reshape(-1, 3))


def add_patches(xsec, ax=None, out=None):
    """Generate patches to show the canal, etc."""
    gr = xsec.gr
    
    # === Add ARK water body and sheet pilings as patches
    pARK   = extent_patch(gr.ARK_extent, fc='blue', ec='none', alpha=1)
    pDamwL = extent_patch(gr.damwL_extent, fc='k', ec='k', alpha=1)
    pDamwR = extent_patch(gr.damwR_extent, fc='k', ec='k', alpha=1)

    # --- Add patches to our current axes.
    ax.add_patch(pARK)
    ax.add_patch(pDamwL)
    ax.add_patch(pDamwR)
    
    # --- Set the zorder to render on top of contour lines
    pARK.set_zorder(3)
    pDamwL.set_zorder(3)
    pDamwR.set_zorder(3)

    # --- Add the thickness of the resistance layer patch at the bottom of the ARK
    cbp = c_bot_patch(xsec)
    cbp.set_zorder(4)
    ax.add_patch(cbp)


def add_legends(xsec, ax=None, out=None):
    # === Legend for the geologic background
    # --- Generate from xsec.geo_units
    # --- First the handles    
    handles = [
        Patch(facecolor=gu['color'], edgecolor='none',
            label=f"{k}: kh={gu['kh']:.1f}, kv={gu['kv']:.2f}, rhow={gu['rhow']:.0f}")
        for k, gu in xsec.geo_units.items()
    ]
    # --- Generate legend using these handgles and put the bbox to anchor
    leg_geo = ax.legend(handles=handles, loc='center left',
            bbox_to_anchor=(1, 0.5),
            title='Soil properties')
    ax.add_artist(leg_geo)

    # === Generate a legend for the contoured heads and streamlines
    leg_contours = ax.legend(
        handles =[
            Line2D([0], [0], color='k', lw=1.0, linestyle='--', label=r'$\phi$[z=-7m]'),
            Line2D([0], [0], color='b', lw=0.5, label='stijghoogte'),
            Line2D([0], [0], color='r', lw=0.5, label='stroomlijn')
        ],
        bbox_to_anchor=(1, 1),
        loc='upper left',
        title=r"Contours en $\phi$[z=-7]"
    )
    ax.add_artist(leg_contours)
    
    # === Legend for the ratio contours
    leg_ratio = ax.legend(handles=[
        Line2D([0], [0], ls='-.', color='red',    lw=1, label=r"$\sigma_{wat}/\sigma{tot}$=0.8"),
        Line2D([0], [0], ls='-.', color='orange', lw=1, label=r"$\sigma_{wat}/\sigma{tot}$=1.0"),
        Line2D([0], [0], ls='-.', color='green',  lw=1, label=r"$\sigma_{wat}/\sigma{tot}$=1.2"),
        ],
        bbox_to_anchor=(1, 0),
        loc='lower left',
        title='Water overdruk')
    ax.add_artist(leg_ratio)


def add_overpressure(xsec, ax=None, out=None):
    """Add the water over pressure to potentially lift overlaying layers."""
    
    # --- gravity and water densityWet density
    gr = xsec.gr    
    rho_wat = 1000.

    # --- Compute head at the bottom of the cells
    Phi_bot = out['Phi'].copy()  # --- at cell centers
    Qz = out['Qz'].copy()        # --- flow in upward z-direction (at bottom of cells)
    kv = gr.kv.copy()            # --- Vertical conductivity

    # --- Special cells
    Qz[kv[1:] == 0] = 0                  # --- Vertical flow at horizontal cell boundaries
    kv += 1e-6                           # --- Cells inside sheet piling have kv=0
    dPhi = Qz * gr.DZ[1:] / 2 / kv[1:]   # --- Correction for shift to bottom of the cells
    Phi_bot[:-1] += dPhi                 # --- Head ad cell bottoms

    # --- Water-pressure head at the bottom of cells
    h_wat = (Phi_bot - gr.Z[1:])

    # --- Total pressure caused by overlying wet density at bottom of cells gr.Z[1:].
    h_wet = np.cumsum(gr.rhow / rho_wat * gr.DZ, axis=0)

    # --- Pressure head difference (not used)
    delta_h = h_wet - h_wat
    delta_h = delta_h.clip(-2, 2)
    delta_h[xsec.gr.mask_ARK] = np.nan

    # --- Pressure head ratio (used)
    h_ratio = h_wat / h_wet
    h_ratio[xsec.gr.mask_ARK] = np.nan # --- Use nan inside the ARK water body
    out['h_ratio'] = h_ratio

    # --- Contour the 3 criteria levels hr
    cs = plt.contour(gr.XM[:,0,:], gr.Z[1:,0,:], h_ratio[:, 0, :],
                levels=[0.8, 1.0, 1.2],
                colors=['g','orange','r'],
                linestyles='-.',
                linewidths=1)


def phi_and_hrat(xsec, out=None, zlevels=[-3, -8, -13, -20], xlim=None):
    """Show the head and ratio sigma_water/sigma_total for fixed depths.
    
    Parmeters
    ---------
    xsec: x-section object
        x-section object with all x-section info.
    out: dictionary
        Result of fdm3 simulation with hrad added.
        add_overpressure must have been run before calling this function.
    zlevels: iterable
        Elevations for which the head and h_ratio should be plotted.
    
    """
    gr=xsec.gr
    Phi = out['Phi']
    zlevels = np.array(zlevels)
    
    try:
        hrat = out['h_ratio']
    except KeyError:
        raise KeyError("'h_ration missing: add_overpressure must be called first.")
    
    Iz = np.searchsorted(-gr.zm, -zlevels)
    if len(Iz) == 0:
        raise ValueError("zlevels must fall within the z-range of the model.")
    
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 7))
    
    cb = xsec.c_ARK_bot_dict
    xARK, yARK = xsec.xy_hart_ARK[0]
    ttl = (fr"xy hartlijn ARK in doorsnede : {xARK:.0f},{yARK:.0f}" + "\n" +
           fr"ARK_bodem: $D_0$={cb['D0']:.1f}m, erosie_fractie={cb['fd']:.2f}, c={cb['c']:.0f} d/m"
    )

    fig.suptitle(f'file: "{xsec.name}"' + "\n" + ttl, fontsize=8)
    
    for iz in Iz:
        zlev = np.round(gr.zm[iz], 1)
        ax1.plot(gr.xm, Phi[ iz, 0, :], label=f'z={zlev} m')
        ax2.plot(gr.xm, hrat[iz, 0, :], label=f'z={zlev} m')
        
    ax2.axhline(1.0, xmin=gr.xm[0], xmax=gr.xm[-1],
                c='r', ls='dashed', lw=2, label=r'$\sigma_{water}==\sigma_{total}$')
        
    ax1.vlines([-xsec.b, xsec.b], *ax1.get_ylim(),
               colors='k', linestyles='dashed', lw=1.5, label='ARK-damwand')
    ax2.vlines([-xsec.b, xsec.b], *ax2.get_ylim(),
               colors='k', linestyles='dashed', lw=1.5, label='ARK-damwand')
               
    ax1.set_title("Stijghoogte op vaste diepte")
    ax2.set_title(r"$\sigma_{water}/\sigma_{total}$ op vaste diepte")
    ax2.set_xlabel("x t.o.v. hartlijn ARK")
    ax1.set_ylabel(r'Stijghoogte $\phi$ [m NAP]')
    ax2.set_ylabel(r'ratio $\sigma_{water}/\sigma_{total}$')
    ax1.grid()
    ax2.grid()
    ax1.legend()
    ax2.legend()
    ax2.set_xlim(xlim)
    
    logo(fig, NOTEBOOK_NAME)
    return fig


# --- Steps comprising the workflow generating the final image.
steps = [
    ("background", add_geologic_background),
    ("contours", add_head_contours_and_stream_lines),
    ("patches", add_patches),
    ("legends", add_legends),
    ("overpressure", add_overpressure)
]


def snapshot(ax, label):
    """Print info when snapshot is called."""
    print(f"\n--- {label} ---")
    print(f"lies       : {len(ax.lines)}")
    print(f"patches    : {len(ax.patches)}")
    print(f"collections: {len(ax.collections)}")
    
# %% === Run all

# --- Add all required modeling and presentation info to xsec object
xsec = set_ARK_xsec(xsec)

# --- Parameters for canal bottom variants
# --- D0 = initial resistance layer thickness,
# --- fd = degree of erosion of the resistance layer,
# --- c1 = specific resistance [d/m] of the resistance layer,
# --- L  = width of cross section to plot (to get desired level of detail in plot)
cbot_vars = {
    0: dict(D0=0.01, fd=0.0, c=1,   L=250.),    
    1: dict(D0=1.00, fd=0.0, c=100, L=250.),    
    2: dict(D0=1.00, fd=0.5, c=100, L=250.),    
    3: dict(D0=1.00, fd=1.0, c=100, L=250.),
    4: dict(D0=1.00, fd=2.0, c=100, L=250.),    
    5: dict(D0=1.00, fd=2.0, c=100, L=1000.),
}

# --- Choose case
i = 1

# --- Choose width of cross section polot
L = cbot_vars[i]['L']

# --- Get parameters of canal bottom resistance
xsec.c_ARK_bot_dict = cbot_vars[i]

# --- Apply
set_c_ARK_bot(xsec)

# --- Simulate using fdm3 (finite difference model, steady state)
out = simulate_fdm3(xsec)

# --- Present results graphically
fig, ax = plot_setup(xsec, xlim=(-L,L))

# --- Steps to construct the image from the results
for i, (name, func) in enumerate(steps, 1):
    func(xsec, ax=ax, out=out)
    snapshot(ax, name)
    fig.savefig(os.path.join(dirs.images, 'debug', f"debug_{i:02d}_{name}.png"))
    
# --- Individual worksflow steps (now automated above)
# 1: add_geologic_background(xsec, ax=ax, out=None)
# 2: add_head_contours_and_stream_lines(xsec, ax=ax, out=out)
# 3: add_patches(xsec, ax=ax, out=None)
# 4: add_overpressure(xsec, ax=ax, out=out)
# 5: add_legends(xsec, ax=ax, out=None)

ax.set_xlim(-L, L)

# --- for file_name
xlim = f"_xlim{L:.00f}"

# --- Save the figure for reporting.
cd = xsec.c_ARK_bot_dict
fig.savefig(os.path.join(dirs.images,
        f"ARK_x{xsec.xRD:.0f}_y{xsec.yRD:.0f}_D{cd['D0']:.1f}_fd{cd['fd']:.2f}_c{cd['c']:.0f}{xlim}.pdf"))

fig = phi_and_hrat(xsec, out, zlevels=[-3, -8, -13, -20], xlim=(-L, L))
ax1, ax2 = fig.axes
ax1.set_xlim(-L, L)
ax2.set_xlim(-L, L)
fname = f"ARK_verloop_op_vaste_z_x{xsec.xRD:.0f}_y{xsec.yRD:.0f}_D{cd['D0']:.1f}_fd{cd['fd']:.2f}_c{cd['c']:.0f}{xlim}.pdf"
fig.savefig(os.path.join(dirs.images, fname))

print(f"Total infilration scenario {i} is out['Qfh']['Q'].sum()={out['Qfh']['Q'].sum():.3f} m2/d")

plt.show()

# %% === Heads along path, effect of viscosity

def rmu(T):
    """Return viscosity at TC / vicosity at 10C."""
    return 30 / (T + 20)

def temp_effect_on_heads(qinf=None, phi=None, temps=None, n=1):
    """Return the heads for temps in current scenario.
    
    Parameters
    ----------
    qinf: float
        average canal infiltration in m/d for temp of 10C.
    phi: np.array(np)
        heads at the selected points in the cross section
    temps: np.array(nt)
        temperatures for which the heads should be computed
    n: int
        Number of compartments for which temperature applies
        (compartments is the stretch between consecutive points)
        
    Returns
    -------
    phiT [nt, np]
    """    
    assert phi.ndim==1, "phi must be a 1D vector np long."
    temps = np.array(temps)
    assert temps.ndim==1, "temps must be a 1D vector, nt long"
    
    npnt = len(phi)
    ntmp = len(temps)

    r10 = -np.diff(phi) / qinf
    R10 = np.ones((ntmp, 1)) * r10[None, :]
    
    # --- Compartment resistance at other temperature
    RT = R10.copy()
    for i in range(n):
        RT[:, i] *= rmu(temps)
    
    # --- Infitlration at different temperature
    qinfT = qinf * R10.sum(axis=1) / RT.sum(axis=1)
    
    # --- Heads at different temperatures
    phiT = np.ones((ntmp, 1)) * phi[None, :]
    for j in range(len(r10)):
        phiT[:, j+1] = phiT[:, j] - qinfT * RT[:, j]
        
    assert np.all(phiT.shape == (ntmp, npnt)), f"phiT.shape must be ({ntmp},{npnt}) not {phiT.shape}"
    
    return phiT, qinfT

# %% --- Parameter values of the canal bottom in the 5 cases:
cbot_vars = {
    0: dict(D0=0.01, fd=0.0, c=1,   L=250.),    
    1: dict(D0=1.00, fd=0.0, c=100, L=250.),    
    2: dict(D0=1.00, fd=0.5, c=100, L=250.),    
    3: dict(D0=1.00, fd=1.0, c=100, L=250.),
    4: dict(D0=1.00, fd=2.0, c=100, L=250.),    
}

# --- World coordinates of points along the path in the cross section (x, z)
pts = np.array([(0, -2), (0, -7), (50, -13), (100, -7), (100, -2.25)])
# --- Distance along the flow path (approximately)
dist = np.hstack((0, np.cumsum(np.sqrt(np.diff(pts.T[0]) ** 2 + np.diff(pts.T)[1] ** 2))))

# --- Heads computed at these pts in the given scenarios/cases
heads = {'pts': pts,
         'hds': {
                    0: np.array([-0.4, -0.401, -1.084, -1.745, -2.295]),
                    1: np.array([-0.4, -1.887, -2.046, -2.183, -2.299]),
                    2: np.array([-0.4, -1.771, -1.971, -2.148, -2.298]),
                    3: np.array([-0.4, -1.007, -1.469, -1.916, -2.296]),
                    4: np.array([-0.4, -0.515, -1.156, -1.773, -2.295]),
         },
         'Qinf': {0:3.6113, 1:0.760 , 2:0.987, 3:2.505, 4:3.435}
         
}

str(np.array(pts))

# --- Number of scenarios
Nscen = len(heads['Qinf'])

# --- The flux between these points, considered a stream tube is the same and arbitrary
# --- Take scenario 1 (resistance = 100 d) as default
hds = [heads['hds'][k] for k in heads['hds']]
phiScen10 = np.array(hds)

# --- Infiltration flow in each of the scenarios for mean temp T=10C
Qinf10C = np.array([heads['Qinf'][k] for k in heads['Qinf']])

# --- Average infiltration rate over the width of the canal [m/d]
qin10C = Qinf10C / (2 * xsec.b)

# --- Resistance [scen, traject] for each path section m / [m/d] = d
Rscen10C = np.diff(-phiScen10, axis=1) / qin10C[:, None]

# --- Compare the heads at the points for different viscosities
phiScenT = phiScen10.copy()

linestyles = ['-', '--', '-.', ':', ':']

# --- Temperatures ---
temps = np.array([5, 10, 15, 20, 25])
cmap = plt.get_cmap('rainbow')
temp_clrs = cmap(np.linspace(0, 1, len(temps)))

# --- Setup the figure
fig, axs = plt.subplots(3, 2, sharex=True, sharey=True, figsize=(12, 12))

fig.suptitle("Effect van viscositeit op stijghoogte en infiltratie dsn. ARK\n"
             "Compartimenten met afwijkende temperatuur aangegeven met kleur\n"
             "Coordinaten: [" + ', '.join([f"({p[0]:.1f},{p[1]:.1f})" for p in pts]) + "]"             
             )

for ax, scen in zip(axs.flatten(), range(Nscen)): 

    phi  = heads['hds'][scen]
    qinf = heads['Qinf'][scen] / (2 * xsec.b) # mm/d

    n = 3
    
    phiT, qinfT = temp_effect_on_heads(qinf=qinf, phi=phi, temps=temps, n=n)
        
    clrs = cycle(temp_clrs)

    # --- Run over the temperatures for all scenarios at once
    for it, T in enumerate(temps):
        # --- Next temperature color        
        clr=next(clrs)
        
        # --- Only change visosity in first n compartments
            
        ax.plot(dist[n:], phiT[it][n:], '.-', color='k', ms=6, lw=0.5, label='')
        ax.plot(dist[:n+1], phiT[it][:n+1], '.-', color=clr, ms=6, lw=1.0,
                label=f"T={T}C, qinf={1000*qinfT[it]:.2f} mm/d")
      
    ax.grid(True)
    if scen > 3:
        ax.set_label('Afstand langs pad [m]')
    if scen % 2 == 0:
        ax.set_ylabel('stijghoogte [m+NAP]')
    try:
        phiT[it]
    except IndexError:
        break
        
    cb = cbot_vars[scen]
    D0, fd, c = cb['D0'], cb['fd'], cb['c']
    ax.set_title(f'D0={D0:.1f} m, c={c:.0f} d, fd={fd:.2f}', fontsize=8)
    ax.legend()
    
    logo(fig, NOTEBOOK_NAME)
    fig.savefig(os.path.join(dirs.images, f"temperatuur_effect_{n}_comp.pdf"))
    
    
# %% Exercise with wellen along path
phi0 = -0.4
phie = -2.4

c1, c2, c3, c4a, c4b = 100., 10., 10., 1000., 1000.
c4 = (c4a * c4b) /(c4a + c4b)

c = np.array([c1, c2, c2, c4])
phi = np.zeros(len(c) + 1)
x = np.linspace(0, len(phi), len(phi))

q = 1
phi[0] = phi0
for i in range(len(c)):
    phi[i+1] = phi[i] - c[i] * q
# plt.plot(x, phi, 'x-', label=f"q={q}")
q = (phi0 - phie) / (phi[0] - phi[-1])
for i in range(len(c)):
    phi[i+1] = phi[i] - c[i] * q

plt.plot(x, phi, '.-', label=f"q={q:.2f}")
plt.grid()
plt.legend()
plt.show()



# %% === Abrasing of the resistance layer at the canal bottom

# --- Set up plot
fig, ax = plt.subplots(figsize=(10, 6))
ax.set(title="Verloop dikte uitgeschuurde bodemweerstandslaag",
       xlabel="x vanaf hartlijn kanaal",
       ylabel='Dikte bodemweerstandslaag [m]')
ax.grid()

# --- Initial thickness
D0 = 1.0

# --- Abrased fractions
fds = [0.125, 0.25, 0.50, 1, 1.25, 1.50, 2.]

# ---- Computation and presentation
gr = xsec.gr
mask = gr.mask_ARK_bot
for fd in fds:
       xsec.c_ARK_bot_dict['fd'] = fd
       D = D_ARK_bot(xsec).clip(1e-3, None)
       ax.plot(gr.XM[mask], D,
            label=f'Uitschuring {np.round(100*fd)}%')
    
logo(fig, NOTEBOOK_NAME)
ax.legend()

fig.savefig(os.path.join(dirs.images, "uitschuring_ARK_bodem."))  
plt.show()


print(f"{' ':8}{'kh':8}{'kv':8}{'n':8}{'rho':8}{'rho_wet':8}")
print(f"{' ':8}{'m/d':8}{'m/d':8}{'-':8}{'kg/m3':8}{'kg/m3':8}")
print(f"{'--------'}{'--------'}{'--------'}{'--------'}{'--------'}{'--------'}")
for k, v in xsec.geo_units.items():
    print(f"{k:5}{v['kh']:8.3f}{v['kv']:8.3f}{v['n']:8.2f}{v['rho']:8.0f}{v['rhow']:8.0f}")