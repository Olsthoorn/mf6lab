import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.path import Path
import matplotlib.patches as patches


from src.mf6tools import  Dirs
from coords import plotdigitizerXML as dxml
from fdm.mfgrid import Grid
from etc import newfig

def plot_elev(data, x=None, Znew=None, fault_ix=None,
              title=None, xlabel='x [m]', ylabel='elevation [m]', figsize=(10, 8)):
    """Plot the elevation data"""
    ax = newfig(title, xlabel, ylabel, figsize=figsize)
    for elev in data:
        ax.plot(elev['xw'], elev['yw'], 'o-', mfc='none', lw=0.5)

    clrs = 5 * 'rbgmkc'  
    for iz, (znew, cl, cr) in enumerate(zip(Znew, IzL, IzR)):
        ax.plot(xm[:fault_ix], znew[:fault_ix], '--' , color=clrs[cl])
        ax.plot(xm[fault_ix:], znew[fault_ix:], '--',color=clrs[cr])
            
    plt.show()
    return ax

def gen_patch_x(x, zbot, ztop, alpha=None):
    """Return patch for the layer defined by x and ztop and zbot.
    
    Parameters
    ----------
    x: np.ndarray of x values (not xm values)
        x-coordinates of column boundaries
    zbot, ztop: np.ndarrays, corresponding to xm values
        bottom and top elevation of layer
    alpha: float, 0 <= alpha <= 1
        transparency
    """    
    xm = 0.5 * (x[:-1] + x[1:])
    xm[0], xm[-1] = x[0], x[-1]
    xP = np.hstack((xm, xm[::-1], xm[0]))       
    zP = np.hstack((zbot, ztop[::-1], zbot[0]))
    xy = np.vstack((xP, zP)).T
    codes = np.zeros(len(xy), dtype=int) + Path.LINETO
    codes[0] = Path.MOVETO
    codes[-1] = Path.CLOSEPOLY
    pth = Path(xy, codes)
    return patches.PathPatch(pth, alpha=alpha, ec='black', lw=0.25)

def plot_cross_section(gr=None, title=None, xlabel='x [m]', ylabel='elevation [m]',
                       props=None, IzL=None, IzR=None, fault_ix=None, figsize=(10, 8)):
    """Plot the cross section as generated with the fault.
    
    Parameters
    ----------
    gr: Grid object
        The grid (n, ny, nx)
    titles: (str, str, str)
        The title of the plot and the xlabel and ylabel
    IzL: array of float indices (gr.Nz + 1)
        array showing elvation of layer relatieve to the
        originally digitized layer so that it matches elevations of right layers
    IzR: array of float indices (gr.Nz + 1)
        laike IzL but of right side of fault
    fault_ix: int
        index of gr.x where the fault is.
    kwargs: dict
        extra prameters, now only used for figsize
    """

    x, Znew = gr.x, gr.Z[:, 0, :]
    
    ax = newfig(title, xlabel, ylabel, figsize=figsize)
    
    # ax.plot([x[0], x[-1]], [8.5, 8.5], 'darkblue', label='test watertafel')
    
    # Layer id from Excel file
    clrs, codes, names,  = props['Color'], props['Code'], props['Name']
    
    #IxL = np.arange(len(x ), dtype=int)[:fault_ix] # use fault_ix instead
    #IxR = np.arange(len(xm), dtype=int)[fault_ix:] # use fault_ix instead

    codes_done = []

    for iz, iL, iR in zip(np.arange(gr.nz), IzL, IzR):
        pL = gen_patch_x(x[:fault_ix + 1], Znew[iz + 1, :fault_ix], Znew[iz, :fault_ix], alpha=None)
        
        pR = gen_patch_x(x[fault_ix:], Znew[iz + 1, fault_ix:], Znew[iz, fault_ix:], alpha=None)
        
        pL.set_fc(mcolors.CSS4_COLORS[clrs[iL]])
        pR.set_fc(mcolors.CSS4_COLORS[clrs[iR]])
        if iL not in codes_done:
            pL.set_label(codes[iL] + ' ' + names[iL])
            codes_done.append(iL)
        if iR not in codes_done:
            pR.set_label(codes[iR] + ' ' + names[iR])
            codes_done.append(iR)
        
        ax.add_patch(pL)
        ax.add_patch(pR)

    ax.legend(loc='lower right')
    return ax


HOME = '/Users/Theo/GRWMODELS/python/mf6lab/Projects/DrasseDriehoek/'
assert os.path.isdir(HOME), "Can't find the directory {}".format(HOME)

LENGTH_UNITS = 'meters'
TIME_UNITS = 'days'

dirs = Dirs(HOME)

section_name = 'Drasse Driehoek Tilburg dsn NW-ZO'
sim_name = 'DSN_NW_ZO'
dirs = dirs.add_case(sim_name)
xmlfile = os.path.join(dirs.data, sim_name + '.xml')

os.chdir(dirs.case)

params_wbk = os.path.join(dirs.case, sim_name + '.xlsx')
assert os.path.isfile(params_wbk), "Params_wbk not found: {}".format(params_wbk)

lay = pd.read_excel(params_wbk, sheet_name='LAY', header=0, index_col=0)

props = { 'start_date_time': '2025-01-17',
         'L': 12000.0, # [m] width of model
         'minDz': 0.2, # m min layer thickness
          'drain_depth': 1e-5, # m
          'cDrainage': 100., # d
          'rch': 0.001, # m/d
          'strthd': 0.0, # m initial head
          'dz_pinched_out': 0.01, # [m] thickness of pinched out layers"
}

# === Get the digitized coordinates ======= the free plotdigitizder app was used.
assert os.path.isfile(xmlfile), "Can't find file {}".format(xmlfile)

fault_x = 4190.
kwargs = {'fault_x': fault_x, 'fault_width': 200., 'xjump_min': 200.}

data, meta = dxml.fromPlotdigitizerXML(xmlfile, chunkit=True, **kwargs)

x, Z = dxml.get_layertops(data, meta, dx=10.)
xm = 0.5 * (x[:-1] + x[1:])

Znew, IzL, IzR, fault_ix = dxml.bridge_fault(Z, x, fault_x=fault_x, verbose=True)

gr = Grid(x, [-0.5, 0.5], Znew[:, np.newaxis, :], axial=False, min_dz=props['minDz'])

Lk, Lk33, Lsy, Lss = lay['k'], lay['k33'], lay['Sy'], lay['Ss']
k, k33, sy, ss = gr.const(0.), gr.const(0.), gr.const(0.), gr.const(0.)

icelltype = gr.const(1, dtype=int)

for iz, iL, iR in zip(np.arange(len(Znew) - 1), IzL, IzR):
    k[  iz, :, :fault_ix] = Lk[  iL]
    k[  iz, :, fault_ix:] = Lk[  iR]
    k33[iz, :, :fault_ix] = Lk33[iL]
    k33[iz, :, fault_ix:] = Lk33[iR]
    sy[ iz, :, :fault_ix] = Lsy[ iL]
    sy[ iz, :, fault_ix:] = Lsy[ iR]
    ss[ iz, :, :fault_ix] = Lss[ iL]
    ss[ iz, :, fault_ix:] = Lss[ iR]

fig_kwargs = {'title': section_name, 'xlabel': 'x[m]', 'ylabel': 'NAP [m]',
              'figsize': (16, 8)}

if __name__ == '__main__':
    

    ax = plot_elev(data, x=x, Znew=Znew, fault_ix=fault_ix, **fig_kwargs)

    ax = plot_cross_section(gr=gr, props=lay, IzL=IzL, IzR=IzR, fault_ix=fault_ix, **fig_kwargs)
    plt.show()


