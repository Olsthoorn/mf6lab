import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from src.mf6tools import  Dirs
from coords import fromPlotdigitizerXML
from fdm.mfgrid import Grid
from etc import newfig

HOME = '/Users/Theo/GRWMODELS/python/mf6lab/Projects/Wells/'
assert os.path.isdir(HOME), "Can't find the directory {}".format(HOME)

dirs = Dirs(HOME)

section_name = 'Xection Horstermeer'
sim_name = 'Horstermeer'

dirs = dirs.add_case(sim_name)
xmlfile = os.path.join(dirs.data, sim_name + '.xml')

os.chdir(dirs.case)

params_wbk = os.path.join(dirs.case, sim_name + '.xlsx')
assert os.path.isfile(params_wbk), "Params_wbk not found: {}".format(params_wbk)

lay = pd.read_excel(params_wbk, sheet_name='LAY', header=0, index_col=0)

props = {
    'start_date_time': '2024-01-19',  # Start time of Pennink test series 3
    'length_units': 'meters',
    'time_units': 'days',
    'oc_frequency': 10, # Saving frequency for budget heads and concentrations
    'RWellRing': 200, #m  Radius of ring of well in the Horstermeer polder
    'R1': 1500., # Radius Horstermeer
    'R2': 10000., # Radius boundary model
    'mvDefault': 0.0, # NAP
    'mvPolder': -3.0, # NAP
    'min_dz': 1e-6, # minimum thickness of layers
    'por': 0.35,       # [-] porosity
    'minDz': 0.25, # m min layer thickness
    'dDitch': 1.0, # m Depth of ditch
    'wellTop': -100.0, # m
    'wellBot': -150.0, # m
    'Qwell': -2.1, #  -21000.0, # m3/d, total extraction 7.5e6 / 360
    'rch': 0.001, # m/d
    'strthd': 0.0, # m initial head
    'dz_pinched_out': 0.01, # [m] thickness of pinched out layers"
    'hFar': 0.0, # [m] initial head
    'zIface':  np.nan,      # [cm] elavation of initial interface, will be replaced see below
    'rhoFresh': 1.,    # [g/cm3]
    'rhoSalt': 1.03 / 6, # [g/cm3] milk acc. to Pennink
    'cFresh': 0.0,     # [g/cm3]
    'cSalt':  1.0,     # [g/cm3] 
    'ahl': 1.0,  # [m] longitudal dispersivity
    'ath1': 0.1,  # [m] transversal dispersivity
    'ath2': 0.1,  # [m] transversal dispersivity
    'atv': 0.1,  # [m] vertical dispersivity
    'diffc': 6e-5, # [m2/d] diffusion coefficient
    'drainDepth': 0.5, # [m]
    'cDrainage': 100, # [d]
}

# === Get the digitized coordinates ======= the free plotdigitizder app was used.
x = np.hstack((0., np.logspace(1, np.log10(props['R2']), 50), props['R2'] + 1))
xm = 0.5 * (x[:-1] + x[1:])
elev = props['mvDefault'] + np.zeros_like(xm)
elev[xm < props['R1']] = props['mvPolder']

xStd = 3000.0 # std of iface curve
props['zIface'] = -150 + 100 * np.exp(-((xm / xStd) ** 2))

Z = np.zeros((len(lay) + 1, len(xm)))
for ilay in lay.index:
    Z[ilay + 1] = lay['Tot'][ilay]
Z[0] = elev
            
gr = Grid(x, [-0.5, 0.5], Z[:, np.newaxis, :], axial=True, min_dz=props['min_dz'])

if __name__ == '__main__':
    
    # Show the results
    title = section_name
    ax = newfig(title, "Lijnafstand (m)", "mTAW")

    layer_patches = gr.layer_patches_x(row=0) # Get layer patches
    
    ax.plot([0, props['R2']], [-1.5, -1.5], 'darkblue', label='test watertafel')
    
    ax.plot(xm, props['zIface'], 'cyan', lw=2.5, label='interface')
    
    
    for p, clr, code, name in zip(layer_patches, lay['Color'], lay['Code'], lay['Name']):
        p.set_fc(clr)
        p.set_ec('k')
        p.set_alpha(1.0)
        p.set_lw(0.25)
        p.set_label("{:10s} {}".format(code, name))
        ax.add_patch(p)

    ax.legend()
    plt.show()
