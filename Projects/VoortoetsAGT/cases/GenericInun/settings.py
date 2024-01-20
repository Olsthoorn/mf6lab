import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from src.mf6tools import  Dirs
from coords import fromPlotdigitizerXML
from fdm.mfgrid import Grid
from etc import newfig

HOME = '/Users/Theo/GRWMODELS/python/mf6lab/Projects/VoortoetsAGT/'
assert os.path.isdir(HOME), "Can't find the directory {}".format(HOME)

LENGTH_UNITS = 'meters'
TIME_UNITS = 'days'

dirs = Dirs(HOME)

section_name = 'Generic hand-drawn cross section met inundatielaag'
sim_name = 'GenericInun'

dirs = dirs.add_case(sim_name)
xmlfile = os.path.join(dirs.data, sim_name + '.xml')

os.chdir(dirs.case)

params_wbk = os.path.join(dirs.case, sim_name + '.xlsx')
assert os.path.isfile(params_wbk), "Params_wbk not found: {}".format(params_wbk)

lay = pd.read_excel(params_wbk, sheet_name='LAY', header=0, index_col=0)

props = {
    'start_date_time': '2024-01-19',  # Start time of Pennink test series 3
    'oc_frequency': 10, # Saving frequency for budget heads and concentrations
    'L': 8900., # [m] length of cross section
    'por': 0.35,       # [-] porosity
    'minDz': 0.25, # m min layer thickness
    'rch': 0.001, # m/d
    'strthd': 0.0, # m initial head
    'dz_pinched_out': 0.01, # [m] thickness of pinched out layers"
    'hStrt': 0.0, # [m] initial head
    'zIface': -30,      # [cm] elavation of initial interface
    'rhoFresh': 1.,    # [g/cm3]
    'rhoSalt': 1.03, # [g/cm3] milk acc. to Pennink
    'cFresh': 0.0,     # [g/cm3]
    'cSalt':  1.0,     # [g/cm3] 
    'ahl': 1.0,  # [m] longitudal dispersivity
    'ath1': 0.1,  # [m] transversal dispersivity
    'ath2': 0.1,  # [m] transversal dispersivity
    'atv': 0.1,  # [m] vertical dispersivity
    'diffc': 6e-5, # [m2/d] diffusion coefficient
    'drainDepth': 0.5, # [m]
    'cDrainage': 100, # [d]
    'rch': 0.001, # [m/d]
}

# === Get the digitized coordinates ======= the free plotdigitizder app was used.
assert os.path.isfile(xmlfile), "Can't find file {}".format(xmlfile)

data, meta = fromPlotdigitizerXML(xmlfile)

# === split the digitized coordinates upon jump back of x =====
dx = np.diff(data['x'])
L = np.hstack((0, np.array(np.where(dx < -250)[0]) + 1, len(data)))
elev = []
for i1, i2 in zip(L[:-1], L[1:]):
        print(i1, i2)
        print(data[i1:i2])
        elev.append(data[i1:i2])
        
x = np.linspace(0, props['L'], int(props['L'] / 10 + 1))
xm = 0.5 * (x[:-1] + x[1:])
Z = np.zeros((len(elev), len(xm)))
for iz, e in enumerate(elev):
        Z[iz] = np.interp(xm, e['xw'], e['yw'])
        if iz > 0:
            check = Z[iz] >= Z[iz-1] - props['dz_pinched_out']
            Z[iz][check] = Z[iz-1][check] - props['dz_pinched_out']
            
gr = Grid(x, [-0.5, 0.5], Z[:, np.newaxis, :], axial=False, min_dz=1e-6)


if __name__ == '__main__':
    
    # Show the results
    title = section_name
    ax = newfig(title, "Lijnafstand (m)", "mTAW")

    layer_patches = gr.layer_patches_x(row=0) # Get layer patches
    
    ax.plot([0, props['L']], [0, 0], 'darkblue', label='test watertafel')
    
    for p, clr, code, name in zip(layer_patches, lay['Color'], lay['Code'], lay['Name']):
        p.set_fc(clr)
        p.set_ec('k')
        p.set_alpha(1.0)
        p.set_lw(0.25)
        p.set_label(code + ' ' + name)
        ax.add_patch(p)

    ax.legend()
    plt.show()
