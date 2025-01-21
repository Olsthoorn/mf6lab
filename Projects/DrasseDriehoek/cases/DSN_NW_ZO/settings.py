import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from src.mf6tools import  Dirs
from coords import fromPlotdigitizerXML
from fdm.mfgrid import Grid
from etc import newfig

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
         'minDz': 0.25, # m min layer thickness
          'drain_depth': 1e-5, # m
          'cDrainage': 100., # d
          'rch': 0.001, # m/d
          'strthd': 0.0, # m initial head
          'dz_pinched_out': 0.01, # [m] thickness of pinched out layers"
}


# === Get the digitized coordinates ======= the free plotdigitizder app was used.
assert os.path.isfile(xmlfile), "Can't find file {}".format(xmlfile)

data, meta = fromPlotdigitizerXML(xmlfile, chunckit=True)

fig, ax = plt.subplots()

        
deltax = 1 # m
x = np.linspace(meta['aX1'], meta['aX2'], int((meta['aX2'] - meta['aX1']) / deltax) + 1)
xm = 0.5 * (x[:-1] + x[1:])
Z = np.zeros((len(data), len(xm)))
for iz, e in enumerate(data):
        Z[iz] = np.interp(xm, e['xw'], e['yw'])
        ax.plot(xm, Z[iz])
        ax.plot(e['xw'], e['yw'], '.')
        if iz > 0:
            check = Z[iz] >= Z[iz-1] - props['dz_pinched_out']
            Z[iz][check] = Z[iz-1][check] - props['dz_pinched_out']
            
gr = Grid(x, [-0.5, 0.5], Z[:, np.newaxis, :], axial=False, min_dz=1e-6)
plt.show()

if __name__ == '__main__':
    
    # Show the results
    title = section_name
    ax = newfig(title, "Lijnafstand (m)", "mTAW")

    layer_patches = gr.layer_patches_x(row=0) # Get layer patches
    
    ax.plot([0, props['L']], [0, 0], 'darkblue', label='test watertafel')
    
    for p, clr, code, name in zip(layer_patches, lay['Color'], lay['Code'], lay['Name']):
        p.set_fc(mcolors.CSS4_COLORS[clr])
        p.set_ec('k')
        p.set_alpha(1.0)
        p.set_lw(0.25)
        p.set_label(code + ' ' + name)
        ax.add_patch(p)

    ax.legend()
    plt.show()

