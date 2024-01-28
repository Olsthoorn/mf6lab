import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Path

from src.mf6tools import  Dirs
from fdm.mfgrid import Grid
from etc import newfig

HOME = '/Users/Theo/GRWMODELS/python/mf6lab/Projects/Wells/'
assert os.path.isdir(HOME), "Can't find the directory {}".format(HOME)

LENGTH_UNITS = 'meters'
TIME_UNITS = 'days'

dirs = Dirs(HOME)

section_name = 'Theis Hantush Delayed yield'
sim_name = 'TheisHantush'
dirs = dirs.add_case(sim_name)
xmlfile = os.path.join(dirs.data, sim_name + '.xml')

os.chdir(dirs.case)

params_wbk = os.path.join(dirs.case, sim_name + '.xlsx')
assert os.path.isfile(params_wbk), "Params_wbk not found: {}".format(params_wbk)

lay = pd.read_excel(params_wbk, sheet_name='LAY', header=0, index_col=0)

props = { 'start_date_time': '2024-01-24',
         'oc_frequency': 1,
         'L': 10000.0, # [m] width of model
         'rw': 0.1,
         'minDz': 0.25, # m min layer thickness
          'drain_depth': 1e-5, # m
          'cDrainage': 100., # d
          'rch': 0.001, # m/d
          'strthd': 0.0, # m initial head
          'dz_pinched_out': 0.01, # [m] thickness of pinched out layers"
          'Q' : -12000.0, # m3/d'
          'ztop': 1.0,
}

props['R'] = props['L']
r = np.hstack((0.0, np.logspace(np.log10(props['rw']), np.log10(props['R']), 81), props['R'] + 0.01))
z = props['ztop'] - np.cumsum([0, *lay['D']], axis=0)
         
gr = Grid(r, [], z, axial=True, min_dz = props['dz_pinched_out'])


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
        
    # well screen patches:
    for ilay, screened in lay['Screened'].items():
        if screened:
            zt, zb = gr.z[ilay], gr.z[ilay + 1]            
            p = Rectangle([0.0, zb], props['rw'], zt - zb)
            p.set_fc('black')
            p.set_ec('black')
            p.set_alpha(0.8)
            p.set_lw(0.25)
            p.set_label('screen layer {}'.format(ilay))
            ax.add_patch(p)
            
    ax.set_xlim(0.0, 5)

    ax.legend()
    
    print("Done")
    
    plt.show()

