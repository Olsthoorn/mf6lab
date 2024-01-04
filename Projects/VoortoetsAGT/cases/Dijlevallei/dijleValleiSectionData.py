# Section data from for (see settings.sim_name) in .xml file exported from plotDigitizer app.
##
# TO 20231230

# %%

import settings
import os
import numpy as np
from coords import fromPlotdigitizerXML
import matplotlib.pyplot as plt
from etc import newfig
from fdm.mfgrid import Grid

dirs = settings.dirs
sim_name = settings.sim_name
lay = settings.lay


#%% Row data, from plotdigitzer web app

xmlfile = os.path.join(dirs.data, settings.sim_name + '.xml')
assert os.path.isfile(xmlfile), "Can't find file {}".format(xmlfile)

data, meta = fromPlotdigitizerXML(xmlfile)

# split the data upon jump back of x
dx = np.diff(data['x'])
I = np.hstack((0, np.array(np.where(dx < -250)[0]) + 1, len(data)))
elev = []
for i1, i2 in zip(I[:-1], I[1:]):
        print(i1, i2)
        print(data[i1:i2])
        elev.append(data[i1:i2])

L = 11800 # Width of cross section
x = np.linspace(0, L, int(L / 10 + 1))
xm = 0.5 * (x[:-1] + x[1:])
Z = np.zeros((len(elev), len(xm)))
ztol = 0.01 # m
for iz, e in enumerate(elev):
        Z[iz] = np.interp(xm, e['xw'], e['yw'])
        if iz  > 0:
            check = Z[iz] >= Z[iz-1] - ztol
            Z[iz][check] = Z[iz-1][check] - ztol
            
gr = Grid(x, [-0.5, 0.5], Z[:, np.newaxis, :], axial=False, min_dz=1e-6)

# layer can be exported

if __name__ == '__main__':

    # Show the results
    title = settings.section_name
    ax = newfig(title, "Lijnafstand [m]", "mTWA")

    layer_patches = gr.layer_patches_x(row=0) # Get layer patches
        
    ax.plot([0, L], [0, 0], 'darkblue', label='test watertafel')

    for p, clr, code, name in zip(layer_patches, lay['Color'], lay['Code'], lay['Name']):
        p.set_fc(clr)
        p.set_ec('k')
        p.set_alpha(1.0)
        p.set_lw(0.25)
        p.set_label(code + ' ' + name)
        ax.add_patch(p)

    #for e, clr in zip(elev, lay['Color']):
    #    ax.plot(e['xw'], e['yw'], clr, marker='.') #, label=f'digitized layer {i}')

    ax.legend()
    plt.show()
    
    # %%

