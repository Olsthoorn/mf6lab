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
x = np.linspace(0, L, int(L + 1))
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
    ax = newfig("Cross section {}".format(settings.sim_name), "Lijnafstand [m]", "mTWA [m]")

    layer_patches = gr.layer_patches_x(row=0) # Get layer patches
    
    colors = settings.lay['Color']
    codes  = settings.lay['Code']
    names  = settings.lay['Name']
    
    ax.plot(gr.x[[0, -1]], [36, 36], 'darkblue', label='test water table')

        
    for p, clr, code, name in zip(layer_patches, colors, codes, names):
        p.set_fc(clr)
        p.set_ec('k')
        p.set_alpha(1.0)
        p.set_lw(0.25)
        p.set_label(code + ' ' + name)
        ax.add_patch(p)

    #for i, (e, clr) in enumerate(zip(elev, color_cycler(colors))):
    #    ax.plot(e['xw'], e['yw'], 'k', label=f'digitized layer {i}')

    ax.legend()
    plt.show()
    
    # %%

