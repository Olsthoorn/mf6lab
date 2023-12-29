# Section data from MoervaarDepressieDekzandrugMaldegemStekene.xml
#
# The data have been digitized from the image
#
# MoervaarDepressieDekzandrugMaldegemStekene.png
#
# with the free app plotdigitizer and then exported to xml (only possible with the free app).
#
# TO 20231223

# %%

import settings
import os
import numpy as np
from coords import fromPlotdigitizerXML
import matplotlib.pyplot as plt
from etc import newfig, color_cycler
from fdm.mfgrid import Grid

dirs = settings.dirs

#%% Row data, from plotdigitzer web app

xmlfile = os.path.join(dirs.data, 'MoervaarDpr.xml')
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
        
x = np.linspace(0, 8400, 8401)
xm = 0.5 * (x[:-1] + x[1:])
Z = np.zeros((len(elev), len(xm)))
ztol = 0.01 # m
for iz, e in enumerate(elev):
        Z[iz] = np.interp(xm, e['xw'], e['yw'])
        if iz > 0:
                Z[iz][Z[iz] >= Z[iz-1]] = Z[iz-1][Z[iz] >= Z[iz - 1]] - ztol
            
gr = Grid(x, [-0.5, 0.5], Z[:, np.newaxis, :], axial=False, min_dz=1e-6)

# layer can be exported

if __name__ == '__main__':

    # Show the results
    ax = newfig("The cross section lines", "x [m]", "z [m]")

    layer_patches = gr.layer_patches_x(row=0) # Get layer patches
    
    colors = ['yellow', 'orange', 'gold', 'lightskyblue',
              'lightsalmon', 'violet', 'chocolate', 'yellowgreen']
    
    for p, clr in zip(layer_patches, color_cycler(colors)):
        p.set_fc(clr)
        p.set_ec('k')
        p.set_alpha(1.0)
        p.set_lw(0.25)
        ax.add_patch(p)

    #for i, (e, clr) in enumerate(zip(elev, color_cycler(colors))):
    #    ax.plot(e['xw'], e['yw'], 'k', label=f'digitized layer {i}')

    ax.plot([0, 8000], [0, 0], 'c', label='test line')
    #ax.legend()
    plt.show()
    
    # %%

