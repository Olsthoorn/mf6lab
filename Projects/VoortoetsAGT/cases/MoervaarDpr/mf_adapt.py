# %% Geeneric cross section model
"""
A section model for Moervaar Depressie (Dekzandrug Maldegem-Stekene) is made, based on the digitized elevation of the relevant layers.

The image and the xml file are in the immages and data repectively.

The digitizing has been done from an image of the cross section using plotdigitizer web app.

The idea is to show the stream lines and the water table in the cross section,
which is fed by recharge and drained where the water table reaches ground surface. 
"""

# TO 20231223

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from mf6lab import mf6tools
from moervaarSectionData import gr
import pandas as pd
import mf_adapt

HOME = '/Users/Theo/GRWMODELS/python/VoortoetsAGT/'

LENGTH_UNITS = 'meters'
TIME_UNITS = 'days'

dirs = mf6tools.Dirs(HOME)

sim_name = 'MoervaarDpr'
dirs = dirs.add_case(sim_name)

# Parameters workbook
params_wbk = os.path.join(dirs.case, sim_name + '.xlsx')
assert os.path.isfile(params_wbk), "Params_wbk not found: {}".format(params_wbk)

# chdir to case directory
os.chdir(dirs.case)
sys.path.insert(0, dirs.case)

## Get section data

# %% === tdis ==========  Period Data:
start_date_time = '2023-12-22' # Must be a string.

perDF = mf6tools.get_periodata_from_excel(params_wbk, sheet_name='PER')
period_data = [tuple(sp) for sp in perDF[['PERLEN', 'NSTP', 'TSMULT']].values]

# Stress period start times
sp_start_times = np.datetime64(start_date_time) + np.cumsum(
                  np.hstack((0., np.array([p[0] for p in period_data])[:-1]))
                  ) * np.timedelta64(1, 'D')

Simtdis = {'perioddata': period_data,
           'nper': len(period_data),
           'start_date_time': start_date_time,
           'time_units': TIME_UNITS,
           }

# Conductivities
lay = pd.read_excel(params_wbk, sheet_name='LAY', header=0, index_col=0)

# %% === DIS ========== Grid

IDOMAIN = gr.const(1, dtype=int)
IDOMAIN[gr.DZ < 0.25] = -1

Gwfdis = {'gr': gr,
          'idomain': IDOMAIN,
          'length_units': LENGTH_UNITS}


# %% ==== STO ===============
Gwfsto = {'sy': gr.const(lay['Ss'].values),
          'ss': gr.const(lay['Ss'].values)
          }

# %% === Gwfnpf =========== Horizontal and vertical conductivity

Gwfnpf = {  'k':   gr.const(lay['k'].values),
            'k33': gr.const(lay['k33'].values)
            }

# %% === Gwfic ======== Initial condictions (head)

strthd = gr.const(0.)

Gwfic = {'strt': strthd}

# %% === CHD, Fixed head period data (Only specify the first period)
# Location of wells
# Rwvp1, Rwvp2
# Extraction by wells


# %% === WEL, required to lower head to measured  value
# Location of wells
# Rwvp1, Rwvp2
# Extraction by wells

# %% === DRN,
# Drain Cr [m2/d] is established from dh [m] head loss for c resistance [d].
# q = dh /c [m/d] --> Q = dh dx dy  /c [m2/s] --> Cdr = Q /dh = dx dy / c
c = 100 # d # drainage resistance [d]

hDr = gr.Z[0] - 0.5
Cdr = gr.Area / c
Iz = np.zeros(gr.shape[1:], dtype=int)

for i in range(1, gr.nz):
      Iz[hDr < gr.Z[i]] += 1
      
DRN = [((iz, 0, i), h_, C_) for i, (iz, h_, C_) in enumerate(zip(Iz[0], hDr[0], Cdr[0]))]

Gwfdrn = {'stress_period_data': {0: DRN}}

# %% === Rch

rch = 0.001 # m/d

RCH = [((iz, 0, i), rch) for i, iz in enumerate(Iz[0])]

Gwfrch = {'stress_period_data': {0: RCH}}

# %% === OC ====

Gwfoc = {'head_filerecord':   os.path.join(dirs.SIM, "{}Gwf.hds".format(sim_name)),
         'budget_filerecord': os.path.join(dirs.SIM, "{}Gwf.cbc".format(sim_name)),
         'saverecord': [("HEAD", "ALL"), ("BUDGET", "ALL")],
}

print('Done mf_adapt')

if __name__ == '__main__':
   
   print(dirs)     