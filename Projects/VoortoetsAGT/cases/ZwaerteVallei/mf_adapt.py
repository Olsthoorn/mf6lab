# %% Geeneric cross section model
"""
A section model for Dijle Vallei bij Leuven is made, based on the digitized elevation of the relevant layers.

The image and the xml file are in the immages and data repectively.

The digitizing has been done from an image of the cross section using plotdigitizer web app.

The idea is to show the stream lines and the water table in the cross section,
which is fed by recharge and drained where the water table reaches ground surface. 
"""

# TO 20231223

import os
import numpy as np
from src import mf6tools
from zwaerteBeekvalleiSectionData import gr
import settings

dirs = settings.dirs
sim_name = settings.sim_name
section_name = settings.section_name

# Parameters workbook
params_wbk = settings.params_wbk

## Get section data

# %% === tdis ==========  Period Data:
start_date_time = '2024-01-01' # Must be a string.

perDF = mf6tools.get_periodata_from_excel(params_wbk, sheet_name='PER')
period_data = [tuple(sp) for sp in perDF[['PERLEN', 'NSTP', 'TSMULT']].values]

# Stress period start times
sp_start_times = np.datetime64(start_date_time) + np.cumsum(
                  np.hstack((0., np.array([p[0] for p in period_data])[:-1]))
                  ) * np.timedelta64(1, 'D')

Simtdis = {'perioddata': period_data,
           'nper': len(period_data),
           'start_date_time': start_date_time,
           'time_units': settings.TIME_UNITS,
           }

# Conductivities
lay = settings.lay

# %% === DIS ========== Grid

IDOMAIN = gr.const(1, dtype=int)
IDOMAIN[gr.DZ < 0.25] = -1

Gwfdis = {'gr': gr,
          'idomain': IDOMAIN,
          'length_units': settings.LENGTH_UNITS}


# %% ==== STO ===============
Gwfsto = {'sy': gr.const(lay['Ss'].values),
          'ss': gr.const(lay['Ss'].values)
          }

# %% === Gwfnpf =========== Horizontal and vertical conductivity

Gwfnpf = {  'k':   gr.const(lay['k'].values),
            'k33': gr.const(lay['k33'].values),
            'icelltype': gr.const(lay['ICELLTYPE'].values),
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

hDr = gr.Z[0, 0] - 0.5
Cdr = gr.Area[0] / c
Iz = np.zeros(gr.nx, dtype=int)

J = np.arange(gr.nx, dtype=int)
for i in range(gr.nz):
      z = gr.Z[i + 1, 0]
      J = J[np.logical_or(hDr[J] < z[J], IDOMAIN[i, 0, J] < 0)]
      # print(i, len(J))
      if len(J) > 0:
            Iz[J] += 1
      
DRN = [((iz, 0, i), h_, C_) for i, (iz, h_, C_) in enumerate(zip(Iz, hDr, Cdr))]

Gwfdrn = {'stress_period_data': {0: DRN}}

# %% === Rch

rch = 0.001 # m/d

RCH = [((iz, 0, i), rch) for i, iz in enumerate(Iz)]

Gwfrch = {'stress_period_data': {0: RCH}}

# %% === OC ====

Gwfoc = {'head_filerecord':   os.path.join(dirs.SIM, "{}Gwf.hds".format(sim_name)),
         'budget_filerecord': os.path.join(dirs.SIM, "{}Gwf.cbc".format(sim_name)),
         'saverecord': [("HEAD", "ALL"), ("BUDGET", "ALL")],
}

print('Done mf_adapt')

if __name__ == '__main__':
   
   print(dirs)     