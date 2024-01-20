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
import numpy as np
from src import mf6tools
from genericSectionData import gr
import settings

INACTIVE =[-1, 0]

dirs = settings.dirs
sim_name = settings.sim_name
section_name = settings.section_name
pr = settings.props

# Parameters workbook
params_wbk = settings.params_wbk

# %% === tdis ==========  Period Data:

start_date_time = pr['start_date_time'] # Must be a string.

perDF = mf6tools.get_periodata_from_excel(params_wbk, sheet_name='PER')
period_data = [tuple(sp) for sp in perDF[['PERLEN', 'NSTP', 'TSMULT']].values]

# Stress period start times
sp_start_times = np.datetime64(start_date_time) + np.cumsum(
                  np.hstack((0., np.array([p[0] for p in period_data])[:-1]))
                  ) * np.timedelta64(1, 'D')

Simtdis = {'perioddata': period_data,
           'nper': len(period_data),
           'start_date_time': pr['start_date_time'],
           'time_units': settings.TIME_UNITS,
           }

# Conductivities
lay = settings.lay

# %% === Gwfdis ========== Grid definition

IDOMAIN = gr.const(1, dtype=int)
IDOMAIN[gr.DZ < pr['minDz']] = -1 # Make vertical flow-through cell.

Gwfdis = {'gr': gr,
          'idomain': IDOMAIN,
          'length_units': settings.LENGTH_UNITS}

# %% ==== Gwfsto ====== Storage (transient)
Gwfsto = {'sy': gr.const(lay['Ss'].values),
          'ss': gr.const(lay['Ss'].values)
          }

# %% === Gwfnpf =========== Horizontal and vertical conductivity

Gwfnpf = {  'k':   gr.const(lay['k'].values),
            'k33': gr.const(lay['k33'].values),
            'icelltype': gr.const(lay['ICELLTYPE'].values, dtype=int),
            }


# %% === Gwfic ======== Initial condictions (head)

strthd = gr.const(0.)

Gwfic = {'strt': pr['strthd']}

# %% === Gwfchd ====== Given head cells
# Location of cells with given head


# %% === Gwfwel ===== Wells
# Location of wells

# %% === Gwfdrn ==== Drains

# Drain elevation
hDr = gr.Z[0, 0] - pr['drain_depth']

# Drdain conductance
Cdr = gr.Area[0] / pr['cDrainage']

# Layer of top active cells
Iz = gr.top_active_cells(IDOMAIN)

DRN = [((iz, 0, i), h_, C_) for i, (iz, h_, C_) in enumerate(zip(Iz, hDr, Cdr))]

Gwfdrn = {'stress_period_data': {0: DRN}}

# %% === Rch ===== Recharge

# Recharge into top active cells
RCH = [((iz, 0, i), pr['rch']) for i, iz in enumerate(Iz)]

Gwfrch = {'stress_period_data': {0: RCH}}

# %% === Gwfoc ==== output control

Gwfoc = {'head_filerecord':   os.path.join(dirs.SIM, "{}Gwf.hds".format(sim_name)),
         'budget_filerecord': os.path.join(dirs.SIM, "{}Gwf.cbc".format(sim_name)),
         'saverecord': [("HEAD", "ALL"), ("BUDGET", "ALL")],
}

print('Done mf_adapt')

if __name__ == '__main__':
   
   print(dirs)     