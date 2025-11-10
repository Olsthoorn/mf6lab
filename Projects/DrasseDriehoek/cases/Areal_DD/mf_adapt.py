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
import settings

INACTIVE =[-1, 0]

dirs = settings.dirs
sim_name = settings.sim_name
section_name = settings.section_name
props = settings.props
gr = settings.gr
params_wbk = settings.params_wbk # Parameters workbook

# === tdis ===== period data:
start_date_time = props['start_date_time'] # Must be a string.

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

# === DIS ===== Grid, structured

IDOMAIN = gr.const(1, dtype=int)
IDOMAIN[gr.DZ <= gr.min_dz] = -1

Gwfdis = {'gr': gr,
          'idomain': IDOMAIN,
          'length_units': settings.LENGTH_UNITS}

# ==== Gwfsto ===== Storage, transient

Gwfsto = {'sy': settings.sy,
          'ss': settings.ss,
          'iconvert': settings.icelltype,
          }

# === Gwfnpf ===== Horizontal and vertical conductivity

Gwfnpf = {  'k':   settings.k,
            'k33': settings.k33,
            'icelltype': settings.icelltype,
            }

#  === Gwfic ===== Initial condictions (head)

strthd = gr.const(props['strthd'])

Gwfic = {'strt': strthd}


# === Gwfwel ===== wells

# === Gwfdrn ===== drains
hDr = gr.Z[0, 0] - props['drain_depth']
drn_xyz = np.vstack((gr.xm, np.zeros(gr.nx), hDr)).T
Iz = np.fmax(gr.lrc_from_xyz(drn_xyz)['ic'][:, 0], gr.top_active_cells(IDOMAIN))
gr.top_active_cells(IDOMAIN, Iz)

Cdr = gr.Area[0] / props['cDrainage']
      
DRN = [((iz, 0, i), h_, C_) for i, (iz, h_, C_) in enumerate(zip(Iz, hDr, Cdr))]

Gwfdrn = {'stress_period_data': {0: DRN}}

# === Gwfchd ===== Fixed head period data
hL, hR = 5.5, 11.0
# CHD = ([(lrc, hL) for lrc in gr.I2LRC(gr.NOD[2:10, 0, 0])] + 
#        [(lrc, hR) for lrc in gr.I2LRC(gr.NOD[2:10, 0, -1])]
# )
CHD = [(lrc, hL) for lrc in gr.I2LRC(gr.NOD[2:10, 0, 0])]

Gwfchd = {'stress_period_data': {0: CHD}}

# === Gwfrch ===== recharge
RCH = [((iz, 0, i), props['rch']) for i, iz in enumerate(Iz)]

Gwfrch = {'stress_period_data': {0: RCH}}

# === Gwfoc ====== Output control for flow model

Gwfoc = {'head_filerecord':   os.path.join(dirs.SIM, "{}Gwf.hds".format(sim_name)),
         'budget_filerecord': os.path.join(dirs.SIM, "{}Gwf.cbc".format(sim_name)),
         'saverecord': [("HEAD", "ALL"), ("BUDGET", "ALL")],
}

print('Done mf_adapt')
