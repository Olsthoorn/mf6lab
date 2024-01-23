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

dirs = settings.dirs
sim_name = settings.sim_name
section_name = settings.section_name
pr = settings.props
lay = settings.lay
gr = settings.gr
params_wbk = settings.params_wbk

# === tdis ===== period data:
start_date_time = pr['start_date_time'] # Must be a string.

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
IDOMAIN[gr.DZ < pr['minDz']] = -1

Gwfdis = {'gr': gr,
          'idomain': IDOMAIN,
          'length_units': settings.LENGTH_UNITS}

# ==== Gwfsto ===== Storage, transient

Gwfsto = {'sy': gr.const(lay['Sy'].values),
          'ss': gr.const(lay['Ss'].values),
          'iconvert': gr.const(lay['ICELLTYPE'].values),
          }

# === Gwfnpf ===== Horizontal and vertical conductivity

Gwfnpf = {  'k':   gr.const(lay['k'].values),
            'k33': gr.const(lay['k33'].values),
            'icelltype': gr.const(lay['ICELLTYPE'].values),
            }

#  === Gwfic ===== Initial condictions (head)

strthd = gr.const(pr['strthd'])

Gwfic = {'strt': strthd}

# === Gwfchd ===== Fixed head period data

# === Gwfwel ===== wells

# === Gwfdrn ===== drains
hDr = gr.Z[0, 0] - pr['drain_depth']
drn_xyz = np.vstack((gr.xm, np.zeros(gr.nx), hDr)).T
Iz = gr.lrc_from_xyz(drn_xyz)['ic'][:, 0]
gr.top_active_cells(IDOMAIN, Iz)

Cdr = gr.Area[0] / pr['cDrainage']
      
DRN = [((iz, 0, i), h_, C_) for i, (iz, h_, C_) in enumerate(zip(Iz, hDr, Cdr))]

Gwfdrn = {'stress_period_data': {0: DRN}}

# === Gwfrch ===== recharge
rch = perDF['rch']

stress_period_data = dict()
for iper in rch.index:
    stress_period_data[iper] = [((iz, 0, i), rch[iper]) for i, iz in enumerate(Iz)]

Gwfrch = {'stress_period_data': stress_period_data}

# === Gwfoc ====== Output control for flow model

Gwfoc = {'head_filerecord':   os.path.join(dirs.SIM, "{}Gwf.hds".format(sim_name)),
         'budget_filerecord': os.path.join(dirs.SIM, "{}Gwf.cbc".format(sim_name)),
         'saverecord': [("HEAD", "FREQUENCY", pr['oc_frequency']),
                        ("BUDGET", "FREQUENCY", pr['oc_frequency'])],
}

print('Done mf_adapt')
