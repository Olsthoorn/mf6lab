# Analytic solutions for exatraction on linear cross section (Edelman, 1940 and Carslaw & Jaeger, 1959)
"""
Extration at x=0 from a half-infitte cross section is the perfect copanion of the Theis solutiion in the radial case.
Various cases have been worked out by Edelman (1940), Carslaw & Jaeger > 1928 (1959) and can also be found in the
book of Bruggeman (1999). Th elatter is probably the best reference as all relevant  cases are dealt with by Bruggeman (199).

The idea is to show the stream lines and the water table in the cross section,
"""

# TO 20240124

import os
import numpy as np
from src import mf6tools
import settings

INACTIVE =[-1, 0]

dirs = settings.dirs
sim_name = settings.sim_name
section_name = settings.section_name
lay = settings.lay
pr = settings.props
gr = settings.gr
params_wbk = settings.params_wbk # Parameters workbook

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

Gwfsto = {'sy': gr.const(lay['Sy'].values, axial=False),
          'ss': gr.const(lay['Ss'].values, axial=False),
          'iconvert': gr.const(lay['ICELLTYPE'].values),
          }

# === Gwfnpf ===== Horizontal and vertical conductivity

Gwfnpf = {  'k':   gr.const(lay['k'].values,   axial=False),
            'k33': gr.const(lay['k33'].values, axial=False),
            'icelltype': gr.const(lay['ICELLTYPE'].values),
            }

#  === Gwfic ===== Initial condictions (head)

strthd = gr.const(pr['strthd'])

Gwfic = {'strt': strthd}

# === Gwfwel ===== wells

# Well only where screened
Iglob_well = gr.NOD[:, 0, 0][lay['Screened'] > 0]
stress_period_data = [(lrc, pr['Q']) for lrc in gr.lrc_from_iglob(Iglob_well)]
Gwfwel = {'stress_period_data': stress_period_data,
          'auto_flow_reduce': 0.1,
}

# Wells everywhere in the lowest layer to generate upward seepage
# Qseep = gr.dx * 0.1
# stress_period_data = [(lrc, q) for lrc, q in zip(gr.lrc_from_iglob(gr.NOD[-1].ravel()), Qseep)]
# Gwfwel = {'stress_period_data': stress_period_data,
#           'auto_flow_reduce': 0.1,
# }

# === Gwfdrn ===== drains
hDr = gr.Z[0, 0] - pr['drain_depth']
drn_xyz = np.vstack((gr.xm, np.zeros(gr.nx), hDr)).T
Iz = gr.lrc_from_xyz(drn_xyz)['ic'][:, 0]
gr.top_active_cells(IDOMAIN, Iz)

Cdr = gr.Area[0] / pr['cDrainage']
      
DRN = [((iz, 0, i), h_, C_) for i, (iz, h_, C_) in enumerate(zip(Iz, hDr, Cdr))]

Gwfdrn = {'stress_period_data': {0: DRN}}

# === Gwfchd ===== Fixed head period data

# chd instead of drains
CHD = [((iz, 0, i), pr['strthd']) for i, iz in enumerate(Iz)]
Gwfchd = {'stress_period_data': CHD}

# === Gwfrch ===== recharge

RCH = [((iz, 0, i), pr['rch']) for i, iz in enumerate(Iz)]

Gwfrch = {'stress_period_data': {0: RCH}}

# === Gwfoc ====== Output control for flow model

Gwfoc = {'head_filerecord':   os.path.join(dirs.SIM, "{}Gwf.hds".format(sim_name)),
         'budget_filerecord': os.path.join(dirs.SIM, "{}Gwf.cbc".format(sim_name)),
         'saverecord': [("HEAD", "FREQUENCY", pr['oc_frequency']),
                        ("BUDGET", "FREQUENCY", pr['oc_frequency'])],
}

print('Done mf_adapt')
