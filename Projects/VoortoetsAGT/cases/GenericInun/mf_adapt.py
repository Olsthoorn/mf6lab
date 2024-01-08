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
from genericInunSectionData import gr
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

# IDOMAIN[gr.DZ < 0.25] = -1 # Does not work with Transport, We have to avoid pass-through cells.
params_dict = {'sy':  gr.const(lay['Ss'].values),
               'ss':  gr.const(lay['Ss'].values),
               'k':   gr.const(lay['k'].values),
               'k33': gr.const(lay['k33'].values),
               'icelltype': gr.const(lay['ICELLTYPE'].values),
               'idomain': gr.const(1, dtype=int),
            }

layers = np.array([(i, n) for i, n in enumerate(lay['Split'].values)])
gr_new, new_params = gr.refine_vertically(layers=layers, params_dict=params_dict)


# %% === DIS ========== Grid

Gwfdis = {'gr': gr_new,
          'idomain': new_params['idomain'],
          'length_units': settings.LENGTH_UNITS}


# %% ==== STO ===============
Gwfsto = {'sy': new_params['sy'],
          'ss': new_params['ss'],
          }

# %% === Gwfnpf =========== Horizontal and vertical conductivity

Gwfnpf = {  'k':   new_params['k'],
            'k33': new_params['k33'],
            'icelltype': new_params['icelltype'],
            }

# %% === Gwfic ======== Initial condictions (head)

strthd = gr_new.const(0.)

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

hDr = gr_new.Z[0, 0] - 0.5
Cdr = gr_new.Area[0] / c
# Put the drains and the recharge in the second layer or below if inactive
k = 0
Iz = np.zeros(gr_new.nx, dtype=int) + k

J = np.arange(gr_new.nx, dtype=int)
for i in range(k, gr_new.nz):
      z = gr_new.Z[i + 1, 0]
      J = J[np.logical_or(hDr[J] < z[J], new_params['idomain'][i, 0, J] < 0)]
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

# %% === Gwfbuy (boyancy) ====
FRESH, SALT, rhomin, rhomax = 0.0, 18000.0, 1000, 1025

irhospec = 0,
drhdc = (rhomax - rhomin) / SALT # kg/m3 / (mf/l) (conc FRESH to SALT)
crhoref = FRESH
modelname = sim_name + 'GWT'
auxspeciesname = None

Gwfbuy = {'nrhospecies': 1,
          'packagedata': [irhospec, drhdc, crhoref, modelname, auxspeciesname],
 }

Gwfbuy = {}

# %% ============ T R A N S P O R T ====================

# %% === Gwtfmi ===== Flow model interface
pd = [("GWFHEAD",   os.path.join(dirs.SIM, "{}Gwf.hds".format(sim_name))),
      ("GWFBUDGET", os.path.join(dirs.SIM, "{}Gwf.cbc".format(sim_name))),
]
Gwtfmi = {'packagedata': pd}

# %% === Gwtmst ===== Mobile storage and transfer

Gwtmst = {'porosity': 0.25}
 
# %% === Gwtadv === advection =========

Gwtadv = {'scheme' : 'TVD'} # choose from: upstream, central, TVD

# %% === Gwtdsp === dispersion ========
# diffc = 1e-10 m2/s 1e4 cm2/m2 60 s/min = 6e-5 m2/min
ahl, ath1, ath2, atv, diffc = 1.0, 0.1, 0.1, 0.1, 6e-5
# ahl, ath1, ath2, atv, diffc = 10.0, 1.0, 1.0, 0.1, 6e-5

Gwtdsp = {'alh': ahl, 'ath1': ath1, 'ath2': ath2, 'atv': atv, 'diffc': diffc}

# %% Gwtic === initial concentration ===

STRT = gr_new.const(FRESH)
STRT[gr_new.ZM < -30] = SALT
Gwtic = {'strt': STRT}

# %% Gwtcnc === const conc. ====

# Gwtcnc = {} # Use for sea conc

# %% Gwtoc === output control
Gwtoc = {
      'concentration_filerecord' : os.path.join(dirs.SIM, '{}Gwt.ucn'.format(sim_name)),
      'budget_filerecord':         os.path.join(dirs.SIM, '{}Gwt.cbc'.format(sim_name)),
      'saverecord' : [("CONCENTRATION", "ALL"),
                      ("BUDGET", "ALL")],
}

# %% Gwtssm === source-sink module

Gwtssm = {}

print('Done mf_adapt')

if __name__ == '__main__':
   
   print(dirs)  
   
   # Start only with displacing fresh water (modeling transport, but leaving out the density)
   # Only then add denisty.
   # May be introduces a fixed salinity boundary at the shore.