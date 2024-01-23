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
Grid = settings.Grid

# === Join layers 0 and 1 and remove layer 3 for convergence of density model
lay = lay.iloc[[1, 2, 4]]
gr = Grid(gr.x, gr.y, gr.Z[[0, 2, 4, 5]], min_dz=gr.min_dz)

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

IDOMAIN = gr.const(1, dtype=int)

params_dict = {'sy':  gr.const(lay['Ss'].values),
               'ss':  gr.const(lay['Ss'].values),
               'k':   gr.const(lay['k'].values),
               'k33': gr.const(lay['k33'].values),
               'icelltype': gr.const(lay['ICELLTYPE'].values),
               'idomain': IDOMAIN,
            }

# Vertically refine the model layers
layers = np.array([(i, n) for i, n in enumerate(lay['Split'].values)])
gr_new, new_params = gr.refine_vertically(layers=layers, params_dict=params_dict)

# === DIS ===== Grid, structured
new_params['idomain'][gr_new.DZ <= pr['dz_pinched_out']] = -1

Gwfdis = {'gr': gr_new,
          'idomain': new_params['idomain'],
          'length_units': settings.LENGTH_UNITS}

# ==== Gwfsto ===== Storage, transient
Gwfsto = {'sy': new_params['sy'],
          'ss': new_params['ss'],
          'iconvert': 1,
          }

# === Gwfnpf ===== Horizontal and vertical conductivity

Gwfnpf = {  'k':   new_params['k'],
            'k33': new_params['k33'],
            'icelltype': new_params['icelltype'],
            }

# === Gwfic ======== Initial condictions (head)

strthd = gr_new.const(pr['hStrt'])

Gwfic = {'strt': strthd}

# === Gwfchd ===== fixed head

# === Gwfwel ===== wells

# === Gwfdrn =====
hDr = gr_new.Z[0, 0] - pr['drainDepth']
drn_xyz = np.vstack((gr_new.xm, np.zeros(gr_new.nx), hDr)).T
Iz = gr_new.lrc_from_xyz(drn_xyz)['ic'][:, 0]
gr_new.top_active_cells(IDOMAIN, Iz)

Cdr = gr_new.Area[0] / pr['cDrainage']

Iglob_wt = gr_new.Iglob_from_lrc(np.vstack((Iz, np.zeros(gr_new.nx, dtype=int), gr_new.NOD[0, 0])).T)    
DRN = [((iz, 0, i), h_, C_) for i, (iz, h_, C_) in enumerate(zip(Iz, hDr, Cdr))]

Gwfdrn = {'stress_period_data': {0: DRN}}

# === Gwfrch =====

RCH = [((iz, 0, i), pr['rch']) for i, iz in enumerate(Iz)]

Gwfrch = {'stress_period_data': {0: RCH}}

# === Gwfoc ======

Gwfoc = {'head_filerecord':   os.path.join(dirs.SIM, "{}Gwf.hds".format(sim_name)),
         'budget_filerecord': os.path.join(dirs.SIM, "{}Gwf.cbc".format(sim_name)),
         'saverecord': [("HEAD", "FREQUENCY", pr['oc_frequency']),
                        ("BUDGET", "FREQUENCY", pr['oc_frequency'])],
}

# === Gwfbuy =====  (boyancy)

irhospec = 0
drhdc = (pr['rhoSalt'] - pr['rhoFresh']) / (pr['cSalt'] - pr['cFresh'])
crhoref = pr['cFresh']
modelname = sim_name + 'GWT'
auxspeciesname = "relconc"

Gwfbuy = {'nrhospecies': 1,
          'denseref': pr['rhoFresh'],
          'density_filerecord': os.path.join(dirs.SIM, sim_name + 'Gwf.rho'),
          'packagedata': [irhospec, drhdc, crhoref, modelname, auxspeciesname],
 }

# ============ T R A N S P O R T ====================

# === Gwtfmi ===== Flow model interface
pd = [("GWFHEAD",   os.path.join(dirs.SIM, "{}Gwf.hds".format(sim_name))),
      ("GWFBUDGET", os.path.join(dirs.SIM, "{}Gwf.cbc".format(sim_name))),
]
Gwtfmi = {'packagedata': pd}

# === Gwtmst ===== Mobile storage and transfer

Gwtmst = {'porosity': pr['por']}
 
# === Gwtadv ===== advection =========

Gwtadv = {'scheme' : 'TVD'} # choose from: upstream, central, TVD

# === Gwtdsp ===== dispersion ========

Gwtdsp = {'alh': pr['ahl'],
          'ath1': pr['ath1'],
          'ath2': pr['ath2'],
          'atv': pr['atv'],
          'diffc': pr['diffc']}

# === Gwtic ===== initial concentration ===

cStart= gr_new.const(pr['cFresh'])
cStart[gr_new.ZM < pr['zIface']] = pr['cSalt']
Gwtic = {'strt': cStart}

# === Gwtcnc ===== const conc.


# === Gwtoc ===== output control
Gwtoc = {
      'concentration_filerecord' : os.path.join(dirs.SIM, '{}Gwt.ucn'.format(sim_name)),
      'budget_filerecord':         os.path.join(dirs.SIM, '{}Gwt.cbc'.format(sim_name)),
      'saverecord' : [("CONCENTRATION", "FREQUENCY", pr['oc_frequency']),
                      ("BUDGET", "FREQUENCY", pr['oc_frequency'])],
}

# === Gwtssm ====== source sink mixing

Gwtssm = {}

# =====================================

print('Done mf_adapt')

if __name__ == '__main__':
   
   print(dirs)  
   
   # Start only with displacing fresh water (modeling transport, but leaving out the density)
   # Only then add denisty.
   # May be introduces a fixed salinity boundary at the shore.