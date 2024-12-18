# Axially symmetric Xsection model from center of  Horstermeer plder
"""
A section model for Horsetermeer polder (NL) which is a local depression in the land surface of ca. 3 m.

The idea is to show the stream lines and the water table in the cross section,
which is fed by recharge and drained where the water table reaches ground surface.

@TO 2024-02-07
"""
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
           'time_units': pr['time_units'],
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

##### Below use gr_new instead of gr #######################

# === DIS ===== Grid, structured
new_params['idomain'][gr_new.DZ <= pr['dz_pinched_out']] = -1

Gwfdis = {'gr': gr_new,
          'idomain': new_params['idomain'],
          'length_units': pr['length_units']}

# ==== Gwfsto ===== Storage, transient
Gwfsto = {'sy': new_params['sy'],
          'ss': new_params['ss'],
          'iconvert':new_params['icelltype'],
          }

# === Gwfnpf ===== Horizontal and vertical conductivity
Gwfnpf = {  'k':   new_params['k'],
            'k33': new_params['k33'],
            'icelltype': new_params['icelltype'],
            'alternative_cell_averaging': 'AMT-LMK' if gr_new.axial else None,            
            }

# === Gwfic ======== Initial condictions (head)
strthd = gr_new.const(pr['hFar']) - pr['dDitch']

Gwfic = {'strt': strthd}

# === Gwfchd ===== fixed head at outer boundary
Ichd   = gr_new.NOD[:, 0, -1].flatten()
lrcCHD = gr_new.lrc_from_iglob(Ichd)

zNode = gr_new.ZM.ravel()[Ichd]
cNode = pr['cFresh'] * np.ones_like(zNode)
cNode[zNode <= pr['zIface'][-1]] = pr['cSalt']

stress_period_data = [(lrc, pr['hFar'], c) for lrc, c in zip(lrcCHD, cNode)]
Gwfchd = {'auxiliary': 'relconc',
          'stress_period_data': stress_period_data,
}

# === Gwfwel ===== wells
# Locatiion of well ring:
zW = gr_new.ZM[np.logical_and(gr_new.ZM[:, 0, 0] <= pr['wellTop'], gr_new.ZM[:, 0, 0] >= pr['wellBot']), 0, 0]
xW = np.ones_like(zW) * pr['RWellRing']
yW = np.zeros_like(zW)

Iwell   = gr_new.Iglob_from_xyz(np.vstack((xW, yW, zW)).T)
lrcWell = gr_new.lrc_from_iglob(Iwell)

Qw = pr['Qwell'] / gr_new.DZ.ravel()[Iwell]
stress_period_data = [(lrc, qw, pr['cFresh']) for lrc, qw in zip(lrcWell, pr['Qwell'] / (gr_new.DZ.ravel()[Iwell]))]

Gwfwel = {'auxiliary': 'relconc',
          'stress_period_data': stress_period_data,                
          }

# === Gwfdrn =====  Drains and ditches
xDr = np.arange(0, pr['R2'], 100.0)[1:]
yDr = np.zeros_like(xDr)
zDr = np.interp(xDr, gr_new.xm, gr_new.Z[0, 0, :]) - pr['dDitch']


Idrn = gr_new.Iglob_from_xyz(np.vstack((xDr, yDr, zDr)).T)

Cdr = gr_new.AREA.ravel()[Idrn] / pr['cDrainage']

DRN = [(lrc, h_, C_) for lrc, h_, C_ in zip(gr_new.lrc_from_iglob(Idrn), zDr, Cdr)]

Gwfdrn = {'stress_period_data': {0: DRN}}

# === Gwfrch =====
Iz = gr_new.top_active_cells(new_params['idomain'])
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
Gwtssm = {'sources': [['chd', 'AUX', 'relconc'],
                      ['wel', 'AUX', 'relconc']]}
# =====================================

print('Done mf_adapt')

if __name__ == '__main__':
   
   print(dirs)  
   
   # Start only with displacing fresh water (modeling transport, but leaving out the density)
   # Only then add denisty.
   # May be introduces a fixed salinity boundary at the shore.