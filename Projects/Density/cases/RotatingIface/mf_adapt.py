"""
# The Henry Saltwater Intrucion Problem

The well known Henry problem washes out salt water from a vertical cross section. It's a standard problem allowing analytical verification and a good test of the software.

References:
URL: http://downloads.geo-slope.com/geostudioresources/examples/8/0/CtranW/Henry%20Density%20Dependent.pdf

Henry, H. R. 1964. Effects of dispersion on salt encroachment in coastal aquifers. Sea Water in Coastal Aquifers, U.S. Geol. Surv. Supply Pap., 1613-C, C71-C84.

Simpson, M.J. and Clement, T.B. 2004. Improving the worthiness of the Henry problem as a benchmark for density-dependent groundwater flow models. Water Resources Research, 40 (W01504).

@TO 20240109
"""

import os
import numpy as np
from src import mf6tools
from fdm.mfgrid import Grid
import settings

dirs = settings.dirs
sim_name = settings.sim_name
section_name = settings.section_name
pr = settings.props

# Parameters workbook
params_wbk = settings.params_wbk

## Get section data

# %% === tdis ==========  Period Data:
start_date_time = '2024-01-11' # Must be a string.

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

dx, dz, L, H =  pr['dx'], pr['dz'], pr['L'], pr['H']
x = np.arange(-dx/2, L + dx, dx)
z = np.arange(H, -dz/2, -dz) - H / 2
y = [-0.5, 0.5]
gr = Grid(x, y, z)

Gwfdis = {'gr': gr,
          'idomain':      gr.const(1, dtype=int),
          'length_units': settings.LENGTH_UNITS}

# %% ==== STO ===============
Gwfsto = {'sy': gr.const(lay['Ss'].values),
          'ss': gr.const(lay['Ss'].values),
          }

# %% === Gwfnpf =========== Horizontal and vertical conductivity
Gwfnpf = {  'k':   gr.const(pr['k']),            
            'icelltype': gr.const(lay['ICELLTYPE'].values),
            }

# %% === Gwfic ======== Initial condictions (head)
strthd = gr.const(pr['hStrt'])

Gwfic = {'strt': strthd}

# %% === CHD, Fixed head period data (Only specify the first period)

if True:
      IdxL = gr.NOD[-1, :,  -1].flatten()
      IdxR = np.array([])
else:
      IdxL = gr.NOD[:, :,  0].flatten()
      IdxR = gr.NOD[:, :, -1].flatten()
stress_period_data = [(lrc, pr['hStrt'], pr['cL']) for lrc in gr.LRC(IdxL)] +\
       [(lrc, pr['hStrt'], pr['cR']) for lrc in gr.LRC(IdxR)]
Gwfchd ={'auxiliary': 'relconc',
         'stress_period_data': stress_period_data}

# %% === WEL ====
          
# %% === DRN,

# %% === Rch

# %% === OC ==== output control for Gwf model

Gwfoc = {'head_filerecord':   os.path.join(dirs.SIM, "{}Gwf.hds".format(sim_name)),
         'budget_filerecord': os.path.join(dirs.SIM, "{}Gwf.cbc".format(sim_name)),
         'saverecord': [("HEAD", *pr['oc_frequency']), ("BUDGET", *pr['oc_frequency'])],
}

# %% === Gwfbuy (boyancy) ====
irhospec = 0
modelname = sim_name + 'GWT'
auxspeciesname = "relconc"

Gwfbuy = {'nrhospecies': 1,
          'denseref': pr['rhoref'],
          'density_filerecord': os.path.join(dirs.SIM, sim_name + 'Gwf.rho'),
          'packagedata': [irhospec, pr['drhodc'], pr['crhoref'], modelname, auxspeciesname],
 }

# %% ============ T R A N S P O R T ====================

# %% === Gwtfmi ===== Flow model interface
pd = [("GWFHEAD",   os.path.join(dirs.SIM, "{}Gwf.hds".format(sim_name))),
      ("GWFBUDGET", os.path.join(dirs.SIM, "{}Gwf.cbc".format(sim_name))),
]
Gwtfmi = {'packagedata': pd}

# %% === Gwtmst ===== Mobile storage and transfer

Gwtmst = {'porosity': pr['por']}
 
# %% === Gwtadv === advection =========

Gwtadv = {'scheme' : 'TVD'} # choose from: upstream, central, TVD

# %% === Gwtdsp === dispersion ========

Gwtdsp = {**settings.props['disp']}

# %% Gwtic === initial concentration ===
z0 = (gr.z[0] + gr.z[-1]) / 2
def zIRM(x):
      global z0, pr
      return z0 - pr['IeRM'] / pr['H'] * (x - pr['xmRM'])
def zIML(x):
      global z0, pr
      return z0 - pr['IeML'] / pr['H'] * (x - pr['xmML'])

strtconc = gr.const(pr['cL'])
strtconc[gr.ZM > zIML(gr.XM)] = pr['cM']
strtconc[gr.ZM > zIRM(gr.XM)] = pr['cR']

Gwtic = {'strt': strtconc}

# %% Gwtcnc === const conc. ====

# %% Gwtoc === output control for Gwt model
Gwtoc = {
      'concentration_filerecord' : os.path.join(dirs.SIM, '{}Gwt.ucn'.format(sim_name)),
      'budget_filerecord':         os.path.join(dirs.SIM, '{}Gwt.cbc'.format(sim_name)),
      'saverecord' : [("CONCENTRATION", *pr['oc_frequency']), ("BUDGET", *pr['oc_frequency'])],
}

# %% Gwtssm === source-sink module

Gwtssm = {'sources': [['chd', 'AUX', 'relconc']]}

print('Done mf_adapt')

if __name__ == '__main__':
   
   print(dirs)  
   
   # Start only with displacing fresh water (modeling transport, but leaving out the density)
   # Only then add denisty.
   # May be introduces a fixed salinity boundary at the shore.