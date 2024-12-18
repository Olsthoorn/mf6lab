"""Compuaring MF6 with fdm3t for Hantush's well function.
2024-12-14

A number of Hantush cases is computed and compared with the analytical solution as well as the result from fdm3t.
The cases are defined in the module fdm3t.py under tools/fdm.

@TO 241214 (in mf6lab using USGS's `flopy` and `Modflow 6`)
"""

import os
import numpy as np
from src import mf6tools
import settings

dirs = settings.dirs
sim_name = settings.sim_name
section_name = settings.section_name

# Parameters workbook
params_wbk = settings.params_wbk

# %% Parameters of the model and simulation
pr = settings.props

## Get section data

# %% === tdis ======  time discretization

perDF = mf6tools.get_periodata_from_excel(params_wbk, sheet_name='PER')
period_data = [tuple(sp) for sp in perDF[['PERLEN', 'NSTP', 'TSMULT']].values]

start_date_time = perDF['StartTime'][0] # Must be string

period_data = [(dt_, 1, 1.25) for dt_ in np.diff(settings.t)]

Simtdis = {'perioddata': period_data,
           'nper': len(period_data),
           'start_date_time': start_date_time,
           'time_units': settings.TIME_UNITS,
           }

# %% === Gwfdis ====== spatial discritization of flow model

gr = settings.gr

Gwfdis = {'gr': gr,
          'idomain':      settings.IDOMAIN,
          'length_units': settings.LENGTH_UNITS,
}

# %% === Gwfsto ======= Storage coefficients


Gwfsto = {'sy': settings.ss * 0.,
          'ss': settings.ss,
          'iconvert': 0,
          }

# %% === Gwfnpf ======= Horizontal and vertical conductivity

Gwfnpf = { 'k':   settings.kr,
           'k22': 1e-20,
           'k33': settings.kz,
           'icelltype': 0.,
            'alternative_cell_averaging': 'LOGARITHMIC',
}

# %% === Gwfic ======== Initial condictions (head)

hstrt = gr.const(0.)

Gwfic = {'strt': hstrt}

# %% === Gwfchd ======= fixed head

# Iwell = gr.NOD[-1, :, 0]
# 
# stress_period_data = [(lrc, pr['hCanL'], pr['cCanL']) for lrc in gr.LRC(IcanL, astuples=True)] +\
#                      [(lrc, pr['hCanR'], pr['cCanR']) for lrc in gr.LRC(IcanR, astuples=True)]
#                      
# Gwfchd ={'auxiliary': 'relconc',
#          'stress_period_data': stress_period_data}

# %% === Gwfghb ======= fixed head
# stress_period_data ([cellid, bhead, cond, aux, boundname]) 

c = ((settings.r_ /np.array(settings.rhos)) ** 2  / settings.kD)[:, np.newaxis]
C = (settings.AREA[0, :, :] / c).flatten()

cellid = gr.I2LRC(gr.NOD[0].ravel())
hbead = hstrt[0].flatten()
Gwfghb = {'stress_period_data':
            {0: [(cid, h_, c_) for cid, h_, c_ in zip(cellid, hbead, C)]}
}

# %% === Gwfwel === wells (can also be used for seepage face).
#
Q = 4 * np.pi * settings.kD * np.ones(gr.ny)
cellid = gr.I2LRC(gr.NOD[-1, :, 0])
Gwfwel = {'stress_period_data': 
            {0: [(cid, q) for cid, q in zip(cellid, Q)]}
}

# %% === Gwfoc ==== Output control for flow model
oc_frequency = 1

Gwfoc = {'head_filerecord':   os.path.join(dirs.SIM, f"{sim_name}Gwf.hds"),
         'budget_filerecord': os.path.join(dirs.SIM, f"{sim_name}Gwf.cbc"),
         'saverecord': [("HEAD", "FREQUENCY", settings.oc_frequency), ("BUDGET", "FREQUENCY", settings.oc_frequency)],
}

# %% ============ T R A N S P O R T ====================

#%% === Gwtdis ======= discritization for transport model

# Uses gr above, don't specify


# %% === Gwtfmi ===== Flow model interface

# simultaneous exchange is automatic (see mf_setup.py)

# %% === Gwtmst ===== Mobile storage and transfer

por = 0.35

Gwtmst = {'porosity': por}

# %% === Gwtadv ====== advection

Gwtadv = {'scheme': 'TVD'} # upstream, central, TVD

# %% === Gwtdsp ====== dispersion & diffusion


# %% === Gwtic ====== initial concentration

# Gwtic = {'strt': pr['cFresh']}

# %% === Gwtssm ====== Source-Sink mixing

Gwtssm = {'sources': [['chd', 'AUX', 'relconc']]}
Gwtssm = {}

# %% === Gwtcnc ====== constant concentraction


# Concentration cells [(l, r, c) conc] of ink injection points.
#concOn  = [(tuple(lrc_), pr['cInk']) for lrc_ in lrc]
#concOff = [(tuple(lrc_),  0.0)       for lrc_ in lrc]
#
#stress_period_data = dict()
#for isp, ink in enumerate(perDF['Ink']):
#      if ink == 0:
#            stress_period_data[isp] = concOff
#      else:
#            stress_period_data[isp] = concOn
#
#Gwtcnc = {'stress_period_data': stress_period_data}

# %% === Gwtoc ====== output control

Gwtoc = {
      'concentration_filerecord' : os.path.join(dirs.SIM, '{}Gwt.ucn'.format(sim_name)),
      'budget_filerecord':         os.path.join(dirs.SIM, '{}Gwt.cbc'.format(sim_name)),
      'saverecord' : [("CONCENTRATION", "FREQUENCY", settings.oc_frequency),
                      ("BUDGET", "FREQUENCY", settings.oc_frequency)],
}
print('Done mf_adapt')

if __name__ == '__main__':
   
   print(dirs)     