"""Wims aquifer

Will generate an aquifer with facies to simulate dispersion as Wim de Lange see it.

The Idea is to generate a 1D cross section with a geostatistical structure that mimics
how Wim sees the aquifer, so that dispersion with it can be simulated and the results
can be compared with the theory Wim has developed and published.

@ TO 2025/04/10
"""

# %%
import os
import numpy as np
from src import mf6tools
import settings

from importlib import reload


reload(settings)
dirs = settings.dirs
sim_name = settings.sim_name
section_name = settings.section_name

# Parameters Excel workbook
params_wbk = settings.params_wbk

# Parameters of the current case model and simulation
pr = settings.props

# Get the data for the packages specific to this simulation.
# Note that the params.wbk sheet NAM determines which packages are used
# and so for which packages data must be provided below.
# Default parameters come from params_wbk.

# %% === tdis ======  time discretization

perDF = mf6tools.get_periodata_from_excel(params_wbk, sheet_name='PER')
period_data = [tuple(sp) for sp in perDF[['PERLEN', 'NSTP', 'TSMULT']].values]

start_date_time = pr['start_date_time'] # Must be string

Simtdis = {'perioddata': period_data,
           'nper': len(period_data),
           'start_date_time': start_date_time,
           'time_units': settings.TIME_UNITS,
           }

# %% === Gwfdis ====== spatial discritization of flow model

gr = settings.props['gr']

# IDOMAIN: Start all cells inactive vertically flow-through cells
idomain = gr.const(1, dtype=int) # Inactive cells are vertially pass-through cells

Gwfdis = {'gr': gr,
          'idomain': idomain,
          'length_units': settings.LENGTH_UNITS,
}

# %% === Gwfsto ======= Storage coefficients

Gwfsto = {'sy': gr.const(pr['sy']),
          'ss': gr.const(pr['ss']),
          'iconvert': 1,
          }

# %% === Gwfnpf ======= Horizontal and vertical conductivity
k = pr['k_field']
Gwfnpf = { 'k':   k,            
            'icelltype': pr['icelltype'],
}

# %% === Gwfic ======== Initial condictions (head)

hstrt = gr.const(pr['hStrt'])
Gwfic = {'strt': hstrt}

# %% === Gwfchd ======= fixed head

hL, hR, relconc = pr['hL'], pr['hR'], 0.01
chd_ = ([((k, i, j), hL,  relconc) for (k, i, j) in gr.lrc_from_iglob(gr.NOD[:, 0, 0],  astuples=True)] +
        [((k, i, j), hR, 0.) for (k, i, j) in gr.lrc_from_iglob(gr.NOD[:, 0, -1], astuples=True)]
)

Gwfchd ={'auxiliary': 'relconc',
         'stress_period_data': {0: chd_}
         }

# %% === Gwfdrn === drains (can also be used for seepage face).


# %% === Gwfwel === wel (milk injection point)


# %% === Gwfrch === drains (can also be used for seepage face).

# %% === Gwfbuy (boyancy) ====

# %% === Gwfoc ==== Output control for flow model
Gwfoc = {'head_filerecord':   os.path.join(dirs.SIM, "{}Gwf.hds".format(sim_name)),
         'budget_filerecord': os.path.join(dirs.SIM, "{}Gwf.cbc".format(sim_name)),
         'saverecord': [("HEAD",   "FREQUENCY", pr['oc_frequency']),
                        ("BUDGET", "FREQUENCY", pr['oc_frequency'])],
}

# %% ============ T R A N S P O R T ====================

#%% === Gwtdis ======= discritization for transport model

# Uses gr above, do not specify

# %% === Gwtfmi ===== Flow model interface

# simultaneous exchange is automatic (see mf_setup.py)

# %% === Gwtmst ===== Mobile storage and transfer

por = gr.const(pr['por'])

Gwtmst = {'porosity': por}

# %% === Gwtadv ====== advection

Gwtadv = {'scheme': 'TVD'} # upstream, central, TVD

# %% === Gwtdsp ====== dispersion & diffusion

# Gwtdsp ={**pr['disp']}

# %% === Gwtic ====== initial concentration

# cstrt = gr.const(pr['cFresh'])

# %% === Gwtssm ====== Source-Sink mixing

Gwtssm = {'sources': [['chd', 'AUX', 'relconc'],
                      ['wel', 'AUX', 'relconc']]}

Gwtssm = {'sources': [['chd', 'AUX', 'relconc']]}


# %% === Gwtcnc ====== constant concentraction

# %% === Gwtoc ====== output control

Gwtoc = {
      'concentration_filerecord' : os.path.join(dirs.SIM, '{}Gwt.ucn'.format(sim_name)),
      'budget_filerecord':         os.path.join(dirs.SIM, '{}Gwt.cbc'.format(sim_name)),
      'saverecord' : [("CONCENTRATION", "FREQUENCY", pr['oc_frequency']),
                      ("BUDGET",        "FREQUENCY", pr['oc_frequency'])],
}
print('Done mf_adapt')

if __name__ == '__main__':
   
   print(dirs)     
