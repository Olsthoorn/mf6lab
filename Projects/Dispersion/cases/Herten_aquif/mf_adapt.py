'''Series 1 of MF6 Pathlines exercices

See http://www.citg.tudelft.nl/live/pagina.jsp?id=68e12562-a4d2-489a-b82e-deca5dd32c42&lang=en

https://flopy.readthedocs.io/en/stable/Notebooks/modpath7_create_simulation_example.html


@TO 240323 (using USGS's `flopy` and `Modflow 6` and Modpath 7)
'''

import os
import numpy as np
from src import mf6tools
import settings
from fdm.mfgrid import Grid

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

start_date_time = pr['start_date_time'][0] # Must be string

Simtdis = {'perioddata': period_data,
           'nper': len(period_data),
           'start_date_time': start_date_time,
           'time_units': settings.TIME_UNITS,
           }

# %% === Gwfdis ====== spatial discritization of flow model

x = np.arange(0, pr['L'] + pr['dx'] / 2, +pr['dx'])
z = np.arange(pr['H'], 0 - pr['dz'] / 2, -pr['dz'])
y = [-pr['D'] / 2, pr['D'] / 2]

gr = Grid(x, y, z)

# IDOMAIN
# Start all cells inactive
IDOMAIN = gr.const(1, dtype=int) # Inactive cells are vertially pass-through cells

Gwfdis = {'gr': gr,
          'idomain':      IDOMAIN,
          'length_units': settings.LENGTH_UNITS,
}

# %% === Gwfsto ======= Storage coefficients

Gwfsto = {'sy': gr.const(pr['sy']),
          'ss': gr.const(pr['ss']),
          'iconvert': 1,
          }

# %% === Gwfnpf ======= Horizontal and vertical conductivity
k = gr.const(pr['k'])

Gwfnpf = { 'k':   k,            
            'icelltype': pr['icelltype'],
}

# %% === Gwfic ======== Initial condictions (head)

hstrt = gr.const(0.) # Reasonalble hstrt
Gwfic = {'strt': hstrt}

# %% === Gwfchd ======= fixed head

IL = gr.NOD[:, 0,  0]
IR = gr.NOD[: ,0, -1]

stress_period_data = [(lrc, pr['hL'], pr['cL']) for lrc in gr.LRC(IL, astuples=True)] +\
                     [(lrc, pr['hR'], pr['cR']) for lrc in gr.LRC(IR, astuples=True)]
                     
Gwfchd ={'auxiliary': 'relconc',
         'stress_period_data': stress_period_data}

# %% === Gwfdrn === drains (can also be used for seepage face).
# TODO

# %% === Gwfoc ==== Output control for flow model
Gwfoc = {'head_filerecord':   os.path.join(dirs.SIM, "{}Gwf.hds".format(sim_name)),
         'budget_filerecord': os.path.join(dirs.SIM, "{}Gwf.cbc".format(sim_name)),
         'saverecord': [("HEAD", "FREQUENCY", pr['oc_frequency']), ("BUDGET", "FREQUENCY", pr['oc_frequency'])],
}

# %% ============ T R A N S P O R T ====================

#%% === Gwtdis ======= discritization for transport model

# Uses gr above, don't specify


# %% === Gwtfmi ===== Flow model interface

# simultaneous exchange is automatic (see mf_setup.py)

# %% === Gwtmst ===== Mobile storage and transfer

por = gr.const(pr['por'])
Gwtmst = {'porosity': por}

# %% === Gwtadv ====== advection

Gwtadv = {'scheme': 'TVD'} # upstream, central, TVD

# %% === Gwtdsp ====== dispersion & diffusion

Gwtdsp ={**pr['disp']}

# %% === Gwtic ====== initial concentration

Gwtic = {'strt': pr['cFresh']}

# %% === Gwtssm ====== Source-Sink mixing

Gwtssm = {'sources': [['chd', 'AUX', 'relconc']]}
Gwtssm = {}

# %% === Gwtcnc ====== constant concentraction

# Ink injection points (constant concentration at injection locations)

lrc = gr.lrc_from_iglob(gr.NOD[:, 0, 0])

# Concentration cells [(l, r, c) conc] of ink injection points.
concOn  = [(tuple(lrc_), pr['cSalt'])  for lrc_ in lrc]
concOff = [(tuple(lrc_), pr['cFresh']) for lrc_ in lrc]

stress_period_data = dict()
for isp, On in enumerate(perDF['On']):
      if On == 0:
            stress_period_data[isp] = concOff
      else:
            stress_period_data[isp] = concOn

Gwtcnc = {'stress_period_data': stress_period_data}

# %% === Gwtoc ====== output control

Gwtoc = {
      'concentration_filerecord' : os.path.join(dirs.SIM, '{}Gwt.ucn'.format(sim_name)),
      'budget_filerecord':         os.path.join(dirs.SIM, '{}Gwt.cbc'.format(sim_name)),
      'saverecord' : [("CONCENTRATION", "FREQUENCY", pr['oc_frequency']),
                      ("BUDGET", "FREQUENCY", pr['oc_frequency'])],
}

# %% === MP7 =========

Istart = np.NOD[:, :, 0]
nodes = gr.lrc_from_iglob(Istart)

MP7mp7 = {
      'trackdir': 'forward',
      'nodes': nodes, # list of nodes
      'rowcelldivision': 5,
      'columncelldivision': 5,
      'layercelldivision': 5,
}

# Output see files in directory MP7
# met extensies
# .mppth, .mpend .mp???


print('Done mf_adapt')

if __name__ == '__main__':
   
   print(dirs)     