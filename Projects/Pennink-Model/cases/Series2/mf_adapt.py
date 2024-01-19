'''Series 2 of Pennink(1915)'s 1903 - 1905 sand box model experiments.
Experiments series 2 (Septemner 1904)

See http://www.citg.tudelft.nl/live/pagina.jsp?id=68e12562-a4d2-489a-b82e-deca5dd32c42&lang=en

In this experiment, Pennink(1915) studies freshwater flow between a recharge
canal on the right of the model to an extraction canal at the left.

Pennink then injects ink at different depths near the recharge canal and
shows in a series of photos the movement of the ink.

The flow is from the right to the left canal and the flow is
1.25 L/h (Pennink (1905) p485 with a (water table) gradient of 1:130

To easier match the actual times of the tests, we will use time in minutes
instead of days.

@TO 090312 (in mflab using Seawat)
@TO 240112 (in mf6lab using USGS's `flopy` and `Modflow 6`)
'''

import os
import numpy as np
from src import mf6tools
import settings
from fdm.mfgrid import Grid, index

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
IDOMAIN = gr.const(-1, dtype=int) # Inactive cells are vertially pass-through cells

# Only where sand --> active
IDOMAIN[:, 0, :][gr.inpoly(pr[  'sand'], row=0)] =  1
IDOMAIN[:, 0, :][gr.inpoly(pr['canalL'], row=0)] = pr['IDCL'] # Cells in Canal Left
IDOMAIN[:, 0, :][gr.inpoly(pr['canalR'], row=0)] = pr['IDCR'] # Cells in Canal Right

# Inactive cells are within canals but above their stage level
IDOMAIN[np.logical_and(IDOMAIN == pr['IDCL'], gr.ZM > pr['hCanL'])] = 0
IDOMAIN[np.logical_and(IDOMAIN == pr['IDCR'], gr.ZM > pr['hCanR'])] = 0

Gwfdis = {'gr': gr,
          'idomain':      IDOMAIN,
          'length_units': settings.LENGTH_UNITS,
}

# %% === Gwfsto ======= Storage coefficients

Gwfsto = {'sy': gr.const(pr['sy']),
          'ss': gr.const(pr['ss']),
          }

# %% === Gwfnpf ======= Horizontal and vertical conductivity
k = gr.const(pr['k'])

k[gr.ZM > pr['zCapZone']] /= 5 # Unsat zone cond. above full capillary zone

Gwfnpf = { 'k':   k,            
            'icelltype': pr['icelltype'],
}

# %% === Gwfic ======== Initial condictions (head)

hstrt = pr['hCanL'] + gr.XM / pr['L'] * (pr['hCanR'] - pr['hCanL']) # Reasonalble hstrt
hstrt[IDOMAIN == pr['IDCL']] = pr['hCanL'] # in left canal
hstrt[IDOMAIN == pr['IDCR']] = pr['hCanR'] # in right canal

Gwfic = {'strt': hstrt}

# %% === Gwfchd ======= fixed head

IcanL = gr.NOD[IDOMAIN == pr['IDCL']]
IcanR = gr.NOD[IDOMAIN == pr['IDCR']]

stress_period_data = [(lrc, pr['hCanL'], pr['cCanL']) for lrc in gr.LRC(IcanL, astuples=True)] +\
                     [(lrc, pr['hCanR'], pr['cCanR']) for lrc in gr.LRC(IcanR, astuples=True)]
                     
Gwfchd ={'auxiliary': 'relconc',
         'stress_period_data': stress_period_data}

# %% === Gwfdrn === drains (can also be used for seepage face).
# TODO

# %% === Gwfoc ==== Output control for flow model
Gwfoc = {'head_filerecord':   os.path.join(dirs.SIM, "{}Gwf.hds".format(sim_name)),
         'budget_filerecord': os.path.join(dirs.SIM, "{}Gwf.cbc".format(sim_name)),
         'saverecord': [("HEAD", "FREQUENCY", pr['oc_frequency']),
                        ("BUDGET", "FREQUENCY", pr['oc_frequency'])],
}

# %% ============ T R A N S P O R T ====================

#%% === Gwtdis ======= discritization for transport model

# Uses gr above, don't specify


# %% === Gwtfmi ===== Flow model interface

# simultaneous exchange is automatic (see mf_setup.py)

# %% === Gwtmst ===== Mobile storage and transfer

por = gr.const(pr['por'])
por[gr.ZM > pr['zCapZone']] /= 3 # My estimate for the unsaturated zone

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

lrc  = gr.lrc_from_xyz(pr['xyzInk'])['ic'] # global coords of injection points

IDOMAIN.ravel()[gr.Iglob_from_lrc(lrc)] = pr['iInk'] # mark the cells in IDOMAIN where ink injection takes place

# Concentration cells [(l, r, c) conc] of ink injection points.
concOn  = [(tuple(lrc_), pr['cInk']) for lrc_ in lrc]
concOff = [(tuple(lrc_),  0.0)       for lrc_ in lrc]

stress_period_data = dict()
for isp, ink in enumerate(perDF['Ink']):
      if ink == 0:
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
print('Done mf_adapt')

if __name__ == '__main__':
   
   print(dirs)     