"""Pennink - density flow in Pennnink's (1915) sand box model
 -- SIXTH SERIES OF EXPERIMENTS (MARCH 1905)

see http://www.citg.tudelft.nl/live/pagina.jsp?id=68e12562-a4d2-489a-b82e-deca5dd32c42&lang=en

In his sixth experiment, Pennink (1915) studies simultaneous flow of freshwater and
saltwater flow from a canal at the right and one atthe left corner of the sand-box model
to a well in the center.

The following text is cited from his book, p87:

"C. Experiments like sub B, but central-canal as a means of concentration
substituted by a well-pipe with gauze strainer, placed in the canal axis.
 Sixth series of experiments (March 1905)
 
The water is discharged by means of the well-pipe, by siphoning, on the
left side of the apparatus. Moreover, glass gauging-pipes are place in
the sand in order to determine the pressure of the water just above the
milk, and beside the gause strainer."

The exact construction of the "well-pipe" as not discribed by Pennink, neither  
was its depth and screen length, given it was gauze strainer, it is assumed
that the well is fully screened to a depth decuced from the photos in his book.

We only model the situaton after Pennink pulled up his well screen to its
final position of 32 cm above the original interface.

Erroneous milk-head values may cause the milk to flow fast into or out of the model.
Measuring directly from the photo was difficult and computing was not accurate
enough due to vertical gradients in the sand above the interface.
Therefore, the head in the milk reservoir was fine-tuned in several trial runs.
The applied head values are in the workbook on sheet PER.

@TO 090320 100520
@TO 20220118
"""

import os
import numpy as np
from src import mf6tools
from fdm.mfgrid import Grid, index
import settings

dirs = settings.dirs
sim_name = settings.sim_name
section_name = settings.section_name

# Parameters Excel workbook
params_wbk = settings.params_wbk

# Parameters of the cuuent case model and simulation
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

x = np.arange(0, pr['L'] + pr['dx'] / 2, +pr['dx'])
z = np.arange(pr['H'], 0 - pr['dz'] / 2, -pr['dz'])
y = [-pr['D'] / 2, pr['D'] / 2]

gr = Grid(x, y, z)

# IDOMAIN: Start all cells inactive vertically flow-through cells
IDOMAIN = gr.const(-1, dtype=int) # Inactive cells are vertially pass-through cells

IDOMAIN[:, 0, :][gr.inpoly(pr[  'sand'], row=0)] = pr['IDSD']
IDOMAIN[:, 0, :][gr.inpoly(pr['canalL'], row=0)] = pr['ICNL'] # Cells in Canal Left
IDOMAIN[:, 0, :][gr.inpoly(pr['canalM'], row=0)] = pr['ICNM'] # Cells in Canal Center
IDOMAIN[:, 0, :][gr.inpoly(pr['canalR'], row=0)] = pr['ICNR'] # Cells in Canal Right
IDOMAIN[gr.ZM < pr['zIface']] = pr['IDMK']

lrcMilk  = gr.lrc_from_xyz(pr['milkInjPnt']) # global coords milk injection point
IDOMAIN.ravel()[gr.Iglob(lrcMilk['ic'])] = pr['iMlkInjPnt'] # mark milk injection cells

# Location of well in the grid
ixwell = gr.NOD[0, 0][gr.x[:-1] < pr['L'] / 2][-1] # Center ix coordinate

zScreenPoints = np.arange(pr['zScreen'][0], pr['zScreen'][1], -0.1)
Izwell = np.unique(index(zScreenPoints, gr.Z[:, 0, ixwell])['idx'])
Iywell = np.zeros(len(Izwell), dtype=int)
Ixwell = np.ones( len(Izwell), dtype=int) * ixwell

lrcWell = np.vstack((Izwell, Iywell, Ixwell)).T
IDOMAIN.ravel()[gr.Iglob_from_lrc(lrcWell)] = pr['IWEL']

Gwfdis = {'gr': gr,
          'idomain': IDOMAIN,
          'length_units': settings.LENGTH_UNITS,
}

# %% === Gwfsto ======= Storage coefficients

Gwfsto = {'sy': gr.const(pr['sy']),
          'ss': gr.const(pr['ss']),
          }

# %% === Gwfnpf ======= Horizontal and vertical conductivity
k = gr.const(pr['k'])

Gwfnpf = { 'k':   k,            
            'icelltype': pr['icelltype'],
}

# %% === Gwfic ======== Initial condictions (head)

hstrt = gr.const(pr['hStrt'])

Gwfic = {'strt': hstrt}

# %% === Gwfchd ======= fixed head

Ix = gr.NOD[0, 0]
IzTopActive = gr.top_active_cells(IDOMAIN)[gr.NOD[0, 0]]

zTop = np.interp(gr.xm, gr.x[[0, int((gr.nx - 1)/2), -1]], pr['sand'][:, 1][2:5])
lrcTopActive = np.vstack((IzTopActive,                 # iz
                          np.zeros(gr.nx, dtype=int),  # iy
                          gr.NOD[0, 0])).T             # ix

# Salt water head of milk injection point
hMilkPnt = pr['zIface'] + (pr['hStrt'] - pr['zIface']) * pr['rhoFresh'] / pr['rhoSalt']

milkOn  = [(tuple(lrc_), hMilkPnt, pr['cSalt']) for lrc_ in lrcMilk['ic']]
milkOff = [(tuple(lrc_), hMilkPnt, pr['cFresh']) for lrc_ in lrcMilk['ic']]

lrcCanL = gr.lrc_from_iglob(gr.NOD[IDOMAIN == pr['ICNL']])
lrcWell = gr.lrc_from_iglob(gr.NOD[IDOMAIN == pr['IWEL']])
lrcCanR = gr.lrc_from_iglob(gr.NOD[IDOMAIN == pr['ICNR']])

stress_period_data = {0: [],
                      1: (
                        [(lrc, pr['hCNL'], pr['cFresh']) for lrc in lrcCanL] +
                        [(lrc, pr['hWel'], pr['cFresh']) for lrc in lrcWell] +
                        [(lrc, pr['hCNR'], pr['cFresh']) for lrc in lrcCanR]),
                      2: (
                        [(lrc, pr['hCNL'], pr['cFresh']) for lrc in lrcCanL] +
                        [(lrc, pr['hWel'], pr['cFresh']) for lrc in lrcWell]),
                      3: (
                        [(lrc, pr['hWel'], pr['cFresh']) for lrc in lrcWell] +
                        [(lrc, pr['hCNR'], pr['cFresh']) for lrc in lrcCanR]),
                      4:[],
}
  
if False:
   stress_period_data += milkOn
   
Gwfchd ={'auxiliary': 'relconc',
         'stress_period_data': stress_period_data}

# %% === Gwfdrn === drains (can also be used for seepage face).


# %% === Gwfwel === wel (milk injection point)

# Chd is used instead
# stress_period_data = [(lrc, pr['hWel'], pr['cFresh']) for lrc in lrcWell]

# %% === Gwfrch === drains (can also be used for seepage face).

# %% === Gwfbuy (boyancy) ====

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

Gwtdsp ={**pr['disp']}

# %% === Gwtic ====== initial concentration

cstrt = gr.const(pr['cFresh'])
cstrt[IDOMAIN == pr['IDMK']] = pr['cSalt']

Gwtic = {'strt': cstrt}

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
