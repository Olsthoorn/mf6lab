'''Pennink - density flow in Pennnink's (1915) sand box model

FIFTH SERIES OF EXPERIMENTS (MARCH 1905)
See http://www.citg.tudelft.nl/live/pagina.jsp?id=68e12562-a4d2-489a-b82e-deca5dd32c42&lang=en

In this experiment, Pennink (1915) studies simultaneous fresh and saltwater
flow from a canal in the right and one in the left corner of his sand-box model
to a canal in the center.

The following text is cited from his Pennink (1915), p63:

"To accelerate the motion and to observe the phenomena more quickly and
more clearly, the slopes are now made steeper. For this purpose a canal
has been reserved in the middle of the sand mass, in which a lower level
may be kept up.

The discharge from that canal is brought about by means of a pipe
connected with the back part of the apparatus, the cock of which is
visible halfway on the left side."

Pennink simulates upconing below the central canal. Interestingly he also
studied the effect of lateral groundwater flow has on this upconing. Lateral groundwater flow gives an interesting shape of the cone, demonstrating that it is possible for the groundwater below the canal to remain fresh, while the salt water enters from the side.

@TO 090320 100508 100719
@TO 20220116
'''

import os
import numpy as np
from src import mf6tools
from fdm.mfgrid import Grid
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
lrcCanM = gr.lrc_from_iglob(gr.NOD[IDOMAIN == pr['ICNM']])
lrcCanR = gr.lrc_from_iglob(gr.NOD[IDOMAIN == pr['ICNR']])

stress_period_data = {0: [],
                      1: (
                        [(lrc, pr['hCNL'], pr['cFresh']) for lrc in lrcCanL] +
                        [(lrc, pr['hCNM'], pr['cFresh']) for lrc in lrcCanM] +
                        [(lrc, pr['hCNR'], pr['cFresh']) for lrc in lrcCanR]),
                      2: (
                        [(lrc, pr['hCNL'], pr['cFresh']) for lrc in lrcCanL] +
                        [(lrc, pr['hCNM'], pr['cFresh']) for lrc in lrcCanM]),
                      3: (
                        [(lrc, pr['hCNM'], pr['cFresh']) for lrc in lrcCanM] +
                        [(lrc, pr['hCNR'], pr['cFresh']) for lrc in lrcCanR]),
                      4:[],
}
  
if False:
   stress_period_data += milkOn
   
Gwfchd ={'auxiliary': 'relconc',
         'stress_period_data': stress_period_data}

# %% === Gwfdrn === drains (can also be used for seepage face).


# %% === Gwfwel === wel (milk injection point)

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
