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
from fdm.mfgrid import Grid

dirs = settings.dirs
sim_name = settings.sim_name
section_name = settings.section_name

# Parameters workbook
params_wbk = settings.params_wbk

# Parameters of the model and simulation
pr = settings.props

## Get section data

# === tdis ======  time discretization
perDF = mf6tools.get_periodata_from_excel(params_wbk, sheet_name='PER')
period_data = [tuple(sp) for sp in perDF[['PERLEN', 'NSTP', 'TSMULT']].values]


start_date_time = pr['start_date_time'] # Must be string

Simtdis = {'perioddata': period_data,
           'nper': len(period_data),
           'start_date_time': start_date_time,
           'time_units': settings.TIME_UNITS,
           }

# === Gwfdis ====== spatial discritization of flow model

xMin, xMax = np.unique(pr['sand'][:, 0])[[0, -1]]
zMin, zMax = np.unique(pr['sand'][:, 1])[[0, -1]]
zMax = pr['hL'] + pr['dz']
x = np.arange(xMin, xMax + pr['dx'] / 2, +pr['dx'])
z = np.arange(np.ceil(zMax), np.floor(zMin) - pr['dz'] / 2, -pr['dz'])
y = [-0.5, 0.5]

gr = Grid(x, y, z, min_dz=pr['min_dz'])

# IDOMAIN
# Start all cells inactive
IDOMAIN = gr.const(-1, dtype=int) # Inactive cells are vertially pass-through cells

# Only where sand --> active
IDOMAIN[:, 0, :][gr.inpoly(pr[ 'sand'], row=0)] =  1
IDOMAIN[:, 0, :][gr.inpoly(pr['canal'], row=0)] = pr['IDC'] # Cells in Canal Left

Gwfdis = {'gr': gr,
          'idomain':      IDOMAIN,
          'length_units': settings.LENGTH_UNITS,
}

# === Gwfsto ======= Storage coefficients

Gwfsto = {'sy': gr.const(pr['sy']),
          'ss': gr.const(pr['ss']),
          'iconvert': 1,
          }

# === Gwfnpf ======= Horizontal and vertical conductivity
k   = gr.const(pr['k'])
k33 = gr.const(pr['k33'])

Gwfnpf = { 'k':   k,
          'k33': k33,          
          'icelltype': pr['icelltype'],
}

# === Gwfic ======== Initial condictions (head)

hstrt = pr['hL']
Gwfic = {'strt': hstrt}

# === Gwfchd ======= fixed head
active = IDOMAIN > 0
Izwt = np.arange(gr.nz)[gr.Z[1:, 0, 0] < pr['hL']]
lrcL = gr.lrc_from_iglob(gr.NOD[Izwt, 0,  0][active[Izwt, 0,  0]])
lrcR = gr.lrc_from_iglob(gr.NOD[Izwt, 0, -1][active[Izwt, 0, -1]])
lrcC = gr.lrc_from_iglob(gr.NOD[IDOMAIN == pr['IDC']])
stress_period_data = [(lrc, pr['hL']) for lrc in lrcL] +\
                     [(lrc, pr['hR']) for lrc in lrcR] +\
                     [(lrc, pr['hC']) for lrc in lrcC]
                     
Gwfchd ={ #'auxiliary': 'relconc',
         'stress_period_data': stress_period_data}

# === Gwfdrn === drains (can also be used for seepage face).

# Here just to get the top active cells
hDr = gr.Z[0, 0] - pr['drain_depth']
drn_xyz = np.vstack((gr.xm, np.zeros(gr.nx), hDr)).T
Iz = gr.lrc_from_xyz(drn_xyz)['ic'][:, 0]
Iz = gr.top_active_cells(IDOMAIN, Iz)

# Do nothing with drains for now

# === Gwfoc ==== Output control for flow model
Gwfoc = {'head_filerecord':   os.path.join(dirs.SIM, "{}Gwf.hds".format(sim_name)),
         'budget_filerecord': os.path.join(dirs.SIM, "{}Gwf.cbc".format(sim_name)),
         'saverecord': [("HEAD", "FREQUENCY", pr['oc_frequency']), ("BUDGET", "FREQUENCY", pr['oc_frequency'])],
}

# === T R A N S P O R T =====

# === Gwtdis ======= discritization for transport model

# Uses gr above, don't specify


# === Gwtfmi ===== Flow model interface

# simultaneous exchange is automatic (see mf_setup.py)

# === Gwtmst ===== Mobile storage and transfer

# === Gwtadv ====== advection

# === Gwtdsp ====== dispersion & diffusion

# === Gwtic ====== initial concentration

# === Gwtssm ====== Source-Sink mixing

# === Gwtcnc ====== constant concentraction

# === Gwtoc ====== output control

Gwtoc = {
      'concentration_filerecord' : os.path.join(dirs.SIM, '{}Gwt.ucn'.format(sim_name)),
      'budget_filerecord':         os.path.join(dirs.SIM, '{}Gwt.cbc'.format(sim_name)),
      'saverecord' : [("CONCENTRATION", "FREQUENCY", pr['oc_frequency']),
                      ("BUDGET", "FREQUENCY", pr['oc_frequency'])],
}
print('Done mf_adapt')

if __name__ == '__main__':
   
   print(dirs)     