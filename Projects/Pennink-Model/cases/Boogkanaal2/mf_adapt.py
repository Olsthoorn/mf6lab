"""Cross section through the Boogcanaal taken from the book Pennink(1915)'s.

See http://www.citg.tudelft.nl/live/pagina.jsp?id=68e12562-a4d2-489a-b82e-deca5dd32c42&lang=en

This cross section was concstructed by De Hogestegeer, engeneer of Pennink. De Hogesteger also caried out as series of sand-box model tests simulating a cross section through the Dunes south of Zandvoort perpendicular to the coast, so show upconing under permanent lateral flow towards the flanks of the dune area. These experiments are in the last chapter of  Penninks (1915) book.

Next ot this the theory of groundwater flow in coss sections betwen recharge and extraction canals was discussed and visualized by manually contructing the head and stream lines in a drawing, which was not physically tested. One drawing shows the flow in a cross section between two half canals at the sides of the model to one in the center. The side canals are half due to the impermeable boundaries that were assumed there to construct teh drawding. The drawing only served at theoretical excercise, constructing the flow, as it was believed to be at the time after all the other physical experiments by Pennink but we can nevertheless model it.

@TO 090312 (in mflab using Seawat)
@TO 240112 (in mf6lab using USGS's `flopy` and `Modflow 6`)
@TO 240201
"""
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

# Generate the grid with refinement around the canals
zMax, zMin = pr['hL'], pr['frame'][:, 1].min()
xMin, xMax = pr['frame'][:, 0].min(), pr['frame'][:, 0].max()
xC = 0.5 * (xMin + xMax)

x = np.hstack((np.arange(xMin, xMin + 1, pr['dx1']),
               np.arange(xMin + 1, xC - 1, pr['dx']),
               np.arange(xC - 1, xC + 1, pr['dx1']),               
               np.arange(xC + 1, xMax - 1, pr['dx']),
               np.arange(xMax - 1, xMax + pr['dx1'] / 2, pr['dx1'])               
               ))
z = np.hstack((np.arange(zMax, zMax - 2.0, -pr['dz1']),
               np.arange(zMax - 2, zMin - pr['dz1'], -pr['dz'])
            ))

y = [-0.5, 0.5]

gr = Grid(x, y, z, min_dz=pr['min_dz'])

# IDOMAIN
# Start all cells inactive
IDOMAIN = gr.const(-1, dtype=int) # Inactive cells are vertially pass-through cells

# Only where sand --> active
IDOMAIN[:, 0, :][gr.inpoly(pr[ 'sand'], row=0)] =  1
IDOMAIN[:, 0, :][gr.inpoly(pr['canalL'], row=0)] = pr['IDL'] # Cells in Canal Left
IDOMAIN[:, 0, :][gr.inpoly(pr['canalC'], row=0)] = pr['IDC'] # Cells in Canal Left
IDOMAIN[:, 0, :][gr.inpoly(pr['canalR'], row=0)] = pr['IDR'] # Cells in Canal Left
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

# Adapt cells inside the canals to actually being in the water, k = infinit.
k_water = 1000.
k[  IDOMAIN == pr['IDL']] = k_water
k[  IDOMAIN == pr['IDC']] = k_water
k[  IDOMAIN == pr['IDR']] = k_water
k33[IDOMAIN == pr['IDL']] = k_water
k33[IDOMAIN == pr['IDC']] = k_water
k33[IDOMAIN == pr['IDR']] = k_water

Gwfnpf = { 'k':   k,
          'k33': k33,          
          'icelltype': pr['icelltype'],
}

# === Gwfic ======== Initial condictions (head)

hstrt = pr['hL']
Gwfic = {'strt': hstrt}

# === Gwfchd ======= fixed head
lrcL = gr.lrc_from_iglob(gr.NOD[IDOMAIN == pr['IDL']])
lrcC = gr.lrc_from_iglob(gr.NOD[IDOMAIN == pr['IDC']])
lrcR = gr.lrc_from_iglob(gr.NOD[IDOMAIN == pr['IDR']])

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