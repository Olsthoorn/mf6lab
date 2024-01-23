'''Pennink - density flow in Pennnink's (1915) sand box model

FOURTH SERIES OF EXPERIMENTS (MARCH 1905)
See http://www.citg.tudelft.nl/live/pagina.jsp?id=68e12562-a4d2-489a-b82e-deca5dd32c42&lang=en

In this experiment, Pennink (1915) studies simultaneous fresh and saltwater
flow from both a canal at the right and one at the left corner of hise sand-box model
to an axtracting canal in the center. The following text is cited from his book, p63:

"B. Phenomena of movement of the liquids of different specific gravity,
viz water and milk in case the discharge canal is placed in the centre
and the water symmetrically supplied on both sides."

Height of the apparatus   = 0.96 m
Breath of the apparatus   = 0.97 m
Thickness of the sandmass = 0.018 m

The milk (i.e. the "saltwater"), with a density ca. 1.03 kg/L,is applied
from a reservoir that can be verticallly adjusted, is applied on the right
side near the bottom of of the model.

The initial level of the milk is 25 cm above the bottom of the model.

The reservoir kept at a suitable elevation to create the initial
density distribution, i.e. horizontal at 25 cm above the bottom of the model).

Fresh water is supplied in the two upper corners. The surface of the sand
has a V-shape with the lowest point in the center. This shape can is used
to fix a head gradient at the top of the model:
Initially the water level is above the sand surface, and, therefore,
horizontal, with no flow in the model.

Next, the water flow is initiated and maintained by adding water to the
left and right edges on top of the sand surface and extraction it at the sand
surface in the center of the model. This creates a head gradient equal to
the inclination of the sand surface, 1:16 according to Pennink on P64.

Therefore, the flow is driven by the head, which makes the given water
flow added to the model, 28 L/h, useless information.

Because there is no flexibility in this model, there is only one photo in
Penninks book, showing the situation after some time, which can be
regarded as the final equilibrium.

@TO 20090312 20100523 20100719
@TO 20220115
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

# IDOMAIN
# Start all cells inactive
IDOMAIN = gr.const(-1, dtype=int) # Inactive cells are vertially pass-through cells

# Only where sand --> active
IDOMAIN[:, 0, :][gr.inpoly(pr[  'sand'], row=0)] = pr['IDSD']

IDOMAIN[gr.ZM < pr['zIface']] = pr['IDMK']

lrcMilk  = gr.lrc_from_xyz(pr['milkInjPnt'])['ic'] # global coords milk injection point

IDOMAIN.ravel()[gr.Iglob_from_lrc(lrcMilk)] = pr['iMlkInjPnt'] # mark milk injection cells

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

hstrt = gr.const(pr['hStrt'])

Gwfic = {'strt': hstrt}

# %% === Gwfchd ======= fixed head

Ix = gr.NOD[0, 0]
IzTopActive = gr.top_active_cells(IDOMAIN)[gr.NOD[0, 0]]

zTop = np.interp(gr.xm, gr.x[[0, int((gr.nx - 1)/2), -1]], pr['sand'][:, 1][2:5])
lrcTopActive = np.vstack((IzTopActive, np.zeros(gr.nx, dtype=int), gr.NOD[0, 0])).T

hMilkPnt = pr['zIface'] + (pr['hStrt'] - pr['zIface']) * pr['rhoFresh'] / pr['rhoSalt']

stress_period_data = {0: [(lrcMilk[0], hMilkPnt, pr['cSalt'])],
                      1: [(lrc, zt, pr['cFresh']) for lrc, zt in zip(lrcTopActive, zTop)] +\
                         [(lrcMilk[0], hMilkPnt, pr['cSalt'])],
                      2: [(lrcMilk[0], hMilkPnt, pr['cSalt'])]}

# stress_period_data = {0: [(lrcMilk[0], hMilkPnt, pr['cSalt'])],
#                       1: [(lrcMilk[0], hMilkPnt, pr['cSalt'])],
#                       2: [(lrcMilk[0], hMilkPnt, pr['cSalt'])]}

stress_period_data = {0: [],
                      1: [(lrc, zt, pr['cFresh']) for lrc, zt in zip(lrcTopActive, zTop)],
                      2: []}
                    
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
