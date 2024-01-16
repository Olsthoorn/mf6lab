'''Series 2 of Pennink(1915)'s 1903 - 1905 sand box model experiments.
Experiments series 3 (Septemner 1904)

See http://www.citg.tudelft.nl/live/pagina.jsp?id=68e12562-a4d2-489a-b82e-deca5dd32c42&lang=en

In this experiment, Pennink(1915) studies salt water upconing due to extraction
by a canal under uniform recharge.

To easier match the actual times of the tests, we will use time in minutes
instead of days.

In this experiment, Pennink studies simultaneous cFresh and saltwater
flow in a cross section fed by rain and discharging to a canal at the
westsided of the section.

He creates an initial saltwater zone by letting milk enter slowly at the
bottom of the model. The flow is regulated by a fixed-head reservoir.
After the milk volume has rose to 13 cm height in the model, he started the
recharge by letting 4.2 L/h dripping on the model over a width of 45 cm, which caused to the interface to finally reach an equilibrium position.

During this phase where the milk intefrace is upconing, milk is injected
continuously at a rate of 2.14 ml/min. At this rate the model filled up to 13 cm in 5 hours. Keeping the same rate also works during the entire simulation.

Pennink mentions that the milk was added drop by drop. Thus no fixed saltwater head was applied, so we simulate this milk entry by using a WEL in MODFLOW6.

After a day, Pennink stopped adding milk. Then after another day, he
took last photo in this test series.

As it took about two hours to establish a virtual equilibrium of the
saltwater interface, I shortened the simulation before equilibrium to 2 hours. Then I continu for 24 hours without the milk injection to simulate the washing out of the
milk from the model.

This approach should allow for calibration of the dispersivity in Penninks original model. However dispersion is expected to work differently at the interface between water and fatty milk. Therefore, such calibration was not attempted. Nevertheless the situation in Modflow model after a full day without milk addition, was
similar enough to that on the last photo in Pennink's series 3 that the comparison is considered satisfatory.

@TO 090312 100523 100721 Using seawat through mflab
@TO 240115 (in mf6lab using USGS's `flopy` and `Modflow 6`)
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
IDOMAIN[:, 0, :][gr.inpoly(pr['canalL'], row=0)] = pr['IDCL'] # Cells in Canal Left
IDOMAIN[:, 0, :][gr.inpoly(pr['canalR'], row=0)] = pr['IDCR'] # Cells in Canal Right

# Inactive cells are within canals but above their stage level
IDOMAIN[np.logical_and(IDOMAIN == pr['IDCL'], gr.ZM > pr['hCanL'])] = 0
IDOMAIN[np.logical_and(IDOMAIN == pr['IDCR'], gr.ZM > pr['hCanR'])] = 0

IDOMAIN[gr.ZM < pr['zIface']] = pr['IDMK']

lrcInk   = gr.lrc_recarray(*np.array(pr['xyzInk']).T) # global coords of injection points
lrcMilk  = gr.lrc_recarray(*np.array(pr['milkInjPnt']).T) # global coords milk injection point

IDOMAIN.ravel()[gr.Iglob(lrcInk)]  = pr['iInk'] # mark ink injection cells
IDOMAIN.ravel()[gr.Iglob(lrcMilk)] = pr['iMlkInjPnt'] # mark milk injection cells

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

hwt = pr['hCanL'] + gr.XM / pr['L'] * (pr['hCanR'] - pr['hCanL']) # Reasonalble hstrt
zTopCapZone = hwt[0, 0] + pr['dCapZone']

k[gr.ZM > zTopCapZone[np.newaxis, np.newaxis]] /= 5 # Unsat zone cond. above full capillary zone

Gwfnpf = { 'k':   k,            
            'icelltype': pr['icelltype'],
}

# %% === Gwfic ======== Initial condictions (head)

hstrt = hwt # Reasonalble hstrt
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

# %% === Gwfwel === wel (milk injection point)

milkOn  = [(tuple(lrc_), pr['Qmilk'], pr['cSalt']) for lrc_ in lrcMilk['ic']]
milkOff = [(tuple(lrc_),  0.0,        pr['cSalt']) for lrc_ in lrcMilk['ic']]

stress_period_data = dict()
for isp, im in enumerate(perDF['Milk']):
      if im == 0:
            stress_period_data[isp] = milkOff
      else:
            stress_period_data[isp] = milkOn
            
Gwfwel ={'auxiliary': 'relconc',
         'stress_period_data': stress_period_data}

# %% === Gwfrch === drains (can also be used for seepage face).
IxRch = gr.NOD[0, 0][np.logical_and(gr.xm >= pr['xCrch'] - pr['Wrch'] / 2,
                                    gr.xm <= pr['xCrch'] + pr['Wrch'] / 2)] #
IzRch = gr.top_active_cells(IDOMAIN)[IxRch]

rch = pr['Qrch'] / pr['Wrch'] / pr['D'] # m3/min / cm / cm = cm / min

RCH = [((iz, 0, ix), rch) for ix, iz in zip(IxRch, IzRch)]

Gwfrch = {'stress_period_data': {0: RCH}}

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
por[gr.ZM > zTopCapZone[np.newaxis, np.newaxis]] /= 3 # Estimate for the unsaturated zone

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

# %% === Gwtcnc ====== constant concentraction

# Ink injection points (constant concentration at injection locations)

# Concentration cells [(l, r, c) conc] of ink injection points.
concOn  = [(tuple(lrc_), pr['cInk'])  for lrc_ in lrcInk['ic']]
concOff = [(tuple(lrc_),  0.0)        for lrc_ in lrcInk['ic']]

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
                      ("BUDGET",        "FREQUENCY", pr['oc_frequency'])],
}
print('Done mf_adapt')

if __name__ == '__main__':
   
   print(dirs)     
