import os
import pandas as pd

from src.mf6tools import  Dirs

HOME = '/Users/Theo/GRWMODELS/python/mf6lab/Projects/VoortoetsAGT/'
assert os.path.isdir(HOME), "Can't find the directory {}".format(HOME)

LENGTH_UNITS = 'meters'
TIME_UNITS = 'days'

dirs = Dirs(HOME)

section_name = 'Generic hand-drawn cross section met inundatielaag'
sim_name = 'GenericInun'

dirs = dirs.add_case(sim_name)
os.chdir(dirs.case)

params_wbk = os.path.join(dirs.case, sim_name + '.xlsx')
assert os.path.isfile(params_wbk), "Params_wbk not found: {}".format(params_wbk)

lay = pd.read_excel(params_wbk, sheet_name='LAY', header=0, index_col=0)

props = {
    'start_date_time': '2024-01-19',  # Start time of Pennink test series 3
    'oc_frequency': 10, # Saving frequency for budget heads and concentrations
    'por': 0.35,       # [-] porosity
    'hStrt': 0.0, # [m] initial head
    'zIface': -30,      # [cm] elavation of initial interface
    'rhoFresh': 1.,    # [g/cm3]
    'rhoSalt': 1.03, # [g/cm3] milk acc. to Pennink
    'cFresh': 0.0,     # [g/cm3]
    'cSalt':  1.0,     # [g/cm3] 
    'ahl': 1.0,  # [m] longitudal dispersivity
    'ath1': 0.1,  # [m] transversal dispersivity
    'ath2': 0.1,  # [m] transversal dispersivity
    'atv': 0.1,  # [m] vertical dispersivity
    'diffc': 6e-5, # [m2/d] diffusion coefficient
    'drainDepth': 0.5, # [m]
    'cDrainage': 100, # [d]
    'rch': 0.001, # [m/d]
}
