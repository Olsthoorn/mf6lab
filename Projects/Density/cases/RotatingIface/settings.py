import os
import pandas as pd

from src.mf6tools import  Dirs

# Project name (above cases)
HOME = '/Users/Theo/GRWMODELS/python/mf6lab/Projects/Density/'
assert os.path.isdir(HOME), "Can't find the directory {}".format(HOME)

LENGTH_UNITS = 'meters'
TIME_UNITS = 'days'

dirs = Dirs(HOME)

section_name = 'Rotating Interface Problem (Modflow Development Team (2023), exmaple 53)'
sim_name = 'RotatingIface'
dirs = dirs.add_case(sim_name)
os.chdir(dirs.case)

params_wbk = os.path.join(dirs.case, sim_name + '.xlsx')
assert os.path.isfile(params_wbk), "Params_wbk not found: {}".format(params_wbk)

lay = pd.read_excel(params_wbk, sheet_name='LAY', header=0, index_col=0)

props = {
    'NPER':  1,
    'NSTP': 500,
    'oc_frequency': ['FREQUENCY', 100], # Saving frequency for budget heads and concentrations
    'tmax': 0.5, # [d]
    'L': 300, # [m]
    'H': 40, # [m]
    'dx': 1.0, #[m]
    'dz': 0.5, #[m]
    'Lsys': 150, # [m]
    'k': 2.0, # [m/d]
    'rhoref': 1000., # [kg/m3]
    'crhoref': 0.,
    'drhodc': 0.7,
    'rhoR': 1000.,
    'rhoM': 1012.5,
    'rhoL': 1025.,
    'cR': 0.0, # [kg/m3]
    'cM': 17.5,
    'cL': 35.,
    'IeRM': 40.0, # [m] Interface extent between zone 1 and zone 2
    'IeML': 40.0, # [m] Interface extent between zone 2 and zone 3
    'xmRM': 170.0, # [m] X-midpoing for zone 1 and 2 interface
    'xmML': 130.0, # [m] X-midpoing for zone 2 and zone 3
    'hStrt': 0.0, # [m]
    'por': 0.2,
    'disp' : {'alh': .1, 'ath1': .01, 'ath2': .0, 'atv': .01, 'diffc': 1.89e-5 * 86400}, # [m] and [m2/d]
}
