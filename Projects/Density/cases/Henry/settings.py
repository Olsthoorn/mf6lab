import os
import pandas as pd

from src.mf6tools import  Dirs

# Project name (above cases)
HOME = '/Users/Theo/GRWMODELS/python/mf6lab/Projects/Density/'
assert os.path.isdir(HOME), "Can't find the directory {}".format(HOME)

LENGTH_UNITS = 'meters'
TIME_UNITS = 'days'

dirs = Dirs(HOME)

sim_name = 'Henry'
section_name = 'MF6 Dev. Team (2023), ex 51: {}.'.format(sim_name)
dirs = dirs.add_case(sim_name)
os.chdir(dirs.case)

params_wbk = os.path.join(dirs.case, sim_name + '.xlsx')
assert os.path.isfile(params_wbk), "Params_wbk not found: {}".format(params_wbk)

lay = pd.read_excel(params_wbk, sheet_name='LAY', header=0, index_col=0)

props = {
    'start_date_time': '2024-01-09',
    'NPER':  1,
    'NSTP': 500,
    'oc_frequency': 5, # Saving frequency
    'tmax': 0.5, # [d]
    'L': 2.0, # [m]
    'H': 1.0, # [m]
    'dx': 0.025, #[m]
    'dz': 0.025, #[m]
    'Lsys': 2.0, # [m]
    'k': 864., # [m/d]
    'rhoFresh': 1000., # [kg/m3]
    'rhoSalt': 1024.5,
    'cFresh': 0.0, # [kg/m3]
    'cSalt': 35.,
    'hStrt': 1.0, # [m]
    'por': 0.35,
    'icelltype': 0.,
    'sy': 0.2,
    'ss': 1e-4,
    'disp' : {'alh': None, 'ath1': None, 'ath2': None, 'atv': None, 'diffc': 0.57024},
    'Qin_a': 5.7024, # [m3/d] = (m2/d here)
    'Qin_b': 2.8512, # [m3/d] = (m2/d here)
}
