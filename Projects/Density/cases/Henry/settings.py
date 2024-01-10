import os
import pandas as pd

from src.mf6tools import  Dirs

# Project name (above cases)
HOME = '/Users/Theo/GRWMODELS/python/mf6lab/Projects/Density/'
assert os.path.isdir(HOME), "Can't find the directory {}".format(HOME)

LENGTH_UNITS = 'meters'
TIME_UNITS = 'seconds'

dirs = Dirs(HOME)

section_name = 'Henry (1964) Salt Water Intruston Problem'
sim_name = 'Henry'
dirs = dirs.add_case(sim_name)
os.chdir(dirs.case)

params_wbk = os.path.join(dirs.case, sim_name + '.xlsx')
assert os.path.isfile(params_wbk), "Params_wbk not found: {}".format(params_wbk)

lay = pd.read_excel(params_wbk, sheet_name='LAY', header=0, index_col=0)

props = { # From pdf file see references. GEO-SLOPE Intl Ltd Calgary, Alberta, Canada www.geo-slope.com (no official reference).
    'disp' : {'alh': 1e-20, 'ath1': 1e-20, 'ath2': 1e-20, 'atv': 1e-20, 'diffc': 1.89e-5}, # [m] and [m2/s]
    'k'    : 1e-2, # m/s
    'por'  : 0.35, # [-] porosity
    'qL'    : 6.6e-6, # m/s (unit flux) from left to right
    'qL1'   : 3.3e-6, # m/s (unit flux adapted) from left to right
    'hR'   : 1.0, # Hydrostatic pressure distribution along the right side
    'FRESH': 0.0, # Fresh water concentration
    'SALT' : 1.0, # Constant concentration along right boundary
    'rhomin': 1000, # kg/m3,
    'rhomax': 1025 # kg/m3,
}
