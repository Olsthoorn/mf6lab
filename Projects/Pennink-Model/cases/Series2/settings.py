import os
import pandas as pd
import numpy as np

from src.mf6tools import  Dirs

# Project name (above cases)
HOME = '/Users/Theo/GRWMODELS/python/mf6lab/Projects/Pennink-Model/'
assert os.path.isdir(HOME), "Can't find the directory {}".format(HOME)

LENGTH_UNITS = 'cm'
TIME_UNITS = 'minutes'

dirs = Dirs(HOME)

sim_name = 'Series2'
section_name = 'Pennink (1915). Case {}'.format(sim_name)
dirs = dirs.add_case(sim_name)
os.chdir(dirs.case)

params_wbk = os.path.join(dirs.case, sim_name + '.xlsx')
assert os.path.isfile(params_wbk), "Params_wbk not found: {}".format(params_wbk)

# lay = pd.read_excel(params_wbk, sheet_name='LAY', header=0, index_col=0)

props = {
    'oc_frequency': ['FREQUENCY', 25], # Saving frequency for budget heads and concentrations
    'L': 65, # [cm]
    'H': 65, # [cm]
    'D': 1.8, # [cm]
    'dx': 1.0, #[cm]
    'dz': 1.0, #[cm]
    'zCapZone': 51.0,  # [cm] Top of full capillary zone (see descripition)
    'icelltype': 0,
    'k_mpd': 650., # [m/d] calbrated from data in Pennink's series 1 experiments 
    'k': 65000 / (24 * 60), # [cm/min] 
    'sy': 0.2,
    'ss': 1e-4,
    'disp' : {'alh': 0.1, 'ath1': 0.01, 'ath2': 0.01, 'atv': 0.01, 'diffc': 3.96e-4}, # Diffc in cm2/min
    'hCanL': 45.2, # [cm] head in left canal
    'hCanR': 45.7, # [cm] head in right canal
    'por': 0.35, # [-] porosity
    'cCanL': 0.0, # concentration in canal left side
    'cCanR': 0.0, # concentration in canal right side
    'rhoFresh': 1., # [g/cm3]     not used in this model
    'rhoSalt': 1.0245, # [g/cm3]  not used in this model
    'cFresh': 0.0, # [g/cm3]      not used in this model
    'cSalt': 0.035,  # [g/cm3]    not used in this model
    'cInk': 1.0, # ink concentration (relative)
    'IDCL': 2, # IDOMAIN value for cells in Canal Left (For easy finding of cells.)
    'IDCR': 3, # IDOMAIN value for cells in Canal Right (For easy finding of cells.)
    'iInk': 4, # IDOMAIN value for ink injection points (For easy finding of cells.)
    'xyzInk': [[50.6, 0, 50.0],
               [52.6, 0, 44.9],
               [53.8, 0, 33.7],
               [56.1, 0, 23.0]],
    'sand' :np.array([[ 0.0,  0.0], # Contour sand mass cm vs LL of model
                      [66.0,  0.0],
                      [65.0, 44.3],
                      [62.6, 44.7],
                      [62. , 47.5],
                      [61.5, 53. ],
                      [61. , 62.9],
                      [55. , 62.1],
                      [48.8, 62.2],
                      [39.8, 61.7],
                      [34.5, 62. ],
                      [26.5, 60.8],
                      [19.2, 58.9],
                      [13.4, 56.7],
                      [ 8.5, 53. ],
                      [ 7.1, 51.9],
                      [ 6.5, 47. ],
                      [ 6. , 44. ],
                      [ 5.2, 42. ],
                      [ 3.9, 40.7],
                      [ 2.2, 40.6],
                      [ 0.5, 40.8],
                      [ 0. ,  0.0]]),
    'canalL': np.array([[ 0.0, 65.00], # Contour of canal on the left side
                       [ 0.0, 41.99],
                       [1.80, 40.72],
                       [3.95, 40.83],
                       [5.61, 42.71],
                       [6.68, 47.54],
                       [7.17, 52.26],
                       [8.33, 57.87],
                       [9.01, 65.00],
                       [ 0.0, 65.00]]),
    'canalR': np.array([[60.79, 65.00], # Contour canal on right side
                       [61.39, 57.70],
                       [61.59, 52.29],
                       [61.99, 48.46],
                       [62.29, 45.32],
                       [62.98, 44.53],
                       [65.00, 44.74],
                       [65.00, 65.00],
                       [60.89, 65.00]]),
}
