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

sim_name = 'Series3'
section_name = 'Pennink (1915) {}.'.format(sim_name)
dirs = dirs.add_case(sim_name)
os.chdir(dirs.case)

params_wbk = os.path.join(dirs.case, sim_name + '.xlsx')
assert os.path.isfile(params_wbk), "Params_wbk not found: {}".format(params_wbk)

# lay = pd.read_excel(params_wbk, sheet_name='LAY', header=0, index_col=0)

props = {
    'start_date_time': '1904-10-05T09:35:00',  # Start time of Pennink test series 3
    'oc_frequency': ['FREQUENCY', 1], # Saving frequency for budget heads and concentrations
    'Qdim': 'cm3/min',
    'L': 65, # [cm]
    'H': 65, # [cm]
    'D': 1.8, # [cm]
    'dx': 1.0, #[cm]
    'dz': 1.0, #[cm]
    'Wrch': 45.0, # Width over which recharge is supplied
    'xCrch': 32.5, # cm x centre of recharge application on top of model
    'dCapZone': 0.04,  # [cm] thickness of the full capillary zone (see Photos)
    'icelltype': 0,
    'k_mpd': 650., # [m/d] calbrated from data in Pennink's series 1 experiments 
    'k': 86500/(24*60), # [cm/min] 
    'sy': 0.2,
    'ss': 1e-4,
    'disp' : {'alh': 0.2, 'ath1': 0.02,
              'ath2': 0.02, 'atv': 0.02,
              'diffc': 3.96e-4}, # Diffc in cm2/min
    'hCanL': 45.0, # [cm] head in left canal
    'hCanR': 46.0, # [cm] head in right canal
    'Qrch': 4200 / 60, # [cm3/min]
    'Qmilk': 2.14, # [cm3/min]
    'por': 0.38, # [-] porosity
    'zIface': 13, # [cm] elavation of initial interface
    'cCanL': 0.0, # concentration in canal left side
    'cCanR': 0.0, # concentration in canal right side
    'rhoFresh': 1., # [g/cm3]     not used in this model
    'rhoSalt': 1.0245, # [g/cm3]  not used in this model
    'cFresh': 0.0, # [g/cm3]      not used in this model
    'cSalt': 35.0,  # [g/cm3]    not used in this model
    'cNoInk' : 0.0,
    'cInk': 1.0, # ink concentration (relative)
    'IDSD': 2, # IDOMAIN value for initial area with sand.
    'IDMK': 3, # IDOMAIN value for intial area with milk. 
    'IDCL': 4, # IDOMAIN value for cells in Canal Left (For easy finding of cells.)
    'IDCR': 5, # IDOMAIN value for cells in Canal Right (For easy finding of cells.)
    'iInk': 6, # IDOMAIN value for ink injection points (For easy finding of cells.)
    'iMlkInjPnt': 7, # IDOMAIN value for milk injection point.
    'xyzInk': [[50.6, 0, 50.0],
               [52.6, 0, 44.9],
               [53.8, 0, 33.7],
               [56.1, 0, 23.0]],
    'milkInjPnt': [[0.1, 0.0, 63.5]],
    'sand' :np.array([[-2.2878, -2.3540], # Contour sand mass cm vs LL of model
                      [67.6080, -2.5445],
                      [67.7985, 60.3133],
                      [43.8015, 60.1229],
                      [30.0890, 59.1705],
                      [23.2327, 57.8371],
                      [17.7096, 55.9323],
                      [11.4247, 52.6942],
                      [10.4725, 52.6942],
                      [ 8.5679, 50.4085],
                      [-2.2878, 49.2656]]),
    'canalL': np.array([[-2.2878, 39.7417], # Contour of canal on the left side
                         [3.2353, 38.7893],
                         [5.7112, 41.0750],
                         [7.4252, 45.2656],
                         [7.8061, 50.7894],
                         [10.4725, 65.8372],
                         [-2.2878, 65.6467]]),
    'canalR': np.array([[68.3698, 63.3610], # Contour canal on right side
                        [66.2748, 66.2182],
                        [61.1326, 65.2658],
                        [61.3231, 51.7418],
                        [61.8944, 46.9799],
                        [63.0372, 42.4084],
                        [65.5130, 41.2655],
                        [69.1316, 41.4560],
                        [68.7507, 64.3134]])
}
