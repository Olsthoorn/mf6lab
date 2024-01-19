import os
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

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
    'photo': 'Series3_01_p50.jpg', # background phohto
    'extent': (-20.5, 84.0, -9, 83), # extent of photo
    'oc_frequency': 2, # Saving frequency for budget heads and concentrations
    'Qdim': 'cm3/min',
    'L': 65,           # [cm]
    'H': 65,           # [cm]
    'D': 1.8,          # [cm]
    'dx': 1.0,         # [cm]
    'dz': 1.0,         # [cm]
    'Wrch': 45.0,      # Width over which recharge is supplied
    'xCrch': 32.5,     # cm x centre of recharge application on top of model
    'dCapZone': 0.04,  # [cm] thickness of the full capillary zone (see Photos)
    'icelltype': 0,    # [-] 0 = not convertable (confined) 1 is convertable
    'k_mpd': 650.,     # [m/d] calbrated from data in Pennink's series 1 experiments 
    'k': 86500/(24*60), # [cm/min] 
    'sy': 0.2,         # [ -  ]
    'ss': 1e-4,        # [1/cm]
    'disp' : {'alh':  0.2, 'ath1': 0.02, 'ath2': 0.02, 'atv': 0.02, 'diffc': 3.96e-4}, # Diffc in cm2/min
    'hCanL': 45.0,     # [cm] head in left canal
    'hCanR': 46.0,     # [cm] head in right canal
    'Qrch': 4200 / 60, # [cm3/min]
    'Qmilk': 2.14,     # [cm3/min]
    'por': 0.38,       # [-] porosity
    'zIface': 13,      # [cm] elavation of initial interface
    'cCanL': 0.0,      # concentration in canal left side
    'cCanR': 0.0,      # concentration in canal right side
    'rhoFresh': 1.,    # [g/cm3]
    'rhoSalt': 1.0245, # [g/cm3]
    'cFresh': 0.0,     # [g/cm3]
    'cSalt': 35.0,     # [g/cm3] 
    'cNoInk' : 0.0,    # [g/cm3]
    'cInk': 35.0,       # ink concentration (relative)
    'IDSD': 2,         # IDOMAIN value for initial area with sand.
    'IDMK': 3,         # IDOMAIN value for intial area with milk. 
    'IDCL': 4,         # IDOMAIN value for cells in Canal Left (For easy finding of cells.)
    'IDCR': 5,         # IDOMAIN value for cells in Canal Right (For easy finding of cells.)
    'iInk': 6,         # IDOMAIN value for ink injection points (For easy finding of cells.)
    'iMlkInjPnt': 7, # IDOMAIN value for milk injection point.
    'xyzInk': np.array([[50.6, 0, 50.0],
                        [52.6, 0, 44.9],
                        [53.8, 0, 33.7],
                        [56.1, 0, 23.0]]),          # xyz-Coords of ink injection points (not used here)
    'milkInjPnt': np.array([[63.5, 0.0, 0.1]]),     # xyz-Coords of milk injection point
    'sand' :np.array([[ 0.0,  0.0], # Contour sand mass cm vs LL of model
                      [65.0,  0.0],
                      [65.0, 44.7],
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
                      [1.80, 40.7],
                      [ 0.0, 42.0],
                      [ 0.0,  0.0]]),
    'canalL': np.array([[0.0, 65.00], # Contour of canal on the left side
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


if __name__ == '__main__':
    
    pr = props
    
    print('Props')
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_title("Matching the size of the model to the photo.")
    ax.set_xlabel('x [cm]')
    ax.set_ylabel('z [cm]')
    
    foto = Image.open(os.path.join(dirs.photos, pr['photo']))
    ax.imshow(foto, extent=pr['extent'])

    ax.plot([0, pr['L'], pr['L'], 0,  0], [0, 0, pr['H'], pr['H'], 0], 'b', label='bbox around model')
    ax.plot(*pr['sand'].T, 'brown', label='sand')
    ax.plot(*pr['canalL'].T, 'black', label='canalL')
    ax.plot(*pr['canalR'].T, 'black', label='canalR')
    ax.plot(pr['xyzInk'][:, 0], pr['xyzInk'][:, 2], 'ro', label='Inkk injection points')
    plt.show()
