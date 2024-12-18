import os
import numpy as np

from src.mf6tools import  Dirs

# Project name (above cases)
HOME = '/Users/Theo/GRWMODELS/python/mf6lab/Projects/Pennink-Model/'
assert os.path.isdir(HOME), "Can't find the directory {}".format(HOME)

LENGTH_UNITS = 'cm'
TIME_UNITS = 'minutes'

dirs = Dirs(HOME)

sim_name = 'Series6'
section_name = 'Pennink (1915) {}.'.format(sim_name)
dirs = dirs.add_case(sim_name)
os.chdir(dirs.case)

params_wbk = os.path.join(dirs.case, sim_name + '.xlsx')
assert os.path.isfile(params_wbk), "Params_wbk not found: {}".format(params_wbk)

props = {
    'start_date_time': '1904-10-05T09:35:00',  # Start time of Pennink test series 3
    'photo': 'Series5_01_p70.jpg', # background photo
    'extent': (-28.7, 156.3, -18.9, 134.8), # Extent of the Photo background.
    'oc_frequency': 5, # Saving frequency for budget heads and concentrations
    'Qdim': 'cm3/min', # Flow dimension
    'L': 97,           # [cm]
    'H': 96,           # [cm]
    'D': 1.8,          # [cm]
    'dx': 1.0,         # [cm]
    'dz': 1.0,         # [cm]
    'zScreen': (75, 55), # [cm] top and bottom of screen
    'icelltype': 0,    # [-] 0 = not convertable (confined) 1 is convertable
    'k_mpd': 650.,     # [m/d] calbrated from data in Pennink's series 1 experiments 
    'k': 86500/(24*60), # [cm/min] 
    'sy': 0.2,         # [ -  ]
    'ss': 1e-4,        # [1/cm]
    'disp' : {'alh':  0.2, 'ath1': 0.02, 'ath2': 0.02, 'atv': 0.02, 'diffc': 3.96e-4}, # Diffc in cm2/min
    'por': 0.38,       # [-] porosity
    'zIface': 23,      # [cm] elavation of initial interface
    'rhoFresh': 1.,    # [g/cm3]
    'rhoSalt': 1.03, # [g/cm3] milk acc. to Pennink
    'cFresh': 0.0,     # [g/cm3]
    'cSalt':  1.0,     # [g/cm3]     
    'IDSD': 2,         # IDOMAIN value for initial area with sand.
    'IDMK': 3,         # IDOMAIN value for intial area with milk. 
    'ICNL': 4,         # IDOMAIN value for canal left. 
    'ICNM': 5,         # IDOMAIN value for canal in center. 
    'ICNR': 6,         # IDOMAIN value for canal right. 
    'iMlkInjPnt': 7,   # IDOMAIN value for milk injection point.
    'IUNS' : 8,        # IDOMAIN value for the unsaturated zone
    'IWEL': 9,         # IDOMAIN value for the well cells 
    'milkInjPnt': [(91.5, 0., 1.5)],     # xyz-Coords of milk injection point
    'hCNL': 91.0, # head in the center top of model (drain point)
    'hCNM': 87.5,   # [cm] head in right hand canal
    'hWel': 87.5,   # [wel] head
    'hCNR': 91.0, # head in the center top of model (drain point)
    'hMilk': 87.5 - 6.5, # [cm] saline had as mentioned by Pennink
    'Qout' : 60e3 / 60, # [cm/min, Pennink mentions 60 L/h
    'hStrt' : 90.0, # Uniform start head.
    'sand' : np.array([ # [cm] contour of sand body.               
                [ 0.0,  0.0],
                [97.0,  0.0],
                [97.0, 92.1],
                [48.5, 90.7],
                [ 0.0, 92.1],
                [ 0.0,  0.0]]),
    'unsa': np.array([
                [ 0.0, 88.0],
                [48.5, 85.0],
                [97.0, 88.0],
                [97.0, 89.0],
                [48.5, 86.0],
                [ 0.0, 89.0],
                [ 0.0, 88.0]]),
    'canalR': np.array([ # [cm] contour around right canal
                [ 97.00, 96.00],
                [ 91.80, 96.00],
                [ 93.45, 89.20],
                [ 94.59, 82.70],
                [ 97.00, 82.00],
                [ 97.00, 96.00]]),
    'canalL': np.array([ # [cm] contour around left canal
                [ 0.00, 96.00],
                [ 0.00, 82.00],
                [ 2.41, 82.70],
                [ 3.55, 89.20],
                [ 5.20, 96.00],
                [ 0.00, 96.00]]),
    'canalM': np.array([ # [cm] contour around center canal
                [43.30, 96.00],
                [44.50, 80.97],
                [46.50, 76.22],                
                [47.00, 75.48],
                [50.00, 75.48],                
                [50.50, 76.22],
                [52.50, 80.97],
                [53.70, 96.00],                
                [43.30, 96.00]]),
    'wt': np.array([
                [ 0.00, 87.30],
                [ 3.20, 87.30],
                [44.48, 77.40],
                [50.24, 77.40],
                [93.80, 87.30],
                [97.00, 87.30]]),
}


if __name__ == '__main__':
    
    import matplotlib.pyplot as plt
    from PIL import Image
    
    pr = props
    
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_title("Matching the size of the model to the photo.")
    ax.set_xlabel('x [cm]')
    ax.set_ylabel('z [cm]')
    
    foto = Image.open(os.path.join(dirs.photos, pr['photo']))
    
    if True:
        extent = pr['extent']
        ax.imshow(foto, extent=extent)
        plt.plot([0, pr['L'], pr['L'], 0., 0.], [0., 0., pr['H'], pr['H'], 0.], label='model outer size')
        
        plt.plot(*pr['sand'].T, 'brown', lw=1, label='Zandlichaam')
        plt.plot(*pr['CanL'].T, 'c', lw=1, label='Canal left')
        plt.plot(*pr['CanM'].T, 'b.-', lw=1, label='Canal center')
        plt.plot(*pr['CanR'].T, 'c', lw=1, label='Canal right')
        plt.plot(*pr['milkInjPnt'][0][:3:2], 'ro', label='milkInjPnt')
        plt.plot(*pr['wt'].T, 'b', lw=2, label='Water table')
        plt.plot([0, 97], [pr['zIface'], pr['zIface']], 'w', lw=3, label="Interface Initial")
        
        plt.legend()
    else:
        ax.imshow(foto, extent=None)

    plt.show()
    
    
