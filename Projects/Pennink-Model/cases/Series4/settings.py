import os
import numpy as np

from src.mf6tools import  Dirs

# Project name (above cases)
HOME = '/Users/Theo/GRWMODELS/python/mf6lab/Projects/Pennink-Model/'
assert os.path.isdir(HOME), "Can't find the directory {}".format(HOME)

LENGTH_UNITS = 'cm'
TIME_UNITS = 'minutes'

dirs = Dirs(HOME)

sim_name = 'Series4'
section_name = 'Pennink (1915) {}.'.format(sim_name)
dirs = dirs.add_case(sim_name)
os.chdir(dirs.case)

params_wbk = os.path.join(dirs.case, sim_name + '.xlsx')
assert os.path.isfile(params_wbk), "Params_wbk not found: {}".format(params_wbk)

# lay = pd.read_excel(params_wbk, sheet_name='LAY', header=0, index_col=0)

props = {
    'start_date_time': '1904-10-05T09:35:00',  # Start time of Pennink test series 3
    'photo': 'Series4_01_p64.jpg', # background photo
    'extent': (-24.3, 155.2, -22.8, 128.8), # Extent of the Photo background.
    'oc_frequency': 2, # Saving frequency for budget heads and concentrations
    'Qdim': 'cm3/min', # Flow dimension
    'L': 97,           # [cm]
    'H': 96,           # [cm]
    'D': 1.8,          # [cm]
    'dx': 1.0,         # [cm]
    'dz': 1.0,         # [cm]
    'icelltype': 0,    # [-] 0 = not convertable (confined) 1 is convertable
    'k_mpd': 650.,     # [m/d] calbrated from data in Pennink's series 1 experiments 
    'k': 86500/(24*60), # [cm/min] 
    'sy': 0.2,         # [ -  ]
    'ss': 1e-4,        # [1/cm]
    'disp' : {'alh':  0.2, 'ath1': 0.02, 'ath2': 0.02, 'atv': 0.02, 'diffc': 3.96e-4}, # Diffc in cm2/min
    'por': 0.38,       # [-] porosity
    'zIface': 25,      # [cm] elavation of initial interface
    'rhoFresh': 1.,    # [g/cm3]
    'rhoSalt': 1.03, # [g/cm3] milk acc. to Pennink
    'cFresh': 0.0,     # [g/cm3]
    'cSalt':  1.0,     # [g/cm3] 
    'IDSD': 2,         # IDOMAIN value for initial area with sand.
    'IDMK': 3,         # IDOMAIN value for intial area with milk. 
    'iMlkInjPnt': 4, # IDOMAIN value for milk injection point.
    'milkInjPnt': [(91.5, 0., 1.5)],     # xyz-Coords of milk injection point
    'hMid': 91.0, # head in the center top of model (drain point)
    'hStrt' : 94.0, # Uniform start head.
    'hfMilk': 92.0, # cm, head in milk injection point.
    'Qmilk': 2.14, # cm3/min (Series 3, not Series 4). # TODO
    'sand' :np.array([[ 0.0,  0.0],
                      [97.0,  0.0],
                      [97.0, 94.0],
                      [48.0, 91.0],
                      [ 0.0, 94.0],
                      [ 0.0,  0.0]]),
}


if __name__ == '__main__':
    
    import matplotlib.pyplot as plt
    from PIL import Image
    
    pr = props
    
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_title("Matching the size of the model to the photo.")
    ax.set_xlabel('x [cm]')
    ax.set_ylabel('z [cm]')
    
    foto = Image.open(os.path.join(dirs.photos, 'Series4_01_p64.jpg'))
    
    if True:
        extent = (-24.3, 155.2, -22.8, 128.8)
        ax.imshow(foto, extent=extent)
        plt.plot([0, pr['L'], pr['L'], 0., 0.], [0., 0., pr['H'], pr['H'], 0.], label='model outer size')
        
        plt.plot(*pr['sand'].T, 'r', lw=1, label='Zandlichaam')
        plt.plot(*pr['milkInjPnt'][0][:3:2], 'ro', label='milkInjPnt')
        plt.plot([0, 97], [pr['zIface'], pr['zIface']], 'w', lw=3, label="Interface Initial")
        
        plt.legend()
    else:
        ax.imshow(foto, extent=None)

    plt.show()
    
    
