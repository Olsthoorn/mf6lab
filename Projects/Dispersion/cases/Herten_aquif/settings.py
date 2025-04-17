import os
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

from src.mf6tools import  Dirs

# Project name (above cases)
HOME = '/Users/Theo/GRWMODELS/python/mf6lab/Projects/Pathlines/'
assert os.path.isdir(HOME), "Can't find the directory {}".format(HOME)

LENGTH_UNITS = 'cm'
TIME_UNITS = 'minutes'

dirs = Dirs(HOME)

sim_name = 'Herten_aquif'
section_name = 'Exercise1 {}.'.format(sim_name)
dirs = dirs.add_case(sim_name)
os.chdir(dirs.case)

params_wbk = os.path.join(dirs.case, sim_name + '.xlsx')
assert os.path.isfile(params_wbk), "Params_wbk not found: {}".format(params_wbk)

date_times = [
    "1904-05-01T09:22:00",   # ink added    year,month,1,10,19;   % opname p12
]

props = {
    'start_date_time': '1905-05-01T09:22:00',
    'photo': 'herten_outcrop.tif', # background phohto
    'extent': (0, 16, 0, 7), # extent of photo
    'oc_frequency': 3, # Saving frequency for budget heads and concentrations
    'Qdim': 'cm3/min',
    'L': 16, # [cm]
    'H':  7, # [cm]
    'D':  1, # [cm]
    'dx': 0.05, #[cm]
    'dz': 0.05, #[cm]
    'icelltype': 0,
    'k': 650., # [m/d] get it from the facies
    'sy': 0.2,
    'ss': 1e-5,
    'disp' : {'alh': 0.2, 'ath1': 0.02,
              'ath2': 0.02, 'atv': 0.02,
              'diffc': 3.96e-4}, # Diffc in cm2/min
    'hL': 1.0, # [cm] head in left canal
    'hR': 1.0, # [cm] head in right canal
    'por': 0.38, # [-] from herten facies
    'cL': 1.0, # concentration in canal left side
    'cR': 0.0, # concentration in canal right side
    'cFresh': 0.0,
    'cSalt': 1.0
}

if __name__ == '__main__':
    
    pr = props
    
    print('Props')
    fig, axs = plt.subplots(4,1, figsize=(10, 10), sharex=True, sharey=True)
    fig.suptitle("{}\nMatching the size of the model to the photo.".format(sim_name))
    
    axs[0].set_title("{} on the model grid".format(pr['photo']))
    axs[1].set_title("{} facies".format(pr['photo']))
    axs[2].set_title("{} conductivities".format(pr['photo']))
    axs[3].set_title("{} porosity".format(pr['photo']))
    
    axs[0].set_ylabel('z [m]')
    axs[3].set_xlabel('x [m]')
    
    foto = Image.open(os.path.join(dirs.photos, pr['photo']))
    axs[0].imshow(foto, extent=pr['extent'])
    axs[1].imshow(foto, extent=pr['extent'])
    axs[2].imshow(foto, extent=pr['extent'])
    axs[3].imshow(foto, extent=pr['extent'])

    #TODO: add herten facies and herten conductivity and herten porosity to the images directory
    
    for ax in axs:
        ax.plot([0, pr['L'], pr['L'], 0,  0], [0, 0, pr['H'], pr['H'], 0], 'b', label='bbox around model')
    
    plt.show()
    
    
