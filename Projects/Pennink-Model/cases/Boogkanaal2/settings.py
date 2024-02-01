import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from PIL import Image

from src.mf6tools import  Dirs

# Project name (above cases)
HOME = '/Users/Theo/GRWMODELS/python/mf6lab/Projects/Pennink-Model/'
assert os.path.isdir(HOME), "Can't find the directory {}".format(HOME)

LENGTH_UNITS = 'm'
TIME_UNITS = 'days'

dirs = Dirs(HOME)

sim_name = 'Boogkanaal2'
section_name = 'Pennink (1915) {}.'.format(sim_name)
dirs = dirs.add_case(sim_name)
os.chdir(dirs.case)

params_wbk = os.path.join(dirs.case, sim_name + '.xlsx')
assert os.path.isfile(params_wbk), "Params_wbk not found: {}".format(params_wbk)

props = {
    'start_date_time': '1905-01-01',
    'photo': 'PenninkBoogkanaal2_correct.jpg', # background phohto
    'extent': np.array([-0.35, 10.35, -9.31, 2.00]), # extent of photo [L, R, B, T]
    
    'oc_frequency': 1, # Saving frequency for budget heads and concentrations
    'Qdim': 'm2/d',
    'dx':   0.25, # [m]
    'dz':   0.20, # [m]
    'dx1': 0.05, # [m] dx where grid is refined
    'dz1': 0.05, # [m] dz where grid is refined
    'min_dz': 0.001, # [m]
    'icelltype': 1,
    'k':   10.0, # [m/d]
    'k33': 10.0, #  [m/d]
    'sy': 0.2,  # [-1]
    'ss': 1e-4, # [1/m]
    'drain_depth': 0.1, # Used to determine top active cells
    'hL': 1.00, # [m] head left boundary
    'hC': 0.00, # [m] head canal
    'hR': 1.00, # [m] head ricgt boundary
    'IDL': 2, # IDOMAIN value for fixed head cells left boundary
    'IDC': 3, # IDOMAIN value for fixed head in canal
    'IDR': 4, # IDOMAIN value for fixed head cells right boundary
    'hTop': np.array([
                [ 0.00, 1.00],
                [ 0.42, 1.00],
                [ 0.47, 1.52],
                [ 0.57, 1.52],
                [ 0.72, 1.50],
                [ 0.98, 1.46],
                [ 1.51, 1.42],
                [ 1.89, 1.39],
                [ 2.40, 1.40],
                [ 2.66, 1.39],
                [ 2.98, 1.32],
                [ 3.44, 1.27],
                [ 4.12, 1.28],
                [ 4.32, 1.23],
                [ 4.52, 1.25],
                [ 4.63, 0.00],
                [ 5.37, 0.00],
                [ 5.48, 1.21],
                [ 5.63, 1.28],
                [ 5.87, 1.30],
                [ 6.43, 1.33],
                [ 7.03, 1.35],
                [ 7.41, 1.40],
                [ 7.90, 1.51],
                [ 8.45, 1.56],
                [ 8.92, 1.57],
                [ 9.22, 1.62],
                [ 9.49, 1.60],
                [ 9.55, 1.00],
                [10.00, 1.00]]),
    'hBot': np.array([
                [ 0.0,-8.58],
                [10.0,-8.58]]),
    'sand': np.array([
                [ 0.00, 1.00],
                [ 0.42, 1.00],
                [ 0.47, 1.52],
                [ 0.57, 1.52],
                [ 0.72, 1.50],
                [ 0.98, 1.46],
                [ 1.51, 1.42],
                [ 1.89, 1.39],
                [ 2.40, 1.40],
                [ 2.66, 1.39],
                [ 2.98, 1.32],
                [ 3.44, 1.27],
                [ 4.12, 1.28],
                [ 4.32, 1.23],
                [ 4.52, 1.25],
                [ 4.63, 0.00],
                [ 5.37, 0.00],
                [ 5.48, 1.21],
                [ 5.63, 1.28],
                [ 5.87, 1.30],
                [ 6.43, 1.33],
                [ 7.03, 1.35],
                [ 7.41, 1.40],
                [ 7.90, 1.51],
                [ 8.45, 1.56],
                [ 8.92, 1.57],
                [ 9.22, 1.62],
                [ 9.49, 1.60],
                [ 9.55, 1.00],
                [10.00, 1.00],
                [10.00,-8.62],
                [ 0.00,-8.62],
                [ 0.00, 1.00]]),
    'canalL': np.array([
                [0.00,1.00],
                [0.00,0.35],
                [0.07,0.35],
                [0.18,0.40],
                [0.28,0.48],
                [0.36,0.59],
                [0.39,0.70],
                [0.41,0.85],
                [0.42,1.00],
                [0.00,1.00]]),
    'canalC': np.array([
                [4.63, 0.00],
                [4.64,-0.05],
                [4.68,-0.15],
                [4.74,-0.24],
                [4.82,-0.30],
                [4.91,-0.34],
                [4.99,-0.35],
                [5.11,-0.33],
                [5.21,-0.27],
                [5.28,-0.21],
                [5.33,-0.13],
                [5.36,-0.06],
                [5.37, 0.00],
                [4.63, 0.00]]),
    'canalR': np.array([
                [ 9.56,1.00],
                [ 9.57,0.96],
                [ 9.59,0.80],
                [ 9.63,0.70],
                [ 9.68,0.62],
                [ 9.75,0.55],
                [ 9.83,0.51],
                [ 9.92,0.48],
                [10.00,0.47],
                [10.00,1.00],
                [ 9.56,1.00]]),
    'frame': np.array([
                [ 0.00, 1.62],
                [ 0.00,-8.58],
                [10.00,-8.58],
                [10.00, 1.62],
                [ 0.00, 1.62]]),
    'photoExtent': np.array([
                [-0.35, 2.00],
                [-0.35,-9.31],
                [10.35,-9.31],
                [10.35, 2.00],
                [-0.35, 2.00]])
    }


def canal_patches(props, alpha=None):
    """Return a list of patches for the canals."""
    # Build patches
    ptchs = []
    # Set left and  right most points equal to gr.x[0] and gr.x[-1]
    for canal in ['canalL', 'canalC', 'canalR']:
        xy = props[canal]
        codes = np.zeros(len(xy), dtype=int) + Path.LINETO
        codes[0] = Path.MOVETO
        codes[-1] = Path.CLOSEPOLY
        pth = Path(xy, codes)
        ptchs.append(patches.PathPatch(pth, fc='blue', ec='k', lw=0.25, alpha=alpha, zorder=5))        
    return ptchs

props['canalPatches'] = canal_patches(props)


if __name__ == '__main__':
    
    pr = props
    
    print('Props')
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_title("Matching the size of the model to the photo.")
    ax.set_xlabel('x [cm]')
    ax.set_ylabel('z [cm]')
    
    foto = Image.open(os.path.join(dirs.photos, pr['photo']))
    ax.imshow(foto, extent=pr['extent'])
    
    ax.plot(*pr['frame'].T,  'red',  label='frame')
    ax.plot(*pr['photoExtent'].T, 'blue', label='photoExtent')
    ax.plot(*pr['sand'].T, 'brown', label='sand')
    ax.plot(*pr['canalL'].T, 'black', label='canalL')
    ax.plot(*pr['canalC'].T, 'black', label='canalC')
    ax.plot(*pr['canalR'].T, 'black', label='canalR')
    
    ax.legend(loc='best')
    
    plt.show()
    
    
