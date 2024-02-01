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

sim_name = 'Boogkanaal1'
section_name = 'Pennink (1915) {}.'.format(sim_name)
dirs = dirs.add_case(sim_name)
os.chdir(dirs.case)

params_wbk = os.path.join(dirs.case, sim_name + '.xlsx')
assert os.path.isfile(params_wbk), "Params_wbk not found: {}".format(params_wbk)

props = {
    'start_date_time': '1905-01-01',
    'photo': 'PenninkBoogkan1_correct.png', # background phohto
    'extent': (-19.1, 19.5, -7.90, 12.50), # extent of photo [L, R, B, T]
    'oc_frequency': 1, # Saving frequency for budget heads and concentrations
    'Qdim': 'm2/d',
    'dx':   0.10, # [m]
    'dz':   0.10, # [m]
    'dx1':  0.05, # [m] in refined area
    'dz1':  0.05, # [m] in refined area
    'min_dz': 0.001, # [m]
    'icelltype': 1,
    'k':  10.0, # [m/d]
    'k33': 7.5, #  [m/d] Optimal value (calibrated)
    'sy': 0.2,  # [-1]
    'ss': 1e-4, # [1/m]
    'drain_depth': 0.1, # Used to determine top active cells
    'hL': 1.52, # [m] head left boundary
    'hC': 1.26, # [m] head canal
    'hR': 1.52, # [m] head ricgt boundary
    'IDL': 2, # IDOMAIN value for fixed head cells left boundary
    'IDR': 3, # IDOMAIN value for fixed head cells right boundary
    'IDC': 4, # IDOMAIN value for fixed head in canal
    'hTop': np.array([    
                [-17.400, 7.4000],
                [-15.021, 7.4016],
                [-10.910, 5.3379],
                [ -9.944, 5.3402],
                [ -1.803, 1.3175],
                [ -0.502, 0.6799],
                [  0.498, 0.6473],
                [  1.969, 1.3379],
                [ 10.758, 5.4817],
                [ 11.701, 5.4606],
                [ 16.296, 7.6843],
                [ 17.60, 7.6875]]),
    'hBot': np.array([
                [-17.40,-7.0000],
                [ 17.60,-7.0000]]),
    'sand': np.array([  
                [-17.400, 7.4000], # help point leave (to start grid at -17.4)                  
                [-17.300, 7.4000],
                [-15.021, 7.4016],
                [-10.910, 5.3379],
                [ -9.944, 5.3402],
                [ -1.803, 1.3175],
                [ -0.502, 0.6799],
                [  0.498, 0.6473],
                [  1.969, 1.3379],
                [ 10.758, 5.4817],
                [ 11.701, 5.4606],
                [ 16.296, 7.6843],
                [ 17.600, 7.6875],                
                [ 17.600,-7.0000],
                [-17.400,-7.0000],
                [-17.400, 7.4000],
                ]), # Leave 17.60 rounded value
    'canal': np.array([    
                [ -1.80, 1.26],
                [ -0.50, 0.66],
                [  0.50, 0.66],
                [  1.95, 1.26],
                [ -1.80, 1.26]]),
    'filters_15': np.array([
                [-15.00, 0.55, 1.48],
                [-15.00,-1.08, 1.49],
                [-15.00,-2.60, 1.49],
                [-15.00,-4.00, 1.50]]),
    'filters_10': np.array([
                [-9.95,  0.50, 1.41],
                [-9.95, -1.46, 1.42],
                [-9.95, -2.94, 1.43],
                [-9.95, -4.75, 1.46]]),
    'filters_05': np.array([
                [-5.00, 0.60, 1.33],
                [-5.00,-1.35, 1.35],
                [-5.00,-3.30, 1.36],
                [-5.00,-5.10, 1.41]]),
    'filters_00': np.array([
                [0.05, -0.60, 1.32],
                [0.05, -2.10, 1.33],
                [0.05, -3.60, 1.36],
                [0.05, -5.05, 1.39],
                [0.05, -6.60, 1.40]]),
    }

def canal_patches(props, alpha=None):
    """Return a list of patches for the canals."""
    # Build patches
    ptchs = []
    # Set left and  right most points equal to gr.x[0] and gr.x[-1]
    for canal in ['canal']:
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

    ax.plot(*pr['sand'].T, 'brown', label='sand')
    ax.plot(*pr['canal'].T, 'black', label='canal')
    for f in ['filters_15', 'filters_10', 'filters_05', 'filters_00']:
        ax.plot(props[f][:, 0], props[f][:, 1], 'o', label=f)
    plt.show()
    
    
