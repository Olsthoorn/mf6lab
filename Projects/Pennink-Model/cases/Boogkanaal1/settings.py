import os
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

from src.mf6tools import  Dirs

# Project name (above cases)
HOME = '/Users/Theo/GRWMODELS/python/mf6lab/Projects/Pennink-Model/'
assert os.path.isdir(HOME), "Can't find the directory {}".format(HOME)

LENGTH_UNITS = 'm'
TIME_UNITS = 'days'

dirs = Dirs(HOME)

sim_name = 'Boogkanaal1'
section_name = 'Pennink Boogkanaal (1915) {}.'.format(sim_name)
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
    'dx':   0.25, # [m]
    'dz':   0.20, # [m]
    'min_dz': 0.001, # [m]
    'icelltype': 1,
    'k':  10.0, # [m/d]
    'k33': 5.0, #  [m/d]
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
                [-20.000, 7.4000],
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
                [ 20.000, 7.6875]]),
    'hBot': np.array([
                [ 20.000,-7.0000],
                [ 20.000,-7.0000]]),
    'sand': np.array([    
                [-20.000, 7.4000],
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
                [ 17.630, 7.6875],
                [ 20.000, 7.6875],
                [ 20.000,-7.0000],
                [-17.300,-7.0000],
                [-20.000, -7.000]]),
    'canal': np.array([    
                [ -1.80, 1.26],
                [ -0.50, 0.66],
                [  0.50, 0.66],
                [  1.95, 1.26],
                [ -1.80, 1.26]]),
    'filters_15': np.array([
                [-15.00, 0.55],
                [-15.00,-1.08],
                [-15.00,-2.60],
                [-15.00,-4.00]]),
    'filters_10': np.array([
                [-9.95,  0.50],
                [-9.95, -1.46],
                [-9.95, -2.94],
                [-9.95, -4.75]]),
    'filters_05': np.array([
                [-5.00, 0.60],
                [-5.00,-1.35],
                [-5.00,-3.30],
                [-5.00,-5.10]]),
    'filters_00': np.array([
                [0.05, -0.60],
                [0.05, -2.10],
                [0.05, -3.60],
                [0.05, -5.05],
                [0.05, -6.60]]),
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

    ax.plot(*pr['sand'].T, 'brown', label='sand')
    ax.plot(*pr['canal'].T, 'black', label='canal')
    for f in ['filters_15', 'filters_10', 'filters_05', 'filters_00']:
        ax.plot(props[f][:, 0], props[f][:, 1], 'o', label=f)
    plt.show()
    
    
