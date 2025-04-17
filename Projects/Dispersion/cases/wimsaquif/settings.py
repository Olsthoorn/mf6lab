# %%
import os
import numpy as np
import matplotlib.pyplot as plt
from src.mf6tools import  Dirs
from importlib import reload
import fdm.mfgrid as mfgrid

dtypes = {
    'chd': np.dtype([('cellId', '<i8', 3), ('head', '<f4')]),
    'ghb': np.dtype([('k', '<i8'), ('i', '<i8'), ('j', '<i8'), ('bhead', '<f4'), ('cond', '<f4')]),
    'drn': np.dtype([('k', '<i8'), ('i', '<i8'), ('j', '<i8'), ('elev', '<f4'), ('cond', '<f4')]),
    'riv': np.dtype([('k', '<i8'), ('i', '<i8'), ('j', '<i8'), ('stage', '<f4'), ('cond', '<f4'), ('rbot', '<f4')]),
    'wel': np.dtype([('k', '<i8'), ('i', '<i8'), ('j', '<i8'), ('flux', '<f4')]),
}

# Project name (above cases)
# %%
HOME = '/Users/Theo/GRWMODELS/python/mf6lab/Projects/Dispersion/'
assert os.path.isdir(HOME), "Can't find the directory {}".format(HOME)

LENGTH_UNITS = 'meter'
TIME_UNITS = 'days'

# Sets up the directory structure and serves as a namespace for the directories in the project
dirs = Dirs(HOME)

# Get casename and add the directories for working with this case
sim_name = 'wimsaquif'
section_name = sim_name
dirs = dirs.add_case(sim_name)
os.chdir(dirs.case)

# %%
params_wbk = os.path.join(dirs.case, sim_name + '.xlsx')
assert os.path.isfile(params_wbk), "Params_wbk not found: {}".format(params_wbk)

# section_name = 'geotopVeluwe3'

props = {
    'sim_name': sim_name,
    'section_name': section_name,
    'oc_frequency': 1, # Saving frequency for budget heads and concentrations
    'icelltype': 0,    # [-] 0 = not convertable (confined) 1 is convertable
    # 'k_field': section_name + '.npy',          # [m/d]
    'start_date_time': str(np.datetime64('now'))[:10],
    'Qdim': 'm2/d',
    'sy': 0.2,
    'ss': 1e-5,
    'por':0.35,    
    'hStrt' : 0.0,
    'hL' : 3.6,
    'hR' : 0.0,
    'cFresh': 0.,
}

# from kfieldGeneration import k_field, gr
from MIM_k_field import k_field, k_field_pars, gr 

props.update(gr=gr, k_field=k_field, k_field_pars=k_field_pars)
props['por'] = k_field_pars['theta']
props['hL'] = props['hR'] + k_field_pars['J'] * (gr.xm[-1] - gr.xm[0])

extent = (gr.x[0], gr.x[-1], gr.z[-1], gr.z[0])
props['extent']=extent

if __name__ == '__main__':
            
    pr = props
    
    if False:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.set(title=f"{section_name}., vertical cross section, shape={gr.shape}", xlabel='x [m]', ylabel='z [m]')
        gr.plot_grid(ax=ax, row=0)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set(title=f"{section_name}, vertical cross section, shape={gr.shape}", xlabel='x [m]', ylabel='z [m]')
    kplot = ax.imshow(np.log10(k_field), aspect='auto', extent=props['extent'])
    plt.colorbar(kplot, label="log10(k)", orientation='horizontal', ax=ax)

    plt.show()
    
    
