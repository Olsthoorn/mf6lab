import os
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image

from src.mf6tools import  Dirs

HOME = '/Users/Theo/GRWMODELS/python/mf6lab/Projects/VoortoetsAGT/'
assert os.path.isdir(HOME), "Can't find the directory {}".format(HOME)

LENGTH_UNITS = 'meters'
TIME_UNITS = 'days'

dirs = Dirs(HOME)

section_name = 'Generic hand-drawn cross section (test)'
sim_name = 'Generic'
dirs = dirs.add_case(sim_name)
os.chdir(dirs.case)

params_wbk = os.path.join(dirs.case, sim_name + '.xlsx')
assert os.path.isfile(params_wbk), "Params_wbk not found: {}".format(params_wbk)

lay = pd.read_excel(params_wbk, sheet_name='LAY', header=0, index_col=0)

props = { 'start_date_time': '2024-01-01',
         'minDz': 0.25, # m min layer thickness
          'drain_depth': 1e-5, # m
          'cDrainage': 100., # d
          'rch': 0.001, # m/d
          'strthd': 0.0, # m initial head
}

if __name__ == '__main__':
    
    print("Run genericSectionData.py to show the cross section as used.")
