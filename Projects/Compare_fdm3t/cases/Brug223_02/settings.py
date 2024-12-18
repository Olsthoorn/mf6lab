import os
import numpy as np
import matplotlib.pyplot as plt
from pprint import pprint

from fdm.fdm3t import fdm3t, cases
from fdm.mfgrid import Grid

from src.mf6tools import  Dirs

# Project name (above cases)
HOME = '/Users/Theo/GRWMODELS/python/mf6lab/Projects/Compare_fdm3t/'
assert os.path.isdir(HOME), "Can't find the directory {}".format(HOME)

LENGTH_UNITS = 'm'
TIME_UNITS = 'days'

dirs = Dirs(HOME)

sim_name = 'Brug223_02'
section_name = 'Bruggeman(1999) problem 223_02. Sim_name={}.'.format(sim_name)
dirs = dirs.add_case(sim_name)
os.chdir(dirs.case)

params_wbk = os.path.join(dirs.case, sim_name + '.xlsx')
assert os.path.isfile(params_wbk), "Params_wbk not found: {}".format(params_wbk)

case_name = 'Brug223_02'
props = cases[case_name]

kD = np.sum(props['kr'] * props['D'])
S  = np.sum(props['ss'] * props['D'])
t = props['t']

gr   = Grid(props['r'], [-0.5, 0.5], -np.cumsum(np.hstack((0, props['D']))), axial=False)
grAx = Grid(props['r'], [-0.5, 0.5], -np.cumsum(np.hstack((0, props['D']))), axial=True)

IDOMAIN = gr.const(1, dtype=int)

AREA = np.pi * (gr.x[1:] ** 2 - gr.x[:-1] ** 2)[np.newaxis, np.newaxis, :] * gr.const(1.)
CIRCUMF = 2 * np.pi * gr.xm[np.newaxis, np.newaxis, :] * gr.const(1.0)

kr = props['kr'][:, np.newaxis, np.newaxis] * CIRCUMF
kz = props['kz'][:, np.newaxis, np.newaxis] * CIRCUMF
ss = props['ss'][:, np.newaxis, np.newaxis] * CIRCUMF

hb = 1. # Boundary head at r=R_cylinder

oc_frequency = 1

title = props['title'] + f", nt={len(t) -1}, grMF.shape={gr.shape}, grFDM3t.shape={grAx.shape} epsilon={props['epsilon']:3g}"

if __name__ == '__main__':

    pprint(props)
    
    
