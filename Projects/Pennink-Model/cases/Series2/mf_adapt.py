# %% Pennink1 - density flow in Pennnink's (1905) sand box model
# Experiments series 2 (Septemner 1904)
#
# Run takes about 75 seconds on a 2GHz Mac
#
#see http://www.citg.tudelft.nl/live/pagina.jsp?id=68e12562-a4d2-489a-b82e-deca5dd32c42&lang=en
#
# in this experiment, Pennink studies freshwater flow between a recharge
# canal on the right of the model to an extraction canal at the left.
# He then injects ink at different depths near the recharge canal and
# shows in a series of photos the movement of the ink.
#
# The flow is from the right to the left canal and the flow is
# 1.25 L/h (Pennink (1905) p485 with a (water table) gradient of 1:130
#
#
# This model can be run by mt3dms as well as by seawat.
# To easier match the actual times of the tests, we will use time in minutes
# instead of days.
#
# You may want to try VDF 1 instead of -1 in the NAM worksheet
# which cause the Ink to be of seawater density
# You may also swicth off VDF in that sheet and run mt3dms to do this
# switch off VDF adn GCG and switch on FTL and LMT in the NAM worksheet
#
# TO 090312

import os
import sys
import numpy as np
from fdm import Grid
from mf6lab import mf6tools
import pandas as pd

HOME = '/Users/Theo/GRWMODELS/python/Pennink-Model/'

LENGTH_UNITS = 'centimeters'
TIME_UNITS = 'minutes'

dirs = mf6tools.Dirs(HOME)

sim_name = 'Series2' # Should match with the excel workbook with the parameters
dirs = dirs.add_case(sim_name)


params_wbk = os.path.join(dirs.case, sim_name + '.xlsx')

# chdir to case directory
os.chdir(dirs.case)
sys.path.insert(0, dirs.case)

# %% Parameters
por = 0.38       # calibrated
k = 65000/24/60  # [cm/min] Conductivity calbrated from data in Pennink's series 1 experiments 
FRESHWATER = 0   # Freshwater concentration
SEAWATER   = 1   # not used in this simulation

# %% Times of Photos taken as printed in the book Pennink(1915)

Times=['1904-09-01T09:00',  # start of experiment
       '1904-09-01T09:40',  # ink added to injection points
       '1904-09-01T09:55',
       '1904-09-01T10:20',
       '1904-09-01T10:55',
       '1904-09-01T11:25',
       '1904-09-01T13:20',  # addition of ink stopped at this time
       '1904-09-01T14:15',
       '1904-09-02T05:50',  # last photo
]

start_date_time = Times[0] # Must be a string.

datetime_photo = pd.read_excel(params_wbk, 'PER',
                     index_col=0, header=1)[['time1', 'time2', 'Foto', 'C_Ink']]
datetime_photo.index = datetime_photo.index.astype(int)


# Contour of sand mass (obtained with ginput form photo)

sand = np.array([[ 0.0,  0.0],
              [66.0,  0.0],
              [65.0, 44.3],
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
              [ 0.5, 40.8],
              [ 0. ,  0.0]])

# Contour of canal on the left side
canalL = np.array([[ 0.0, 65.00],
                   [ 0.0, 41.99],
                   [1.80, 40.72],
                   [3.95, 40.83],
                   [5.61, 42.71],
                   [6.68, 47.54],
                   [7.17, 52.26],
                   [8.33, 57.87],
                   [9.01, 65.00],
                   [ 0.0, 65.00]])

# Contour canal on right side
canalR = np.array([[60.79, 65.00],
                   [61.39, 57.70],
                   [61.59, 52.29],
                   [61.99, 48.46],
                   [62.29, 45.32],
                   [62.98, 44.53],
                   [65.00, 44.74],
                   [65.00, 65.00],
                   [60.89, 65.00]])


# %% The grid, the box is 96x96 cm and the sand-box model is 1.8 cm thick

MW    = 65.0  # [cm] Width of sand-box model, see Pennink p6
MH    = 65.0  # [cm] Top of model
D     =  1.8  # [cm] Thickness of model

zCapZone  = 51.0  # [cm] Top of full capillary zone (see descripition)
hCanL = 45.2  # [cm] Stage canal left side, see photo p32
hCanR = 45.7  # [cm] Stage right-hand canal, from given gradient and left canal stage

# %% Ink injection coordinates [cm]
xyzInk=[
   [50.6, 0, 50.0],
   [52.6, 0, 44.9],
   [53.8, 0, 33.7],
   [56.1, 0, 23.0],
]

iInk = 4   # zone number for ink injection points
cInk = 1.0 # Ink concentration

t_inkOn  = np.datetime64('1904-09-01T09:40')
t_inkOff = np.datetime64('1904-09-01T13:20')


# %% === tdis ==========  Period Data:
"""The period data are determined by the fixed period length and the moments on which
a photo was taken. The period numbers are computed from the photo times. It is at these times that one or more the boundary conditions may change, i.e. switching on or off of the ink to
show the stream lines.
"""
perDF = mf6tools.get_periodata_from_excel(params_wbk, sheet_name='PER')
nper = len(perDF)
period_data = [tuple(sp) for sp in perDF[['PERLEN', 'NSTP', 'TSMULT']].values]

# Stress period start times
sp_start_times = np.datetime64(Times[0]) + np.cumsum(
                  np.hstack((0., np.array([p[0] for p in period_data])[:-1]))
                  ) * np.timedelta64(1, 'm')

# %% === Gwfdis ========== Grid everything in cm
dx = 1.0 # [cm] model grid cell width
dy =   D # [cm] model grid cell in row direciton
dz = 1.0 # [cm] model grid cell heigth

xGr = np.arange(0, MW + dx, dx) # mesh grid lines
yGr = [-D/2, D/2]               # mesh grid lines, however the real model was 1.8 cm thick
zGr = np.arange(MH, -dz, -dz)   # mesh grid lines

gr = Grid(xGr,yGr,zGr) # mesh houskeeping

# %% ==== IDOMAIN ===============
# IDOMAIN cell values of the canal left and right
IDCL, IDCR = 2, 3 

# Start all cells inactive
IDOMAIN = gr.const(0, dtype=int)

# Only where sand --> active
IDOMAIN[:, 0, :][gr.inpoly(sand, row=0)] =  1

# Boundary cells are indacted for later recognition
IDOMAIN[:, 0, :][gr.inpoly(canalL, row=0)] = IDCL
IDOMAIN[:, 0, :][gr.inpoly(canalR, row=0)] = IDCR

# Inactive cells are within canals but above their stage level
IDOMAIN[np.logical_and(IDOMAIN == IDCL, gr.ZM > hCanL)] = 0
IDOMAIN[np.logical_and(IDOMAIN == IDCR, gr.ZM > hCanR)] = 0

# %% === Gwfnpf =========== Horizontal and vertical conductivity
HK = gr.const(k)
VK = gr.const(k)

# === Unsaturated zone conductivity (above full capillary zone)
# HK = 0 above full capillary zone, but VK > 0 to let recharge seep through
# We could also use flowthrough cells in mf6
HK[np.logical_and(HK > 0, gr.ZM > zCapZone)] = k / 10

# %% Unsaturated zone
PEFF = gr.const(por)
PEFF[gr.ZM > zCapZone] = por / 3 # My estimate for the unsaturated zone

# %% === Gwfic ======== Initial condictions (head)
STRTHD = gr.const(1) * (hCanL + hCanR) / 2. + 5 # Average between the two canal levels, extra to avoid useless first contour
STRTHD[IDOMAIN == IDCL] = hCanL # in left canal
STRTHD[IDOMAIN == IDCR] = hCanR # in right canal

# %% === CHD, Fixed head period data (Only specify the first period)
IcanL = gr.NOD.ravel()[IDOMAIN.ravel() == IDCL] # Left canal
IcanR = gr.NOD.ravel()[IDOMAIN.ravel() == IDCR] # Right canal
chd = [(lrc, hCanL) for lrc in gr.LRC(IcanL, astuples=True)] +\
      [(lrc, hCanR) for lrc in gr.LRC(IcanR, astuples=True)]
CHD = {0: chd}

# %% === DRN === for seepage face.
# TODO

# %% === OC ====
Gwfoc = {'head_filerecord':   os.path.join(dirs.SIM, "{}Gwf.hds".format(sim_name)),
         'budget_filerecord': os.path.join(dirs.SIM, "{}Gwf.cbc".format(sim_name)),
         'saverecord': [("HEAD", "ALL"), ("BUDGET", "ALL")],
}

# %% ============ T R A N S P O R T ====================

# %% === Gwtfmi ===== Flow model interface
pd = [("GWFHEAD",   os.path.join(dirs.SIM, "{}Gwf.hds".format(sim_name))),
      ("GWFBUDGET", os.path.join(dirs.SIM, "{}Gwf.cbc".format(sim_name)))
]
Gwtfmi = {'packagedata': pd}

# %% === Gwtmst ===== Mobile storage and transfer
Gwtmst = {'porosity': 0.2}
 
# %% === Gwtadv === advection =========
SCHEME = 'TVD' # upstream, central, TVD

# %% === Gwtdsp === dispersion ========
# diffc = 1e-10 m2/s 1e4 cm2/m2 60 s/min = 6e-5 m2/min
diffc = 6e-5
#diffc = 0.0 # test
DISPERSIVITIES = {'alh': 1.0, 'ath1': 0.01, 'ath2': 0.01, 'atv': 0.1, 'diffc': diffc}

# %% Gwtic === initial concentration ===

STRTC = gr.const(FRESHWATER)

# %% Gwtcnc === const conc. ====

# Ink injection points (constant concentration at injection locations)

lrc  = gr.lrc(*np.array(xyzInk).T) # global coords of injection points
IDOMAIN.ravel()[gr.I(lrc)] = iInk # For convenence, mark the cells with ink injection

# Concentration cells [(l, r, c) conc] of ink injection points.
concOn  = [(lrc_, cInk) for lrc_ in lrc]
concOff = [(lrc_,  0.0) for lrc_ in lrc]

CONSTCONC = {     perDF.index[sp_start_times <=  t_inkOn][-1]: concOn,
                  perDF.index[sp_start_times <= t_inkOff][-1]: concOff
                  }

# %% Gwtoc === output control
GWTOC = {
      'concentration_filerecord' : os.path.join(dirs.SIM, '{}Gwt.ucn'.format(sim_name)),
      'budget_filerecord':         os.path.join(dirs.SIM, '{}Gwt.cbc'.format(sim_name)),
      'saverecord' : [("CONCENTRATION", "ALL"),
                      ("BUDGET", "ALL")],
}
print('Done mf_adapt')

if __name__ == '__main__':
   
   print(dirs)     