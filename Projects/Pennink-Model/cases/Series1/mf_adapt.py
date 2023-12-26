""" Pennink1 - density flow in Pennnink's (1905) sand box model
   Experiments series 1
   The computation takes about 50 sec on a 2 GHz mac
   ee http://www.citg.tudelft.nl/live/pagina.jsp?id=68e12562-a4d2-489a-b82e-deca5dd32c42&lang=en
   in his sand-box experiment, Pennink studies freshwater flow between from
   rainwater recharge to aa extraction canal at the left.
   He then injects ink at 2 points at the top of the model
   and shows in a series of photos the movement of the ink.
   This model can be run by mt3dms as well as by seawat.
   To easier match the actual times of the tests, we will use time in minutes
   instead of days.

TO 090312 100523 100722
"""

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

sim_name = 'Series1' # Should match with the excel workbook with the parameters
dirs = dirs.add_case(sim_name)


params_wbk = os.path.join(dirs.case, sim_name + '.xlsx')

# chdir to case directory
os.chdir(dirs.case)
sys.path.insert(0, dirs.case)

print('Pennink (1915), series 1, 2 Sep 1904')


# %% Countour of sand mass obtained form photo by ginput
xSand =[
   -2.8568
   68.3082
   68.6354
   55.7112
   21.8465
   20.3741
   13.8302
    7.6135
   -1.2208
]

zSand =[
   -2.6501
   -3.0078
   65.3106
   58.6934
   57.6204
   56.3684
   50.8243
   44.7436
   44.2071
   ]

xCanL =[
   -2.3660
    5.6503
    6.7955
    8.4314
    9.2494
    9.2494
   -2.8568
   ];

zCanL =[
   37.4110
   37.5898
   43.1340
   53.1493
   59.4088
   63.8799
   64.5953
   ]

xCanR=[
   68.6354
   60.7827
   61.6007
   62.2551
   62.7459
   63.4003
   68.4718
]

zCanR =[
   65.3106
   62.4491
   51.1820
   46.1743
   43.3128
   41.5244
   41.3456
]

# %% Times of the photos in Pennink (1915) 2 sept 1904
#  these times are used in workbook, sheet PER
year = 1904; month=5;  # day unknown, first day in known month assumed
TIMES = datenum([
    year,month,1, 9,22,0;   % ink added    year,month,1,10,19;   % opname p12
    year,month,1,10,39,0;   % p14
    year,month,1,11,09,0;   % p16
    year,month,1,11,12,0;   % some ink drops added to show flow path
    year,month,1,11,39,0;   % p18
])


# %% Parameters
por = 0.38       # calibrated
k = 65000/24/60  # [cm/min] Conductivity calbrated from data in Pennink's series 1 experiments 
FRESHWATER = 0   # Freshwater concentration
SEAWATER   = 1   # not used in this simulation

k = 86500 / (24 * 60) # [cm/min] calibrated
peff = 0.38           # [-] effective porosity, calibrated

# %% The grid, the box is 96x96 cm and the sand-box model is 1.8 cm thick

MW    = 65.0  # [cm] Width of sand-box model, see Pennink p6
MH    = 65.0  # [cm] Top of model
D     =  1.8  # [cm] Thickness of model

zCapZone  = 51.0  # [cm] Top of full capillary zone (see descripition)
hCanL = 45.2  # [cm] Stage canal left side, see photo p32
hCanR = 45.7  # [cm] Stage right-hand canal, from given gradient and left canal stage


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
%% Point sources locations
xyzInk=[
    54.5 0 55
    38.5 0 55 ]

idx = xyzindex(xyzInk,gr)


ink1 = 5; IBOUND(idx(1)) = ink1
ink2 = 6; IBOUND(idx(2)) = ink2 

zoneVals = {ink1 0 0;ink2 0 0};
zoneConc = {'C_Ink1';'C_ink2'};

[~,PNTSRC] = bcnZone(basename,'CCC',IBOUND,zoneVals,zoneConc);

%% RCH is 4.2 L/h, converted to cm3/min

W=45; % [cm] width of rain added to top of model
M=37; % [cm] center of rain added to top of model (mid between 2 screws see photo)
N=4.2*1e3/60/W/dy;  % [cm/min]


RECH=zeros(gr.Ny,gr.Nx);                  % specify recharge for all stress periods
RECH(:,gr.Xm>=M-0.5*W & gr.Xm<=M+0.5*W,:)=N;   % only where rain is applied in model
RECH = bsxfun(@times,RECH,YS(1:NPER));
