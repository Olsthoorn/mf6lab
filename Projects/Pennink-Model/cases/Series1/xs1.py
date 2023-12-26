"""Cross section model 1 by Pennink (1915)"""

# Pennink (1915) - Grondwater Stroombanen
# Experiments series 1
#
# In his sand-box experiment, Pennink studies freshwater flow between from
# rainwater recharge to an extraction canal at the left.
# He then injects ink at 2 points at the top of the model
# and shows the developments of the stream lines in a series of photos
# that show the movement of the injected black ink.
#
# This model can be run by mt3dms as well as by seawat.
# To easier match the actual times of the tests, we will use time in minutes
# instead of days.
#
# TO 090312 100523 100722

# clear variables; close all;

import numpy as np
import matplotlib.pyplot as plt
import flopy
from fdm import Grid

basename='PenninkSeries1'

print('Tests Pennink (1915), series 1, 2 Sep 1904\n')

Scale = 65 / 113  # [cm/mm] foto p32
MW    = 65        # [cm] model width. Pennink p6
MH    = 65        # [cm] model height
zCap  = 52        # [cm] Elevation of top of full capillary zone (see description)

# Countour of sand mass obtained form photo by ginput
xSand =[-2.8568, 68.3082, 68.6354, 55.7112, 21.8465, 20.3741, 13.8302, 7.6135, -1.2208]
zSand =[-2.6501, -3.0078, 65.3106, 58.6934, 57.6204, 56.3684, 50.8243, 44.7436, 44.2071]
xCanL =[-2.3660, 5.6503, 6.7955, 8.4314, 9.2494, 9.2494, -2.8568]
zCanL =[37.4110, 37.5898, 43.1340, 53.1493, 59.4088, 63.8799, 64.5953]
xCanR =[68.6354, 60.7827, 61.6007, 62.2551, 62.7459, 63.4003, 68.4718]
zCanR =[65.3106, 62.4491, 51.1820, 46.1743, 43.3128, 41.5244, 41.3456]

## Point sources locations
xyzInk = [[54.5, 0., 55.], [38.5, 0., 55.]]

idx = xyzindex(xyzInk , gr)

# initial head

# The grid, the box is 96x96 cm and 1.8 cm thick
FRESHWATER = 0    # Relative minimum concentration
SEAWATER   = 1    # Relative maximum concentration

k = 86500 / (24 * 60) # [cm/min] hydraulic conductivity, calibrated
peff = 0.38       # [-] effective porosity, calibrated



## Grid is always 3D
dx, dy, dz = 1.0, 1.8, 1.8  # [cm] grid cell width
xGr = np.linspace(0, MW, int(MW / dx) + 1) # [cm] grid
yGr = [-dy / 2, dy / 2]      # [cm] grid
zGr = np.linspae(MH, 0, int(MH / dz) + 1)    # [cm] grid

gr = Grid(xGr,yGr,zGr)

## Times of the photos in Pennink (1915) 2 sept 1904
#  these times are used in workbook, sheet PER
year, month = 1904, 5  # day unknown, first day in known month assumed
TIMES = [np.datetime64(year, month, 1, 9, 22, 0),
         np.datetime64(year, month, 1, 9, 22, 0),  # ink added    year,month,1,10,19   # opname p12
         np.datetime64(year, month, 1, 10, 39, 0), # p14
         np.datetime64(year, month, 1, 11, 09, 0), # p16
         np.datetime64(year, month, 1, 11, 12, 0), # some ink drops added to show flow path
         np.datetime64(year, month, 1, 11, 39, 0), # p18
]

## Model arrays
idomain = gr.const(0)

# use idomain as zone array
idomain[gr.inpoly(gr.XM, gr.ZM, xSand, zSand)] =  1
idomain[gr.inpoly(gr.XM, gr.ZM, xCanL, zCanL)] = -2
idomain[gr.inpoly(gr.XM, gr.ZM, xCanR, zCanR)] =  3  # not a fixed head

# dry above canal
idomain[np.logical_and(idomain == -2, gr.ZM > hCanL)] = 0
idomain[np.logical_and(idomain ==  3, gr.ZM > hCanL)] = 0

HK = gr.const(k)
VK = gr.const(k)

HK[np.logical_and(HK > 0, gr.ZM > zCap)] = k / 10. # HK=0 above full capillary zone

# peff in saturated and unsaturation zone
PEFF = gr.const(peff)
PEFF[gr.ZM > zCap] = peff / 3  # unsaturated zone

ICBUND = idomain

# Initial head
hCanL = 45.2      # [cm] photo on p32 (canal stage on left side of model)
hCanR = ...  # [cm] photo on p32 (canal stage on right side of model)

STRTHD = gr.const(hCanL)
STCONC = gr.const(FRESHWATER)

# use unique ints to mark cells in idomain
ink1, ink2 = 5, 6
idomain[idx[1]] = ink1
idomain[idx[2]] = ink2

zoneVals = [ink1, 0, 0, ink2, 0, 0]
zoneConc = {'C_Ink1', 'C_ink2'}

-, PNTSRC = bcnZone(basename, 'CCC', idomain, zoneVals, zoneConc)

## RCH is 4.2 L/h, converted to cm3/min
W = 45. # [cm] width of rain added to top of model
M = 37. # [cm] center of rain added to top of model (mid between 2 screws see photo)
N = 4.2e3 / 60 / W / dy  # [cm/min]

NPER = getPeriods(basename)

RECH = np.zeros((gr.Ny, gr.Nx))                  # specify recharge for all stress periods
RECH[:, np.logical_and(gr.Xm >= M-0.5 * W, gr.Xm <= M + 0.5 * W),:] = N   # only where rain is applied in model
RECH = bsxfun(@times,RECH, YS[0:NPER])

# save underneath xSand zSand xCanL zCanL # needed to for PSIMASK in mf_analyze

sim = flopy.mf6.MFSimulation(..., ..., ...)
tdis = flopy.mf6....(sim, ..., ...)

model = flopy.mf6.Modflow(sim, ...)
dis = flopy.mf6.MflowDis(model, ...,)
ic = flopy.mf6.MFlowIc(model, ...)
chd = flopy.mf6.MFlowChd(model, ...)
