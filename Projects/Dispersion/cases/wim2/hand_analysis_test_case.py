# Hand analysis test case. This allowed verification of the modflow and modpath outcomes.
# Both succeeded.

# TO 2025-04-22

# The test case has the same cross section of 100 rows by 1000 column.
# This test is layered with no heterogeneities.
# Each block of 10 layers has the same conductivity.
# The conductivity of successive blocks increase from [1, 2, 3, 4, ... 10]
# The gradient and porosity are the same as in all other cases. J=0.0036 and por=0.31.

# Hand computation is now straightforward. The hand computed discharge, times and distances
# are compared with the values from modflow and modpath.

# %%
import numpy as np
import matplotlib.pyplot as plt
from pass_times_analysis import *

# %%
k_field = props['k_field']
kx = k_field[:, 0]
J = props['k_field_pars']['J']
gr = props['gr']
por = props['por']

Q = kx.mean() * gr.dz.sum() * J # Klopt met modflow
print(f"Q_model = {Q} m2/d")

# %% Travel times
xObs = pass_times['xObs'][0]
ku = np.unique(kx)
times = xObs[None, :] / (ku[: ,None] * J / por)
print("Hand computed:")
print(np.sort(np.round(times, 0), axis=0))

print("Modpath pathlines interpolated")
print(np.unique(np.round(pass_times['time'], 0), axis=0))

print("Modflow pass times / hand pass times")
print(np.round(np.sort(times, axis=0) / np.unique(np.round(pass_times['time'], 0), axis=0), 2))
# %%
tObs = pass_xvals['tObs'][0]
print("Hand computed  x values at the observation times:")
xVals = np.round(tObs[None, :] * ku[:, None] * J / por, 2)
print(xVals)
print("Modpath computed:")
xMP = np.unique(np.round(pass_xvals['x'], 2), axis=0)
print(xMP)
print("Ratio hand computed / Modpath computed:")
print(np.round(xVals / xMP, 2))

# %%
