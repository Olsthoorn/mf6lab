# Generate a profile

import numpy as np
import matplotlib.pyplot as plt

# 1. Borehole positions (x) and depths (z)
borehole_x = [0, 50, 100, 150, 200]  # meters
depth = np.linspace(0, 30, 100)      # 30 meters deep

# 2. Synthetic lithology values for each borehole (values = codes for materials)
lithologies = [
    np.piecewise(depth, [depth < 10, (depth >= 10) & (depth < 20), depth >= 20], [1, 2, 3]),
    np.piecewise(depth, [depth < 12, (depth >= 12) & (depth < 22), depth >= 22], [1, 2, 3]),
    np.piecewise(depth, [depth < 8,  (depth >= 8) & (depth < 18), depth >= 18], [1, 2, 3]),
    np.piecewise(depth, [depth < 15, (depth >= 15) & (depth < 25), depth >= 25], [1, 2, 3]),
    np.piecewise(depth, [depth < 10, (depth >= 10) & (depth < 20), depth >= 20], [1, 2, 3]),
]

# 3. Interpolate lithology in x-z space
X, Z = np.meshgrid(borehole_x, depth)
Lith = np.array(lithologies).T  # shape = (len(depth), len(boreholes))

# 4. Plot
fig, ax = plt.subplots(figsize=(10, 5))
cmap = plt.get_cmap('Set3', 3)  # 3 discrete colors
cs = ax.pcolormesh(X, Z, Lith, shading='auto', cmap=cmap)
cbar = plt.colorbar(cs, ax=ax, ticks=[1.33, 2, 2.66])
cbar.ax.set_yticklabels(['Sand', 'Silt', 'Clay'])  # Example labels

ax.invert_yaxis()
ax.set_xlabel('Distance (m)')
ax.set_ylabel('Depth (m)')
ax.set_title('Synthetic Cross Section from Boreholes')
plt.show()
