# profile2 

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Borehole positions
x_bh = np.array([0, 50, 100, 150, 200])

# Shared "key" layer: base and top elevations
top_key = np.array([5, 6, 4, 7, 6])
bot_key = top_key + 4  # constant thickness

# Minor layer, only in middle boreholes (will pinch out)
top_minor = np.array([np.nan, 10, 9, 11, np.nan])
bot_minor = top_minor + 1.5

# Interpolate key layer surfaces
x_dense = np.linspace(x_bh[0], x_bh[-1], 300)
top_key_i = interp1d(x_bh, top_key, kind='cubic')(x_dense)
bot_key_i = interp1d(x_bh, bot_key, kind='cubic')(x_dense)

# Interpolate minor only where not nan
mask = ~np.isnan(top_minor)
top_minor_i = interp1d(x_bh[mask], top_minor[mask], kind='linear', bounds_error=False, fill_value=np.nan)(x_dense)
bot_minor_i = interp1d(x_bh[mask], bot_minor[mask], kind='linear', bounds_error=False, fill_value=np.nan)(x_dense)

# Plot
fig, ax = plt.subplots(figsize=(10, 5))
ax.fill_between(x_dense, top_key_i, bot_key_i, color='tan', label='Key sand layer')
ax.fill_between(x_dense, top_minor_i, bot_minor_i, color='lightblue', label='Pinching silt lens')
ax.set_ylim(15, 0)
ax.set_xlabel("Distance (m)")
ax.set_ylabel("Depth (m)")
ax.legend()
ax.set_title("Cross Section with Key Layer and Pinching Lens")
plt.show()
