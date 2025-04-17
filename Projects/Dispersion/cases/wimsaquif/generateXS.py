# Generate a k field to generate a cross section
# %%
import numpy as np
import matplotlib.pyplot as plt
from fdm.mfgrid import Grid

# %% One block
nxb, nzb = 125, 25
tiles_x, tiles_z = 8, 3
# Aquifer:
nx, nz = tiles_x * nxb, tiles_z * nzb

kb, kc = 1., 100.

k = np.zeros((nzb, nxb), dtype=float) + kb
k[10:15, 25:50] = kc

k_tiled = np.tile(k, (tiles_z, tiles_x))

x = np.arange(nx + 1, dtype=float)
z = np.arange(nz + 1, dtype=float)[::-1]

gr = Grid(x, None, z, axial=False)

if __name__ == "__main__":
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.imshow(k_tiled)
    plt.show()


# %%
