# Generate a k field to generate a cross section

# The cross section is filled with overlapping blocks each of which
# has a distinct randomly selected width and height and condutivity.
# The height is selected from a small list [1, 2, 3, 4 5]
# The width is selected using and lognormal distribution.
# The conductivity is also selected from a lognormal disribution.
# First a label_fld is fille with the blocks, using integer labels.
# Due to overlap the set of labels is not contiuous.
# Finally each zone with a unique label get its unique k-value.

# %%

# %%
import numpy as np
import matplotlib.pyplot as plt
from fdm.mfgrid import Grid

# %% === size of the cross section. Assume dx=dz = 1 m ====
nz, nx = 100, 1000

x = np.arange(nx + 1, dtype=float)
z = np.arange(nz + 1, dtype=float)[::-1]

gr = Grid(x, None, z, axial=False)

label_fld = np.zeros((nz, nx), dtype=int) - 1
k_field   = np.zeros((nz, nx), dtype=float)

# %% == elect and fill subzones
ixlist = np.arange(nx, dtype=int)
izlist = np.arange(nz, dtype=int)
dzlist  = np.arange(1, 5)
wxlist = np.random.lognormal(mean=1.2, sigma=1, size=1000).astype(int) + 1

# We'll put an overlapping block on each randomly chosen cell of the cross section
NOD = list(np.random.choice(np.arange(nx * nz, dtype=int), size=nx * nz))

ix  = np.arange(nx, dtype=int)
iz  = np.arange(nz, dtype=int)
IX, IZ = np.meshgrid(ix, iz)

ilabel = 0
while len(NOD) > 0:    
    inod = NOD.pop() # Select a position in the cross section at random
    
    # Select width and height of the block and make
    # sure it's within the section limits
    ix1 = IX.ravel()[inod]
    iz1 = IZ.ravel()[inod]
    ix2 = ix1 + np.random.choice(wxlist)
    iz2 = iz1 + np.random.choice(dzlist)
    ix1, ix2 = max(ix1, 0), min(ix2, nx)
    iz1, iz2 = max(iz1, 0), min(iz2, nz)
    
    # Assign label
    label_fld[iz1:iz2, ix1:ix2] = ilabel
    ilabel += 1
 
 # See what and how many labels we have   
labels = np.unique(label_fld)
nlabel = len(labels)  

# %% === k_field ===============================
k = np.random.lognormal(mean=np.log10(0.89), sigma=np.sqrt(6.6), size=len(labels))

for i, label in enumerate(labels):
    k_field.ravel()[label_fld.ravel() == label] = k[i]
k_field = k_field.reshape(nz, nx)

# %% === show the k_field

if __name__ == "__main__":
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set(title="Conductivity field", xlabel='x [m]', ylabel='z [m]')
    
    c = ax.imshow(np.log10(k_field), aspect='auto', extent=(gr.xm[0], gr.xm[-1], gr.zm[-1], gr.zm[0]))
    cbar = plt.colorbar(c, cmap='viridis', orientation='horizontal', ax=ax)
    cbar.set_label("log10(k)", labelpad=10)
    ax.set_xlim(gr.x[0], gr.x[-1])
    ax.set_ylim(gr.z[-1], gr.z[0])
    plt.show()
