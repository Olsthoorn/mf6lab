import numpy as np
import matplotlib.pyplot as plt

# Parameters
nx, nz = 300, 100           # Grid size
length, depth, ampl, nwave, zmin = 300, 100, 2., 3, -10.    # Physical dimensions (m)
dx, dz = length/nx, depth/nz

x = np.linspace(0, length, nx)
z = np.linspace(0, depth, nz)
X, Z = np.meshgrid(x, z)

# Begin met een gewelfd basisvlak (bijv. sinusvormig)
base = zmin + ampl * np.sin(2 * np.pi * x / length * nwave)  # 3 golven
surface = base.copy()

# Veld voor lenzen
facies = np.zeros((nz, nx))

# Functie om een lens toe te voegen
def add_lens(facies, surface, center_x, height, width_x, thickness):
    global X, Z
    x0 = center_x
    z0 = np.interp(center_x, x, surface)  # hoogte van het basisvlak op x0

    # Maak een 2D dikteprofiel (elliptisch)
    mask = ((X - x0) / width_x)**2 + ((Z - z0) / thickness)**2 <= 1
    new_lens_top = z0 + thickness * np.sqrt(1 - ((X - x0)/width_x)**2)
    new_lens_top[~mask] = -np.inf

    # Alleen toevoegen als binnen bereik en nog geen lens
    for i in range(nz):
        for j in range(nx):
            if mask[i, j] and Z[i, j] <= new_lens_top[i, j] and Z[i, j] >= surface[j]:
                facies[i, j] = 1
                surface[j] = max(surface[j], Z[i, j])

    return facies, surface

# Voeg lenzen toe tot oppervlak gevuld is
np.random.seed(0)
while np.max(surface) < depth:
    cx = np.random.uniform(0.1 * length, 0.9 * length)
    thk = np.random.uniform(2, 8)     # lensdikte
    wdx = np.random.uniform(10, 40)   # lensbreedte
    facies, surface = add_lens(facies, surface, cx, height=0, width_x=wdx, thickness=thk)

# Plot
plt.figure(figsize=(12, 6))
plt.contourf(X, Z, facies, levels=[-0.1, 0.5, 1.1], colors=['lightgrey', 'goldenrod'])
plt.plot(x, base, 'k--', lw=1, label='basisvlak')
plt.plot(x, surface, 'k-', lw=1, label='gevuld oppervlak')
plt.gca().invert_yaxis()
plt.title("Gesimuleerde lenzen in een aquifer")
plt.xlabel("x [m]")
plt.ylabel("z [m]")
plt.legend()
plt.show()
