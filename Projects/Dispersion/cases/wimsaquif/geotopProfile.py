# Generate an array of conductivities from an image of a cross section
# The colors of the image are clustered. Then each pixel gets an integer
# cluster label. The number of distinct labels is equal to the number of
# clusters. Then the colors are mapped to conductivities. Finally
# the array of cluster labels is converted to an array of conductivities
# that can thereafter be used as the conductivties of the cross section to be modeled.
# TO20251011 with help of ChatGPT

import os
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from matplotlib.colors import ListedColormap
from matplotlib.patches import Path, PathPatch

home = '/Users/Theo/GRWMODELS/python/mf6lab/Projects/Dispersion/cases/wimsaquif'

os.chdir(home)

section_name = 'geotopVeluwe3'

# Stap 1: Afbeelding inladen
img = Image.open(os.path.join(home, 'images', section_name + '.png')).convert("RGB")  # Zorg dat je afbeelding in dezelfde folder staat
data = np.array(img)

# Plot afbeelding
fig, ax1 = plt.subplots()
plt.imshow(data)
ax1.set_title("Oorspronkelijke afbeelding")
plt.axis('off')

# Maak de kleurdata plat: (rows * cols, 3)
pixels = data.reshape(-1, 3)

# Aantal verwachte geologische facies (kleurklassen)
n_colors = 4

kmeans = KMeans(n_clusters=n_colors, random_state=0, n_init='auto')
kmeans.fit(pixels)
labels = kmeans.labels_
colors = kmeans.cluster_centers_.astype(np.uint8)


# Maak plaatje van de clusterkleuren en hun labels
fig, ax =plt.subplots()
ax.set_title("Cluster-kleuren met labels", fontsize=14)
for i, clr in enumerate(colors):  

    # Maak een vierkant voor elke cluster
    vertices = np.array([(0, i - 0.5), (1, i - 0.5), (1, i + 0.5), (0, i + 0.5), (0, i - 0.5)], dtype=float)
    codes = [1, 2, 2, 2, 79]  # Path codes voor gesloten pad
    pth = Path(vertices=vertices, codes=codes)

    p = PathPatch(path=pth, ec='k', fc=clr / 255, label=f"Cluster {i}")
    ax.add_patch(p)

    ax.plot(vertices[:, 0], vertices[:, 1], 'k', label="") # Rand van het vakje
    
    # Voeg label toe aan het vakje
    ax.text(0.5, i, f"Cluster {i}", ha="center", va="center", color="white", fontsize=12)

ax.invert_yaxis()
ax.axis('off')
ax.legend(loc='upper right')


# Maak een 2D veld met labels (gefacetteerde afbeelding)
labeled_img = labels.reshape(data.shape[:2])


# Visualiseer de clusters
custom_cmap = ListedColormap(colors / 255.0)

fig, ax2 = plt.subplots()

plt.imshow(labeled_img, cmap=custom_cmap)
ax2.set_title("Geclusterde afbeelding (facies)")
plt.colorbar(label="Cluster label", orientation='horizontal')
plt.axis("off")

# Stel zelf een mapping op:
k_values = {
    0: 1e1,  # bijvoorbeeld klei
    1: 1e-2,  # silt
    2: 1e2,  # fijn zand
    3: 1e3,  # grof zand
    4: 1e-3,  # grind
    5: 1e-1,  # tussenlaag
}

# Genereer het doorlatendheidsveld
# Dit is de array waarin elk pixel zijn doorlatendheid heeft gekregen
k_field = np.vectorize(k_values.get)(labeled_img)

# Save the array
np.save(os.path.join(home, 'data', section_name + '.npy'), k_field)

extent = [10500., 17250., -50., 0.]

fig, ax3 = plt.subplots()

plt.imshow(np.log10(k_field), cmap="viridis", origin="upper", extent=extent, aspect='auto')
plt.colorbar(label="log10(k) [m/dag]", orientation='horizontal')
ax3.set_title("Doorlatendheidsveld afgeleid uit afbeelding")
ax3.set_xlabel("Horizontale afstand (px)")
ax3.set_ylabel("Diepte (px)")
plt.show()

