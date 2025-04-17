# Criteria om te beoordelen of het aantal clusters dat bij Kmeans clustering
# moet worden gekozen optimaal is:
# TO 2025/04/11 met ChatGPT

import os
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

os.chdir('/Users/Theo/GRWMODELS/python/mf6lab/Projects/Dispersion/cases/wimsaquif')

# Stap 1: Afbeelding inladen
img = Image.open("geotopVeluwe3.png").convert("RGB")  # Zorg dat je afbeelding in dezelfde folder staat
data = np.array(img)

# Maak de kleurdata plat: (rows * cols, 3)
pixels = data.reshape(-1, 3)

# Elbow-methode
# Bij de elbow-methode plot je de inertia (de som van kwadratische afstanden tussen
# punten en hun dichtstbijzijnde clustercentrum) tegen het aantal clusters:

# Interpretatie: Kies het aantal clusters bij de "knik" (elbow) waar de afname
# in inertia plots minder snel gaat. Daar voeg je nog maar weinig toe met extra clusters.
 
inertias = []
k_range = range(1, 11)
for k in k_range:
    kmeans = KMeans(n_clusters=k, random_state=0).fit(pixels)
    inertias.append(kmeans.inertia_)

fig, ax = plt.subplots(figsize=(10, 6))
ax.set(title="Elbow methode", xlabel="Number of clusters", ylabel='Inertia (totale kwadratische fout)')
ax.plot(k_range, inertias, 'bo-')
ax.grid()

# Silhouette-score
# De silhouette-score meet hoe goed een datapunt past binnen zijn cluster vergeleken met andere
# clusters. De score ligt tussen -1 (slecht) en 1 (goed). Een hogere gemiddelde score wijst
# op betere clustering.

# Interpretatie: Kies het aantal clusters met de hoogste silhouette-score. Dit betekent dat
# je clusters goed gescheiden én compact zijn.

from sklearn.metrics import silhouette_score

scores = []
for k in range(2, 11):
    kmeans = KMeans(n_clusters=k, random_state=0).fit(pixels[:1000, :])
    score = silhouette_score(pixels[:1000, :], kmeans.labels_)
    scores.append(score)

fig, ax = plt.subplots(figsize=(10, 6))
ax.set(title='Silhouette-analyse', xlabel='Aantal clusters', ylabel='Silhouette score')
ax.plot(range(2, 11), scores, 'ro-')
ax.grid()
plt.show()

