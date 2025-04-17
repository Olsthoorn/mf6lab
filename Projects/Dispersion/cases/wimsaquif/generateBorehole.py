import numpy as np
import matplotlib.pyplot as plt

facies = ['klei', 'silt', 'fzand', 'zand', 'gzand', 'grind']
transitions = {
    'klei':  [0.6, 0.4, 0.3, 0.2, 0.1, 0.05],
    'silt':  [0.4, 0.5, 0.4, 0.3, 0.15, 0.1],
    'fzand': [0.3, 0.4, 0.5, 0.4, 0.2, 0.1],
    'zand':  [0.2, 0.3, 0.4, 0.5, 0.2, 0.1],
    'gzand': [0.1, 0.15, 0.2, 0.3, 0.5, 0.25],
    'grind': [0.05, 0.1, 0.1, 0.25, 0.4, 0.5],
}

facies_colors= {
    'klei':'darkgray',
    'silt': 'lightgray',
    'fzand': 'yellow',
    'zand': 'gold',
    'gzand': 'orange',
    'grind':'brown'}

for k, item in transitions.items():
    item = np.array(item, dtype=float) / np.array(item, dtype=float).sum()
    transitions[k] = [float(v) for v in item]
    

def generate_borehole(n_layers=10, facies=None, probs=None):
    borehole = []
    current = np.random.choice(facies)
    for _ in range(n_layers):
        borehole.append(current)
        probs = transitions[current]
        current = np.random.choice(facies, p=probs)
    return borehole

def generate_layer_thicknesses(n_layers=10, mean_thickness=1.5):
    return np.random.exponential(scale=mean_thickness, size=n_layers)

def generate_synthetic_borehole(n_layers=None, facies=None, probs=None, mean_thickness=None):
    facies_seq = generate_borehole(n_layers=n_layers, facies=facies, probs=probs)
    thicknesses = generate_layer_thicknesses(len(facies_seq), mean_thickness=mean_thickness)
    return list(zip(facies_seq, thicknesses))


def plot_borehole(borehole, location=0):
    depths = np.cumsum([0] + [t for _, t in borehole])
    facnames = set()
    for i, (facies, thickness) in enumerate(borehole):        
        plt.fill_betweenx([depths[i], depths[i+1]], location, location + 1, label=facies if facies not in facnames else "", color=facies_colors[facies])
        facnames.add(facies)
    plt.gca().invert_yaxis()


# Voorbeeld
bh = generate_synthetic_borehole(n_layers=100, facies=facies, probs=transitions, mean_thickness=1.0)
plot_borehole(bh)
plt.legend()
plt.show()
