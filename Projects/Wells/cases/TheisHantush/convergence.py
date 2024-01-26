# Read and analyze convergence
import os
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import exp1

from settings import dirs

os.chdir(dirs.case) #os.chdir(os.path.dirname(__file__))
logging.info("Running from {}".format(os.getcwd()))

inner = pd.read_csv(os.path.join(dirs.SIM, 'inner.csv'), header=0)
outer = pd.read_csv(os.path.join(dirs.SIM, 'outer.csv'), header=0)


# Compare with Threis:
Q, kD, Sy, r = 1200, 50 * 20, 0.1, 1e4
hstart = 0
times = np.logspace(0, np.log10(10958), 21)
r = np.logspace(0, np.log10(10000), 51)

fig, ax = plt.subplots()
ax.set_title("Theis")
ax.set_xlabel("r [m]")
ax.set_ylabel("h [m]")
ax.set_xscale('log')
ax.grid(True)

for t in times:
    s = Q / (4 * np.pi * kD) * exp1(r ** 2 * Sy / (4 * kD * t))
    ax.plot(r, hstart - s)
    
plt.show()


