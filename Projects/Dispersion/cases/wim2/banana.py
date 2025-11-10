# Banana shape for lenses
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import filtfilt

def banana(xc=0., b=100., A=1., n=2, x=None):
    """Return banana thickness.
    
    Parameters
    ----------
    xc: float
        Center of banana.
    b: float
        Half-length of banana.
    A: float
        Thickness at center of banana.
    sigma_x: float
        Sharpness of banana nodes (Don't know if required for the shape).
    """    
    y = np.zeros_like(x)
    mask = np.abs(x - xc) < b
    y[mask] = A * (1 - (np.abs(x[mask] - xc) / b) ** n)
    return y

bbase = 200.
Abase = 1.0
nbase = 1.

A, b, n = Abase, bbase, nbase
xbase = np.linspace(-4 * bbase, 4 * bbase, 201)
x     = xbase[(xbase >= - 2 * bbase) & (xbase <= 2 * bbase)]

# Begin met een gewelfd basisvlak (bijv. sinusvormig)
zbase, ampl, wlen, nwave = -10., 0.3, 2 * b, 3

base = zbase + ampl * np.random.rand(len(xbase))

fig, ax = plt.subplots(figsize=(10, 6))
title="base"
m = 10
b = np.ones(m) / m 
for i in range(1):
    base = filtfilt(b, 1., base)
    ax.plot(xbase, base, label=f'i = {i}')

base = base[ (xbase >= -2 * bbase) & (xbase <= 2 * bbase)]

fig, ax = plt.subplots(figsize=(10, 6))
title=f"Generated lenses using A={A}, b={b}, n{n}"

Y = np.zeros((100, len(x)))
Y[0] = base

y = np.zeros_like(x)

for i in range(len(Y) - 1, 0, -1):
    xc = x[np.argmin(Y[i] + zbase)]    
    A = (np.random.random(1) * 0.5 + 0.2) * Abase
    b = (np.random.random(1) * 0.5 + 0.2) * bbase
    n = np.random.random(1) * 8  + 2
    dy = banana(xc=xc, A=A, b=b, n=n, x=x)
    y[:] = Y[i]
    Y[i -1] = Y[i] + 0.5 * dy
    Y[i   ] = y - 0.5 * dy
    
 
for y1, y2 in zip(Y[:-1], Y[1:]):
    mask = y1 < y2
    y2[mask] = y1[mask]
 
for y in Y:   
    ax.plot(x, y)

ax.grid()
ax.legend()
plt.show()




