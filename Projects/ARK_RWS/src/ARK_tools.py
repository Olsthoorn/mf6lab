import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd
from PIL import Image
from glob import glob
from tools.fdm.src import fdm3Blom, mfgrid
from tools.etc import color_cycler, logo, descr, pickleto, picklefrom
from tools.analytic.nsec1 import Nsec1
from importlib import reload

print("=== Current Session ===")
print("os.getcwd():", os.getcwd())
print("sys.executable: ", sys.executable)
print("sys.path:")
for p in sys.path:
    print(p)
print("=======================")


NOTEBOOK_NAME = "explore_data.ipynb"
print(f"NOTEBOOK_NAME = '{NOTEBOOK_NAME}'")

# --- project directory namespace
class Dirs:
    def __init__(self):
        self.rws = "/Users/Theo/Development/python/mf6_tools/mf6lab/Projects/ARK_RWS/"
        self.home   = os.path.join(self.rws, "src")
        self.images = os.path.normpath(os.path.join(self.home, '../images'))
        self.data   = os.path.normpath(os.path.join(self.home, '../data'))
        self.pb_data = os.path.normpath(os.path.join(self.data, 'peilbuizen'))
        self.pb_imag = os.path.normpath(os.path.join(self.images, 'peilbuizen'))
        
dirs = Dirs()

# --- dtypes for boundary conditions of Fdm3
dtypes = fdm3Blom.Fdm3.dtypes

print("Project directory namespace:")
descr(dirs)


raai_colors = [
    '#1f77b4',  # blue
    '#ff7f0e',  # orange
    '#2ca02c',  # green
    '#d62728',  # red
    '#9467bd',  # purple
    '#8c564b',  # brown
    '#e377c2',  # pink
    '#7f7f7f',  # gray
    '#bcbd22',  # olive
    '#17becf',  # cyan
    '#aec7e8',  # light blue
    '#ffbb78',  # light orange
    '#98df8a',  # light green    
]

leg_colors = [
    "#0000FF",  # blue
    "#FF0000",  # red
    "#00FF00",  # green
    "#000033",  # dark navy
    "#FF00B6",  # magenta
    "#005300",  # dark green
    "#FFD300",  # yellow
    "#009FFF",  # sky blue
    "#9A4D42",  # brown
    "#00FFBE",  # turquoise
    "#783FC1",  # purple
    "#1F9698",  # teal
    "#FFACFD",  # light pink
    "#B1CC71",  # olive
    "#F1085C",  # crimson
    "#FE8F42",  # orange
    "#DD00FF",  # violet
    "#201A01",  # dark brown
    "#720055",  # dark magenta
    "#766C95"   # slate gray
]

def muT_mu10C(temp):
    """Return ratio visc(TC) / visc(10C).
    
    Good approximation (from ChatGPT, 2026)
    """
    return np.exp(-0.0264 * (temp - 10))

class Peilbuizen():
    def __init__(self, workbook_name, sheet_name=None):
        self.data = pd.read_excel(os.path.join(dirs.pb_imag, workbook_name), sheet_name=sheet_name, index_col='name')
        
    def plot(self, ax, **kwargs):
        """Plot all observation well locations at once."""
        ax.plot(self.data['x'], self.data['y'], 'o', **kwargs)
        
    def put_raai_names(self, ax, xoffset=100., verbose=False, **kwargs):
        """Put the piezometer 'raai' (transect) names to the right of the last in the 'raai' + xoffset."""
        # We have transects 100, 200 ..1300 with wells 101, 102, ... 203, 204, ...1304,  etc'
        transects = np.arange(100, 1500, 100)
        for r1, r2 in zip(transects[:-1], transects[1:]):
            idx = (pbxy.data.index >= r1) & (pbxy.data.index < r2)
            
            # -- y is mean of wells in transect
            y = pbxy.data.loc[idx, 'y'].mean().round()
            
            # ---  x is rightmost well in transect
            x = pbxy.data.loc[idx, 'x'].values[-1].round() + xoffset
            
            ax.text(x, y, f'raai{r1/100:.0f}', **kwargs)
            if verbose:
                print(f"raai{r1/100:.0f}, x={x:.0f}, y={y:.0f}") 


# --- Example:
pbxy = Peilbuizen(os.path.join(dirs.pb_imag, 'pb_tauw_XY.xlsx'), sheet_name='pbTauw')

# --- Plot the data (could be on a map).
# fig, ax = plt.subplots(figsize=(8, 12))
# pbxy.plot(ax, mec='r', mfc='none', ms=6)
# pbxy.put_raai_names(ax, ha='left', va='center')
# ax.plot(ARK_center_line_GE[:,0], ARK_center_line_GE[:,1], 'b', label='ARK-hartlijn')
# ax.set_aspect(1)
# 
# plt.show()

def point_line_distance(P, A, B):
    """
    Distance from points P to infinite line through A and B.

    Parameters
    ----------
    P : (..., 2) array_like
        One or more points.
    A, B : (2,) array_like
        Two points defining the line.

    Returns
    -------
    d : ndarray
        Distances.
    """

    P = np.asarray(P)
    A = np.asarray(A)
    B = np.asarray(B)

    AB = B - A
    AP = P - A

    cross = AB[0] * AP[..., 1] - AB[1] * AP[..., 0]

    return np.round(np.abs(cross) / np.linalg.norm(AB), 2), np.sign(cross)




def point_segment_distance(P, A, B):
    P = np.asarray(P)
    A = np.asarray(A)
    B = np.asarray(B)

    AB = B - A
    AP = P - A

    L2 = np.dot(AB, AB)

    t = np.sum(AP * AB, axis=-1) / L2
    t = np.clip(t, 0.0, 1.0)

    Q = A + t[..., None] * AB

    return np.linalg.norm(P - Q, axis=-1)


def point_polyline_distance(P, polyline):
    """
    Distance from points P to a polyline.

    Parameters
    ----------
    P : (..., 2) array
        One or more points.
    polyline : (n, 2) array
        Vertices of the polyline.

    Returns
    -------
    dmin : (...)
        Minimum distance.
    Qmin : (..., 2)
        Closest points on the polyline.
    """

    P = np.asarray(P)
    polyline = np.asarray(polyline)

    dmin = np.inf * np.ones(P.shape[:-1])
    Qmin = np.zeros(P.shape)

    for A, B in zip(polyline[:-1], polyline[1:]):

        AB = B - A
        AP = P - A

        L2 = np.dot(AB, AB)

        t = np.sum(AP * AB, axis=-1) / L2
        t = np.clip(t, 0.0, 1.0)

        Q = A + t[..., None] * AB

        d = np.linalg.norm(P - Q, axis=-1)

        mask = d < dmin

        dmin[mask] = d[mask]
        Qmin[mask] = Q[mask]

    return dmin, Qmin

