# Getting the coordinates of layer boundaries using digitized points (https://automeris.io/wpd/)

# %%
import numpy as np
import matplotlib.pyplot as plt

# %%
def get_raw_coordinates(fname, xmin, xmax, xbreak):
    """Return the raw digitized coordinates using webplotdigitizer website.
    
    Parameters
    ----------
    fname: str
        text file name  or path with the raw
        coordinates in two columns without header (x, y)
        Each layer elevation starts with low x-value and ends with a high x-value.
    xmin, xmax: float
        minium and maximum x-value
        Each layer should start with an x < xmin and end with an x > xmax
    xbreak: float
        Location of vertical fault(s)
    """
    xy = np.load(fname)
    with open(fname, 'r') as fp:
        xy = fp.readlines
    
    planes = chunkit(xy, xmin, xmax)

# %%
def chunkit(fname):
    xy = np.load(fname)
    return xy

fname = 
xy = chunkit(fname)