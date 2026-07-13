# %%
"""Intersecting layers with a canal cross section.

Class Ditch_profile set the profile coordinates.

Using the class object and running method layer_path generates the path
of the polygon surrounding the layer.

The object path returns the path of the canal profile
The method patch returns the patch of the canal profile.

By running layer_path for a series of layers, the paths
around all these layers are obtaind. Then using
mpatch.PathPatch(path, close=True, fc=fc, ec=ec, **kwargs)
then generates patches to show the layers.
"""

# %%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatch
from itertools import cycle

# %%
class Ditch_profile:
    def __init__(self, xPr, zPr):        
        self.x = xPr
        self.z = zPr

    @property
    def vertices(self):
        x = np.r_[self.x, self.x[ 0], self.x[0]]
        z = np.r_[self.z, self.z[-1], self.z[0]]
        return np.c_[x, z]

    @property
    def path(self):
        return mpath.Path(self.vertices, closed=True)
    
    def patch(self, fc='blue', ec='k', **kwargs):
        return mpatch.PathPatch(self.path, fc=fc, ec=ec, **kwargs)

    def layer_path(self, zTop, zBot, xright):
        """Return path around layer by intersecting layer with ditch_profile
        """
        if zTop <= zBot:
            raise ValueError('zTop <= zBot')
        
        iTop = np.searchsorted(self.z, zTop)
        iBot = np.searchsorted(self.z, zBot)
        
        if iTop >= len(self.z) - 1:
            xn = self.x[-1]
        elif iTop == 0:
            xn = self.x[0]
        else:            
            x1, x2 = self.x[iTop-1:iTop+1]
            z1, z2 = self.z[iTop-1:iTop+1]

            t = (zTop-z1)/(z2-z1)
            xn = x1 + t*(x2-x1)
            
        if iBot >= len(self.z) -1:
            x0 = self.x[-1]
        elif iBot == 0:
            if self.z[1] == zBot:
                x0 = self.x[1]
                iBot = 1
            else:
                x0 = self.x[0]
        else:
            x1, x2 = self.x[iBot-1:iBot+1]
            z1, z2 = self.z[iBot-1:iBot+1]

            t = (zBot-z1)/(z2-z1)
            x0 = x1 + t*(x2-x1)
        
        # --- Entirely above profile
        xPoly = np.r_[x0,   self.x[iBot:iTop],   xn, xright, xright, x0]
        zPoly = np.r_[zBot, self.z[iBot:iTop], zTop, zTop,   zBot, zBot]
                
        vertices = np.c_[xPoly, zPoly]
        return mpath.Path(vertices, closed=True)

    def plot(self, ax=None, **kwargs):
        if ax is None:            
            ax = plt.gca()

        ax.plot(
            np.r_[self.x, self.x[[ 0, 0]]],
            np.r_[self.z, self.z[[-1, 0]]],
            **kwargs
        )
   
if __name__ == '__main__':   
    xPr = np.array([0.0, 20.0, 35.0, 45.0, 50.0, 50.0])  # --- x of profile
    zPr = np.array([-6.4, -6.4, -5.0, -3.3, -2.5,  -0.40]) # --- z of profile
    zTops = [-7.0, -6.4, -5.5, -2.1, -1.0,  0.0, 0.6] # --- top of layers
    zBots = [-7.5, -7.0, -6.4, -5.5, -2.1, -1.0, 0.0] # --- bottom of layers

    xmax = 2000. # --- extent of layers to the right

    # --- Layer colors
    clrs = cycle(["brown", "green", "khaki", "gold", "orange", 'gray'])

    # --- Set the profile of the right half of the disk
    dprof = Ditch_profile(xPr, zPr)

    fig, ax = plt.subplots()
    ax.add_patch(dprof.patch(fc='blue', ec='k'))

    # --- Plot of the ditch right half profile
    dprof.plot(color='b', lw=0.5)

    # --- Plot each layer as a patch
    for zTop, zBot in zip(zTops, zBots):
        clr = next(clrs)
        pth = dprof.layer_path(zTop, zBot, xmax)
        P = mpatch.PathPatch(pth, fc=clr, ec='k')
        ax.add_patch(P)

    # ax.plot(*P._path._vertices.T, 'bo')

    ax.set_xlim(0, 100)
    ax.set_ylim(-8, 1)

    ax.set_title("Layers intersecting a canal profile")
    ax.set(xlabel='x from canal center [m]', ylabel='z [m]')
    plt.show()

    #mask = poly.path.contains_points(
    #    np.c_[X.ravel(), Z.ravel()]
    #).reshape(X.shape)

