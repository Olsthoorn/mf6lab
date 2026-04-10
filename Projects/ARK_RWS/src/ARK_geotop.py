# ARK_fdm
# %%
"""Set up some classes to facilitate importing cross-section images
    of the Geotop subsurface model from Dinoloket.nl.

    Importing a geotop pdf x-section yields two pages
    
        1. The mage of the cross section itself
        2. The image with the legend and a small map showing
            where the x-section lies on the map.

    The ImagePicker below is used to
    1. to sample the colors from the legend boxes
    2. to sample the coordinates of the x-section in the maps
    3. to sample the extent of the x-section.
    4. to sample the tick marks of the axes of the x-section

    The ticmarks are used to  compute the x of the end of each
    x-section.

    The various sampled and computed data sets are pickled
    for later use. They are in dirs.data directory.
    
    All manipulation is done here on the voxels defined by
    the dx and dy of the geotop grid. Normally dx=100 m, dy=0.5 m.
    This is true as long as the x-section is parallel to the x
    or y axis on the maps. For other directions, the dx may
    be adapted. If not everything will work, but the result may
    be slightly less accurate as the sampling point will not
    coincide with each Geotop voxel.

    The CrossSectionDigitizer class will sample the cross section
    in the image when all other data are present (pxl_extent,
    world_extent etc.). This sampling is done aoutmatically.
    The result is an array with legend indices idx_arr. The
    value 0 corresponds to empty cells (with legend 'none') the
    higher values correspond to the repective legend color and label.

    The Getop_xsec class carries all relevant data for
    each section while its methods allow to visualize
    the x-secxtions and compute the parameter values
    for each voxel. These values are obtained from
    the geo_units dictionary pertaining to each x-section.
    
    The x-sections will initially show some gaps, for instance due to
    an incicion canal. These can be filled by their nearest horizontal
    neighbor. After that the x-section will be completely filled which
    may facilitate generated model property cross sections, in which
    the model grid differs from the geotop grid by another choice of
    layer and column size.
    
    For filling layers, see ARK_fdm.py / ARK_fdm.ipynb
       
    @TO 2026-03-24, 2026-04-10
    """
# %%
import os
import re
from glob import glob
from pathlib import Path
import pdf2image

import numpy as np
import matplotlib.pyplot as plt

from tools.fdm.src.mfgrid import Grid


# %%
class Dirs:
    """Local project directory namespace.
    
    To facilitate location of resources for the project
    """
    def __init__(self):
        parts = Path(os.getcwd()).parts
        idx = parts.index('ARK_RWS')
        self.home = os.path.join(*parts[:idx + 1])

        self.data   = os.path.join(self.home, 'data/')
        self.dino   = os.path.join(self.data, 'dinoloket/')
        self.gis    = os.path.join(self.data, 'gis/')
        self.doc    = os.path.join(self.home, 'doc/')
        self.images = os.path.join(self.home, 'images/')
        self.videos = os.path.join(self.home, 'videos/')
        self.src    = os.path.join(self.home, 'src/')
        self.notebooks = os.path.join(self.home, 'notebooks/')
 
 
def parse_geotop_filename(geotop_pdf_name):
    """Return xsec_type, x and y form geotop_pdf file name."""
    
    # --- Everything after "doorsnede" in the file name
    after = geotop_pdf_name.split("doorsnede", 1)[1]

    # --- Extract xsec_type, x and y
    match = re.search(r' (\D+) (\d+\.?\d*),(-?\d+\.?\d*)', after)
    
    if match:
        xsec_type = (match.group(1))
        x = float(match.group(2))
        y = float(match.group(3))
    else:
        xsec_type, x, y = None, None, None

    return xsec_type, x, y

# %%
class ImagePicker:
    def __init__(self, image):
        self.image = image

    def pick_points(self, n=1, zoom=False, title="Title of X-section"):
        """
        Click n points in the image.
        Returns list of (x, y) pixel coordinates.
        """
        plt.close('all') # temp test.
        self.fig, self.ax = plt.subplots()
        self.ax.imshow(self.image) # Default arguments
        self.ax.set_title(title + "\n" +"Click to select points, to stop press ENTER")

        if zoom:
            plt.axis('on')
        else:
            plt.axis('off')

        pts = plt.ginput(n, timeout=0)
        plt.close(self.fig)

        # Convert to integer pixel coordinates
        pts = [(int(px), int(py)) for px, py in pts]
        return pts
    
    def get_colors(self, n=-1, size=5):
        """
        Click n points and return sampled RGB colors.
        Use right-click to remove point and Enter to finish.
        """
        title="Pick the colors in sequence from the legend-boxes."
        pts = self.pick_points(n, title=title)

        colors = []
        for px, py in pts:
            half = size // 2
            patch = self.image[
                py-half:py+half+1,
                px-half:px+half+1
            ]

            # --- Median turns to float array of shap (1, 3)
            color = np.median(patch, axis=(0,1)) / 255.

            # -- List of (1,3) RGB arrays with values 0.0-1.0
            colors.append(color)
                
        return np.array(colors)
    
    def get_pxl_bbox(self, n=-1):
        """Return bbox. Zoom in and press corners. Enter when done."""
        pts = self.pick_points(n=n, title="Pick points for the bounding box.")
        print("Points picked in pixels: ", pts)
        
        # --- To avoid wrong points due to zooming, just use the
        # --- min and max of the x and y of the array of points
        px, py = np.array(pts).T
        
        pxmin, pxmax = np.min(px), np.max(px)
        pymin, pymax = np.min(py), np.max(py)
        pxl_extent = (pxmin, pxmax, pymin, pymax)
        
        return pxl_extent

    
    def show_click(self, px, py):
        """Show where you click."""
        fig, ax = plt.subplots()
        ax.imshow(self.image)
        ax.plot(px, py, 'ro')
        plt.show()
        
    def snap_to_grid(self, x, z):
        """
        Snap world coordinates to nearest voxel center.
        """
        x_min, _, z_min, _ = self.world_bbox

        col = int((x - x_min) / self.dx)
        row = int((z - z_min) / self.dz)

        # --- center of voxel
        x_snap = x_min + (col + 0.5) * self.dx
        z_snap = z_min + (row + 0.5) * self.dz

        return x_snap, z_snap
    
    def snap_pixel(self, px, py):
        """Full click -> snap -> pixel pipeline."""
        # --- pixel → world
        x, z = self.pixel_to_world(px, py)

        # --- snap
        x_snap, z_snap = self.snap_to_grid(x, z)

        # --- back to pixel
        px_snap, py_snap = self.world_to_pixel(x_snap, z_snap)

        return px_snap, py_snap
    
    def sample_snapped(self, px, py, size=3):
        px, py = self.snap_pixel(px, py)

        half = size // 2
        patch = self.image[
            py-half:py+half+1,
            px-half:px+half+1
        ]

        return patch.mean(axis=(0,1)), (px, py)
    
    def show_snap(self, px, py):
        px_s, py_s = self.snap_pixel(px, py)
        
        fig, ax = plt.subplots()
        ax.imshow(self.image)

        ax.plot(px, py, 'ro', label='click')
        ax.plot(px_s, py_s, 'gx', label='snapped')

        ax.legend()
        plt.show()


class CrossSectionDigitizer:
    """
    I want to generate a cross section as an array of soil-type indices form an image of a cross
    section in which each color refers to a soil type. The image has a legend with the soil type
    color in a small box with the soil type name next to it. The colors of these boxes link
    the colors of the image to the legend and to the soil types. The soil type indices will be
    the number of the boxes as they occur in the legend.
    
    The actual cross section is on an axes within the downloaded image. The legend is a small box
    next to the cross section image but both are on the same downloaded image.
    The actual cross section consists of rows and columns froming rectangular voxels of unifrom size
    (0.5 m vertical and 100 m horizontal). Because each voxel has its own distinct soil type and, therefore,
    color, matching the legend, the pattern of voxels clearly shows up in the cross section.
    
    The downloaded image as pixel coordinates and so has the actual cross section and the legend on it.
    The actual cross section extent on the image can be given in both pixel coordinates and in
    real world x and z coordinates which then maps each pixel to a voxel when the voxel grid is
    defined as well in real-world x-z coordinates.
    
    The color of the center of each voxel in the actual cross section (or the median of a few pixels
    around this center to make sure on always get a unique result) should be matched with the colors of
    the small boxes in the legend to get the soil-type index.
    
    This allows to fill a voxel array with soil-type indices which can then be use
    to generate a X-section groundwater flow finite difference model.
    
    What to do:
    1. Extract legend colors -> soil type indices
    2. Map image pixels -> real world coordinates
    3. Sample voxel centers -> get color
    4. Math color -> sooil index
    5. Fill array

    """
    def __init__(self, image):
        """Initialize by the Geotop image."""
        self.image = image  # numpy array (H, W, 3)

        self.legend_colors = None   # (n_types, 3)
        self.legend_labels = None   # optional

        self.pxl_bbox   = None    # (pxmin, pxmax, pymin, pymax)
        self.world_bbox = None    # (x_min, x_max, z_min, z_max)

        self.grid = None          # (nx, nz)
        self.dx = None
        self.dz = None
        
    def set_pxl_bbox(self, pxmin, pxmax, pymin, pymax):
        """The bbox of the cross section in pixels."""
        self.pxl_bbox = (pxmin, pxmax, pymin, pymax)

    def set_world_bbox(self, x_min, x_max, z_min, z_max):
        """The actual cross sections in real-world coordinates."""
        self.world_bbox = (x_min, x_max, z_min, z_max)
        
    def set_grid(self, dx, dz):
        """Define voxel grid."""
        self.dx = dx
        self.dz = dz

        x_min, x_max, z_min, z_max = self.world_bbox

        # --- Number of cells in x and z direction.
        ncols = int((x_max - x_min) / dx)
        nrows = int((z_max - z_min) / dz)

        self.grid = (ncols, nrows)
        
    def set_legend_colors(self, colors, labels=None):
        """
        colors: list of RGB tuples
        
        Manually click them.
        Or define small rectangles around legend boxes.
        """
        self.legend_colors = np.array(colors)
        self.legend_labels = labels
        
    def pixel_to_world(self, px, py):
        pxmin, pxmax, pymin, pymax = self.pxl_bbox
        x_min, x_max, z_min, z_max = self.world_bbox

        x = x_min + (px - pxmin) / (pxmax - pxmin) * (x_max - x_min)

        # --- note: image y increases downward → z increases upward
        z = z_max - (py - pymin) / (pymax - pymin) * (z_max - z_min)

        return x, z
    
    def world_to_pixel(self, x, z):
        pxmin, pxmax, pymin, pymax = self.pxl_bbox
        x_min, x_max, z_min, z_max = self.world_bbox

        px = pxmin + (x - x_min) / (x_max - x_min) * (pxmax - pxmin)
        py = pymin + (z - z_max) / (z_min - z_max) * (pymax - pymin)

        return int(px), int(py)
        
    def sample_color(self, px, py, size=5, dark_thresh=50):
        half = size // 2

        # safe slicing
        y0 = max(py - half, 0)
        y1 = min(py + half + 1, self.image.shape[0])
        x0 = max(px - half, 0)
        x1 = min(px + half + 1, self.image.shape[1])

        patch = self.image[y0:y1, x0:x1]

        # reshape to list of pixels
        pixels = patch.reshape(-1, 3)

        # compute brightness (Euclidean norm or simple sum)
        brightness = np.linalg.norm(pixels, axis=1)

        # filter out dark pixels
        mask = brightness > dark_thresh

        if np.any(mask):
            filtered = pixels[mask]
        else:
            # fallback: use all pixels if everything was filtered out
            filtered = pixels

        # robust representative color
        return np.median(filtered, axis=0)
    
    def match_color(self, color255):
        """Use nearest color in RGB space.
        
        Make sure that the color's range 0-255 matches that
        of the legend_colors which are in the range 0.0-1.0
        
        Later improvement: Use LAB color space (much better perceptually).
        
        """
        # --- color from image is 0-255 convert to range 0-1 of legend_colors
        diffs = self.legend_colors - np.array(color255) / 255.
        dist = np.sqrt((diffs**2).sum(axis=1))
        return np.argmin(dist) # is an int
    
    def build_idx_array(self):
        """Return array of size (nz,nx) filled with legend indices (dtype int).        
        """
        assert np.all(self.legend_colors[0] > 0.95), (
            "leg_colors must have wite [1., 1., 1.] as first color."
        )
        
        nx, nz = self.grid
        x_min, x_max, z_min, z_max = self.world_bbox
        x = np.linspace(x_min, x_max, nx + 1) 
        z = np.linspace(z_max, z_min, nz + 1)
        xm = 0.5 * (x[:-1] + x[1:])
        zm = 0.5 * (z[:-1] + z[1:])

        # --- Notice dtype int
        arr = np.zeros((nz, nx), dtype=int)

        for ix in range(nx):
            for iz in range(nz):
                px, py = self.world_to_pixel(xm[ix], zm[iz])
                
                # --- A bus-stop to be used for debugging
                if (ix == nx - 50) and (iz == nz - 50):
                    pass

                color = self.sample_color(px, py)
                soil_idx = self.match_color(color) # an int

                arr[iz, ix] = soil_idx
                
        # --- Very ugly and arbitrarily, cancel index 0 in line 1
        # --- Reason: due to vertical extent mismatch, the horizontal black line
        #     in the image interfers with the color of the top voxels.
        arr[0, arr[0, :] == 1] = 0
        return arr    # dtype int


def plot_result(arr, world_extent=None):
    """Plot the cross section array with voxels now the legend index.
    """
    fig, ax = plt.subplots(figsize=(12, 5))
    mappable = ax.imshow(arr, origin='upper', extent=world_extent)

    fig.colorbar(mappable, ax=ax, label='Soil type index', location='bottom')
    
    # --- Add grid lines surrounding the voxels
    nz, nx = arr.shape
    for x in np.linspace(world_extent[0], world_extent[1], nx+1):
        ax.axvline(x, color='k', lw=0.2)
    for y in np.linspace(world_extent[2], world_extent[3], nz+1):
        ax.axhline(y, color='k', lw=0.2)
     
    ax.set_xlabel('x langs doorsnede [m]')
    ax.set_ylabel('NAP [m]')
    
    ax.set_aspect(50)
    plt.show()
    

def show_leg_index_array(xsec):
    """Show the array arr with the colors given."""
    
    fname = xsec['fname']
    arr = xsec['arr']
    leg_colors = xsec['leg_colors']
    leg_labels = xsec['leg_labels']
    
    print("Plotting:", fname)
    
    arr_RGB = np.zeros((*arr.shape, 3), dtype=int)
    arr_RGB[:,:,0] = legend_colors[arr.ravel(), 0].reshape(arr.shape)
    arr_RGB[:,:,1] = legend_colors[arr.ravel(), 1].reshape(arr.shape)
    arr_RGB[:,:,2] = legend_colors[arr.ravel(), 2].reshape(arr.shape)
    
    fig, ax = plt.subplots(figsize=(10,6))
    
    fig.suptitle(fname)

    # --- Plot the legend labels in their correct color    
    nl = len(leg_labels)
    fxs = np.linspace(0.05, 0.095, nl + 1)[1:]
    fy = 0.9
    for fx, color, label in zip(fxs, leg_colors, leg_labels):
        fig.text(fx, fy, label, color=color, fontsize=10, transform=transFigure)
    
    ax.imshow(arr_RGB, origin='upper', extent=world_extent, )
    ax.set_aspect(50)
    
    xmin, xmax, zmin, zmax = world_extent
    nx, nz = arr.shape
    x = np.linspace(xmin, xmax, nx + 1)
    z = np.linspace(zmin, zmax, nz + 1)
    
    ax.vlines(x, ymin=zmin, ymax=zmax, color='k', lw=0.2)
    ax.hlines(z, xmin=xmin, xmax=xmax, color='k', lw=0.2)
    return ax
    
# --- Lithoclasses from legend of geotop X-sections
LITHO_CLASSES = {
    'a' :  {'kh':   2., 'kv': 0.2, "n":0.35, 'rho': 2600, "descr": "antrop."},
    'v' :  {'kh':   2., 'kv': 0.2, "n":0.70, 'rho': 1400, "descr": "veen"},
    'k' :  {'kh':  0.1, 'kv':0.01, "n":0.50, 'rho': 2600, "descr": "klei"},
    'kz':  {'kh':   1., 'kv': 0.1, "n":0.40, 'rho': 2600, "descr": "klei-zand"},
    'zf':  {'kh':   5., 'kv': 0.5, "n":0.35, 'rho': 2600, "descr": "fijn zand"},
    'zm':  {'kh':  15., 'kv': 1.5, "n":0.35, 'rho': 2600, "descr": "m.f. zand"},
    'zg':  {'kh':  30., 'kv': 3.0, "n":0.35, 'rho': 2600, "descr": "grof zand"},
    'g':   {'kh': 100., 'kv': 10., "n":0.30, 'rho': 2600, "descr": "grind"},
    'she': {'kh':  10., 'kv': 1.0, "n":0.45, 'rho': 2450, "descr": "schelpen"},
}

# --- Geological unis from legend of geotop X-sections
GEO_UNITS = {
    "NUAAOP"         : {"kh": 5.,   "kv": 5.,	 "n":0.40, "rho": 2600, "descr": "Anthro. Opgebr."},
    "NUECga"         : {"kh": 2.,   "kv": 0.2,	 "n":0.38, "rho": 2600, "descr": "F.v. Echteld"},
    "NUECgb"         : {"kh": 2.,   "kv": 0.2,	 "n":0.38, "rho": 2600, "descr": "F.v. Echteld"},
    "NUEC1"          : {"kh": 2.,   "kv": 0.2,	 "n":0.38, "rho": 2600, "descr": "F.v. Echteld"},
    "NUNIHO"         : {"kh": 1.,   "kv": 0.1,	 "n":0.50, "rho": 1200, "descr": "F.v.Nieuwkoop Hollandveen"},
    "NUNIBA"         : {"kh": 0.05, "kv": 0.005, "n":0.50, "rho": 1400, "descr": "F.v.Nieuwkoop Basisveenlaag"},
    "NUNBXWI-SI-KO"  : {"kh": 8.,   "kv": 0.8,   "n":0.38, "rho": 2600, "descr": "F.v.Boxtel-laagpakketten van Wierden-Singraven-Kootwijk"},
    "NUBX"           : {"kh": 5.,   "kv": 0.5,   "n":0.38, "rho": 2600, "descr": "F.v.Boxtel"},
    "NUKR-BXDE"      : {"kh": 25.,  "kv": 5.,    "n":0.35, "rho": 2600, "descr": "F.v.Krefenheye-Boxtel laagpakket Terlijnen"},
    "NUDR" 	         : {"kh": 30.,  "kv": 5.,	 "n":0.35, "rho": 2600, "descr": "F.v.Drenthe"},
    "NUgs" 	         : {"kh": 10.,  "kv": 1.,	 "n":0.35, "rho": 2600, "descr": "???"},
    "NUUR2"          : {"kh": 20.,  "kv": 2.,    "n":0.35, "rho": 2600, "descr": "F.v.Urk"},
    "NUST"           : {"kh": 40.,  "kv": 4.,    "n":0.35, "rho": 2600, "descr": "F.v.Sterksel"},
}

class Geotop_xsec:
    """Class to store and manipulate getop x-sections"""
    
    def __init__(self, fname,
                 xsec_type,
                xRD, yRD,
                pxl_extent,
                world_extent,
                xy_map_pxl,
                leg_colors,
                leg_labels,
                geo_units,
                idx_arr):
        self.name = fname
        self.xsec_type = xsec_type
        self.xRD = xRD
        self.yRD =yRD
        self.pxl_extent = pxl_extent
        self.world_extent = world_extent
        self.xy_map_pxl = xy_map_pxl
        self.leg_colors = np.array(leg_colors)
        self.leg_labels = leg_labels
        self.geo_units = geo_units
        self.idx_arr = idx_arr
        self.shape = idx_arr.shape
        
        # Set some useful properties
        self.nx, self.ny = self.shape
        
        xmin, xmax, zmin, zmax, = self.world_extent
        self.dx = (xmax - xmin) / self.nx
        self.dz = (zmax - zmin) / self.nz
        
        # --- verify leg_labels with geo_units.keys()
        s = set(self.leg_labels).difference(self.geo_units.keys()).difference(['none'])
        if not len(s) == 0:
            print("Missing leg_labels in geo_units.keys():")
            print(s)
            raise ValueError("One or more leg_labels not in geo_units.keys()")

        # --- Verify leg_labels and leg_colors:
        if not len(self.leg_labels) == len(self.leg_colors):
            raise ValueError("len(leg_colors) != len(leg_labels)")
        
        self.map_xy()
        
        return None
    
    @property
    def dx(self):
        xmin, xmax, _, _ = self.world_extent
        nx = self.shape[1]
        return (xmax - xmin) / nx
    @property
    def dz(self):
        _, _, zmin, zmax = self.world_extent
        nz = self.shape[0]
        return (zmax - zmin) / nz
    @property
    def x(self):
        nx = self.shape[1]
        xmin, xmax, _, _ = self.world_extent
        return np.linspace(xmin, xmax, nx + 1)
    @property
    def z(self):
        nz = self.shape[0]
        _, _, zmin, zmax = self.world_extent
        return np.linspace(zmax, zmin, nz + 1)
    @property
    def xm(self):
        x_ = self.x
        return 0.5 * (x_[:-1] + x_[1:])
    @property
    def zm(self):
        z_ = self.z
        return 0.5 * (z_[:-1] + z_[1:])
    @property
    def X(self):
        nz, nx = self.shape
        return np.broadcast_to(self.x[None, :], (nz + 1, nx + 1))
    @property
    def Z(self):
        nz, nx = self.shape
        return np.broadcast_to(self.z[:, None], (nz + 1, nx + 1))
    @property
    def XM(self):
        nz, nx = self.shape
        return np.broadcast_to(self.xm[None, :], (nz, nx))
    @property
    def ZM(self):
        nz, nx = self.shape
        return np.broadcast_to(self.z[:, None], (nz, nx))
    @property
    def Area(self):
        return -np.diff(self.z)[:, None] * np.diff(self.x)[None, :]

    
    def map_xy(self):
        """Add the coordinates of the X-section to self.
        
        The coordinates are obtained from the pixel coordinates
        of the map on the legend page of the geotop pdf,
        the x,y in the name of the geotop_pdf file,
        and the length of the x_section, which is the
        the xmax of the world_extent.
        
        This function is invoked only at the instantiation of this class.
        """
        # --- pixel coordinates of xsec on the small map
        # --- of the geotop pdf, p2
        dxy_pxl = np.diff(self.xy_map_pxl, axis=0)
        ds_pxl = np.sqrt((dxy_pxl ** 2).sum(axis=1))
        
        ex =  dxy_pxl.T[0] / ds_pxl # cos
        ey = -dxy_pxl.T[1] / ds_pxl # sin
        
        # --- world_dist  / pxl_dist
        L = self.world_extent[1] # --- xmax
        scale = L / ds_pxl.sum()

        # --- world length of line pieces
        ds = ds_pxl * scale

        # --- First point
        start_point = np.array([self.xRD, self.yRD])
         
        # --- Points along the cross section (bending points and end points)
        points = np.zeros_like(self.xy_map_pxl, dtype=float) + np.array([start_point])
        for i, (_ds, _ex, _ey) in enumerate(zip(ds, ex, ey)):
            points[i + 1] = points[i] + _ds * np.array([[_ex, _ey]])
        
        # --- Distance along xsec in m
        self.dist_m = np.round(np.cumsum(np.hstack((0, ds))), 0)
        
        # --- XY coordinates of Xsec points (RD-coordinaten m)
        self.xyRD = np.round(points, 0)
        return None

                      
    def get_props(self, idx_arr=None):
        """Return property arrays for all properties in geo_units."""
        arrays = {'kh': np.zeros(self.shape, dtype=float),
                  'kv': np.zeros(self.shape, dtype=float),
                  'n' : np.zeros(self.shape, dtype=float),
                  'rho': np.zeros(self.shape, dtype=float),                  
                  }
        
        if idx_arr is None:
            idx_arr = self.idx_arr
        else:
            assert np.issubtype(idx_arr.dtype, np.integer), (
                "idx_arr must be of integer dtype (is index into legend)"
            ) 
        
        for idx, label in enumerate(self.leg_labels):
            if label == 'none':
                continue
            mask = idx_arr == idx
            for variable in arrays.keys():
                arrays[variable][mask] = self.geo_units[label][variable]
                
        arrays['rho_wet'] = arrays['n'] * 1000. + (1 - arrays['n']) * arrays['rho']
        return arrays


    def fill_horizontal(self, skip_valids=0, arr=None):
        """Return horizontally filled legend index array.
        
        Fill gaps (idx=0) of the idx_arr with the nearest
        nonzero value in the same row.
        
        Empty rows remain empty.
        
        Parameters:
        -----------
        skip_valids: int (default=0)
            rows to be let empty if numver of non zeros values is <= skip_valids
        arr: array to be filled | None
            if None, then self.idx_arr will be used.
            Behavior if arr.shape != self.idx_arr.shape is uncertain.
        """
        # Use the xsec's own idx_arr if None
        if arr is None:
            arr = self.idx_arr.copy()
            
        nrows, ncols = arr.shape

        filled = arr.copy()

        # --- cell positions in the row
        x = np.arange(ncols)

        for i in range(nrows):
            row = arr[i]
            
            # --- The nonzero positions in the row
            valid = np.where(row != 0)[0]

            # --- Rows with less than skip_valids non zero indices will be empty
            if len(valid) <= skip_valids:
                row[:] = 0
                continue

            # --- Compute distance to all valid points (shape=(len(x), len(valid))
            dist = np.abs(x[:, None] - valid[None, :])

            # ---- Find nearest valid index for each position using np.argmin along x
            nearest_idx = valid[np.argmin(dist, axis=1)]

            # --- fill the line with the correct legend index values
            filled[i] = row[nearest_idx]

        return filled    
    
    
    def patch_array(self, arr, patch_extent, value):
        """Return patched array.
        
        Parameters
        ----------        
        Holds the modflow5 type grid.
        arr: np.ndarray        
        patch_extent: 4 floats
        xmin, xmax, zmin, zmax of the patch
        value: float
        value to patch
        """
        assert np.all(arr.shape == self.shape), (
            f"Your arr.shape ({arr.shape}) does not match self.shape ({self.shape})."
        )

        xpmin, xpmax, zpmin, zpmax = patch_extent
        xmin, xmax, zmin, zmax = self.world_extent
        nz, nx = self.shape
        
        x = np.linspace(xmin, xmax, nx + 1)
        z = np.linspace(zmax, zmin, nz + 1)
        xm = 0.5 * (x[:-1] + x[1:])
        zm = 0.5 * (z[:-1] + z[1:])
        XM, ZM = np.meshgrid(xm, zm) 
        
        mask = np.logical_and.reduce(
            XM > xpmin, XM < xpmax,
            ZM > zpmin, ZM < zpmax
        )        
        arr[mask] = value
        return arr
    
    def overlap(self, gr):
        """Return legend index array of gr.shape given self.idx_array of self.shape.
        
        Parameters
        ==========
        gr: mfgrid.Grid object
            grid object holding the structured fdm grid and all its properties
        
        Returns
        =======
        idx_arr: int array of gr.shape
            legend indices of the cross-section cells at the gr cell centers
        """
        INVALID = -999
        
        # --- Find the xsec bins
        Ix = np.searchsorted( self.x,  gr.xm, side='right') - 1
        Iz = np.searchsorted(-self.z, -gr.zm, side='right') - 1
        
        # --- Find points outside the xsec's world_extent
        valid_x = (Ix >= 0) & (Ix < self.idx_arr.shape[1])
        valid_z = (Iz >= 0) & (Iz < self.idx_arr.shape[0])

        # --- Only use the valid points (inside the xsec)
        Ix = Ix[valid_x]
        Iz = Iz[valid_z]

        # --- Fill the int array of gr.shape
        gr_idx_array = np.full(gr.shape, INVALID)

        gr_idx_array[np.ix_(valid_z, valid_x)] = self.idx_arr[
            Iz[:, None],
            Ix[None, :]
        ]
        return gr_idx_array

    
    def show_leg_index_array(self):
        """Show the index array with legend colors."""

        print("Plotting:", self.name)
        
        # --- Convenience shorthands
        leg_colors = np.asarray(self.leg_colors * 255, dtype=int)
        idx_arr = self.idx_arr
        
        # --- Map each of the RGB colors to its sheet
        arr_RGB = np.zeros((*self.shape, 3), dtype=int)
        arr_RGB[:,:,0] = leg_colors[idx_arr.ravel(), 0].reshape(self.shape)
        arr_RGB[:,:,1] = leg_colors[idx_arr.ravel(), 1].reshape(self.shape)
        arr_RGB[:,:,2] = leg_colors[idx_arr.ravel(), 2].reshape(self.shape)
        
        # --- Build the plot
        fig, ax = plt.subplots(figsize=(10,6))
        
        fig.suptitle(self.name)

        # --- Plot the legend labels in their correct color
        # --- First get optimal start for lbl boxes
        lx = [0]
        for lbl in self.leg_labels[1:]:
            lx.append(len(lbl) + 2)
        lx = np.array(lx)
        fxs = 0.10 + np.cumsum(lx)/sum(lx) * 0.8                        
        fy = 0.9
        
        # --- Plot the labels
        for fx, color, label in zip(fxs[:-1], self.leg_colors[1:], self.leg_labels[1:]):
            fig.text(fx, fy, label, fontsize=9, transform=fig.transFigure,
                     bbox=dict(ec='k', fc=color))
        
        # --- Show the contents (X-sec in its original Geotop colors)
        ax.imshow(arr_RGB, origin='upper', extent=self.world_extent)
        ax.set_aspect(50)
        
        # --- Plot voxel boundaries
        xmin, xmax, zmin, zmax = self.world_extent
        nx, nz = self.shape
        x = np.linspace(xmin, xmax, nx + 1)
        z = np.linspace(zmin, zmax, nz + 1)
        
        ax.vlines(x, ymin=zmin, ymax=zmax, color='k', lw=0.2)
        ax.hlines(z, xmin=xmin, xmax=xmax, color='k', lw=0.2)
        return ax
    
    def plot_array(self, arr=None, par_name='parameter?'):
        """Plot the cross section array with  values in array.
        """
        if arr is None:
            arr = self.idx_arr
            par_name = "legend-index"
        else:
            assert np.all(arr.shape == self.shape), (
                f"arr.shape {arr.shape} not equal to slf.shape {self.shape}"
            )

        fig, ax = plt.subplots(figsize=(12, 5))
        fig.suptitle(self.name)
        ax.set_title(par_name)
        ax.set(xlabel='x [m]', ylabel='z [m NAP')
        
        mappable = ax.imshow(arr, origin='upper', extent=self.world_extent)

        fig.colorbar(mappable, ax=ax, label=par_name, location='bottom')
        
        # --- Add grid lines surrounding the voxels
        nz, nx = arr.shape
        xmin, xmax, zmin, zmax = self.world_extent
        for x in np.linspace(xmin, xmax, nx+1):
            ax.axvline(x, color='k', lw=0.2)
        for y in np.linspace(zmin, zmax, nz+1):
            ax.axhline(y, color='k', lw=0.2)
        
        ax.set_xlabel('x langs doorsnede [m]')
        ax.set_ylabel('NAP [m]')
        
        ax.set_aspect(50)
        return ax
    
# %%
if __name__ == '__main__':
    dirs = Dirs()
    
    # --- It is crucial to get the geoCodes correct from the Dino-loket cross section image legend    
    geoCodes = ['NUECga', 'NUECgb', 'NUEC1', 'NUNIHO', 'NUNIBA', 'NUBXWI-SI-KO', 'NUBX', 'NUDR', 'NUgs']

    # --- It is also crucial to set proper world extent coordinates for the Dino-loket X-sec image.
    world_extent=(0, 8625, -48.5, 0)
    
    geotop_pdf = glob(dirs.dino + '*.pdf')[-1]
    geotop1, geotop2 = pdf2image.convert_from_path(geotop_pdf, dpi=300)
    geotop1 = np.asarray(geotop1.convert("RGB"))  
    geotop2 = np.asarray(geotop2.convert("RGB"))   
    
    # ---Instantiate the picker with the image to pick from (image with the legend)
    picker = ImagePicker(geotop2)
    
    # --- Zoom before clicking (zoom into legend)
    # --- Step 1: Pick legend colors
    legend_colors = picker.get_colors(n=-1)
    print(legend_colors)
        
    assert len(legend_colors) == len(geoCodes), (
        f"Len(colors) != len(geoCodes): {len(legend_colors)} != {len(geoCodes)}")
 
 
    # --- Initiate a new picker, now with the cross section   
    # --- To get the pixel bounding box
    picker = ImagePicker(geotop1)
    # --- Again, zoom in
    # --- Then pick 2 opposite corners, return to finish
    pxl_extent = picker.get_pxl_bbox(n=-1)
    print(pxl_extent)
    
    # --- Fill an array with soil indices where each index is the number of the
    # --- legend color boxes in that order (the order clicked before)
    
    # --- Instantiate the CrossSectionDigitizer using the cross section.
    digitizer = CrossSectionDigitizer(image=geotop1)
    
    # --- Set pixed extent (see extent obtained above)
    digitizer.set_pxl_bbox(*pxl_extent)
    
    # --- Set world extent (given above in world coordinates)
    digitizer.set_world_bbox(*world_extent)
    
    # --- Set the grid specifying the voxel width and the voxel height
    digitizer.set_grid(dx=100, dz=0.5)
    
    # --- Internally set the legend_colors (obtained from clicking the legend)
    digitizer.legend_colors = legend_colors
    
    # === Fill the grid array with the soil-indices (color box or legend item index)
    x = np.linspace(world_extent[0], world_extent[1], 201)
    z = np.linspace(world_extent[2], world_extent[3],  101)
    
    # === Fill in an array of a cross section using  dx and dz and world_extent
    arr = digitizer.build_idx_array()
    
    # --- Show the cross section using imshow, which fills the voxels
    plot_result(arr, world_extent)

    # --- Add title and save
    fig = plt.gcf()
    fig.suptitle(f"""{os.path.basename(geotop_pdf)}
                 with colors converted to soil-indices
                 """)
    fig.savefig(os.path.join(dirs.images, f"{os.path.basename(geotop_pdf)}"))
    
    plt.show()
    

# %%
