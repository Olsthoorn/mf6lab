# ARK_fdm
# %%
"""Set up some classes to facilitate importing cross-section images like
    Geotop from Dinoloket.nl.

    When downloading such an image and subsequently importing them
    you get two images from the same download. The first is the
    image of the cross section itself, the second is the image
    with the legend, a small map showing where the image is
    and some additional information.

    The picker below allows first importing the image with the legend.
    You can pick the colorboxes of the legend to pick the colors and
    capture them in a RGB array. The last color will always be white
    (255., 255., 255.) for layer use with the sampling of the actual
    cross section.

    Then the picker is instantiated with the image of the actual cross
    section and the ll and ur corners are picked giving the
    image's extent in pixels. Together with the separtely given
    world_extent of the image's X-sec and the horizontal and vertical
    size of the Geotop voxels in the image, a voxel array is then
    automatically filled by sampling its voxel colors. The array
    values are the index of the legend boxes in that order,
    where the last corresponds to pure whte, meaning an empty
    voxel.

    It's best to stick with the size of the Geotop voxels when
    sampling because using a finer grid may result in background
    lines in the image being interpreted as legend color, which is
    not what you want.

    Clearly, the legend index can be converted to anything else that
    corresponds to the legend index, for instance layer names, layer
    types, conductivity etc. A convenient way is to link such
    properties with the legend in a pandas DataFrame.


    The second subject is to fill a cross section model grid
    with properties. The voxels of this model grid may not
    correspond with that of the Geotop cross section used above.
    
    Given the Geotop world_extent and its dx and dz voxel size, the
    voxel of each coordinate pair is uniquely difined. So given a
    normal model x-section grid, the cell value can be sampled in the
    Geotop X-section uniquely. Moreover one can specify a slice of the
    Getop X-section to match a slice of the actual model grid and fill that.
       
    @TO 2026-03-24, 04-03
    """
# %%
import os
from glob import glob
from pathlib import Path
import pdf2image

import numpy as np
import matplotlib.pyplot as plt

from tools.fdm.src.mfgrid import Grid


# %%
class PropSection:
    """Class specifying extended cross section properties.
    
    The PropSection consists of a set of PropBlocks that
    together constitute a complete cross section.
    
    It's purpose is mainly to manage the set of PropSec objects.
    """
    def __init__(self, propsecs: list | tuple):
        self.secs = propsecs
        
    def plot(self):
        for sec in self.secs:
            sec.plot()

    def fill(self, gr, prop_name):
        A = gr.const(0.)
        for sec in self.secs:
            A = self.fill(A, prop_name)
        return A
    
    

class PropsSec:
    """Class specifying the properties of part of a vertical X-section independent of the model grid.
    
    The FDM model grid can be filled with values from these property sections.
    The FDM grid for a property will be completely filled if the total set of
    property sections cover the entire grid X-section.
    
    Property section may overlap to overwrite parts that have already been filled.
    This helps refining the model, adding things like sheet piling and exacavations
    as well a specifying different scenarios.
    
    The PropSec is instantiated with two inputs:
    
    1) The extent, i.e. the coordinates of the ends of the section.
    2) The data, which are provided as a pandas DataFrame.
        
    """
    
    def __init__(self, extent, data):
        """Instantiate a block.
        
        Parameters
        ----------
        extent: tuple
            The spatial exent of the section.
        data: pd.DataFrame
            The layer data. The data DataFrame has at least the following columns:
            id z1 z2 k1 k3 S n name type color
            z2 < z1, will be verified.
        """
        self.extent=extent
        self.props = data
        
    def plot(self):
        """Plot the prop section using patches"""
        pass
    
    
    def fill(self, gr, A, prop_name='kx'):
        """Return grid array A with property of this section filled in.
        
        Parameters
        ----------
        A is an array of gr.shape. Row 1 of this array will be
        overwritten by the property values of the currenct parameter within the
        extent of the current PropSec object.
        """
        pass
        

class ImagePicker:
    def __init__(self, image):
        self.image = image

    def pick_points(self, n=1, zoom=False):
        """
        Click n points in the image.
        Returns list of (x, y) pixel coordinates.
        """
        fig, ax = plt.subplots()
        ax.imshow(self.image)
        ax.set_title(f"Click {n} point(s), then press ENTER")

        if zoom:
            plt.axis('on')
        else:
            plt.axis('off')

        pts = plt.ginput(n, timeout=15)
        plt.close(fig)

        # Convert to integer pixel coordinates
        pts = [(int(x), int(y)) for x, y in pts]
        return pts
    
    def get_colors(self, n=-1, size=5):
        """
        Click n points and return sampled RGB colors.
        Use right-click to remove point and Enter to finish.
        """
        pts = self.pick_points(n)
        
        # --- Remove points caused by zooming. They have color [255., 255., 255.]
        pts = [p for p in pts if not np.all(np.isclose(p, 255.))]

        colors = []
        for px, py in pts:
            half = size // 2
            patch = self.image[
                py-half:py+half+1,
                px-half:px+half+1
            ]
            color = np.median(patch, axis=(0,1))
            colors.append(color)

        colors = [clr for clr in colors if not np.all(np.isclose(clr, 255.))]
        colors.append([255., 255., 255.])
        return np.array(colors)
    
    def get_bbox(self, n=-1):
        """Return bbox. Zoom in and press corners. Enter when done."""
        pts = picker.pick_points(n=n)
        print(pts)
        # --- To avoid wrong points due to zooming, just use the last two points
        pts = pts[-2:]
        print(pts)
        
        # --- Make extent
        ((x1, y1), (x2, y2)) = pts
        xmin, xmax = sorted([x1, x2])
        ymin, ymax = sorted([y1, y2])
        extent = (xmin, xmax, ymin, ymax)
        
        print(extent)
        return extent

    
    def show_click(self, px, py):
        fig, ax = plt.subplots()
        ax.imshow(self.image)
        ax.plot(px, py, 'ro')
        plt.show()
        
    def snap_to_grid(self, x, z):
        """
        Snap world coordinates to nearest voxel center.
        """
        x_min, _, z_min, _ = self.world_bbox

        ix = int((x - x_min) / self.dx)
        iz = int((z - z_min) / self.dz)

        # --- center of voxel
        x_snap = x_min + (ix + 0.5) * self.dx
        z_snap = z_min + (iz + 0.5) * self.dz

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

        import matplotlib.pyplot as plt
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

        self.section_bbox = None    # (xmin, xmax, ymin, ymax) in pixels
        self.world_bbox = None      # (x_min, x_max, z_min, z_max)

        self.grid = None            # (nx, nz)
        self.dx = None
        self.dz = None
        
    def set_section_bbox(self, xmin, xmax, ymin, ymax):
        """The bbox of the cross section in pixels."""
        self.section_bbox = (xmin, xmax, ymin, ymax)

    def set_world_bbox(self, x_min, x_max, z_min, z_max):
        """The actual cross sections in real-world coordinates."""
        self.world_bbox = (x_min, x_max, z_min, z_max)
        
    def set_grid(self, dx, dz):
        """Define voxel grid."""
        self.dx = dx
        self.dz = dz

        x_min, x_max, z_min, z_max = self.world_bbox

        # --- Number of cells in x and z direction.
        nx = int((x_max - x_min) / dx)
        nz = int((z_max - z_min) / dz)

        self.grid = (nx, nz)
        
    def set_legend_colors(self, colors, labels=None):
        """
        colors: list of RGB tuples
        
        Manually click them.
        Or define small rectangles around legend boxes.
        """
        self.legend_colors = np.array(colors)
        self.legend_labels = labels
        
    def pixel_to_world(self, px, py):
        pxmin, pxmax, pymin, pymax = self.section_bbox
        x_min, x_max, z_min, z_max = self.world_bbox

        x = x_min + (px - pxmin) / (pxmax - pxmin) * (x_max - x_min)

        # --- note: image y increases downward → z increases upward
        z = z_max - (py - pymin) / (pymax - pymin) * (z_max - z_min)

        return x, z
    
    def world_to_pixel(self, x, z):
        xmin, xmax, ymin, ymax = self.section_bbox
        x_min, x_max, z_min, z_max = self.world_bbox

        px = xmin + (x - x_min) / (x_max - x_min) * (xmax - xmin)
        py = ymin + (z_max - z) / (z_max - z_min) * (ymax - ymin)

        return int(px), int(py)
    
    # def sample_color(self, px, py, size=3):
    #     """Robust sampling. Don't use a single pixel, use a small window."""
    #     half = size // 2
    #     patch = self.image[
    #         py-half:py+half+1,
    #         px-half:px+half+1
    #     ]
    #     return patch.mean(axis=(0,1))
    
    def sample_color(self, px, py, size=5, dark_thresh=40):
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
    
    def match_color(self, color):
        """Use nearest color in RGB space.
        Later improvement: Use LAB color space (much better perceptually).
        
        """
        diffs = self.legend_colors - np.array(color)
        dist = np.sqrt((diffs**2).sum(axis=1))
        return np.argmin(dist)
    
    def build_array(self, gr=None):
        """Return array of size (nz,nx) filled with legend index.
        
        Parameters
        ----------
        gr: Grid object | None
            grid object used with fdmr (in tools/fdm/src/mfgrid)
            Allows filling an arbitrary rectangular grid.
        """
        
        if gr is None:
            nx, nz = self.grid
            x_min, x_max, z_min, z_max = self.world_bbox
            x = np.linspace(x_min, x_max, nx + 1) 
            z = np.linspace(z_max, z_min, nz + 1)
            xm = 0.5 * (x[:-1] + x[1:])
            zm = 0.5 * (z[:-1] + z[1:])
        else:
            self.world_bbox = (gr.x[0], gr.x[-1], gr.z[-1], gr.z[0])
            nx, nz = gr.nx, gr.nz
            x = gr.x
            z = gr.z
            xm = gr.xm
            zm = gr.zm

        arr = np.zeros((nz, nx), dtype=int)

        for ix in range(nx):
            for iz in range(nz):
                px, py = self.world_to_pixel(xm[ix], zm[iz])

                color = self.sample_color(px, py)
                soil_idx = self.match_color(color)

                arr[iz, ix] = soil_idx
        return arr    


def plot_result(arr, world_extent=None):
    """Plot the cross section array with voxels now the legend index.
    """
    fig, ax = plt.subplots(figsize=(12, 5))
    mappable = ax.imshow(arr[::-1], origin='lower', extent=world_extent)

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
    # plt.show()
    
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
    
# --- It is crucial to get the geoCodes correct from the Dino-loket cross section image legend    
geoCodes = ['NUECga', 'NUECgb', 'NUEC1', 'NUNIHO', 'NUNIBA', 'NUBXWI-SI-KO', 'NUBX', 'NUDR', 'NUgs']

# --- It is also crucial to set proper world extent coordinates for the Dino-loket X-sec image.
world_extent=(0, 3630, -20.2, -0.80)

    
# %%
if __name__ == '__main__':
    dirs = Dirs()
    
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
    
    # --- The last legend_color is always [255., 255., 255.] to indicate empty cells later on
    # --- Therefore we add a legend index 'none' for this code.
    geoCodes.append('none')
    
    assert len(legend_colors) == len(geoCodes), (
        f"Len(colors) != len(geoCodes): {len(legend_colors)} != {len(geoCodes)}")
 
    # --- Initiate a new picker, now with the cross section   
    # --- To get the pixel bounding box
    picker = ImagePicker(geotop1)
    # --- Again, zoom in
    # --- Then pick 2 opposite corners, return to finish
    extent = picker.get_bbox(n=-1)
    print(extent)
    
    # --- Fill an array with soil indices where each index is the number of the
    # --- legend color boxes in that order (the order clicked before)
    
    # --- Instantiate the CrossSectionDigitizer using the cross section.
    digitizer = CrossSectionDigitizer(image=geotop1)
    
    # --- Set pixed extent (see extent obtained above)
    digitizer.set_section_bbox(*extent)
    
    # --- Set world extent (given above in world coordinates)
    digitizer.set_world_bbox(*world_extent)
    
    # --- Set the grid specifying the voxel width and the voxel height
    digitizer.set_grid(dx=100, dz=0.5)
    
    # --- Internally set the legend_colors (obtained from clicking the legend)
    digitizer.legend_colors = legend_colors
    
    # === Fill the grid array with the soil-indices (color box or legend item index)
    x = np.linspace(world_extent[0], world_extent[1], 201)
    z = np.linspace(world_extent[2], world_extent[3],  101)
    
    # -- Fill in an array of a cross section according to the grid object
    gr = Grid(x, None, z, axial=False)
    arr = digitizer.build_array(gr=gr)
    
    # === Fill in an array of a cross section using  dx and dz and world_extent
    arr = digitizer.build_array(gr=None)
    
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
