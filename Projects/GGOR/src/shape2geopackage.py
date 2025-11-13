"""
ggor main shapefile / database
===============

Support for reading and showing the basic ggor database provided as a shapefile

This file should be run only once (and has been bun already). Adapt for new ggor cases.

Main Components
---------------
- Read the shapefile into a GeoDataFrame
- Inspect the columns
- Save it as a geopackage for later direct retrieval
- Visualize it.

Example
-------
>>> import mf6lab.regional as reg
>>> reg.setup_model("case01")
>>> reg.run_model("case01")
>>> reg.plot_results("case01")
"""
import os
import geopandas as gpd
from pathlib import Path
import matplotlib.pyplot as plt
import etc

import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S"
)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
DEFAULT_GRID_SIZE = 100
CASES_DIR = Path(__file__).parent.parent  / "cases"

# -----
if __name__ == '__main__':
    
    logging.info("Starting process...")
    
    cases = ['AAN_GZK', 'Noorderpark']
    crs = "EPSG:28992"
    
    for case in cases:
        # --- make sure directories and files exist
        logging.info(f"handling {case}")

        shpfile = os.path.join(CASES_DIR, case, 'data', case + '.shp')
        assert os.path.isfile(shpfile), f"No such file {shpfile}"
        outdir = os.path.join(CASES_DIR, case, 'data')
        assert os.path.isdir(outdir), f"No such directory {outdir}"
    
        # --- get the shape file and add missing crs
        dbase = gpd.read_file(shpfile, engine='fiona')
        dbase.crs = crs

        # --- prepare and save GeoDataFrame to geopackage
        gpkg_file = os.path.join(outdir, case + '.gpkg')
        dbase.to_file(gpkg_file, driver='GPKG')
        
        logging.info(f"{case} saved to {gpkg_file}")
        
        ax = etc.newfig(f"Percelen {case}, crs={crs} Amersfoort / RD New",
                        "xRD [m]", "yRD [m]", figsize=(8, 7))
        ax.set_aspect(1)
        dbase.plot(ax=ax)
        
        outfile = os.path.join(CASES_DIR, case, 'images', case + '.png')
        ax.figure.savefig(outfile)
        
        logging.info(f"{case} figure saved to {outfile}")
                
    
    logging.debug("Shapefiles succesfully converted to geopackages")
    print("Quit figures to finish.")
    plt.show()
    
