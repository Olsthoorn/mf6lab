# -*- coding: utf-8 -*-
"""_summary_

Settings (this file) can be run on its own, and will be imported
by mf_adapt to start running Modflow.

It specifies:
    the directory structure.
    properties of the simulation
    model properties not already in the origonal data (shpaefile)

"""

# %% --- Imports
import os
import sys
from pathlib import Path

import time
from contextlib import contextmanager

import pandas as pd
import KNMI  # To convert KNMI data to a pandas DataFrame
from etc import newfigs

import ggor_tools as ggt
from mf6tools import Dirs  # Name space for this project

# --- setting up the logger
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S"
    )

# --- Setting up timing
@contextmanager
def log_timed(logger, msg):
    start = time.perf_counter()
    yield
    logger.info(f"{msg} in {time.perf_counter() - start:.2f} seconds")

# --- make sure ggor_src is in sys.path to find ggor_tools
mf6_src  = str(Path.cwd() / '..' / '..' / 'src')
ggor_src = str(Path.cwd() / 'src')

for path in (ggor_src, mf6_src):
    if path not in sys.path:
        sys.path.insert(0, path)
        logger.info(f"Added to sys.path: {path}")


print("\nSys.path:\n--------")
for p in sys.path:
    print(p)
print('--------')

start_script = time.perf_counter()

# %% --- Set the home directory and name space for this project
HOME = '/Users/Theo/GRWMODELS/python/mf6lab/Projects/GGOR/' # /cases/case
    
assert os.path.isdir(HOME), "Can't find the directory {}".format(HOME)

sim_name = 'AAN_GZK'

# Generate the name space for the directories of this project
dirs = Dirs(HOME)
dirs.add_case(sim_name) # Add case to this directory namespace
dirs.add_to_home('src')
dirs.add_to_home('data/meteo')
dirs.add_to_home('data/bofek')

# Move to this case directory
os.chdir(dirs.case)
logger.info("cwd = {}".format(os.getcwd()))

# %% [markdown]
# # Properties for the simulation
# These properties are used in the GGOR_modflow class and in the GGOR_data class
# and in the GGOR_modflow class.
# They are also used in the GGOR_modflow class to create the MODFLOW 6 input files.
# The properties are used to create the MODFLOW 6 input files and to run the simulation.
# The properties are also used to create the GGOR_modflow object and to run the simulation.
# The properties are also used to create the GGOR_data object and to run the simulation.
# %%
props = {
        'nper_test': 1000, # max number of stress periods when test is True
        'test': False,
        'length_units': 'meters',
        'time_units': 'days',
        'use_w_not_c': False,     # Use anal. form. instead of entry resistance and circumference.
        'start_date_time': '2024-01-27',
         'oc_frequency': 1,
         'icelltype': 1,
         'dx':   1.0,    # [m] cell width
         'minDz': 0.01, # m min layer thickness (also for pinched-out layers)
         'drain_depth': 0.15, # m
         'cDrainage':  100.0, # d
         'rch':        0.001, # m/d
         'strthd':       0.0, # m initial head
}

def get_tdata(stn=240, start='20100101', end='20191231', folder=None):
    """Return the meteo data, handle it and return it as a pd.DataFrame
    
    Hanling implies adding columns "hyear" and "hand
    
    Parameters
    ----------
    stn: int
        KNMI station (240 = De Bilt)
    start, end: 'yyyymmdd'
        start and end of the index of tdata
    folder: path
        folder with the KNMI meteo data (dirs.data.meteo)
    """
    meteo_data = KNMI.knmi.get_weather(
            stn=240,  # KNMI station number for De Bilt
            start=start,
            end=end,
            folder=folder # dirs.data_meteo
    )
    tdata = ggt.handle_meteo_data(meteo_data,
                                  summer_start=4,
                                  summer_end=10)
    return tdata

def get_parcel_data(defaults=None, BMINMAX=(5, 250)):
    """Return the geopandas.GeoDataFrame with the parcel data.
    
    Bofek data, converts from Bofek codes to
    the soil properties kh, kv and sy.
    The BOFEK column represents a Dutch standardized soil type.
    The corresponding values for 'kh', 'kv' and 'Sy' are currently read from an Excel worksheet into a pandas DataFrame (thus becoming a table)
    
    Parameters
    ----------
    defaults: ggr.defaults a dict
        defaults from ggor_tools
    BMINMAX: (float, float)
        minimum and maximum values of parcel half-width
    """    
    bofek = pd.read_excel(
        os.path.join(dirs.data_bofek, "BOFEK eenheden.xlsx"),
        sheet_name = 'bofek',
        index_col=0,
        engine="openpyxl")
    
    # --- Create a GGOR_modflow object and
    # --- get the upgraded parcel_data from it,
    # --- excluding parcels that are too small BMIN or too wide BMAX    
    parcel_data = ggt.GGOR_data(dirs=dirs,
                            defaults=defaults, # ggt.defaults,
                            bofek=bofek,
                            BMINMAX=BMINMAX,
                            ).data
    return parcel_data

def get_grid(parcel_data=None, dx=None):
    gr = ggt.grid_from_parcel_data(parcel_data=parcel_data,
                                   dx=props['dx'])
    return gr


# %% 
if __name__ == '__main__':
    pass
    
    # %%    
    xlims, ylims = None, None
        
    axs = newfigs(titles=['title1', 'title2'],
                xlabel='time',
                ylabels=['heads [m]', 'flows [m2/d]'],
                xscale='linear',
                yscales=['linear', 'linear'],
                sharex=True,
                sharey=False,
                xlim=None,
                ylims=None,
                figsize=(12, 6),
    )

    test=False
    
    tdata = get_tdata(stn=240, start='20100101', end='20191231',
                      folder=dirs.data_meteo)
    
    if props['test']:
        with log_timed(logger, "Test data generated"):
            # --- Limits len(tdata) for testing          
            tdata = tdata.iloc[:props['nper_test']]
            
            # --- convert the actual data into test data
            tdata = ggt.gen_testdata(tdata=tdata,
                                        RH  =(270, 0.0, 0.002, 0.004),
                                        EV24=(180, 0.0, 0.001, 0.002),
                                )
            
            # --- get parcel data for testing
            parcel_data = ggt.get_test_parcels(os.path.join(
                                    dirs.case, 'pdata_test.xlsx'), 'parcel_tests1')
            # --- Special test
            parcel_data = parcel_data.iloc[0:4]
    else:
        with log_timed(logger, "Parcel data obtaind"):
            parcel_data = get_parcel_data(defaults=ggt.defaults,
                                          BMINMAX=(5, 250))    

    logger.info(f"Settings run in {time.perf_counter() - start_script:.2f} seconds")
        
# %%
