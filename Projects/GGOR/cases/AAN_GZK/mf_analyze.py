# -*- coding: utf-8 -*-

""" This file mf_analyze.py is used to analyze the results of the mf6 simulation. It is
run after the simulation is complete to visualize and interpret the results.

Note that mf6 does not support frf, fff and flf, not even for the structured grid.
Therefore a separate python file is used to compute those old-fashioned arrays
from the cell-by-cell flows (cbc) file, as they are very transparent and easy to understand
and used.

@ TO 2025-07
"""
import os
import sys
sys.path.insert(0, "/Users/Theo/GRWMODELS/python/mf6lab/src")
sys.path.insert(0, "/Users/Theo/GRWMODELS/python/mf6lab/projects/ggor/src")

import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import flopy

import time
from contextlib import contextmanager

import etc
import ggor_tools as ggt
from settings import sim_name, dirs, props, get_grid

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

# --------------------
# --- rest of the code
# --------------------
start_script = time.perf_counter()

with log_timed(logger, 'Load gr, tdata, parcel_data'):
    # --- unpickle tdata
    tdata_file = os.path.join(dirs.data, 'tdata.pkl')
    tdata = pd.read_pickle(tdata_file)

    # --- unpickle parcel_data
    pdata_file = os.path.join(dirs.data, 'parcel_data.pkl')
    parcel_data = pd.read_pickle(pdata_file)

    # --- regenerate the grid
    gr = get_grid(parcel_data=parcel_data, dx=props['dx'])

with log_timed(logger, "loading MFsimulation"):
    sim = flopy.mf6.MFSimulation.load(sim_name=sim_name,
                                  version='mf6',
                                  sim_ws=dirs.SIM,
                                  lazy_io=True)

with log_timed(logger, "sim.get_model"):
    gwf = sim.get_model('{}Gwf'.format(sim.name).lower()) # list(sim.model_names)[0])

with log_timed(logger, "heads_obj loaded"):
    heads_obj = ggt.Heads_obj(sim=sim, tdata=tdata, gr=gr)

with log_timed(logger, "Watbal_obj loaded"):
    watbal = ggt.Watbal_obj(sim=sim, dirs=dirs, gr=gr)

parcels = [0]
with log_timed(logger, f"Plotting parcels [{', '.join([str(p) for p in parcels])}]"):    
    title=['Parcel heads']
    ax = heads_obj.plot(tdata=tdata,
            parcel_data=parcel_data,
            parcels=parcels,
            plotGXG=True,
            figsize=(14, 8))

    ax = watbal.plot(parcel_data=parcel_data,
                        tdata=tdata,
                        parcels=parcels,   # over all parcels
                        sharey=True)
    plt.gcf().savefig(dirs.SIM + '/watbal.png', dpi=300)
    
with log_timed(logger, "Plotting Watbal"):
    ax = watbal.plot(parcel_data=parcel_data,
                        tdata=tdata,
                        parcels=None,   # over all parcels
                        sharey=True)
    plt.gcf().savefig(dirs.SIM + '/watbal_all_parcels.png', dpi=300)

with log_timed(logger, "GXG added, pickled and plotted"):
    parcel_data['GHG'] = heads_obj.GXG['GHG']
    parcel_data['GVG'] = heads_obj.GXG['GVG']
    parcel_data['GLG'] = heads_obj.GXG['GLG']
    
    # --- overwite old parcel_data.pkl file
    parcel_data.to_pickle(pdata_file)

    # --- map GHG
    for what in ['GHG', 'GVG', 'GLG']:
        ax = etc.newfig(f"{what} Gooise Aanzichts Kade", "xRD", "yRD",
                        figsize=(10, 8))
        parcel_data.plot(what,
                     cmap='viridis',
                     legend=True,
                     edgecolor='black',
                     figsize=(10, 8),
                     ax=ax)
    plt.gcf().savefig(dirs.SIM + f"/map_{what}.png", dpi=300)
                     

plt.show()

# --- total elapsed time in script
logger.info(f"mf_analyse finished in {time.perf_counter()-start_script:.2f} seconds")
