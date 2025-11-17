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
import mf6_bootstrap # noqa: F401
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import flopy

import time
from timing import log_timed

import etc
import ggor_tools as ggt
from settings import case_name, dirs, props, get_grid

# --- setting up the logger
import logging
import logging_setup  # noqa: F401
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

    
# --------------------
# --- rest of the code
# --------------------
start_script = time.perf_counter()

with log_timed(logger, 'Load gr, tdata, parcel_data'):
    # --- get saved tdata
    tdata_file = os.path.join(dirs.data, 'tdata.pkl')
    tdata = pd.read_pickle(tdata_file)

    # --- get saved parcel_data
    pdata_file = os.path.join(dirs.data, 'parcel_data.pkl')
    parcel_data = pd.read_pickle(pdata_file)

    # --- regenerate the grid object
    gr = get_grid(parcel_data=parcel_data, dx=props['dx'])

# --- Read the modflow simulation
with log_timed(logger, "loading MFsimulation"):
    sim = flopy.mf6.MFSimulation.load(sim_name=case_name,
                                  version='mf6',
                                  sim_ws=dirs.SIM,
                                  lazy_io=True)

# --- read the groundwater flow model
# with log_timed(logger, "sim.get_model"):
#    gwf = sim.get_model('{}Gwf'.format(sim.name).lower()) # list(sim.model_names)[0])

# --- Load the heads from Modflow and plot a parcel
with log_timed(logger, "heads_obj loaded"):
    heads_obj = ggt.Heads_obj(sim=sim, tdata=tdata, gr=gr)
    
    # --- plot time line for parcels
    # --- choose 1 or 2 ok, otherwise timeline plot becomes messy
    parcels = [0]
    title=['Parcel heads']
    axs = heads_obj.plot(tdata=tdata,                      
            parcel_data=parcel_data,
            parcels=parcels,
            plotGXG=True,
            figsize=(14, 8))
    axs[0].figure.suptitle(case_name)


# --- loading Modflows CBC output is expensive (takes about 50 seconds)
with log_timed(logger, "Watbal_obj loaded"):
    
    # --- Loading water budget components for all cells from CBC file
    print("Loading Watbal ... may take several minutes ...")
    
    watbal = ggt.Watbal_obj(sim=sim, dirs=dirs, gr=gr)

    # --- plot running water budget for all parcels
    axs = watbal.plot(parcel_data=parcel_data,
                        tdata=tdata,
                        parcels=None,   # over all parcels
                        sharey=True)
    axs[0].figure.suptitle(case_name)
    plt.gcf().savefig(dirs.images + '/watbal_all_parcels.png', dpi=300)

    # --- plot running water budget for all selected parcels
    axs = watbal.plot(parcel_data=parcel_data,
                        tdata=tdata,
                        parcels=parcels,   # id's of selected parcels
                        sharey=True)
    axs[0].figure.suptitle(case_name)
    plt.gcf().savefig(dirs.images + '/watbal.png', dpi=300)
    
# --- Plot 3 maps of the parcels GXG
with log_timed(logger, "GXG added, pickled and plotted"):

    vmin, vmax = np.inf, -np.inf
    for gxg in ['GHG', 'GVG', 'GLG']:
        parcel_data[gxg] = heads_obj.GXG[gxg]
        vmin = np.fmin(vmin, parcel_data[gxg].min())
        vmax = np.fmax(vmax, parcel_data[gxg].max())
        
    # --- overwite old parcel_data.pkl file
    parcel_data.to_pickle(pdata_file)

    # --- map GHG, GVG, GLG
    for what in ['GHG', 'GVG', 'GLG']:
        ax = etc.newfig(f"{what} {case_name}", "xRD", "yRD",
                        figsize=(10, 8))
        # --- plot map of G?G of parcels
        parcel_data.plot(what,
                     cmap='viridis',
                     vmin=vmin,
                     vmax=vmax,
                     legend=True,
                     edgecolor='black',
                     figsize=(10, 8),
                     ax=ax)
    ax.figure.suptitle(case_name)
    plt.gcf().savefig(dirs.images + f"/map_{what}.png", dpi=300)
                     

plt.show()

# --- total elapsed time in script
logger.info(f"mf_analyse finished in {time.perf_counter()-start_script:.2f} seconds")
