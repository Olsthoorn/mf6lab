# -*- coding: utf-8 -*-

# Note that the PYTHONPATH=/Users/Theo/GRWMODELS/python/tools/ was set in ~/zshrc to allow including modules in tools.
# Use VENV flopy by setting it


"""Simultate GGOR for nparcels using MODFLOW.

The parcel (area) data are in a shape file.
The attributes are in the dbf
The .dbd file is read into pd.DataFrame.
The meteo is reada from existing file or obtained from KNMI site.
Modflow is used to simultaneously simulate all the parcels dynamically.
The results are shown for selected paracels (hds, GXG)
The running water budget is shown for all the parcels combined.

Different scenarios can be dealt with as cases. Scenarios are used simulation of
a regular case or for testing the behavior of the model in given test-circumstances.
A scenario of testing has (for instance) max 5 parcels. The reason is that the number of head
time series plotted is limited to 5. Hence a DataFrame needs to
be generated from the data that define the properties os these parcels. This DataFrame
can be read from an Excel workbook with the name of the case. This gives full control
over the parcel properties of the test case.
The next data for a case is the meteo information. This to can be read from
an excel sheet from the same workbook. The required fields are in the
example workbook named "test_basic.xlsx".

@ TO 2020-09-06

TO 2023-08-18, system was rerun using mf2005.
Next, needs verification and specifying which answers are to be obtained.
Question whether to converto to mf6 with strucured or unstructured grid.
Need demo using the portion if __name__ == __main__: converted to Jupyter for ease of access and understanding.
There each step should be demonstrated and verified, graphically if possible.

"""
import os
import sys
import time
from contextlib import contextmanager
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
from KNMI import knmi
from fdm.src.mfgrid import Grid
from fdm.src.mf6_face_flows import get_structured_flows_as_dict
import etc

import mf6tools as mf6tools

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

NOT = np.logical_not
AND = np.logical_and

#% Names in the databbase and their replacements used in GGOR
# Names not mentioned are not replaced.
colDict = {
        'Bofek'     : 'bofek',
        'Greppels'  : 'n_trench',
        'Med_Cdek'  : 'c_CB',
        'Med_Ddek'  : 'D1',
        'Med_Kwel'  : 'q_up',
        'Med_Phi2'  : 'phi',
        'Med_mAHN3' : 'AHN',
        'Shape_Area': 'A_parcel',
        'Shape_Leng': 'O_parcel',
        'Winterpeil': 'h_winter',
        'Zomerpeil' : 'h_summer',
}

#% Defaults used for required parameters that, however, are not in the databse.
# They will be used when not in the database.
defaults = {
        'd_drain': 0,    # [m] Tile drai depth below ground ele.; 0 where no drains.
        'd_trench': 0.3, # [m] Trench depth in case present.
        'c_drain': 5.,   # [d] Area tile drainage resistance. Also used for trenches.
        'wi_ditch' : 2., # [d] Ditch resist, used for flow from ditch to ground. (analytical)
        'wi_ditch2' : 2., # [d] Ditch resis. reg. aquif. for flow to ground. (analytical)
        'wo_ditch' : 1.,  # [d] Ditch resis. when flow is to ditch. (anlaytical)
        'wo_ditch2' : 1., # [d] Ditch resis. reg. aquif. for flow to ditch. (analytical)
        'ci_ditch' : 2.,  # [d] Ditch bottom and side entry resistance (applied to Omega)
        'co_ditch' : 1.,  # [d] Ditch bottom and side entry resistance (applied to Omega)
        'd_ditch' : 1.0,  # [m] depth of ditch below ground surface
        'b_ditch' : 0.75, # [m] half-width of the ditch
        'D_CB' : 0.1,     # [m] (dummy thickness) van basisveenlaag (CB=confining bed)
        'D2' : 40.,       # [m] thickness of regional aquifer
        'S2' : 1e-3,      # [-] total (elastic) storage coeff. of regional aquifer
        'kh2': 30.,       # [m/d]  horizontal condutivity of regional aquifer
        'kv2':  6.,       #[ m/d]  vertical conductivity of regional aquifer
        'ET_surfd': 1.0,  # [m] depth of surf in ET below ground surface.
        'ET_exdp': 2.5,  # [m] Modflow's extinction depth (see ET package)
}


# --- Modflow cell-by-cell flow labels:
# --- Translates short labels to those in the CBC file.
cbc_labels = {
        'STO': 'STORAGE',
        'FLF': 'FLOW LOWER FACE ',
        'WEL': 'WELLS',         # used for seepage. Are in second layer.
        'EVT': 'ET',
        'GHB': 'HEAD DEP BOUNDS',
        'RIV': 'RIVER LEAKAGE',
        'DRN': 'DRAINS',        # to capture surf. runoff, or actual drains and trenches.
        'RCH': 'RECHARGE',
}

# --- For legend when of running water budget plot
# --- leg is legend for this label in the graph
# --- clr is the color of the filled graph
watbal_label = {
        'RCH': {'leg': 'RCH', 'clr': 'green'},
        'EVT': {'leg': 'EVT', 'clr': 'gold'},
        'WEL': {'leg': 'WEL(in wvp2)' , 'clr': 'blue'},
        'DRN': {'leg': 'DRN(drn/trenches/runoff)', 'clr': 'lavender'},
        'RIV': {'leg': 'RIV(ditch out)', 'clr': 'magenta'},
        'GHB': {'leg': 'GHB(ditch in+out)', 'clr': 'indigo'},
        'FLF': {'leg': 'FLF(leakage)', 'clr': 'gray'},
        'STO': {'leg': 'STO', 'clr': 'cyan'}
        }

def selection_check(parcels=None, n=None):
    """Returns selected parcels after verification of the selection.
    
    Parameters
    ----------
    parcels: sequence of slice telling which parcels to pick
        selection of parcels from the database (DataFrame)
    n: int or None
        Maximum number of parcels in case parcels is a slice
    
    Returns
    -------
    sequence of ints telling which parcels to select from the DataFrame.
    Slice is turned into a sequence, which requires n as int <= number of parcels
    in the parcel DataFrame.
    """
    if parcels is None:
        parcels = slice(0, n, 1)
    elif isinstance(parcels, int):
        parcels = slice(0, parcels, 1)
    elif isinstance(parcels, (tuple, list, np.ndarray, range, slice)):
        pass
    else:
        raise ValueError("parcels must be None (for all), an int or a sequence")
    if isinstance(parcels, slice):
        pass
    else:
        parcels = np.array(parcels)
    return parcels


def gen_testdata(tdata, **kwargs):
    """Return copy of tdata with altered or added input columns for testing.

    The tdata tuples consist of:
        (nday, value1, value2, value3 ..)
        The nday is the duration of each succcessive value.
        The values will be implied one after the other, each during the
        nday peirod. After the last value is used, the sequence is repeated.
        Therefore, the number of values is immaterial.
        Note that columns are replaced or added to the tdata DataFrame.
        A new DataFrame is returnd. Theo old one is left untouched.

    Parameters
    ----------
    tdata: pd.DataFrame witih datetime index
        input for modelling
    kwargs: a dict of tuples
        each kwarg has a key and a tuple as argument. The tuple consists
        of a number_of_days followd by values. Like (120, 0.5, 0.7, -03).
        The first value is the duration during which each of the following
        values will be used in turn, repeated after the last value has been
        consumed.
        Due to how Python works, the call can be like so
        gen_testdata(tdata, RH=(200, 0.002, 0.003, 0.02), hLR=(150, 0.7, -0.3),
                     h1=(300, -0.4, 0.2, -0.1))
        The arguments will then be available in kwargs like so
            kwargs={'RH': (200, 0.002, 0.003, 0.02),
                    'hLR': (150, 0.7, -0.3),
                    'h1': (300, -0.4, 0.2, 0.02)}
        This makes this function extremely flexable to generate testdata.
    Example
    -------
    tdatanew = gen_testdata(tdata,
                        RH=(200, 0.02, 0.01, 0.03, 0.01)),
                        EV24=(365, 0., -0.001)
                        q_up=(150, -0.001, 0.001, 0, -0.003)
                       )
    tdatanew is then a pd.DataFrame with columns 'RH', 'EV24', 'q_up', next
        to other columns that were already in tdata. Existing columns may be
        oerwritten if they already exist.
    """
    tdata = tdata.copy() # leave tdata intact
    for key in kwargs:
        # --- index array telling which of the tuple values to pick
        # --- Day number since start of tdata.index
        daynum = (tdata.index - tdata.index[0]) / np.timedelta64(1, 'D')
        period = int(kwargs[key][0]) # days
        values = np.array(kwargs[key][1:])
        Id = np.asarray((daynum // period) % len(values), dtype=int)
        
        # --- Add or replace column in tdata copy
        tdata[key] = np.array(values[Id])
    return tdata


def handle_meteo_data(meteo=None, summer_start=4, summer_end=9):
    """Add and fill columns []'summer', 'hyear' 'hand'].

    Added columns are
        summer: bool
        hyear: hydrological year. They atart at April 1 and end March 31 to get proper GVG. Do not change this !
            GVG will be the mean of the values on 14/3, 28/3 and 14/4 of each hydological year'
        hand: official groundwater measurement dates: 14th and 28th of every month)

    Parameters
    ----------
    meteo: pd.DataFrame
        meteo data, with timestamp as index and columns 'RH' and 'EVT24'
        for precipitation and (Makkink) evaporation respectively in m/d.
    summer_start: int
        month coniciding with start of hydological summer. Default 4 (April).
    summer_end: int
        month coinciding with end of hydrologiacl summer. Default 10 (October = start hydrological winter).
    """
     # --- Check for existance of both columns in meteo.columns.
    dcol = {'RH', 'EVT24'}.difference(meteo.columns)
    
    if dcol:
        KeyError("Missing column[s] [{}] in meteo DataFrame".format(', '.join(dcol)))

    # --- verify, data are in m/d
    if not meteo['RH'].median() < 0.01 and meteo['EVT24'].median() < 0.01:
        AssertionError("Median of Precipitration = {:5g} and median evapotranspiration = {:4g}\n"
                          .format(meteo['RH'].median(), meteo['EVT24'].median()) +
                       "Percipication and or evapotranspiration likely not in m/d!")

    # --- Add boolean column indicating summer (needed to set summer and winter ditch levels)
    meteo.loc[:, 'summer'] = [True if t.month in range(summer_start, summer_end) 
                               else False
                               for t in meteo.index]

    # --- hydrological year column 'hyear'
    hyear_start_month = 3   # Don't change! It's needed in the GXG class
    hyear_start_day   = 14  # Don't change! It's needed in the GXG class
    meteo.loc[:, 'hyear'] = [t.year
        if t.month >= hyear_start_month and t.day >= hyear_start_day
        else t.year - 1 for t in meteo.index]

    meteo.loc[:, 'hand'] = [t.day % 14 == 0 for t in meteo.index]

    return meteo


def grid_from_parcel_data(parcel_data=None, dx=None, laycbd=(1, 0)):
    """Get gridobject for GGOR simulation.

    Parameters
    ----------
    parcel_data: pd.DataFrame
        table with all parcel property data
    dx: float
        cell width (all cells have the same width in GGOR, e.g. 1 m)
    LAYCBD: list (or sequence) of length 2 (nlayer)
        value of 1 for each layer with a confining unit below it, else value of 0
    """
    bmax = parcel_data['b'].max() # width on widest parcel in the database
    nx = int(np.ceil(bmax) // dx) # width of the grid in cells
    ny = len(parcel_data)
    nz = 2 # GGOR always has a cover layer and a regional layer with a  confining bed in between.

    Dx  = np.ones(nx) * dx
    xGr = np.hstack((0, np.cumsum(Dx)))
    yGr = np.arange(ny, -1, -1, dtype=float)

    Z = np.zeros((nz + 1, ny, nx))

    Z[0] = parcel_data['AHN'].values[:, np.newaxis] * np.ones((1, nx))
    Z[1] = Z[0] - parcel_data['D1'  ].values[:, np.newaxis] * np.ones((1, nx))
    Z[2] = Z[1] - parcel_data['D2'  ].values[:, np.newaxis] * np.ones((1, nx))

    return Grid(xGr, yGr, Z)


def set3D(layvals, shape=None):
    """Return array of the desired shape from parcel_data column[s].
    
    Use this to generate arrays such as the initial heads and conductivities from the
    given parcel values.

    The shape of the output array can be (nz, ny, nx), or (nz, ncpl).
    layvals is a pd.DataFrame of a pd.Series object. If it is a series or a DataFrame with
    a single column, then all layers will have the same values.
    Otherwise, the number of layval columns must equal the number of layeres, i.e. shape[0].
    
    The rows in the layvals correspond to the parcels. Therefore, the row values are
    generally different. The columns refer to the layers.
    
    The output array will get the row values for all cells in each layer of the rows.
    In the structured Modflow grid, each row has ncol layers.
    
    In our case the rows of the datafram correspond to irow in the structured Modflow grid.
    
    The columns will be the layer values in sequence. The index
    of the layvals DataFrame is the row number of the grid. The first column fills
    the rows of layer iz=0 (all ix values will be the same (uniform layer)).
    The second column of the DataFrame will fill the rows of the layer 1, etc.
    If layvals has one column or is a pd.Series, the 3D array with shape will have
    the same values in each layer, such as for STRTHD.
    

    Parameters
    ----------
    layvals: pd.DataFrame or pd.Series
        The model layer values to be used.
        The number of columns in the DataFrame must be 1 or equal to the number of model layeres shape[0].
    shape: tuple (nz, ny, nx) or (nz, ncpl)
        the shape of the array to be generated.
    """
    ny = len(layvals.index)
    nz = shape[0]
    nColumns = 1 if isinstance(layvals, pd.Series) else layvals.shape[1]
    assert nColumns == 1 or nColumns == nz,\
        f'The number of columns in layvals=({nColumns}) but must be 1 or equal to nlay=({nz}).'
            
    ncpl = np.prod(shape[1:] if len(shape) == 3 else shape[-1])
    nx   = shape[-1] if len(shape) == 3 else ncpl // ny
        
    vals = layvals.values
    if nColumns == 1:
        vals = vals[np.newaxis, :, np.newaxis] * np.ones((nz, ny, nx))
    else:        
        vals = vals.T[:,        :, np.newaxis] * np.ones((nz, ny, nx)) 
        
    return vals if len(shape)==3 else vals.reshape((nz, ny, nx))


def get_drain_elev_with_trenches(pdata=None, gr=None, d_drn=None):
    """Retrun drain elevations corrected for trenches.

    Parcels with 1 or more trenches are assumed to have no tile drainage.
    Their drainage depth will be set to zero except for the trench locations,
    where the drainage depth is set to 'dtrench'.

    Paramters
    ---------
    parcel_data: pd.DataFrame
        parcel properties
    gr: mfgrid.Grid object
        object holding the mesh information
    Returns
    -------
    elev: np.ndarray of shape (nparcel, ncol)
        Updated drainage elevations. I.e. ground surface in parcels that have
        one or more trenches except on the trench locations, where the
        elevation is set to groundsurface minus the trench depth.
        The trench locations follows from the number of trenches in the parcel ntr.
    """
    
    elev = ((pdata['AHN'] - pdata['d_drain']).values[:, np.newaxis]
                                                    * np.ones((1, gr.nx)))
    Id = np.arange(gr.nx, dtype=int)
    xLeft = gr.X[0][:,:-1] # left side of cells in top layer for all rows)
    
    for iparcel, (b, ntr, ahn, dtr) in enumerate(
            zip(pdata['b'], pdata['n_trench'], pdata['AHN'], pdata['d_trench'])):
        
        # --- default drain depth (simulates surface runoff)
        elev[iparcel] = ahn - d_drn
        
        if ntr:
            # --- only if parcel has one or more trenches
            atr = b / ntr # distance between trenches or trench and ditch
            xtr = np.cumsum(np.ones(ntr) * atr) # locations of trenches
            
            # --- cell row indices of trenches
            Itr = [Id[xLeft[iparcel] <= _x][-1] for _x in xtr]
            elev[iparcel][Itr] = ahn - dtr
            
    return elev # Ny * Nx array

    for iparcel, (b, ntr, ahn, dtr) in enumerate(
            zip(pdata['b'], pdata['n_trench'], pdata['AHN'], pdata['d_trench'])):
        print(b, ntr, ahn, dtr)


def get_cond_DRN(pdata=None, gr=None):
    """Return conductance for use by DRN (Surface Runoff).
    
    Parameters
    ==========
    pdata: dict
        parcel properties pd.DataFrame
    """
    c_drain = pdata['c_drain'].values[:, np.newaxis] * np.ones((1, gr.nx))
    condDRN = gr.Area / c_drain
    return condDRN # Conductance for all cells in top layer as an Ny * Nx array


def get_cond_GHB(pdata=None, gr=None, use_w_not_c=None):
    """Retrun conductance for use by GHB for two layers first column.
    
    Parameters
    ==========
    pdata: pd.DataFrame
        parcel properties
    gr: structured Grid object
        Class that holds the grid
    """
    if use_w_not_c:
        wi = np.vstack((pdata['wi_ditch'], pdata['wi_ditch2']))
        
    else: # Use the real ditch resistance, dicth circumference
        # TODO use either two or three layer grid it's now inconsistent
        wi = np.vstack((pdata['ci_ditch'] * pdata['D1'] / pdata['ditch_omega1'],
                        pdata['ci_ditch'] * pdata['D2'] / pdata['ditch_omega2']))
        # extra resistance due to partial penetration of ditch
        wi += np.vstack((pdata['wpp1'], pdata['wpp2']))
        
    dy    = np.vstack((gr.Dy[:, 0], gr.Dy[:, 0]))
    cond  = np.vstack((pdata['D1'], pdata['D2'])) / wi * dy
    cond[np.isnan(cond)] = 0. # This is where ditch_omega is zero, hence no ditch present (2nd layer)
    return cond


def get_RIV_Cond(pdata=None, gr=None, use_w_not_c=None):
    """Return the Conductances of ditches for use by RIV package.
    
    parameters
    ==========
    pdata: dict
        parameter data
    """
    if use_w_not_c:
        # Use analytic ditch resistance with layer thickness and no partial penetration
        wi = np.vstack((pdata['wi_ditch'], pdata['wi_ditch2'])) # Same value for both layers
        wo = np.vstack((pdata['wo_ditch'], pdata['wo_ditch2'])) # Same value for both layers
        assert np.all(wi >= wo), "ditch entry resist. wi must be >= to the ditch exit resist. wo!"
        dw = wi - wo; eps=1e-10; dw[dw==0] = eps # prevent (handle) division by zero
        w     = (wo * wi / dw)
        
    else: # Use real ditch resistance with ditch circumference
        ci = np.vstack((pdata['ci_ditch'], pdata['ci_ditch'])) # Same value for both layers
        co = np.vstack((pdata['co_ditch'], pdata['co_ditch'])) # Same value for both layers
        assert np.all(ci >= co), "ditch entry resist. must be >= the ditch exit resist!"
        dc = ci - co; eps=1e-10; dc[dc==0] = eps # prevent (handle) division by zero
        c     = (co * ci / dc)
        # To analytic resistance, using the ditch circumference
        w = c * np.vstack((pdata['D1'] / pdata['ditch_omega1'],
                        pdata['D2'] / pdata['ditch_omega2']))
        # Add partial penetration to resistance
        w += np.vstack((pdata['wpp1'], pdata['wpp2']))
        
    dy    = np.vstack((gr.Dy[:, 0], gr.Dy[:, 0]))
    cond  = np.vstack((pdata['D1'], pdata['D2'])) / w * dy
    cond[np.isnan(cond)] = 0. # When ditch_omega is zero
    return cond


class GGOR_data:
    """Cleaned parcel data object. Only its self.parcel_data will be used (pd.DataFrame)."""

    def __init__(self, dirs=None, bofek=None, BMINMAX=(5., 500.), defaults=None):
        """Get parcel data for use by GGOR.

        Parameters
        ----------
        dirs: mf6tools.Dirs object
            holds paths to current case directories
        bofek: pd.DataFrame
            bofek values for ['kh', 'Sy', 'staring', 'ksat_cmpd'], the index
            is the bofek_id number.
            bofek are Dutch standardized soil parameters, see Internet.
        BMINMAX: tuple of 2 floats
            min and max halfwidth value for parcels to be considered.
        """

        sim_name = os.path.basename(dirs.case)
        
        # --- read dbf file into pd.DataFrame
        
        # --- Obsolete read the .dbf file of the shape file for the case
        # self.data = data_from_dbffile(os.path.join(dirs.data, sim_name + '.dbf'))
        
        # --- Use geopandas read the .pgkg file directly
        case_gpkg_file = os.path.join(dirs.data, sim_name + '.gpkg')
        self.data = gpd.read_file(case_gpkg_file)

        # replace column names to more generic ones
        self.data.columns = [colDict[h] if h in colDict else h
                                         for h in self.data.columns]

        # --- compute parcel width to use in GGOR
        self.compute_parcel_width(BMINMAX=BMINMAX)

        # --- set kh, kv and Sy from bofek
        self.apply_bofek(bofek) # bofek is one of the kwargs a pd.DataFrame

        # --- add required parameters if not in dbf
        self.apply_defaults(defaults)

        self.compute_and_set_omega()

        self.compute_and_set_wpp()


    def compute_and_set_omega(self):
        """Compute and set the half wetted ditch circumference in the two model layers.

        ditch_omega1 [m] is half the width of the ditch plus its wetted sided.
        ditch_omega2 [m] ia the same for the underlying regional aquifer.

        calls normal function to allow using it with test data
        """
        compute_and_set_omega(self.data)


    def compute_and_set_wpp(self):
        """Compute and set the the extra resistance due to parital ditch penetration.

        dwpp1 is the extra resistance in [d] for the top layer.
        dwpp2 is the extra resistance in [d] for the regional aquifer

        calls normal function to allow using it with test data.
        """
        compute_and_set_wpp(self.data)


    def compute_parcel_width(self, BMINMAX=(5., 10000.)):
        """Add computed parcel width to the to dataFrame self.data.
        
        The parcel width is computed fromt the parcel area A and parcel cirumference
        present in in the original database. It's the longest axis of a rectangle area
        deduced from the given area and circumference of each parcel.

        Parameters
        ----------
        BMIN : float
            minimum parcel width
        BMAX : float
            maximum parcel width.
            The final model width will be that of the
            widest parcel in the database, not the initial BMAX.
        """
        A     = np.asarray(self.data['A_parcel'])  # Parcel area
        Om    = np.asarray(self.data['O_parcel'])  # Parcel cifcumference
        det   = Om ** 2 - 16 * A            # determinant
        L     = det>=0                    # determinant>0 ? --> real solution
        PW    = np.nan * np.zeros_like(L)   # init paracel width
        PL    = np.nan * np.zeros_like(L)   # init parcel length
        PW[ L] = (Om[L] - np.sqrt(det[L]))/4  # width, smallest of the two values
        PL[ L] = (Om[L] + np.sqrt(det[L]))/4  # length, largest of the two values
        PW[NOT(L)] = np.sqrt(A[NOT(L)])      # if no real solution --> assume square
        PL[NOT(L)] = np.sqrt(A[NOT(L)])      # same, for both width and length

        PW = np.fmin(PW, max(BMINMAX)) # Arbitrarily limit the width of any parcel to BMAX.
        PW = np.fmax(PW, min(BMINMAX))

        # --- Add column 'b' to Data holding half the parcel widths.
        self.data['b'] = PW/2

        # Use only the parcles that have with > BMIN and that have bofek data
        L=np.where(AND(PW > min(BMINMAX), NOT(self.data['bofek']==0)))

        self.data = self.data.iloc[L]

        # --- Any data left?
        assert len(L) > 0, "Cleaned parcel database has length 0, check this."


    def apply_bofek(self, bofek=None):
        """Add columns 'kh', 'sy', 'st' and 'kv' to self.data.

        Bofek is a standard Dutch database of soil properties linked to soil
        types indentified by the Bofek index.
        
        It uses the bofek DataFram containing a table of bofek data.

        Parameters
        ----------
        bofek: pd.DataFrame
            table of bofek data, with bofek id in index column having
            at least the following columns ['kh', 'Sy', 'staring', 'ksat_cmpd']
        """
        # --- Verify that the required bofek parameters are in bofek columns
        required_cols = {'kh', 'Sy', 'staring', 'ksat_cmpd'}
        dset = set.difference(required_cols, set(bofek.columns))
        if not dset:
            pass
        else:
            raise KeyError("missing columns [{}] in bofek DataFrame".format(','.join(dset)))

        # --- Verify that all self.data['BOFEK'] keys are in bofek.index,
        # --- so that all parcels get their values!
        dindex = set(self.data['bofek'].values).difference(set(bofek.index))
        if not dindex:
            pass
        else:
            raise KeyError('These keys [{}] in data are not in bofek.index'.format(', '.join(dindex)))

        if 'bofek' in self.data.columns:
            self.data['kh'] = np.array([bofek['kh'].loc[i] for i in self.data['bofek']])
            self.data['sy'] = np.array([bofek['Sy'].loc[i] for i in self.data['bofek']])
            #self.data['st'] = np.array([bofek['staring'].loc[i] for i in self.data['bofek']])
            self.data['kv'] = np.array([bofek['ksat_cmpd'].loc[i] / 100. for i in self.data['bofek']])
        else:
            missing = set(['kh', 'sy', 'kv']).difference(self.data.columns)
            if len(missing) > 0:
                raise ValueError(
                    "If 'bofek' not in database, columns [{}] must be prseent."
                                                .format(', '.join(missing)))


    def apply_defaults(self, defaults=None):
        """Add data missing values to self.data using defaults dict.

        There are more parameters necessary than provided in the database. This
        requires another origin. In this case a set of default values is used.
        These defaults are only applied if there is no corresponding column in self.data.

        Parameters
        ----------
        defaults : dict
            The default parametres with their values.
        """
        defcols = set(defaults.keys()).difference(self.data.columns)
        for dc in defcols: # only for the missing columns
            self.data[dc] = defaults[dc]

# -----------------------------------

def compute_and_set_omega(data=None):
    """Return and set half the wetted ditch circumference in both model layers.
    
    The model has two layers. The ditch may penetrate into the second regional
    model layer. In that case the circumference of the dith in the second layer
    is not zero as it normally would be when the ditch does not reach below
    the bottom of the first layer.

    ditch_omega1 [m] is half the width of the ditch plus its wetted sided.
    ditch_omega2 [m] ia the same for the regional aquifer.

    Fields are added in place.

    Parameters
    ----------
    data: pd.DataFrame
        the parcel data
    """
    #Omega for the cover layer
    hLR  =  0.5 * (data['h_winter'] + data['h_winter'])
    zditch_bottom  = data['AHN'] - data['d_ditch']
    zdeklg_bottom  = data['AHN'] - data['D1']
    b_effective =np.fmax(0,
            np.fmin(zditch_bottom - zdeklg_bottom, data['b_ditch']))
    data['ditch_omega1'] = b_effective + (hLR - zditch_bottom)

    # Omega for the regional aquifer
    zaquif_top     = data['AHN'] - data['D1'] - data['D_CB']
    data['ditch_omega2'] = (data['b_ditch'] +
        (zaquif_top - zditch_bottom)) * (zaquif_top - zditch_bottom >= 0)


def compute_and_set_wpp(data=None):
    """Compute and return extra resistance due to contraction of flow lines.

    Partial penetration of the ditch into the layer (first, and or second layer)
    causes additonal ditch resistance. The resistance is computed analytically.
    See theory for its defivation.

    Parameters
    ----------
    data: pd.DataFrame
        parcel properties
    """
    data['wpp1'] = 2 /  (np.pi * np.sqrt(data['kh'] * data['kv'])) * np.log((
        data['D1'] * np.sqrt(data['kh'] / data['kv']))/(0.5 * data['ditch_omega1']))

    data['wpp2'] = 2 /  (np.pi * np.sqrt(data['kh2'] * data['kv2'])) * np.log((
        data['D2'] * np.sqrt(data['kh2'] / data['kv2']))/(0.5 * data['ditch_omega2']))
    data.loc[np.isnan(data['wpp2']), 'wpp2'] = np.inf
    return


def data_from_dbffile(dbfpath):
    """Return parcel info shape.dbf file into pandas.DataFrame.

    ------------------------------------------------------    
    As of Nov 2025, this is obsolete, we now use geopandas
    ------------------------------------------------------
    Also make sure that the data type is transferred from shapefile to DataFrame.
    
    The original databse is a shpafile with fields that hold the parcel parameters and
    shape that gives the parcel circumference coordinates. The .dbf file is read
    into a pd.DataFrame and serves that the database of parcel parameters, one
    record per parcel, from which the Modflw GGOR model is generated.

    Parameters
    ----------
    dbfpath: str
        name of path to file with .dbf extension, holding parcel data.
    """
    import shapefile
    
    try:
        sf   = shapefile.Reader(dbfpath)
    except FileNotFoundError:
        raise FileNotFoundError(f"Unable to open '{dbfpath}")

    # --- Read shapefile data into pd.DataFrame
    records = [y[:] for y in sf.records()] # turns records into list
    columns=[c[0] for c in sf.fields[1:]]
    data = pd.DataFrame(data=records, columns=columns)

    # --- Get the dtype of each column of the shapefile
    tp   = [t[1] for t in sf.fields[1:]]
    tt = []
    for t, in tp:
        if   t=='N': tt.append(int)
        elif t=='F': tt.append(float)
        elif t=='C': tt.append(str)
        else:        tt.append(object)

    return data.astype({h: t for h, t in zip(data.columns, tt)}) # set column types and return DataFrame


def model_parcel_areas(gr=None, IBOUND=None):
    """Return the model parcel area, for all parcels.

    Parameters
    ----------
    gr: mfgrid.Grid object
    IBOUND: ndarray
        modflow's IBOUND array

    Returns
    -------
    Areas: ndarray
        ndarray of the active cells in each row in the model
    """
    return ((IBOUND[0] != 0) * gr.Area).sum(axis=1)


class Heads_obj:
    """Heads object, to store and plot head data."""

    def __init__(self, sim=None, tdata=None, gr=None):
        """Return Heads_obj.

        Parameters
        ----------
        sim: flopy simulation object
            with access to all models and packages.
            model = sim.get_model(sim.get_name())
        gr: fdm.mfgrid.Grid object
            holds the Modflow grid.
        """
        start = time.perf_counter()
        logger.info("Creating the heads obj with GXGs")
        self.model =  sim.get_model(list(sim.model_names)[0])
        self.gr = gr
        self.HDS = self.model.output.head() # Flopy heads object
        heads = self.HDS.get_alldata() # (nper, nlay, nrow, ncol)
        self.tdata = tdata # original GeoDataFrame with extra columns

        # --- Active cells
        active = self.model.dis.idomain.get_data(); active[active!=0] = 1
        Arel = gr.AREA * active / (gr.AREA * active).sum(axis=-1)[:, :, np.newaxis]
        
        # --- Row-average heads (cell width and active cells taken into account)
        self.avgHds = (heads * Arel[np.newaxis, :, :, :]).sum(axis=-1) # (nper, nlay, nrow)
        self.GXG = self.get_GXG()
        logger.info(f"Heads_obj + GXG created, in {time.perf_counter() - start:.2f} seconds.")


    def get_GXG(self):
        """Return GXG_dict with GXG, HG3, VG3, LG3.

        This GXG object holds the GLG, GVG and GHG, i.e. the lowest, hightes and spring
        groundwater head information and their long-term averaged values based on
        the number of hydrological years implied in the given tdata (meteo data).
        
        (A hydrological year runs form Apriol1 through March 31 the next year.
        The GHG and GHG are based on the 14th and 28th of each month in hydrological year. But the GVG is the aveage of March 14, March 28 and April 14 of the current year (beginning of the hydrlogical year))

        self.gxg is a recarray with all the individual records. (nyear * 9, nparcel)
        self.GXG is a recarray with the long-time averaged values (nparcel).

        @TO 2020-08-31
        """
        _start = time.perf_counter()
        
        logger.info("Computing HG3, LG3, VG3 and GXG")
        
        # --- average parcel head if the first layer at the hand-pick dates
        # --- continue with just these measurements
        ahds = self.avgHds[self.tdata['hand'], 0, :]   # (nthand, nparcel)
        nparcel = ahds.shape[-1]
        
        # --- meteo with extra columns at the hand-picked dates
        tdath = self.tdata.loc[self.tdata['hand']]      # just use the hand data
        tdxh  = tdath.index 
        

        # --- GLG and GHG
        nparcel = ahds.shape[-1] # avgHds (nthand, nparcel)
        hyears = np.asarray(np.unique(self.tdata['hyear'].values), dtype=int)

        # --- Cut off incomplete start and end hydrological years
        # --- A robust way is to discard hyears with less than 24 hand values
        # --- Because we have only a point at the 14t adn 28th of each month (24 in total)
        mask = tdath['hyear'] == hyears[0]
        if len(tdath.loc[mask]) < 24: # We have 24 'hand' measurements per year
            hyears = hyears[1:]
            
        mask = tdath['hyear'] == hyears[-1]
        if len(tdath.loc[mask]) < 24: # we have 24 'hand' measurements per year
            hyears = hyears[:-1]
        
        # --- Skip the first hydological year because it's a run-in year for Modflow
        hyears = hyears[1:]
                
        self.VG3 = dict() # per hyear (3, nparcel) heads and (3, nparcel) timestamps
        self.HG3 = dict() # per hyear (3, nparcel) heads and (3, nparcel) timestamps
        self.LG3 = dict() # per hyear (3, nparcel) heads and (3, nparcel) timestamps

        for hyear in hyears:
            ah = ahds[tdath['hyear'] == hyear]
            td = np.array(tdxh[tdath['hyear'] == hyear])

            # --- indices to select which time data (no argsort here)
            Ias = np.arange(len(td), dtype=int).reshape(len(td), 1) @ np.ones((1, nparcel), dtype=int)
            
            # --- to make this work, our hydrological year starts at March 14th, don't change!
            self.VG3[hyear] = {'h': np.take_along_axis(ah,          Ias[:3, :], axis=0),
                               't': np.take_along_axis(td[:, None], Ias[:3, :], axis=0)
            }
            
            # --- argsort along time axes to allow selecting the 3 highest and lowest values in hyear
            Ias = np.argsort(ah, axis=0)  # Indices of (len(ah), nparcel))
            
            self.LG3[hyear] = {'h': np.take_along_axis(ah,          Ias[:3, :], axis=0),
                               't': np.take_along_axis(td[:, None], Ias[:3, :], axis=0)
            }
            
            self.HG3[hyear] = {'h': np.take_along_axis(ah,          Ias[-3:, :], axis=0),
                               't': np.take_along_axis(td[:, None], Ias[-3:, :], axis=0)
            }

        GVG = np.zeros(nparcel)
        for hyear, vg3 in self.VG3.items():
            GVG += vg3['h'].mean(axis=0)
        GVG /= len(self.VG3) # same as len(hyears)

        GLG = np.zeros(nparcel)
        for hyear, lg3 in self.LG3.items():
            GLG += lg3['h'].mean(axis=0)
        GLG /= len(self.LG3) # same as len(hyears)

        GHG = np.zeros(nparcel)
        for hyear, hg3 in self.HG3.items():
            GHG += hg3['h'].mean(axis=0)
        GHG /= len(self.HG3) # same as len(hyears)
        
        self.GXG = {'hyears': hyears,
                    'GHG': GHG,
                    'GVG': GVG,
                    'GLG': GLG
        }
        
        logger.info(f"Computing GXG done in {time.perf_counter() - _start:.2f} seconds")
        return self.GXG
      

    def plot(self, axs=None, tdata=None, parcel_data=None,
                   parcels=[0, 1, 2, 3, 4],
                   titles=['Heads0', 'Heads1'], xlabel='time [d]', ylabels=['head0 [m]', 'head1 [m]'],
                   figsize=(14, 8), loc='best', plotGXG=True, **kwargs):
        """Plot the running heads in both layers.

        Parameters
        ----------
        ax: plt.Axes
            Axes to plot on.
        tdata: pd.DataFrame with columns 'RH', 'EVT24 and 'summer'
            time datacorresponding to the avgHds data
        parcel_data: geopandas.GeoDataFrame with parcel data
            parcel properties data (used to generate labels)
        parcels: sequence of ints (tuple, list)
            None is all parcels.
            The parcel nrs to show in the graph
        titles: str
            The titles of the 2 charts
        xlabel: str
            The xlabel
        ylabels: str
            The ylabels of the 2 charts
        figsize: tuple of two
            Width and height om image in inches if image is generated and ax is None.
        loc: str (default 'best')
            location to put the legend
        GXG: boolean
            whether or not to plot the GXG on axs[0]
        kwargs: Dict
            Extra parameters passed to newfig or newfigs if present

        Returns
        -------
        The one or two plt.Axes`ax
        """
        parcels = selection_check(parcels, n=self.gr.ny)

        if axs is None:
            axs = etc.newfigs(titles, xlabel, ylabels, figsize=figsize, **kwargs)  
            
        clrs = etc.color_cycler()
        lw = 1
        
        for ilay, ax in enumerate(axs):
            for ip in parcels:
                clr = next(clrs)
                if ilay == 0:
                    # --- plot DRN height
                    hDr = (parcel_data['AHN'] - parcel_data['d_drain']
                                           ).loc[ip] * np.ones(len(tdata))
                    
                    # --- get h_winter and h_summer
                    hLR = parcel_data['h_winter' ].loc[ip] * np.ones(len(tdata))                    
                    hLR[tdata['summer']] = parcel_data['h_winter'].loc[ip]

                    # --- plot h_winter and h_summer
                    ax.plot(tdata.index, hLR, clr, ls='dashed', lw=lw,
                           label=f"parcel {ip}, hLR")
                    ax.plot(tdata.index, hDr, clr, ls='dashdot', lw=lw,
                           label=f"parcel {ip}, zdr")

                ax.plot(tdata.index, self.avgHds[:, ilay, ip], clr, ls='solid',
                             lw=lw, label=f"parcel {ip}")
                    
                if plotGXG:
                    ax = axs[0]
                    hyears = self.GXG['hyears']
                    
                    ax.set_title(ax.get_title() + f" {hyears[0]}-{hyears[-1]}")            

                    ax.axhline(self.GXG['GHG'][ip], c='b', label='GHG')
                    ax.axhline(self.GXG['GVG'][ip], c='g', label='GVG')
                    ax.axhline(self.GXG['GLG'][ip], c='r', label='GLG')
                    
                    for iy, hyear in enumerate(self.VG3):
                        if iy == 0:
                            labels = ["HG3", "VG3", "LG3"]
                        else:
                            labels = ["", "", ""]
                        ax.plot(self.HG3[hyear]['t'][:, ip], self.HG3[hyear]['h'][:, ip], 'b^', label=labels[0])                    
                        ax.plot(self.VG3[hyear]['t'][:, ip], self.VG3[hyear]['h'][:, ip], 'go', label=labels[1])                    
                        ax.plot(self.LG3[hyear]['t'][:, ip], self.LG3[hyear]['h'][:, ip], 'rv', label=labels[2])
                    
            ax.legend(loc=loc, fontsize='xx-small')
            
            plot_hydrological_year_boundaries(ax=ax, tindex=self.tdata.index)
                        
        return axs


def plot_hydrological_year_boundaries(ax=None, tindex=None):
    """Plot hydrological year boundaries on a given axis.

    Parameters
    ----------
    ax: plt.Axes
        an existing axes with a datatime x-axis
    tindex: DateTime index
        tindex to use for this graph.
    """
    years = np.unique(np.array([t.year for t in tindex]))

    if isinstance(ax, plt.Axes): ax = [ax]
    for a in ax:
        for yr in years:
            t = np.datetime64(f'{yr}-04-01')
            if t > tindex[0] and t < tindex[-1]:
                a.axvline(t, color='gray', ls=':')



def show_boundary_locations(lbls=None, CBC=None, iper=0, figsize=(10,8.5)):
    """Show the location of the nodes in recarray given CBC data.

    The refers to ['WEL', 'DRN', 'GHB', 'RIV', 'CHD'] for which the data
    from the CB files are returned as recarrays with fields 'node' and 'q'

    Parameters
    ----------
    lbls: list, tuple
        one or more of ['WEL', 'DRN', 'GHB', 'RIV']
    CBC: open file handle
        CBC file handle
    iper: int
        stress period number
    figsize: 2 tuple of floats
        size of the resulting figure.
    """
    if lbls is None:
        lbls = ['WEL', 'DRN', 'GHB', 'RIV']
    elif isinstance(lbls, str):
        lbls = [lbls]

    for lbl in lbls:
        IB = np.zeros((CBC.nlay, CBC.nrow, CBC.ncol), dtype=int)
        nodes = CBC.get_data(text=cbc_labels[lbl])[iper]['node'] - 1
        IB.ravel()[nodes] = 1

        titles=[f"Top layer, lbl={lbl}, iper={iper}",
                f"Bottom layer, lbl={lbl}, iper={iper}"]
        ax = etc.newfigs(titles=titles, xlabel='column', ylabels=['row', 'row'],
                     sharx=True, sharey=True, figsize=figsize)
        ax[0].spy(IB[0], marker='.', markersize=2)
        ax[1].spy(IB[1], marker='.', markersize=2)


class Watbal_obj:
    """Water balance object."""

    def __init__(self, sim=None, dirs=None, gr=None):
        """Return Watbal object carrying the water budget for all cross sections in m/d.

        Parameters
        ----------
        sim: flopy sim object
            simulation object, model = model.get_model(sim.model_name) object
            holds the simulation with access to all its models and packages
        dirs: mf6tools.Dirs object
            Holds the paths to the directories involbed in this simulation.
        gr: fdm_tools.mfgrid.Grid Object
            holds the structured model grid

        Generates
        ---------
        self.W : np.recarray having labels ['RCH', 'EVT', 'GHB etc'] and values that
            are np.ndarrays with shape (1, nlay, nrow, nper)
            Each row has the cross sectional discharge in m/d.

        Note that labels must be adapted if new CBC packages are to be included.
        
        If everything is ok, then the water balance should add up to zero.
        We should check this.

        @TO 20170823 as function
        @TO 20200907 turned Watbal into class Wabal_Obj
        @TO 20231023 update to mf6, changed W to dict for easier understanding
        """
        start_watbal = time.perf_counter()
        
        logger.info("Creating the Watbal_obj")
        
        required_modflow_packages = ['RCH', 'EVT', 'WEL', 'GHB', 'RIV', 'DRN',
                                     'FLF', 'STO-SS', 'STO-SY']
        self.labels = required_modflow_packages #Layer 0 labels

        self.gr = gr
        
        model = sim.get_model(list(sim.model_names)[0])
        
        logger.info("model.output.budget() takes time ...")
        
        with log_timed(logger, 'loading model.output.budget()'):
            self.CBC = model.output.budget()
        
        # Get structured flows because we need the flow_lower_face
        grb_file = os.path.join(dirs.GWF, sim.name + 'Gwf.dis.grb') # Flopy's grid
        
        with log_timed(logger, 'computing stuctured flows'): 
            fflows = get_structured_flows_as_dict(self.CBC, grb_file=grb_file)

        kstpkper = self.CBC.get_kstpkper()
        
        # To get the CBC values in m/d for the entire region represented in the database:
        # Because each modeled parcel is only a cross of the real parcel, the model
        # does not know the parcel's true, only the database knows it.
        # To get the CBC values in m/d for the entire model:
        # First sum the model's CBC values over the parcel, which yields m3/d, and divide
        # by the parcel's active model are, which yields the values in m/d for each parcel.
        # The to get the parcel's contribution to the total CBC  values in m/d for the
        # entire modeled region, divide by the size of the model's region and multiply bye
        # The parcel's true area. The region's true area is the sum over all parcels true area.
        
        # --- Multiply array with active cells. Active is 1 and inactive is 0
        active = model.get_package('DIS').idomain.data
        active[active !=0 ] = 1

        # --- Area of each cross section
        A_xsec = (gr.Area * active[0]).sum(axis=-1)
         
        # --- Contribution of each true parcel area to total of the modeled region                
        self.W = dict() # keys are watbal labels,  values are float (nlay, nrow, nper)
        
        _vals3D = np.zeros(gr.shape).ravel()
        print('Setting up Watbal_Obj ...')
        for lbl in self.labels:
            self.W[lbl[:3]] = np.zeros((self.CBC.nlay, self.CBC.nrow, self.CBC.nper))
            print(lbl, end='')
            
            # For these labels, the packages of which CBC values are obtained
            # in a recarray with 'node' and 'q' with:
            # dtype = [('node', '<i4'), ('node2', '<i4'), ('q', '<f8')])
            # The dimension will be [L3/T]
            if lbl in ['WEL', 'GHB', 'RIV', 'DRN', 'RCH', 'EVT']:
                
                vals = self.CBC.get_data(text=lbl) # as recarray with given dtype
                
                for iper in range(self.CBC.nper):
                    # Vals contains all non-zero node values at once
                    Ivals = vals[iper]['node'] - 1 # zero-based Modflow node numbers
                    Q     = vals[iper]['q']
                    
                    _vals3D[:] = 0. # reset _vals3D to zeros for next loop cycle
                    _vals3D[Ivals] = Q
                    
                    # Sum _vals3D over the columns of the model strutured grid
                    self.W[lbl][:, :, iper] = _vals3D.reshape(gr.shape).sum(axis=-1)
                    
                    # Show progress
                    if iper % 100 == 0: print('.',end='')                
            elif lbl in ['STO-SY', 'STO-SS']:
                vals = self.CBC.get_data(text=lbl) # as a full 3D arrwy over the grid
                for iper in range(self.CBC.nper):
                    self.W[lbl[:3]][:, :, iper] += np.sum(vals[iper], axis=-1) # W['STO']
                    if iper % 100 == 0: print('.', end='') # Show progress                
            elif lbl in ['FLF']:                
                for iper in range(self.CBC.nper):
                    self.W[lbl][-1, :, iper] = +fflows[kstpkper[iper]]['flf'][0][0].sum(axis=-1)
                    self.W[lbl][ 0, :, iper] = -fflows[kstpkper[iper]]['flf'][0][0].sum(axis=-1)
                    if iper % 100 == 0: print('.', end='') # Show progress                
            else:
                raise ValueError(f"Illegal label {lbl}")
            print(f"Last iper for {lbl}: {iper}")
        for lbl in self.W: # prevents clash with STO-SS and STO-SY            
            self.W[lbl] /= A_xsec[np.newaxis, : , np.newaxis] # from m3/d to m/d # not yet summed over all parcels!

        print('Done. See self.CBC and self.W')
        logger.info(f"Watbal_obj created in {time.perf_counter() - start_watbal:.2f} seconds")


    def plot(self, parcel_data=None, tdata=None,
                    parcels=None, sharey=False, ax=None):
        """Plot the running water balance of the GGOR selected parcels or entire area in mm/d.

        Parameters
        ----------
        parcel_data: pd.DataFrame
            parcel spatial and property data
        tdata: pd.DataFrame with time index
            time_index must correspond with MODFLOW's output
            if not specified, then CBB.times is used
        parcels: list or sequence
            the water budget will be taken over the indices in parcels.
            Use None for all.
        sharey: bool
            the y-axis of the two charts will be shared if True
        ax: plt.Axis object or None (axes will the be created)
            axes for figure.
        """
        m2mm = 1000. # from m to mm conversion

        parcels = selection_check(parcels, n=self.gr.ny)

        tindex = self.CBC.times if tdata is None else tdata.index # time index

        # Sum over all (selected) parcels. The watbal values are in mm/d. To sum over all
        # parcels multiply by their share of the regional area [-]
        Arel = (parcel_data['A_parcel'].values[parcels] /
                parcel_data['A_parcel'].values[parcels].sum()
                ) # (nlay, parcels, nper)

        clrs = [watbal_label[L]['clr'] for L in watbal_label]
        lbls = [watbal_label[L]['leg'] for L in watbal_label]

        # Axes title ttl
        if isinstance(parcels, slice):
            ttl = 'Taken over parcels[{}:{}:{}]'.format(
                        parcels.start, parcels.stop, parcels.step)
        else:
            ttl = 'Taken over parcels [{}]'.format(
                ', '.join(['{}'.format(s) for s in parcels]))

        # Two axes, one for the running budget of the top layer
        # and one for the running budget of the bottom layer
        if ax is None:
            ax = etc.newfigs(titles=(
                    'Water balance top layer. '   + ttl,
                    'Water balance botom layer. ' + ttl), xlabel='time', ylabels=['mm/d', 'mm/d'],
                         figsize=(14, 8), sharey=False, sharex=True)

        if not isinstance(ax, (list, tuple, np.ndarray)):
            raise ValueError("Given ax must be an sequence of 2 axes.")

        V0 = np.zeros((len(watbal_label), self.CBC.nper)) # running budget of top layer
        V1 = np.zeros((len(watbal_label), self.CBC.nper)) # running budget of bot layer
        # Accumulate to generate an plt.Axes.stackplot
        for i, lbl in enumerate(watbal_label):
            V0[i] = (self.W[lbl][ 0, parcels, :] * 
                     Arel[np.newaxis, :, np.newaxis]).sum(axis=1) * m2mm
            V1[i] = (self.W[lbl][-1, parcels, :] * 
                     Arel[np.newaxis, :, np.newaxis]).sum(axis=1) * m2mm
            
        # TODO: Verify all budget sums to zero.

        ax[0].stackplot(tindex, V0 * (V0>0), colors=clrs, labels=lbls)
        ax[0].stackplot(tindex, V0 * (V0<0), colors=clrs) # no labels
        ax[1].stackplot(tindex, V1 * (V1>0), colors=clrs, labels=lbls)
        ax[1].stackplot(tindex, V1 * (V1<0), colors=clrs) # no labels

        ax[0].legend(loc='best', fontsize='xx-small')
        ax[1].legend(loc='best', fontsize='xx-small')

        return ax


def get_test_parcels(path, sheet_name, test_id_col='Test'):
    """Return the parcel test data from given workbook and worksheet.

    Parameters
    ----------
    path to excel_workbook: str
        excel workbook file with extension
    sheet_name: str
        sheet name within workbook
    """
    parcel_data = pd.read_excel(path, sheet_name=sheet_name, engine='openpyxl')

    # Fill empty fields with the base valuesl
    other_cols = parcel_data.columns.tolist()
    test_id_col = other_cols.pop(other_cols.index(test_id_col))

    for col in other_cols:
        missing = parcel_data[col].isna().values
        parcel_data.loc[missing, col] = parcel_data[col].iloc[0]

    parcel_data[test_id_col] = parcel_data[test_id_col].fillna(method='ffill')

    compute_and_set_omega(parcel_data)
    compute_and_set_wpp(parcel_data)

    return parcel_data

if __name__ == '__main__':
    start_main = time.perf_counter()
    
    logger.info("Runing ggor_tools __main__ section")
    logger.info("TODO: check this main section")
    
    test = True
    
    HOME = os.getcwd()
    
    logger.warning("cwd = {}".format(os.getcwd()))
    dirs = mf6tools.Dirs(HOME)
    dirs.add_case('AAN_GZK')

    mf_parameters_wbk = os.path.join(dirs.mf_parameters, 'mf_parameters.xlsx')

    #Get the meteo data from an existing file or directly from KNMI
    meteo = knmi.get_weather(stn=240, start='20100101', end='20191231',
                                  folder=dirs.meteo)

    # Add columns "summer' and "hyyear" to it"
    tdata = handle_meteo_data(meteo, summer_start=4, summer_end=10)
    tdata = tdata.iloc[:1000] # Limits data set just for testing
    if test:
        tdata = gen_testdata(tdata=tdata,
                              RH  =(270, 0.0, 0.002, 0.004),
                              EV24=(180, 0.0, 0.001, 0.002),
                              )

        parcel_data = get_test_parcels(os.path.join(
                                dirs.case, 'pdata_test.xlsx'), 'parcel_tests1')

        # Special test
        parcel_data = parcel_data.iloc[0:4]
    else:
        # Bofek data, coverting from Bofek codes to the soil properties kh, kv and sy.
        # The BOFEK column represents a Dutch standardized soil type. It is used.
        # The corresponding values for 'kh', 'kv' and 'Sy' are currently read from
        # an Excel worksheet into a pandas DataFrame (thus becoming a table)
        bofek = pd.read_excel(os.path.join(dirs.bofek, "BOFEK eenheden.xlsx"),
                              sheet_name = 'bofek', index_col=0, engine="openpyxl")

        # Create a GGOR_modflow object and get the upgraded parcel_data from it,
        # excluding parcels that are too small or too wide
        parcel_data = GGOR_data(
            dirs=dirs, defaults=defaults, bofek=bofek, BMINMAX=(5, 250)).data
    
    # MODFLOW grid
    gr = grid_from_parcel_data(parcel_data=parcel_data, dx=1.0) 

    logger.info(f"ggor_tools __main__ section done in {time.perf_counter() - start_main:.2f} seconds")
