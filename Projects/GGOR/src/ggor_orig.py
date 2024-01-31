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
import numpy as np
import pandas as pd
import shapefile
from inspect import signature
import matplotlib.pyplot as plt
from KNMI import knmi
from fdm import mfgrid
import flopy
import flopy.utils.binaryfile as bf
from collections import OrderedDict
import logging

from get_fff import get_flow_lower_face

logging.basicConfig(level=logging.WARNING, format=' %(asctime)s - %(levelname)s - %(message)s')

NOT = np.logical_not
AND = np.logical_and
OR  = np.logical_or

# User settings to change execution behavior at runtime
user_settings = {
    'case_name': 'AAN_GZK',    # Name of read contained in available database.
    'use_w_not_c': False, # Use anal. form. instead of entry resistance and circumf.
    'modflow' : flopy.mf6,
    'd_drn': 0.,  # m, Drain depth for surface runoff.
    'use_packages': ['wel', 'drn',  'ghb', 'riv', 'rcha', 'evta'],
    'silent': False,
    'dx': 1.
}

#% Names in the databbase and their replacements used in GGOR
# Names not mentioned are not replaced.
colDict = {  'Bofek'     : 'bofek',
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
defaults = {'d_drain': 0,  # [m] Tile drainage depth below local ground elevation
                           #, may be zero if no drains are present.
            'd_trench': 0.3, # [m] Trench depth in case present.
            'c_drain': 5., # [d] Tile drainage areal resistance. Also used for trenches.
            'wi_ditch' : 2.,  # [d] Ditch resist. for flow from ditch to ground. (analytical)
            'wi_ditch2' : 2.,  # [d] Ditch resis. reg. aquif. for flow to ground. (analytical)
            'wo_ditch' : 1.,  # [d] Ditch resis. when flow is to ditch. (anlaytical)
            'wo_ditch2' : 1.,  # [d] Ditch resis. reg. aquif. for flow to ditch. (analytical)
            'ci_ditch' : 2.,  # [d] Ditch bottom and side entry resistance (applied to Omega)
            'co_ditch' : 1.,  # [d] Ditch bottom and side entry resistance (applied to Omega)
            'd_ditch' : 1.0, # [m] depth of ditch below ground surface
            'b_ditch' : 0.75, # [m] half-width of the ditch
            'D_CB' : 0.1, # [m] (dummy thickness) van basisveenlaag (CB=confining bed)
            'D2' : 40., # [m] thickness of regional aquifer
            'S2' : 1e-3,# [-] total (elastic) storage coeff. of regional aquifer
            'kh2': 30., # [m/d]  horizontal condutivity of regional aquifer
            'kv2':  6., #[ m/d]  vertical conductivity of regional aquifer
            'ET_surfd': 1.0, # [m] depth of surf in ET below ground surface.
            'ET_exdp': 2.5, # [m] Modflow's extinction depth (see ET package)
            }


#% Modflow cell-by-cell flow labels: Translates short labels to those in the CBC file.
cbc_labels = {
    'STO': 'STORAGE',
    'FLF': 'FLOW LOWER FACE ',
    'WEL': 'WELLS',              # used for seepage. Are in second layer.
    'EVT': 'ET',
    'GHB': 'HEAD DEP BOUNDS',
    'RIV': 'RIVER LEAKAGE',
    'DRN': 'DRAINS',             # to capture surface runoff, also used for actual drains and for trenches.
    'RCH': 'RECHARGE',
    }

# For legend when of running water budget plot
#leg is legend for this label in the graph
#clr is the color of the filled graph
watbal_label = OrderedDict({
                'RCH': {'leg': 'RCH', 'clr': 'green'},
                'EVT': {'leg': 'EVT', 'clr': 'gold'},
                'WEL': {'leg': 'WEL(in wvp2)' , 'clr': 'blue'},
                'DRN': {'leg': 'DRN(drn/trenches/runoff)', 'clr': 'lavender'},
                'RIV': {'leg': 'RIV(ditch out)', 'clr': 'magenta'},
                'GHB': {'leg': 'GHB(ditch in+out)', 'clr': 'indigo'},
                'FLF': {'leg': 'FLF(leakage)', 'clr': 'gray'},
                'STO': {'leg': 'STO', 'clr': 'cyan'}})

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
        parcels = slice(0, None, 1)
    elif isinstance(parcels, int):
        parcels = slice(0, parcels, 1)
    elif isinstance(parcels, (tuple, list, np.ndarray, range, slice)):
        pass
    else:
        raise ValueError("parcels must be None (for all), an int or a sequence")
    if isinstance(parcels, slice):
        parcels = np.arange(n)[parcels]
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
        # index array telling which of the tuple values to pick
        # Daynumber since start of tdata.index
        daynum = (tdata.index - tdata.index[0]) / np.timedelta64(1, 'D')
        period = int(kwargs[key][0]) # days
        values = np.array(kwargs[key][1:])
        I = np.asarray((daynum // period) % len(values), dtype=int)
        # Add or replace column in tdata copy
        tdata[key] = np.array(values[I])
    return tdata


class Dir_struct:
    """GGOR directory structure.

    Expected directory structure

    GGOR_home/ # GGOR_home is the GGOR_home directory specified above.
            bin/
                {modflowname}.exe
                {modflowname}.mac
            src/
                 analytic/
                 numeric/
            data/
                 mf_parameters
                    mf_parameters.xlsx
                 bofek/
                 meteo/
                 spatial/
                         AAN_GZK # Stores QGIS data (shapefiles etc)
                         ....
            doc/
            notebooks/
            cases/ # Stores MODFLOW input and output files
                  AAN_GZK
                  ....
    [tools] i.e. ~/python/tools in the PYTHONPATH, see ~/.zshrc
    """    
    def __init__(self, home='.', case=user_settings['case_name'] ):
        """Generate GGOR directory structure.

        Parameters
        ----------
        home: str
            path to GGOR home directory
        case: str
            the name of the current case
        """
        self.case_name = case
        
        self.home = os.path.abspath(os.path.expanduser(home))        
        self.bin = os.path.join(self.home, 'bin') # Binary folder
        self.src = os.path.join(self.home, 'src')
        self.data = os.path.join(self.home, 'data')
        self.cases = os.path.join(self.home, 'cases') # folder where cases are stored."""
        self.mf_parameters = os.path.join(self.data, 'mf_parameters', 'mf_parameters.xlsx') # mf6 params wbk
        self.meteo = os.path.join(self.data, 'meteo') # meteo data folder
        self.spatial = os.path.join(self.data, 'spatial') # folder with spatial data.
        self.bofek = os.path.join(self.data, 'bofek') # Folder of bofek units.
        self.case = os.path.join(self.spatial, self.case_name) # Folder of current case        
        self.mf6 = os.path.join(self.case, 'mf6')  # Folder where case output goes

        if not os.path.isdir(self.cases):
            os.mkdir(self.cases)
        if not os.path.isdir(self.mf6):
            os.mkdir(self.mf6)
            
        self.wd = os.path.join(self.cases, self.case_name)
        
        if os.name == 'posix':
            # Version 2021 mf6.2.1. Compiled on Mac using pymake.
            self.exe = os.path.join(self.bin, 'mf6.mac')            
        else:
            self.exe = os.path.join(self.bin, 'mf2005.exe') # Version 2021 mf6.2.1
                
        #### verify existance of required folders and files ###################
        for folder_name in [self.home, self.bin, self.src, self.data, self.cases,
                     self.meteo, self.spatial,
                     self.bofek, self.wd, self.case, self.mf6]:
            assert os.path.isdir(folder_name), "Folder not found: {}".format(folder_name)
        
        for file_name in [self.mf_parameters, self.exe]:
            assert os.path.isfile(file_name), "File not found: {}".format(file_name)
        
        os.chdir(self.wd)
        print("Current working directory: {}".format(os.getcwd()))
        return

def get_mf6_params_from_excel(wbk):
    """Read mf6 parameters from workbook (sheetname='mf6').
    
    Read mf6 parameters from Excel workbook as a pd.DataFrame.
    Turn that DataFrame into a dicionary where the keys are the modflow packages
    and the subkeys are the variable names.
    All values are read as strings. The colun "type" is used to convert
    values that should not be strings into the desired type.
    
    The headers in the sheet 'MF6' of the workbook are:
    ['package', 'parameter', 'value', 'type', 'Meaning']
    
    Parameters
    ----------
    wbk: str (path)
        Excel workbook name that holds the mf6 parameters in sheet 'MF6'
        
    Returns
    -------
    params: dictionary
        a dict with package name as key, in which each item is a dictionary
        specifying the values for this key.
        
    @TO 220413
    """
    paramsDf = pd.read_excel(wbk, 'MF6', header=0, usecols="A:D", engine="openpyxl").dropna(axis=0)
    mf6_params = dict()
    for pkg in np.unique(paramsDf['package'].values):
        mf6_params[pkg]=dict()
        
    ### Convert specified strings to correct type
    for i in paramsDf.index:
        pkg, parameter, value, type = paramsDf.loc[i]
        if value == 'None':   continue # Skip (emtpy line)
        elif type == 'None':  value = None
        elif type == 'str':   pass
        elif type == 'float': value = float(value)
        elif type == 'int':   value = int(value)
        elif type == 'bool':  value = True if value == 'True' or value == 1 else False
        elif type == 'list':  value = exec(value)
        else: # Immediately verify unknown types
            raise ValueError("Unknown parameter type pkg={}, {}, {}, {}".format(
                            pkg, parameter, value, type))
        mf6_params[pkg][parameter]=value
    return mf6_params


def handle_meteo_data(meteo_data=None, summer_start=4, summer_end=9):
    """Set and store meteo data and add the columns summer, hyear and hand.

    Added columns are
        summer: bool
        hyear: hyddrological year. They atart at March 14 to get proper GVG. Do not change this !
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
    dcol = {'RH', 'EVT24'}.difference(meteo_data.columns) # Check for existance of both columns in meteo_data.columns.
    if dcol:
        KeyError("Missing column[s] [{}] in meteo DataFrame".format(', '.join(dcol)))

    #verify, data are in m/d
    if not meteo_data['RH'].median() < 0.01 and meteo_data['EVT24'].median() < 0.01:
        AssertionError("Median of Precipitration = {:5g} and median evapotranspiration = {:4g}\n"
                          .format(meteo_data['RH'].median(), meteo_data['EVT24'].median()) +
                       "Percipication and or evapotranspiration likely not in m/d!")

    # Add boolean column indicating summer (needed to set summer and winter ditch levels)
    meteo_data['summer'] = [True if t.month in range(summer_start, summer_end) 
                               else False
                               for t in meteo_data.index]

    # hydrological year column 'hyear'
    hyear_start_month = 3   # Don't change! It's needed in the GXG class
    hyear_start_day   = 14  # Don't change! It's needed in the GXG class
    meteo_data['hyear'] = [t.year
        if t.month >= hyear_start_month and t.day >= hyear_start_day
        else t.year - 1 for t in meteo_data.index]

    meteo_data['hand'] = [t.day % 14 == 0 for t in meteo_data.index]

    return meteo_data


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

    return mfgrid.Grid(xGr, yGr, Z)


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


def get_drain_elev_with_trenches(pdata=None, gr=None):
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
    d_drn = user_settings['d_drn']
    
    elev = ((pdata['AHN'] - pdata['d_drain']).values[:, np.newaxis]
                                                    * np.ones((1, gr.nx)))
    I = np.arange(gr.nx, dtype=int)
    xLeft = gr.X[0][:,:-1] # left side of cells in top layer for all rows)
    for iparcel, (b, ntr, ahn, dtr) in enumerate(
            zip(pdata['b'], pdata['n_trench'], pdata['AHN'], pdata['d_trench'])):
        elev[iparcel] = ahn - d_drn # default drain depth (simulates surface runoff)
        
        if ntr: # only if parcel has one or more trenches
            atr = b / ntr # distance between trenches or trench and ditch
            xtr = np.cumsum(np.ones(ntr) * atr) # locations of trenches
            Itr = [I[xLeft[iparcel] <= _x][-1] for _x in xtr] # cell row indices of trenches
            elev[iparcel][Itr] = ahn - dtr
    return elev # Ny * Nx array

    for iparcel, (b, ntr, ahn, dtr) in enumerate(
            zip(pdata['b'], pdata['n_trench'], pdata['AHN'], pdata['d_trench'])):
        print(b, ntr, ahn, dtr)


def get_cond_DRN(pdata=None, gr=None):
    """Return conducntance for use by DRN (Surface Runoff).
    
    Parameters
    ==========
    pdata: dict
        parcel properties pd.DataFrame
    """
    c_drain = pdata['c_drain'].values[:, np.newaxis] * np.ones((1, gr.nx))
    condDRN = gr.Area / c_drain
    return condDRN # Conductance for all cells in top layer as an Ny * Nx array


def get_cond_GHB(pdata=None, gr=None):
    """Retrun conductance for use by GHB for two layers first column.
    
    Parameters
    ==========
    pdata: pd.DataFrame
        parcel properties
    gr: structured Grid object
        Class that holds the grid
    """
    if user_settings['use_w_not_c']:
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


def get_RIV_Cond(pdata=None, gr=None):
    """Return the Conductances of ditches for use by RIV package.
    
    parameters
    ==========
    pdata: dict
        parameter data
    """
    if user_settings['use_w_not_c']:
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

    def __init__(self, GGOR_home=None, case=None,  bofek=None, BMINMAX=(5., 500.), defaults=None):
        """Get parcel data for use by GGOR.

        Parameters
        ----------
        GGOR_home: str
            path to the GGOR+home directory
        case: str
            name of the case. Must correspond with folder in numeric and with folder in cases
            foldder in cases will be made if necessary.
        bofek: pd.DataFrame
            bofek values for ['kh', 'Sy', 'staring', 'ksat_cmpd'], the index
            is the bofek_id number.
            bofek are Dutch standardized soil parameters, see Internet.
        BMINMAX: tuple of 2 floats
            min and max halfwidth value for parcels to be considered.
        """
        dirs = Dir_struct(home=GGOR_home, case=case)

        # read dbf file into pd.DataFrame
        self.data = data_from_dbffile(os.path.join(dirs.case, case + '.dbf'))

        # replace column names to more generic ones
        self.data.columns = [colDict[h] if h in colDict else h
                                         for h in self.data.columns]

        # compute parcel width to use in GGOR
        self.compute_parcel_width(BMINMAX=BMINMAX)

        # set kh, kv and Sy from bofek
        self.apply_bofek(bofek) # bofek is one of the kwargs a pd.DataFrame

        # add required parameters if not in dbf
        self.apply_defaults(defaults)

        self.compute_and_set_omega()

        self.compute_and_set_wpp()


    def compute_and_set_omega(self):
        """Compute and set the half wetted ditch circumference in the two model layers.

        ditch_omega1  [m] is half the width of the ditch plus its wetted sided.
        ditch_omega2 [m] ia the same for the regional aquifer.

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
        
        The parcel width is computed fromt the parcel area A and parcel cirumference Om
        available in the database. It's the longest axis of a rectangulre area
        deduced from A and Om.

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
        det   = Om ** 2 - 16 * A             # determinant
        I     = (det>=0);                   # determinant>0 ? --> real solution
        B     = np.nan * np.zeros_like(I)   # init paracel width
        L     = np.nan * np.zeros_like(I)   # init parcel length
        B[ I] = (Om[I] - np.sqrt(det[I]))/4  # width, smallest of the two values
        L[ I] = (Om[I] + np.sqrt(det[I]))/4  # length, largest of the two values
        B[NOT(I)] = np.sqrt(A[NOT(I)])      # if no real solution --> assume square
        L[NOT(I)] = np.sqrt(A[NOT(I)])      # same, for both width and length

        B = np.fmin(B, max(BMINMAX)) # Arbitrarily limit the width of any parcel to BMAX.
        B = np.fmax(B, min(BMINMAX))

        # Add column 'b' to Data holding half the parcel widths.
        self.data['b'] = B/2

        # Use only the parcles that have with > BMIN and that have bofek data
        I=np.where(AND(B > min(BMINMAX), NOT(self.data['bofek']==0)))

        self.data = self.data.iloc[I]

        # Any data left?
        assert len(I) > 0, "Cleaned parcel database has length 0, check this."


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
        # Verify that the required bofek parameters are in bofek columns
        required_cols = {'kh', 'Sy', 'staring', 'ksat_cmpd'}
        dset = set.difference(required_cols, set(bofek.columns))
        if not dset:
            pass
        else:
            raise KeyError("missing columns [{}] in bofek DataFrame".format(','.join(dset)))

        # Verify that all self.data['BOFEK'] keys are in bofek.index, so that
        # all parcels get their values!
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
    try:
        sf   = shapefile.Reader(dbfpath)
    except:
        raise FileNotFoundError("Unable to open '{}'.".format(dbfpath))

    # Read shapefile data into pd.DataFrame
    records = [y[:] for y in sf.records()] # turns records into list
    columns=[c[0] for c in sf.fields[1:]]
    data = pd.DataFrame(data=records, columns=columns)

    # Get the dtype of each column of the shapefile
    tp   = [t[1] for t in sf.fields[1:]]
    tt = []
    for t, in tp:
        if   t=='N': tt.append(int)
        elif t=='F': tt.append(float)
        elif t=='C': tt.append(str)
        else:        tt.append(object)

    return data.astype({h: t for h, t in zip(data.columns, tt)}) # set column types and return DataFrame


def newfig2(titles=['title1', 'title2'], xlabel='time',
            ylabels=['heads [m]', 'flows [m2/d]'],
            xscale='linear',
            yscale=['linear', 'linear'],
            sharex=True,
            sharey=False,
            xlims=None,
            ylims=None,
            size_inches=(12, 6),
            **kwargs):
    """Return ax[0], ax[1] for new plot.
    
    This function conveniently sets up a figure and axis for showing results.
    """
    fig, ax = plt.subplots(2, 1, sharex=sharex, sharey=sharey)
    fig.set_size_inches(size_inches)
    for a, title, ylabel in zip(ax, titles, ylabels):
        a.set_title(title)
        a.set_xlabel(xlabel)
        a.set_ylabel(ylabel)
        a.grid()
    ax[0].set_xscale(xscale)
    ax[1].set_xscale(xscale)
    ax[0].set_yscale(yscale[0])
    ax[1].set_yscale(yscale[1])
    if xlims is not None:
        ax[0].set_xlim(xlims[0])
        ax[1].set_xlim(xlims[1])
    if ylims is not None:
        ax[0].set_ylim(ylims[0])
        ax[1].set_ylim(ylims[1])
    return ax


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
        self.model =  sim.get_model(list(sim.model_names)[0])
        self.gr = gr
        self.HDS = self.model.output.head() # Flopy heads object
        heads = self.HDS.get_alldata() # (nper, nlay, nrow, ncol)

        # Active cells
        active = self.model.dis.idomain.get_data(); active[active!=0] = 1
        Arel = gr.AREA * active / (gr.AREA * active).sum(axis=-1)[:, :, np.newaxis]
        
        self.avgHds = (heads * Arel[np.newaxis, :, :, :]).sum(axis=-1) # (nper, nlay, nrow)
        self.GXG    = GXG_object(tdata=tdata, avgHds=self.avgHds)


    def plot(self, ax=None, tdata=None, parcel_data=None,
                   parcels=[0, 1, 2, 3, 4],
                   titles=None, xlabel='time', ylabels=['m', 'm'],
                   size_inches=(14, 8), loc='best', GXG=True,  **kwargs):
        """Plot the running heads in both layers.

        Parameters
        ----------
        ax: plt.Axies
            Axes to plot on.
        tdata: pd.DataFrame with columns 'RH', 'EVT24 and 'summer'
            time datacorresponding to the avgHds data
        parcel_data: pd.DataFrame
            parcel properties data (used to generate labels)
        parcels: sequence of ints (tuple, list).
            None is all parcels.
            The parcel nrs to show in the graph.
        titles: str
            The titles of the charts.
        xlabel: str
            The xlabel
        ylabels: str
            The ylabels of the 2 charts.
        size_inches: tuple of two
            Width and height om image in inches if image is generated and ax is None.
        loc: str (default 'best')
            location to put the legend
        GXG: boolean
            whether or not to plot the GXG also.
        kwargs: Dict
            Extra parameters passed to newfig or newfig2 if present.

        Returns
        -------
        The one or two plt.Axes`ax
        """
        parcels = selection_check(parcels, n=self.gr.ny)

        if ax is None:
            ax = newfig2(titles, xlabel, ylabels, size_inches=size_inches, **kwargs)
            for a in ax:
                plot_hydrological_year_boundaries(a, tdata.index)
        else:
            for a, title, ylabel in ax, titles, ylabels:
                a.grid(True)
                a.set_title(title)
                a.set_xlabel(xlabel)
                a.set_ylabel(ylabel)

        nt, nLay, ny = self.avgHds.shape

        clrs = 'brgkmcy'
        lw = 1
        for ilay, a in zip(range(nLay), ax):
            for iclr, isel in enumerate(parcels):
                clr = clrs[iclr % len(clrs)]
                a.plot(tdata.index, self.avgHds[:, ilay, isel], clr, ls='solid',
                             lw=lw, label="parcel {}".format(isel))
                if ilay == 0:
                    hDr = (parcel_data['AHN'] - parcel_data['d_drain']
                                           ).loc[isel] * np.ones(len(tdata))
                    hLR = parcel_data['h_winter' ].loc[isel] * np.ones(len(tdata))
                    hLR[tdata['summer']] = parcel_data['h_winter'].loc[isel]

                    a.plot(tdata.index, hLR, clr, ls='dashed', lw=lw,
                           label='parcel {}, hLR'.format(isel))
                    a.plot(tdata.index, hDr, clr, ls='dashdot', lw=lw,
                           label='parcel {}, zdr'.format(isel))
            a.legend(loc=loc, fontsize='xx-small')

            if GXG:
                self.GXG.plot(ax[0], parcels=parcels)

        return ax


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
            t = np.datetime64(f'{yr}-03-14')
            if t > tindex[0] and t < tindex[-1]:
                a.axvline(t, color='gray', ls=':')


class GXG_object:
    """Generate GXG object.

    This object hold the GXG (GLG, GVG, GHG)  i.e. the lowest, hightes and spring
    groundwater head information and their long-term averaged values based on
    the number of hydrologi al years implied in the given tdata.
    (A hydrological year runs form March14 through March 13 the next year, but
     the GXG are based on values of the 14th and 28th of each month only.)

    self.gxg is a recarray with all the individual records. (nyear * 9, nparcel)
    self.GXG is a recarray with the long-time averaged values (nparcel).

    @TO 2020-08-31
    """

    def __init__(self, tdata=None, avgHds=None):
        """Initialize GXG object.

        Parameters
        ----------
        tdata: pd.DataFrame
            tdata, we only need its index
        avgHds: np.nd_array shape = (nt, nz, nParcel)
            The xsection-averaged heads for all parcels and all times.
            Heads aveaged along the x-axis taking into account cel width
            and ignoring inactive cells.
        """
        # The 'hand' data are the values at the 14th and 28th of each month.
        ahds = avgHds[tdata['hand'], 0, :]   # (nparcel, nthand)
        tdat = tdata.loc[tdata['hand']]  # just use the hand data

        # GLG and GHG
        nparcel = ahds.shape[-1] # ny  avgHds (nthand, ny)
        hyears = np.unique(tdat['hyear'])

        # Cut off incomplete start and end hydological years
        if tdat.index[ 0].month != 3 or tdat.index[ 0].day != 14: hyears = hyears[ 1:]
        if tdat.index[-1].month != 2 or tdat.index[-1].day != 28: hyears = hyears[:-1]

        # skip the first hydological year because it's a run-in year for Modflow
        hyears = hyears[1:]
        nyear = len(hyears)

        # Format to store the gxg data in a recarray
        # l=low, h=high, v=spring (voorjaar)
        gxg_dtype = [('t', pd.Timestamp), ('hd', float), ('hyear', int),
                 ('l', bool), ('h', bool), ('v', bool)]

        # The gxg recarray has 9 records per hyear (3 glg, 3 ghg, 3 gvg and nparcel layers).
        # We keep them as individual  values to allow showing them in the time-graph.
        # These are nine individual values  that contribute to the GXG in a hydrologial year
        self.gxg = np.zeros((nyear * 9, nparcel), dtype=gxg_dtype)

        T = (True, True, True)    # To tell the 3 vals in hyr are lowest, highest or spring
        F = (False, False, False) # Same to tell they are not lowest, highest or spring

        for iyr, hyear in enumerate(hyears):
            ah = ahds[tdat['hyear'] == hyear]
            td = tdat.loc[tdat['hyear'] == hyear]
            Ias = np.argsort(ah, axis=0)  # Indices of argsort along time axis

            # Make sure hydrological years start at March 14!!
            assert td.index[0].month ==3 and td.index[0].day == 14, "hyears must start at 14th of March"

            hyr = (hyear, hyear, hyear)

            for ip in range(nparcel):
                Iglg = Ias[0:3, ip] # lowest three values in hyr for parcel ip
                Ighg = Ias[-3:, ip] # highest three values in hyr for parcel ip
                Igvg = slice(0, 3, 1) # spring values 14/3, 28/3, 14/4 (hyear starts 14/3)
                # The three lowest values
                self.gxg[iyr * 9 + 0:iyr * 9 + 3, ip] = np.array(
                    [(t, hd, yr, l, h, v) for t, hd, yr, l, h, v in zip(
                        td.index[Iglg], ah[Iglg, ip], hyr, T, F, F)], dtype=gxg_dtype)
                # The three highest values
                self.gxg[iyr * 9 + 3:iyr * 9 + 6, ip] = np.array(
                    [(t, hd, yr, l, h, v) for t, hd, yr, l, h, v in zip(
                        td.index[Ighg], ah[Ighg, ip], hyr, F, T, F)], dtype=gxg_dtype)
                # The three spring values
                self.gxg[iyr * 9 + 6:iyr * 9 + 9, ip] = np.array(
                    [(t, hd, yr, l, h, v) for t, hd, yr, l, h, v in zip(
                        td.index[Igvg], ah[Igvg, ip], hyr, F, F, T)], dtype=gxg_dtype)

        # Compute and store the long-term averaged values, the actual GXG
        dtype = [('id', int), ('GLG', float), ('GHG', float), ('GVG', float)]
        self.GXG = np.ones(nparcel, dtype=dtype)
        for ip in range(nparcel):
            self.GXG[ip] = (
                ip,
                self.gxg[self.gxg[:, ip]['v'], ip]['hd'].mean(),
                self.gxg[self.gxg[:, ip]['l'], ip]['hd'].mean(),
                self.gxg[self.gxg[:, ip]['h'], ip]['hd'].mean())


    def plot(self, ax=None, parcels=[0, 1, 3, 4, 5], **kwargs):
        """Plot GXG.

        Parameters
        ----------
        parcels : np.array of ints
            sequence of indices to select the parcels for plotting
            parcels have already been checked. No need to repeat that here.
        nmax: int
            maximum number of graphs to plot
        """

        clrs = 'brgkmcy'

        for iclr, ip in enumerate(parcels):
            g = self.gxg.T[ip]
            clr = clrs[iclr % len(clrs)]
            ax.plot(g['t'][g['v']], g['hd'][g['v']], clr, marker='o',
                    mfc='none', ls='none', label='vg [{}]'.format(ip))
            ax.plot(g['t'][g['l']], g['hd'][g['l']], clr, marker='v',
                    mfc='none', ls='none', label='lg [{}]'.format(ip))
            ax.plot(g['t'][g['h']], g['hd'][g['h']], clr, marker='v',
                    mfc='none', ls='none', label='hg [{}]'.format(ip))

        hyears = np.unique(self.gxg.T[0]['hyear'])
        t = (pd.Timestamp('{}-{:02d}-{:02d}'.format(hyears[ 0], 3, 14)),
             pd.Timestamp('{}-{:02d}-{:02d}'.format(hyears[-1], 2, 28)))

        lw = 0.5 # line width plot pen
        for iclr, ip in enumerate(parcels):
            clr = clrs[iclr % len(clrs)]
            ax.hlines(self.GXG['GVG'][self.GXG['id']==ip], *t, clr,
                      ls='solid'  , lw=lw, label='GVG parcel {}'.format(ip))
            ax.hlines(self.GXG['GLG'][self.GXG['id']==ip], *t, clr,
                      ls='dashed' , lw=lw, label='GLG parcel {}'.format(ip))
            ax.hlines(self.GXG['GHG'][self.GXG['id']==ip], *t, clr,
                      ls='dashdot', lw=lw, label='GHG parcel {}'.format(ip))

        ax.legend(loc='best')
        return ax
    

def show_boundary_locations(lbls=None, CBC=None, iper=0, size_inches=(10,8.5)):
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
    size_inches: 2 tuple of floats
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

        titles=['Top layer, lbl={}, iper={}'.format(lbl, iper),
                'Bottom layer, lbl={}, iper={}'.format(lbl, iper)]
        ax = newfig2(titles=titles, xlabel='column', ylabels=['row', 'row'],
                     sharx=True, sharey=True, size_inches=size_inches)
        ax[0].spy(IB[0], marker='.', markersize=2)
        ax[1].spy(IB[1], marker='.', markersize=2)


class Watbal_obj:
    """Water balance object."""

    def __init__(self, sim=None, gr=None):
        """Return Watbal object carrying the water budget for all cross sections in m/d.

        Parameters
        ----------
        gr: fdm_tools.mfgrid.Grid Object
            holds the structured model grid
        sim: flopy sim object
            simulation object, model = model.get_model(sim.model_name) object
            holds the simulation with access to all its models and packages

        Generates
        ---------
        self.W : np.recarray having labels ['RCH', 'EVT', 'GHB etc'] and values that
            are np.ndarrays with shape (1, nlay, nrow, nper)
            Each row has the cross sectional discharge in m/d.

        Note that labels must be adapted if new CBC packages are to be included.

        @TO 20170823 as function
        @TO 20200907 turned Watbal into class Wabal_Obj
        @TO 20231023 update to mf6, changed W to dict for easier understanding
        """
        required_modflow_packages = ['RCH', 'EVT', 'WEL', 'GHB', 'RIV', 'DRN',
                                     'FLF', 'STO-SS', 'STO-SY']
        self.labels = required_modflow_packages #Layer 0 labels

        self.gr = gr
        
        model = sim.get_model(list(sim.model_names)[0])
        self.CBC = model.output.budget()
        
        # To get the CBC values in m/d for the entire region represented in the database:
        # Because each modeled parcel is only a cross of the real parcel, the model
        # does not know the parcel's true, only the database knows it.
        # To get the CBC values in m/d for the entire model:
        # First sum the model's CBC values over the parcel, which yields m3/d, and divide
        # by the parcel's active model are, which yields the values in m/d for each parcel.
        # The to get the parcel's contribution to the total CBC  values in m/d for the
        # entire modeled region, divide by the size of the model's region and multiply bye
        # The parcel's true area. The region's true area is the sum over all parcels true area.
        
        # Multiply array with active cells. Active is 1 and inactive is 0
        active = ml.get_package('DIS').idomain.data
        active[active !=0 ] = 1

        # Area of each cross section
        A_xsec = np.sum(gr.Dx * gr.Dy * active[0], axis=-1)
         
        # Contribution of each true parcel area to total of the modeled region
        Arel = (parcel_data['A_parcel'] / parcel_data['A_parcel'].sum()).values
                
        self.W = dict() # keys are watbal labels,  values are float (nlay, nrow, nper)
        _vals3D = np.zeros(gr.shape).ravel()
        print('Setting up Watbal_Obj ...')
        for lbl in self.labels:
            self.W[lbl[:3]] = np.zeros((self.CBC.nlay, self.CBC.nrow, self.CBC.nper))
            print(lbl, end='')
            
            # For these labels, the packages of which CBC values are obtained
            # in a recarray with 'node' and 'q' with:
            # dtype = [('node', '<i4'), ('node2', '<i4'), ('q', '<f8')])
            if lbl in ['WEL', 'GHB', 'RIV', 'DRN', 'RCH', 'EVT']:
                
                vals = self.CBC.get_data(text=lbl) # as recarray with given dtype
                
                for iper in range(self.CBC.nper):
                    # Vals contains all non-zero node values at once
                    Ivals = vals[iper]['node'] - 1 # zero-based Modflow node numbers
                    Q     = vals[iper]['q']
                    
                    _vals3D[:] = 0. # reset _vals3D to zeros for next loop cycle
                    _vals3D[Ivals] = Q
                    
                    # Sum _vals3D over the columns of the model strutured grid
                    self.W[lbl][:, :, iper] = np.sum(_vals3D.reshape(gr.shape), axis=-1)
                    
                    # Show progress
                    if iper % 100 == 0: print('.',end='')                
            elif lbl[:3] == 'STO': # in ['STO-SY', 'STO-SS']:
                vals = self.CBC.get_data(text=lbl) # as a full 3D arrwy over the grid
                for iper in range(self.CBC.nper):
                    self.W[lbl[:3]][:, :, iper] += np.sum(vals[iper], axis=-1) # W['STO']
                    if iper % 100 == 0: print('.', end='') # Show progress                
            elif lbl in ['FLF']:
                flowja=self.CBC.get_data(text='FLOW-JA-FACE')
                grb_file = sim.name + '.dis.grb'
                
                # Get flf as a full 2D array over a grid layer
                flfs = get_flow_lower_face(flowja=flowja, grb_file=grb_file)
                for iper in range(self.CBC.nper):
                    self.W[lbl][-1, :, iper] = +flfs[iper][0].sum(axis=-1)
                    self.W[lbl][ 0, :, iper] = -flfs[iper][0].sum(axis=-1)
            else:
                raise ValueError("Illegal label {}".format(lbl))
            print('Last iper for {}: {}'.format(lbl, iper))
        for lbl in self.W: # prevents clash with STO-SS and STO-SY            
            self.W[lbl] /= A_xsec[np.newaxis, : , np.newaxis] # from m3/d to m/d # not yet summed over all parcels!

        print('Done. See self.CBC and self.W')


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
            ax = newfig2(titles=(
                    'Water balance top layer. '   + ttl,
                    'Water balance botom layer. ' + ttl), xlabel='time', ylabels=['mm/d', 'mm/d'],
                         size_inches=(14, 8), sharey=False, sharex=True)

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

def generate_model_dict(dirs=None, parcel_data=None, tdata=None):
    """Return default and case parameters for all packages."""
    
    gr = grid_from_parcel_data(parcel_data=parcel_data, dx=user_settings['dx'])
    
    always_packages = {'sim', 'gwf', 'dis', 'ic', 'ims', 'npf', 'oc', 'sto', 'tdis'}
    use_packages = [p.lower() for p in user_settings['use_packages']]
    use_packages = always_packages.union(use_packages)

    model_dict = get_mf6_params_from_excel(dirs.mf_parameters)
    
    fp_packages = dict()
    
    ### sim =================================
    if True:    # Simulation
        logging.info("sim")

        model_dict['sim'].update(sim_name=user_settings['case_name'], sim_ws=dirs.wd, exe_name=dirs.exe)
        sim = flopy.mf6.MFSimulation(**model_dict['sim'])
        fp_packages['sim'] = sim


    ### tdis ===============================
    if True:     # Stress period parameters.
        logging.info("tdis")
        start_date_time = str(tdata.index[0])

        nper, nstep, tsmult = len(tdata), 1, 1.0

        if nper == 1:
            dt = np.array([1.])
        else:
            # dt in days, where dt may vary over time.
            dt = np.diff(tdata.index - tdata.index[0]) / np.timedelta64(1, 'D')

            # dt of day 0: assume zero is t[1] - t[0])
            dt = np.hstack((dt[0], dt))

        perioddata = [[sp_time, nstep, tsmult] for sp_time in dt]

        model_dict['tdis'].update(  nper=nper,
                                    perioddata=perioddata,
                                    start_date_time=start_date_time)
        fp_packages['tdis'] = flopy.mf6.ModflowTdis(sim, **model_dict['tdis'])

    ### ims ===============================
    if True: # Solver, geneally use defaults from workbook.
        # model_dict['ims'].update(...=..., )
        logging.info("ims")
        fp_packages['ims'] = flopy.mf6.ModflowIms(sim, **model_dict['ims'])
        
    ### gwf ===============================
    if True: # Groundwater flow model, use name from 'sim').
        logging.info("gwf")
        model_dict['gwf'].update(modelname=model_dict['sim']['sim_name'],
                         exe_name=model_dict['sim']['exe_name'])
        gwf = flopy.mf6.ModflowGwf(sim, **model_dict['gwf'])
        fp_packages['gwf'] = gwf
    

    ### dis ==============================
    if 'dis' in use_packages: # Structured grid like in old MODFLOW
        logging.info("dis")
        idomain = np.ones(gr.shape)
        for iy, b in zip(range(gr.ny), parcel_data['b']):
            idomain[:, iy, gr.Xm[iy] > b] = 0

        model_dict['dis'].update(nlay=gr.nlay, nrow=gr.nrow, ncol=gr.ncol,
                                 delr=gr.dx, delc=gr.dy,
                                top=gr.Z[0], botm=gr.Z[1:],
                                idomain=idomain,
        )
        fp_packages['dis'] = flopy.mf6.ModflowGwfdis(gwf, **model_dict['dis'])

    ### sto ================================
    if 'sto' in use_packages:
        logging.info("sto")
        laytype = np.zeros(gr.nlay); laytype[0] = 1

        # Elastic storage coefficients
        S = set3D(parcel_data[['sy', 'S2']], gr.shape)
        S[[0, 1], :] = 0. # No elastic storage above regional aquifer
        Sy = set3D(parcel_data[['sy', 'sy']], gr.shape)
        # TODO make Sy for regional aquifer explicit

        model_dict['sto'].update(
            save_flows=True,
            storagecoefficient=True, # Interprete as S instead of Ss if confined
            # ss_confined_only=True,   # elastic only when confined
            iconvert=gr.const(laytype), # convertible cell indices # TODO how define ?
            ss=S,
            sy=Sy,
            transient={i: True for i in range(nper)}
           )
        fp_packages['sto'] = flopy.mf6.ModflowGwfsto(gwf, **model_dict['sto'])
    
    ### npf =================================
    if 'npf' in use_packages: # Node properties
        # Add resistance between cover layer and regional aquifer to parcel_data DataFrame
        logging.info("npf")
        parcel_data['kc'] = parcel_data['D_CB'] / parcel_data['c_CB']

        model_dict['npf'].update(
            k=set3D(parcel_data[['kh', 'kh2']], gr.shape),
            k22=1.e-20, # no flow along y=axis
            k33=set3D(parcel_data[['kv', 'kv2']], gr.shape),
        )
        fp_packages['npf'] = flopy.mf6.ModflowGwfnpf(gwf, **model_dict['npf'])

    ### ic ==================================
    if 'ic' in use_packages: # Initial conditions (Initial heads)
        model_dict['ic'].update(strt=set3D(parcel_data['h_winter'], gr.shape))
        fp_packages['ic'] = flopy.mf6.ModflowGwfic(gwf, **model_dict['ic'])

    #### Boundary conditions ####

    ### rcha(recharge)     
    if 'rcha' in use_packages:
        logging.info("Evalutation rcha.")
        print("Evalutation rcha.")
        rch_spd = {isp: tdata['RH'].iloc[isp] for isp in range(len(tdata))}
        model_dict['rcha'].update(readasarrays=True, print_input=True, recharge=rch_spd)
        fp_packages['rcha'] = flopy.mf6.ModflowGwfrcha(gwf, **model_dict['rcha'])
    
    ### evta ================================
    if 'evta' in use_packages: # Evapotranspiration    
        model_dict['evta'].update(readasarrays=True,
            ievt=None,
            surface={0: set3D(parcel_data['AHN'] - parcel_data['ET_surfd'], shape=gr.shape)[0]},
            depth={0: set3D(parcel_data['ET_exdp'], shape=gr.shape)[0]},
            rate={isp: tdata['EV24'].iloc[isp] for isp in range(len(tdata))},
        )
        logging.info("Evalutation evta.")
        print("Evalutation evta ... some patience please")
        fp_packages['evta'] = flopy.mf6.ModflowGwfevta(gwf, **model_dict['evta'])

    active = model_dict['dis']['idomain'] > 0
    
    ### Prepare boundary conditions ==========
    mon = np.array([dt.month for dt in tdata.index])
    day = np.array([dt.day for dt in tdata.index])
    
    Isp_start_summer = np.where(np.logical_and(mon ==  4, day == 1))[0]
    Isp_start_winter = np.where(np.logical_and(mon == 10, day == 1))[0]
        
    pdata = parcel_data

    ### wel ===============Is used to model given seepage from regional aquifer
    if 'wel' in use_packages:
        # (Used for given seepage from second aquifer)
        # qup is still given per parcel and not as a function of time !!
        # The cells are those in the second layer where idomain > 0
        logging.info("wells")
        mask = model_dict['dis']['idomain'][-1] > 0
        LRC_wel = gr.LRC(gr.NOD[-1, :, :][mask], astuples=True)

        Q = (pdata['q_up'].values[:, np.newaxis] * np.ones(gr.shape[1:]) * gr.Area)[mask]
        spd = [(lrc, _Q) for lrc, _Q in zip(LRC_wel, Q)]
        
        # Only the first SP needs data, as long as seepage is constant.
        # TODO: this likely changes in the future, to monthly seepage values.
        # TODO: Specify how monthly seepage values for all parcels are inported.
        # TODO: This import may require an extended data file originating from a regional model.
        # TODO: its dimensions should be (nparcel x months in time series)
        stress_period_data = {0: spd}
        model_dict['wel'].update(stress_period_data=stress_period_data, maxbound=len(LRC_wel))
        fp_packages['wel'] = flopy.mf6.ModflowGwfwel(gwf, **model_dict['wel'])
    

    ### ghb ================Is used for flow from and toward ditches
    if 'ghb' in use_packages:
        # First cell of top and bottom layer (always). Conduction of bottom layer may be zero.
        LRC_ghb = gr.LRC(gr.NOD[0, :, [0, -1]][active[0, :, [0, -1]]], astuples=True) # GHB cells

        # Get conductance for the GHB (connection soil Ditch in top and bottom layer)
        # i.e. one value per parcel in the top layer and one value in the bottom layer.

        # Use analytic ditch resistance with layer thickness and no partial penetration
        logging.info("ghb")
        condGHB = get_cond_GHB(pdata=parcel_data, gr=gr)

        # Get the GHB heads for summer and winter for these ditches
        hw = np.vstack((pdata['h_winter'], pdata['h_winter']))
        hs = np.vstack((pdata['h_summer'], pdata['h_summer']))

        ghb_winter = [(lrc, head, cond) for lrc, head, cond in
                                zip(LRC_ghb, hw.ravel(), condGHB.ravel())]
        ghb_summer = [(lrc, head, cond) for lrc, head, cond in
                                zip(LRC_ghb, hs.ravel(), condGHB.ravel())]

        # Input for the first stress period
        stress_period_data = {0: ghb_summer if tdata.iloc[0]['summer'] == True else ghb_winter}

        # Only generates input when summer changes to wintor or vice versa
        for isp in Isp_start_summer:
            stress_period_data[isp] = ghb_summer
        for isp in Isp_start_winter:
            stress_period_data[isp] = ghb_winter

        model_dict['ghb'].update(stress_period_data=stress_period_data, maxbound=len(LRC_ghb))
        fp_packages['ghb'] = flopy.mf6.ModflowGwfghb(gwf, **model_dict['ghb'])

     
    ### riv ==============Is used for extra flow toward ditch (lower flow to than from)
    if 'riv' in use_packages:
        # TODO used identifyable numers in idomain for layer use in budget file.
        # Riv cells are the same as GHB cells
        logging.info("riv")
        LRC_riv = gr.LRC(gr.NOD[0, :, [0, -1]][active[0, :, [0, -1]]], astuples=True) # RIV cells

        condRIV = get_RIV_Cond(pdata=parcel_data, gr=gr)

        riv_winter = [(lic, stage, cond, rbot) for lic, stage, cond, rbot in
                                zip(LRC_riv, hw.ravel(), condRIV.ravel(), hw.ravel())]
        riv_summer = [(lic, stage, cond, rbot) for lic, stage, cond, rbot in
                                zip(LRC_riv, hs.ravel(), condRIV.ravel(), hs.ravel())]

        # Input for the first stress period
        stress_period_data = {0: riv_summer if tdata.iloc[0]['summer'] == True else riv_winter}

        # Only generates input when summer changes to winter and vice versa.
        for isp in Isp_start_summer:
            stress_period_data[isp] = riv_summer
        for isp in Isp_start_winter:
            stress_period_data[isp] = riv_winter

        model_dict['riv'].update(stress_period_data=stress_period_data, maxbound=len(LRC_riv))
        fp_packages['riv'] = flopy.mf6.ModflowGwfriv(gwf, **model_dict['riv'])


    ### drn ==============Is used for drains, surface runoff and trenches
    if 'drn' in use_packages: # For trenches and surface runoff by means of drains.
        # Drains on top layer, but not in first column and not in inactive cells
        logging.info("drn")
        I  = gr.NOD[0, :, 1:][active[0, :, 1:]]
        LRC_drn = gr.LRC(I, astuples=True)

        # Get drain elevation and drain conductance as Ny * Nx array (top layer)
        elevation = get_drain_elev_with_trenches(pdata=parcel_data, gr=gr).ravel()[I]
        condDRN = get_cond_DRN(pdata=parcel_data, gr=gr).ravel()[I]

        # Stress period data
        #dtype = flopy.modflow.ModflowDrn.get_default_dtype()
        #spd = np.recarray(gr.nod, dtype=dtype)
        spd = [(lrc, elev, cond) for lrc, elev, cond in zip(LRC_drn, elevation, condDRN)]

        # Input is required only for the first SP, because drain data are constant.
        stress_period_data = {0: spd}
        model_dict['drn'].update(stress_period_data=stress_period_data, maxbound=len(LRC_drn))
        fp_packages['drn'] = flopy.mf6.ModflowGwfdrn(gwf, **model_dict['drn'])
    
    ### oc (output control)
    if 'oc' in use_packages:
        logging.info("oc")
        model_dict['oc'].update(head_filerecord="{}.hds".format(model_dict['sim']['sim_name']),
                         budget_filerecord="{}.cbc".format(model_dict['sim']['sim_name']),
                         saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])
        fp_packages['oc'] = flopy.mf6.ModflowGwfoc(gwf, **model_dict['oc'])

    return gr, fp_packages, model_dict, use_packages


def run_modflow(dirs=None, parcel_data=None, tdata=None):
        """Simulate GGOR using MODFLOW.

        Parameters
        ----------
        dirs: DirStruct object
            directory structure object, containing home and case information
        parcel_data: pd.DataFrame
            parcel data. Each row representing one parcel.
        tdata: pd.DataFrame
            time data with columns ['RH', 'EV24', 'summer', 'hyear', 'hand']
            and index datetime.
        dx: float
            uniform cell width of the model (x-direction, perpendicular to long axis of parcel)
        """
        # Build the model dict and flopy packages
        gr, fp_packages, model_dict, use_packages= generate_model_dict(dirs=dirs,
                            parcel_data=parcel_data, tdata=tdata)
    
        sim = fp_packages['sim']
        # Write simulation
        sim.write_simulation(silent=user_settings['silent'])
        
        # Run simulation
        success, buff = sim.run_simulation(silent=user_settings['silent'])
        

        print('Running success = {}'.format(success))
        if not success:
            print(buff)
            logging.critical("Buffer printed because MODFLOW did not terminate normally.")
            raise Exception('MODFLOW did not terminate normally.')

        return gr, fp_packages, model_dict, use_packages


def data_to_excel(dirs=None, data=None, fname=None):
    """Save the current parcel_data DataFrame to an excel file.

    Parameters
    ----------
    dirs: Dir_struct object
        directory structure current case
    data: pd.DataFrame
        data to be written to excel
    fname: str
        name of local excel file to be saved
    """
    fname = os.path.join(dirs.case, os.path.splitext(fname)[0] + '.xlsx')
    parcel_data.to_excel(fname, index=False, engine='openpyxl')
    print("parcel_data saved to '{}'".format(fname))

    return None


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


def save_parcel_data_to_excel(dirs,
                start='pdata_start.xlsx', end='pdata_end.xlsx'):
    """Save the parcel_data to Excel to generate test_data.

    Parameters
    ----------
    dirs: DIr_struct
        GGOR directory structure.
    start: str
        workbook name to store the initial parcel_data
    end: str
        workbook name to store the final/current parcel_dadta
    """
    pdata_start = data_from_dbffile(os.path.join(dirs.case, case + '.dbf'))

    data_to_excel(dirs, pdata_start, start=start)
    data_to_excel(dirs, parcel_data, end=end)

#%% main

if __name__ == "__main__":
    
    logging.warning(sys.version)
    logging.warning(sys.executable)

    test=False

    # Parameters to generate the model. Well use this as **kwargs
    GGOR_home = os.path.expanduser('~/GRWMODELS/python/GGOR') # home directory
    case = 'AAN_GZK'

    #GGOR directory structure
    dirs = Dir_struct(GGOR_home, case=case)
    
    logging.warning("cwd = {}".format(os.getcwd()))

    mf_parameters_wbk = os.path.join(dirs.mf_parameters, 'mf_parameters.xlsx')

    #Get the meteo data from an existing file or directly from KNMI
    meteo_data = knmi.get_weather(stn=240, start='20100101', end='20191231',
                                  folder=dirs.meteo)

    # Add columns "summer' and "hyyear" to it"
    tdata = handle_meteo_data(meteo_data, summer_start=4, summer_end=10)
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
        parcel_data = GGOR_data(defaults=defaults, bofek=bofek, BMINMAX=(5, 250),
                GGOR_home=GGOR_home, case=case).data
    
    # MODFLOW input arrays are int the returned dicts
    if False: # Run Modflow
        gr, fp_packages, model_dict, use_packages = run_modflow(
                    dirs=dirs, parcel_data=parcel_data, tdata=tdata)  
    else:
        gr = grid_from_parcel_data(parcel_data=parcel_data, dx=user_settings['dx']) 

    sim = flopy.mf6.MFSimulation.load(user_settings['case_name'])
    ml = sim.get_model(list(sim.model_names)[0])


    #%% Get the modflow-computed heads and cell by cell flows

    heads_obj = Heads_obj(sim=sim, tdata=tdata, gr=gr)
    
    budObj = ml.output.budget()

    watbal = Watbal_obj(sim=sim, gr=gr)

    #%% Open HDS file and plot heads (with GXG)
    if test:
        for tst in set(parcel_data['Test'].values):
            parcels = list(parcel_data.index[parcel_data['Test'] == tst])
            test_vals_str = '{}'.format(', '.join(
                [str(tv) for tv in parcel_data[tst].iloc[parcels]]))

            titles=['Parcel averaged heads, testvar {} in [{}]'.format(tst, test_vals_str),
                    'Parcel averaged heads, testvar {} in [{}]'.format(tst, test_vals_str)]
            ax = heads_obj.plot(tdata=tdata,
                           parcel_data=parcel_data,
                           parcels=parcels,
                           titles=titles,
                           size_inches=(14, 8))
            if False:
                ax = watbal.plot(parcel_data=parcel_data,
                                 tdata=tdata,
                                 parcels=parcels[0],
                                 sharey=True)
    else:
        parcels = [0, 1, 2, 3]
        titles=['Parcel averaged heads', 'Parcel averaged heads']
        ax = heads_obj.plot(tdata=tdata,
               parcel_data=parcel_data,
               parcels=parcels,
               titles=titles,
               size_inches=(14, 8))
        ax = watbal.plot(parcel_data=parcel_data,
                         tdata=tdata,
                         parcels=parcels,   # over all parcels
                         sharey=True)

#%%
    ax = watbal.plot(parcel_data=parcel_data,
                         tdata=tdata,
                         parcels=None,   # over all parcels
                         sharey=True)

    plt.show()

#%%
    print('---- All done ! ----')
