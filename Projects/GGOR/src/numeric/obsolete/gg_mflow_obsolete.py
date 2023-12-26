# -*- coding: utf-8 -*-
"""Simultate GGOR for nparcels using MODFLOW.

The parcel (area) data are in a shape file.
The attributes are in the dbf
The dbd file is read into pd.DataFrame.
The meteo is reada from existing file or obtained from KNMI site.
Modflow is used to simultaneously simulate all the parcels dynamically.
The results are shown for selected paracels (hds, GXG)
The running water budget is shown for all the parcels combined.

@ TO 2020-09-06
"""
import os
import numpy as np
import pandas as pd
import shapefile
import matplotlib.pyplot as plt
from KNMI import knmi
from fdm import mfgrid
import flopy
import flopy.utils.binaryfile as bf
from collections import OrderedDict
import pdb
NOT = np.logical_not
AND = np.logical_and
OR  = np.logical_or

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


#% Defaults used for required parameters but are not in the databse.
# They will be used when not in the database.
defaults = {'d_drain': 0,  # [m] Tile drainage depth below local ground elevation, may be zero if no drains.
            'd_trench': 0.3, # [m] Trench depth in case present.
            'c_drain': 5., # [d] Tile drainage areal resistance. Also used for trenches.
            'wi_ditch' : 2.,  # [d] Ditch resistance when flow is from ditch to ground. (analytical)
            'wo_ditch' : 1.,  # [d] Ditch resistance when flow is from ground to ditch. (anlaytical)
            'ci_ditch' : 2.,  # [d] Ditch bottom and side entry resistance (applied to Omega)
            'co_ditch' : 1.,  # [d] Ditch bottom and side entry resistance (applied to Omega)
            'd_ditch' : 1.0, #[m] depth of ditch below ground surface
            'b_ditch' : 0.75, # [m] half the width of the ditch
            'D_CB' : 0.1, # [m] (dummy dikte) van basisveenlaag (CB=confining bed)
            'D2' : 40., # [m] dikte regionale aquifer
            'S2' : 1e-3,# [-] total (elastic) storage regional aquifer
            'kh2': 30., # [m/d]  horizontal condutivity regional aquifer
            'kv2':  6., #[ m/d]  vertical conductivity regional aquifer
            'ET_surfd': 1.0, # [m] depth of surf in ET below ground surface.
            'ET_exdp': 2.5, # [m] Modflow's extinction depth (see ET package)
            }


#% Modflow cell-by-cell flow labels: Translates short labels to those in the CBC file.
cbc_labels = {
    'STO': 'STORAGE',
    'FLF': 'FLOW LOWER FACE ',
    'WEL': 'WELLS',              # used for seepage. Are in seconde layer.
    'EVT': 'ET',
    'GHB': 'HEAD DEP BOUNDS',
    'RIV': 'RIVER LEAKAGE',
    'DRN': 'DRAINS',             # also used for trenches.
    'RCH': 'RECHARGE',
    }

# For legend when plotting water balance
#leg is legend for this label in the graph
#clr is the color of the filled graph
watbal_label = OrderedDict()
watbal_label.update(
                {'RCH': {'leg': 'RCH', 'clr': 'green'},
                'EVT': {'leg': 'EVT', 'clr': 'gold'},
                'WEL': {'leg': 'WEL(in wvp2)' , 'clr': 'blue'},
                'DRN': {'leg': 'DRN(drn/trenches/runoff)', 'clr': 'lavender'},
                'RIV': {'leg': 'RIV(ditch out)', 'clr': 'magenta'},
                'GHB': {'leg': 'GHB(ditch in+out)', 'clr': 'indigo'},
                'FLF': {'leg': 'FLF(leakage)', 'clr': 'gray'},
                'STO': {'leg': 'STO', 'clr': 'cyan'}})


def selection_check(selection):
    """Return verified selection."""
    if selection is None:
        selection = slice(0, None, 1)
    elif isinstance(selection, int):
        selection = slice(0, selection, 1)
    elif isinstance(selection, (tuple, list, np.ndarray, range, slice)):
        pass
    else:
        raise ValueError("selection must be None (for all), an int or a sequence")
    return selection


def gen_testdata(tdata, **kwargs):
    """Return copy of tdata with altered or added input columns for testing.

    The tuples consist of:
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
        each kwarg as a key and a tuple as argument. The tuple consists
        of a number_of_days followd by values. Like (120, 0.5, 0.7, -03).
        The first value is the duration during which each of the following
        values will be used in turn, repeated after the last value has been
        consumed.
        Due to Python works, the call can be like so
        gen_testdata(tdata, RH=(200, 0.002, 0.003, 0.02), hLR=(150, 0.7, -0.3),
                     h1=(300, -0.4, 0.2, -0.1))
        The arguments will then be available in kwargs like so
            kwargs={'RH': (200, 0.002, 0.003, 0.02),
                    'hLR': (150, 0.7, -0.3),
                    'h1': (300, -0.4, 0.2, 0.02)}
        This makes this function extremely flexable to generate testdata.
    Example
    -------
    tdatanew = gen_testdata(tdata, RH=(200, 0.02, 0.01, 0.03, 0.01)),
                       EV24 = (365, 0., -0.001)
                       q_up=(150, -0.001, 0.001, 0, -0.003)
                       )
    tdatanew is then a pd.DataFrame with columns 'RH', 'EV24', 'q_up' next
        to other columns that were already in tdata.
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

    GGOR_home/ # GGOR_home is the GGOE_home directory specified above.
            bin/
                mfusg_X64.exe
                mfusg.mac
            src/
                 analytic/
                 numeric/
            data/
                 bofek/
                 meteo/
                 spatial/
                         AAN_GZK
                         ....
            doc/
            notebooks/
            cases/
                  AAN_GZK
                  ....
    ../python/KNMI
    """

    def __init__(self, home='.', case=None ):
        """Generate GGOR directory structure.

        Parameters
        ----------
        home: str
            path to GGOR home directory
        case: str
            the name of the current case
        """
        self.home = os.path.abspath(os.path.expanduser(home))
        self.case_name = case

        #verify existance of required files
        #exe = self.exe_name
        #assert os.path.isfile(exe), "Missing executable '{}'.".format(exe)
        #dbf = os.path.join(self.case, self.case_name + '.dbf')
        #assert os.path.isfile(dbf), "Missing dbase file '{}.'".format(dbf)
        #assert os.path.isdir(self.meteo), "Missing meteodir '{}'.".format(self.meteo)
        #assert os.path.isdir(self.bofek), "Missing bofek folder '{}'.".format(self.bofek)

    # Directory structure
    @property
    def bin(self):
        """Yield bindary folder."""
        return os.path.join(self.home, 'bin')
    @property
    def src(self):
        """Yield source code folder."""
        return os.path.join(self.home, 'src')
    @property
    def data(self):
        """Yield data folder."""
        return os.path.join(self.home, 'data')
    @property
    def cases(self):
        """Yield folder where cases are stored."""
        return os.path.join(self.home, 'cases')
    @property
    def meteo(self):
        """Yield meteo data folder."""
        return os.path.join(self.data, 'meteo')
    @property
    def spatial(self):
        """Yield folder with spatial data.

        Each case corresponds with a folder with the case name.
        """
        return os.path.join(self.data, 'spatial')
    @property
    def bofek(self):
        """Return directory of bofek units."""
        return os.path.join(self.data, 'bofek')
    @property
    def case(self):
        """Return results directory of current case."""
        return os.path.join(self.spatial, self.case_name)
    @property
    def case_results(self):
        """Return folder where case output goes."""
        return self.wd
    @property
    def wd(self):
        """Yield working directory (MODFLOW output) depending on case."""
        if not os.path.isdir(self.cases): os.mkdir(self.cases)
        wd = os.path.join(self.cases, self.case_name)
        if not os.path.isdir(wd): os.mkdir(wd)
        return wd
    def cwd(self):
        """Change to current working directory."""
        os.chdir(self.wd)
    @property
    def exe_name(self):
        """Yield name of code (exe) depending on operating system."""
        if os.name == 'posix':
            exe_path = os.path.join(self.bin, 'mfusg.mac')
        else:
            exe_path = os.path.join(self.bin, 'mfusg_64.exe')
        return exe_path


def handle_meteo_data(meteo_data=None, summer_start=4, summer_end=9):
    """Set and store meteo data and adds columns summer, hyear and hand.

    Added columns are
        summer: bool
        hyear: hyddrological year. They atart at March 14 to get proper GVG!! Don't change
            GVG will be the mean of the values on 14/3, 28/3 and 14/4 of each year'
        hand: groundwater measurement is done at date (14th and 28th of each month)

    Parameters
    ----------
    meteo: pd.DataFrame
        meteo data, with timestamp index and columns 'RH' and 'EVT24'
        for precipitation and (Makkink) evaporation respectively [m/d]
    summer_start: int
        month coniciding with start of summer (hydrologically). Default 4.
    summer_end: int
        month coinciding with end of summer (hydrologically). Default 10.
    """
    dcol = {'RH', 'EVT24'}.difference(meteo_data.columns)
    if not dcol:
        pass
    else:
        KeyError("Missing column[s] [{}] in meteo DataFrame".format(', '.join(dcol)))

    #verify, data are in m/d
    if not meteo_data['RH'].median() < 0.01 and meteo_data['EVT24'].median() < 0.01:
        AssertionError("Median of Precipitration = {:5g} and median evapotranspiration = {:4g}\n"
                          .format(meteo_data['RH'].median(), meteo_data['EVT24'].median()) +
                       "Percipication and or evapotranspiration likely not in m/d!")

    # Add boolean column indicating summer (needed for summer and winter levels)
    meteo_data['summer'] = [True if t.month >= summer_start and t.month < summer_end
                               else False for t in meteo_data.index]

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
        table of parcel property datd
    dx: float
        column width (uniform)
    LAYCBD: list (or sequence) of length 2 (nlayer)
        value of 1 for each layer with confining unit below else value of 0
    """
    bmax = parcel_data['b'].max()
    nx = int(np.ceil(bmax) // dx)
    ny = len(parcel_data)
    nz = 2

    Dx  = np.ones(nx) * dx
    xGr = np.hstack((0, np.cumsum(Dx)))
    yGr = np.arange(ny, -1, -1, dtype=float)

    laycbd = list(laycbd)
    while len(laycbd) < nz: laycbd.append(0.)
    laycbd = np.array(laycbd); laycbd[-1] = 0.

    Z = np.zeros((nz + 1 + np.sum(laycbd), ny, nx))

    Z[0] = parcel_data['AHN'].values[:, np.newaxis] * np.ones((1, nx))
    Z[1] = Z[0] - parcel_data['D1'  ].values[:, np.newaxis] * np.ones((1, nx))
    Z[2] = Z[1] - parcel_data['D_CB'].values[:, np.newaxis] * np.ones((1, nx))
    Z[3] = Z[2] - parcel_data['D2'  ].values[:, np.newaxis] * np.ones((1, nx))

    return mfgrid.Grid(xGr, yGr, Z, LAYCBD=laycbd)


def set3D(layvals, shape=None):
    """Return 3D array with given shape from parcel_data.

    Parameters
    ----------
        layvals: pd.DataFrame (or pd.Series)
            The layer values to be used
        shape: tuple
            the shape of the array to be generated
    """
    vals = layvals.values
    if len(vals.shape) == 1:
        return vals[np.newaxis, :, np.newaxis] * np.ones(shape)
    else:
        return vals.T[:, :, np.newaxis] * np.ones(shape)


def set_spatial_arrays(parcel_data=None, gr=None):
    """Crerate teh spatial arrays for MODFLOW.

    Parameters
    ----------
    parcel_data: pd.DataFrame
        parcel data
    """
    # A string looks up data in the database
    # a number uses this number for all parcels
    shape = gr.shape

    sparr = {}

    sparr['HK']     = set3D(parcel_data[['kh', 'kh2']], shape)
    sparr['VK']     = set3D(parcel_data[['kv', 'kv2']], shape)
    sparr['STRTHD'] = set3D(parcel_data['h_winter'], shape)

    # Use of elastic storage coefficients
    # Because
    sparr['laytyp']=0
    # SY is not used (but generated incase we want to use varying thickness)
    # Using Sy2=S2 makes sure we always have elastic S in the reginal aquifer.
    sparr['SY']     = set3D(parcel_data[['sy', 'S2']], shape) # not used due to LAYTYP=0
    # Ss is used, so we should set ss = sy/D1. However we also set
    sparr['storagecoefficient']=True # in LPF, which interprets ss as s total.
    # Hence, we use 'sy' for 'ss' and 'S2' for 'SS2'
    sparr['SS']     = set3D(parcel_data[['sy', 'S2']], shape)

    # for ET package we need SURF and EXPD
    sparr['SURF']   = set3D(parcel_data['AHN'] - parcel_data['ET_surfd'], shape=(1, *shape[1:]))
    sparr['EXDP']   = set3D(parcel_data['ET_exdp'], shape=(1, *shape[1:]))

    # VKCB for the resistance at bottom of cover layer
    # any resistance inside cover layer stems from kv
    sparr['VKCB'] = np.zeros(gr.shape)
    for ilay, (lcb, icb) in enumerate(zip(gr.LAYCBD, gr.ICBD)):
        if lcb:
            c = set3D(parcel_data['c_CB'], (gr.ncbd, *shape[1:]))
            sparr['VKCB'][ilay] = gr.dz[icb] / c

    # IBOUND:
    # Limit width of rows to width of parcels by making cells beyond parcel width inactive
    # First get index of first cell beyond parcel width
    Ib = np.asarray(np.ceil(np.interp(parcel_data['b'], [gr.x[0], gr.x[-1]], [0, gr.nx + 1])), dtype=int)
    # Then set all cells beyond b[iy] to inactive
    sparr['IBOUND'] = gr.const(1, dtype=int)
    for iy, ix in enumerate(Ib):
        sparr['IBOUND'][:, iy, ix:] = 0 # inactive

    return sparr


def set_stress_period_data(tdata=None):
    """Create stress_period arrays for MODLOW.

    Parameters
    ----------
    tdata: pd.DataFrame
        time data with columns 'RH'', 'EVT24' and 'summer'
        and whose index are timestamps
    """
    spt = dict() # stress period time ddata

    if len(tdata) == 1:
        dt = np.array([1.])
    else:
        dt = np.diff(tdata.index - tdata.index[0]) / np.timedelta64(1, 'D')
        dt = np.hstack((dt[0], dt))

    spt['NPER']   = len(tdata)
    spt['PERLEN'] = dt # days for MODFLOW (assume dt day zero is t[1] - t[0])
    spt['NSTP']   = np.ones(spt['NPER'], dtype=int)
    spt['STEADY'] = np.zeros(spt['NPER'], dtype=bool) if spt['NPER'] > 1\
                    else np.ones(spt['NPER'], dtype=bool) # all False

    # We need only to specify one value for each stress period for the model as a whole?
    spt['RECH'] = {isp: tdata['RH'  ].iloc[isp] for isp in range(spt['NPER'])}
    spt['EVTR'] = {isp: tdata['EV24'].iloc[isp] for isp in range(spt['NPER'])}

    return spt


def set_boundary_data(parcel_data=None, tdata=None, gr=None, IBOUND=None,
                      use_w_not_c=False):
    """Create time-dependent boundary arrays for MODFLOW.

    Parameters
    ----------
    tdata: pd.DataFrame
        time data with columns 'RH', 'EVTR' and 'summer'
        and whose index are timestamps
    """
    spb = dict()
    for what in ['GHB', 'RIV', 'DRN', 'WEL']:
        spb[what]  = set_boundary(what,
                                  parcel_data=parcel_data,
                                  tdata=tdata,
                                  gr=gr,
                                  IBOUND=IBOUND,
                                  use_w_not_c=use_w_not_c)

    spb['OC'] = {(isp, 0): ['save head', 'save budget', 'print budget']
                                             for isp in range(len(tdata))}

    return spb


def get_drain_elev_with_trenches(parcel_data=None, gr=None):
    """Retrun drain elevations corrected for trenches.

    Parcels with 1 or more trenches are assumed to have no tile drainage.
    Their drainage depth will be set to zero except for the trench locations,
    where the drainage depth is set to 'dtrench'.

    Paramters
    ---------
    parcel_data: pd.DataFrame
        parcel_data
    gr: mfgrid.Grid object
        object holding the mesh information
    Returns
    -------
    elev: np.ndarray of shape (nparcel, ncol)
        Updated drainage elevations. I.e. ground surface in parcels that have
        one or more trenches except on the trench locations, where the
        elevation is set to groundsurface minus the trench depth.
    """
    # Distance between two trenches or between ditch and trenches.
    elev = ((parcel_data['AHN'] - parcel_data['d_drain']).values[:, np.newaxis]
            * np.ones((1, gr.nx)))

    y = np.arange(gr.nx + 1)
    x = gr.x
    for ip, ((b, ntr, ahn, dtr), htr) in enumerate(
            zip(parcel_data[['b', 'n_trench', 'AHN', 'd_trench']].values, elev)):
        if ntr: # only if parcel has one or more trenches
            Ltr = 2 * b / (int(ntr) + 1) # dist. betw. trenches or trench and ditch
            xtr = Ltr * np.floor(np.arange(1, ntr / 2 + 1)) # trench locs x <= b
            Itr = np.asarray(np.interp(xtr, y, x), dtype=int) # column indices
            Itr[Itr >= gr.nx] = gr.nx - 1
            htr[:]   = ahn          # default --> ground surface
            htr[Itr] = ahn - dtr    # trench  --> ground surface - trench depth
    return elev


def set_boundary(what=None, parcel_data=None, tdata=None, gr=None, IBOUND=None,
                 use_w_not_c=False):
    """Return dictionary for boundary of given type.

    Parameters
    ----------
    what: str, one of
         'WEL': WEL package used to simulate vertical seepage.
         'DRN': DRN package used to simulate tile drainage (or surface runoff.
         'GHB': GHB package used to simulate ditch in- and outflow.
         'RIV': RIV package used to simulate ditch outflow (together with GHB).
    parcel_dadta: pd.DataFrame
        parcel properties / parcel spatial data
    tdata : pd.DataFrame with time data in columns 'RH', 'EV24', 'hLR', 'summer'
        time data
    gr: gridObject
        the modflow grid
    use_w_not_c: bool
        use the analytical ditch resisance w_ditch insead of the real c_ditch
    """
    boundary_dict = {}
    if what=='WEL':
        # Given seepage in lowest layer.
        L = (IBOUND[-1] != 0).ravel()
        I = gr.NOD[-1].ravel()[L]
        lrc = np.array(gr.I2LRC(I))  # their lrc index tuples

        # Gnerate basic recarray prototype to store the large amount of data efficiently
        dtype = flopy.modflow.ModflowWel.get_default_dtype()
        spd = np.recarray(len(I), dtype=dtype)

        spd['k'] = lrc[:, 0]
        spd['i'] = lrc[:, 1]
        spd['j'] = lrc[:, 2]

        # You can trigger use of monthly seepage values by having fields
        # q_up01, q_up02 .. q_up12 in parcel_data
        monthly_seepage_values = 'q_up01' in parcel_data.columns
        prev_month= tdata.index[0].month - 1
        for isp, t in enumerate(tdata.index):
            if t.month != prev_month:
                # Use monthly values if available in database.
                fld = f'q_up{t.month:02d}' if monthly_seepage_values else 'q_up'
                q_kwel = parcel_data[fld].values[:, np.newaxis] * np.ones((1, gr.nx))
                Q_kwel = q_kwel.ravel()[L] * gr.Area.ravel()[L] # m3/d per cell.
                spd['flux'] = Q_kwel
                boundary_dict[isp] = spd.copy()
            prev_month = t.month

    elif what == 'DRN':
         # Drains in first layer, but not in first column
        elev = get_drain_elev_with_trenches(parcel_data=parcel_data, gr=gr)

        # Get the indices of the drain locations (see IBOUND[0])
        L  = (IBOUND[0, :, 1:-1]!=0).ravel()     # omit inactive cells
        I  = gr.NOD[ 0, :, 1:-1].ravel()[L] # the sought indeces
        lrc  = np.asarray(gr.I2LRC(I), dtype=int)

        c_drain = parcel_data['c_drain'].values[:, np.newaxis] * np.ones((1, gr.nx))
        cond = (gr.Area / c_drain).ravel()[I]
        elev = elev.ravel()[I]

        # Stress period data
        dtype = flopy.modflow.ModflowDrn.get_default_dtype()
        spd = np.recarray(len(I), dtype=dtype)

        spd['k'] = lrc[:, 0]
        spd['i'] = lrc[:, 1]
        spd['j'] = lrc[:, 2]
        spd['elev'] = elev
        spd['cond'] = cond

        # ionly first isp, because constant.
        isp = 0
        boundary_dict  = {isp: spd}

    elif what == 'GHB':
        # Ditch entry resistance, first column both layers.
        pdata = parcel_data
        if use_w_not_c:
            wi = np.vstack((pdata['wi_ditch'], pdata['wi_ditch2']))
        else: # Use the real ditch resistance, dicth circumference
            wi = np.vstack((pdata['ci_ditch'] * pdata['D1'] / pdata['ditch_omega1'],
                            pdata['ci_ditch'] * pdata['D2'] / pdata['ditch_omega2']))
            # extra resistance due to partial penetration of ditch
            wi += np.vstack((pdata['wpp1'], pdata['wpp2']))

        dy    = np.vstack((gr.dy, gr.dy))
        cond  = np.vstack((pdata['D1'], pdata['D2'])) / wi * dy

        cond[np.isnan(cond)] = 0. # This is where ditch_omega is zero

        I = gr.NOD[:, :, 0].ravel()[cond.ravel()>0]
        lrc  = np.array(gr.I2LRC(I.ravel()), dtype=int)

        dtype = flopy.modflow.ModflowGhb.get_default_dtype()
        spd = np.recarray(len(I), dtype=dtype)

        # fixed
        spd['k'] = lrc[:, 0]
        spd['i'] = lrc[:, 1]
        spd['j'] = lrc[:, 2]
        spd['cond'] =  cond[cond > 0]
        sum_prev = not tdata['summer'].iloc[0]
        for isp, summer in enumerate(tdata['summer']):
            if summer != sum_prev:
                fld = 'h_winter' if summer else 'h_winter'
                spd['bhead'] = parcel_data[[fld, fld]].values.T[cond > 0]
                boundary_dict[isp] = spd.copy()
            sum_prev = summer

    elif what == 'RIV':
        pdata=parcel_data
        if use_w_not_c:
            # Use analytic ditch resistance with layer thickness and no partial penetration
            wi = np.vstack((pdata['wi_ditch'], pdata['wi_ditch2'])) # Same value for both layers
            wo = np.vstack((pdata['wo_ditch'], pdata['wo_ditch2'])) # Same value for both layers
            assert np.all(wi >= wo), "ditch entry resist. must be larger or equal to the ditch exit resistance!"

            dw = wi - wo; eps=1e-10; dw[dw==0] = eps # prevent (handle division by zero)
            w     = (wo * wi / dw)
        else: # Use real ditch resistance with ditch circumference
            ci = np.vstack((pdata['ci_ditch'], pdata['ci_ditch'])) # Same value for both layers
            co = np.vstack((pdata['co_ditch'], pdata['co_ditch'])) # Same value for both layers
            assert np.all(ci >= co), "ditch entry resist. must be larger or equal to the ditch exit resistance!"

            dc = ci - co; eps=1e-10; dc[dc==0] = eps # prevent (handle division by zero)
            c     = (co * ci / dc)

            # To analytic resistance, using the ditch circumference
            w = c * np.vstack((pdata['D1'] / pdata['ditch_omega1'],
                               pdata['D2'] / pdata['ditch_omega2']))
            # Add partial penetration to resistance
            w += np.vstack((pdata['wpp1'], pdata['wpp2']))

        dy    = np.vstack((gr.dy, gr.dy))
        cond  = np.vstack((pdata['D1'], pdata['D2'])) / w * dy

        cond[np.isnan(cond)] = 0. # When ditch_omega is zero

        L = cond > 0

        I  = gr.NOD[:, :, 0][L].ravel()
        lrc  = np.asarray(gr.I2LRC(I), dtype=int)

        dtype = flopy.modflow.ModflowRiv.get_default_dtype()
        spd = np.recarray(len(I), dtype=dtype)

        spd['k'] = lrc[:, 0]
        spd['i'] = lrc[:, 1]
        spd['j'] = lrc[:, 2]
        spd['cond'] = cond[L].ravel()
        sum_prev = not tdata['summer'].iloc[0]
        for isp, summer in enumerate(tdata['summer']):
            if summer != sum_prev:
                fld = 'h_winter' if summer else 'h_winter'
                spd['stage'] = parcel_data[[fld, fld]].values.T[L]
                spd['rbot' ] = parcel_data[[fld, fld]].values.T[L]
                boundary_dict[isp] = spd.copy()
            sum_prev = summer
    return boundary_dict

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

        # compute partcel width to use in GGOR
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

        Parameters
        ----------
        BMIN : float
            minimum parcel width
        BMAX : float
            maximum parcel width
            the final maximum parcel width will be that of the
            widest parcel in the database, not the initial BMAX.
            This will also be the width of the modflow model.
        """
        A     = np.asarray(self.data['A_parcel'])  # Parcel area
        O     = np.asarray(self.data['O_parcel'])  # Parcel cifcumference
        det   = O ** 2 - 16 * A             # determinant
        I     = (det>=0);                   # determinant>0 ? --> real solution
        B     = np.nan * np.zeros_like(I)   # init paracel width
        L     = np.nan * np.zeros_like(I)   # init parcel length
        B[ I] = (O[I] - np.sqrt(det[I]))/4  # width, smallest of the two values
        L[ I] = (O[I] + np.sqrt(det[I]))/4  # length, largest of the two values
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

        Defaults are applied only if no corresponding column exists in self.data

        Parameters
        ----------
        defaults : dict
            The default parametres with their values.
        """
        defcols = set(defaults.keys()).difference(self.data.columns)
        for dc in defcols: # only for the missing columns
            self.data[dc] = defaults[dc]


def compute_and_set_omega(data=None):
    """Compute and set the half wetted ditch circumference in the two model layers.

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

    See theory

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
    """Return parcel info shape dpf file  into pandas.DataFrame.

    Also make sure that the data type is transferred from shapefile to DataFrame.

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
    """Return ax[0], ax[1] for new plot."""
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


def model_parcel_areas(gr, IBOUND):
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
    return ((IBOUND[0] > 1) * gr.Area).sum(axis=1)


def get_parcel_average_hds(HDS=None, IBOUND=None, gr=None):
    """Return the parcel arveraged heads for all parcels and times.

    The problem to solve here is handling the inactive cells that make up
    a different part of each parcel. We make use of the Area of the cells
    and of IBOUND to ignore those cells that are inactive.

    Parameters
    ----------
    HDS : headfile object
        obtained by reading the headfile produced by modflow.
    IBOUND: ndarray
        modflow's IBOUND array
    gr: mfgrid.Grid object

    Returns
    -------
    hds: ndarray
        ndarray of heads of shape [nparcels, ntimes]
    """
    #TODO: Check en vereenvoudig de procedure die is nu niet zomaar te doorzoin.
    # Meld dan in de toelichting wat de vorm is van de binnenkomende heads.
    hds = HDS.get_alldata(nodata=-999.99)
    hds[np.isnan(hds)] = 0.
    avgHds = (np.sum(hds * gr.DxLay[np.newaxis, :, :, :], axis=-1) /
              np.sum(      gr.DxLay[np.newaxis, :, :, :], axis=-1))
    return avgHds # shape(nsp, nz, ny)

#55 Meteo data

class Heads_obj:
    """Heads object, to store and plot head data."""

    def __init__(self, dirs, tdata=None, IBOUND=None, gr=None):
        """Return Heads_obj.

        Paeameters
        ----------
        dirs: Dir_struct object
            GGOR directory structure
        tdata: pd.DataFrame
            time data
        IBOUND: ndarray of gr.shape
            Modflow's boundary array
        gr: fdm.mfgrid.Grid object
            holds the Modflow grid.
        """
        self.case = os.path.basename(dirs.case)
        self.gr = gr

        #% Get the modflow-computed heads
        hds_file = os.path.join(dirs.case_results, self.case +'.hds')
        print("\nReading binary head file '{}' ...".format(hds_file))
        self.HDS = bf.HeadFile(hds_file)
        self.avgHds = get_parcel_average_hds(self.HDS, IBOUND, gr=self.gr)
        self.GXG    = GXG_object(tdata=tdata, avgHds=self.avgHds)


    def plot(self, ax=None, tdata=None, parcel_data=None,
                   selection=[0, 1, 2, 3, 4],
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
        selection: sequence of ints (tuple, list).
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
        selection = selection_check(selection)

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
            for iclr, isel in enumerate(selection):
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
                self.GXG.plot(ax[0], selection=selection)

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
        nparcel = ahds.shape[-1] # ny  avgHds (nth, ny)
        hyears = np.unique(tdat['hyear'])

        # Cut off incomplete start and end hydological years
        if tdat.index[ 0].month != 3 or tdat.index[ 0].day != 14: hyears = hyears[ 1:]
        if tdat.index[-1].month != 2 or tdat.index[-1].day != 28: hyears = hyears[:-1]

        # skip the first hydological year
        hyears = hyears[1:]
        nyear = len(hyears)

        # Format to store the gxg data in a recarray
        gxg_dtype = [('t', pd.Timestamp), ('hd', float), ('hyear', int),
                 ('l', bool), ('h', bool), ('v', bool)]

        # The gxg recarray has 9 records per hyear (3 glg, 3 ghg, 3gvg and nparcel layers)
        # These are tine individual values contribution to the GXG
        self. gxg = np.zeros((nyear * 9, nparcel), dtype=gxg_dtype)

        T = (True, True, True)
        F = (False, False, False)

        for iyr, hyear in enumerate(hyears):
            ah = ahds[tdat['hyear'] == hyear]
            td = tdat.loc[tdat['hyear'] == hyear]
            Ias = np.argsort(ah, axis=0)  # Indices of argsort along time axis

            # Make sure hydrological years start at March 14!!
            assert td.index[0].month ==3 and td.index[0].day == 14, "hyears must start at 14th of March"

            hyr = (hyear, hyear, hyear)

            for ip in range(nparcel):
                Iglg = Ias[0:3, ip]
                Ighg = Ias[-3:, ip]
                Igvg = slice(0, 3, 1)
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

        # Comptue and store the long-term averaged values, the actual GXG
        dtype = [('id', int), ('GLG', float), ('GHG', float), ('GVG', float)]
        self.GXG = np.ones(nparcel, dtype=dtype)
        for ip in range(nparcel):
            self.GXG[ip] = (
                ip,
                self.gxg[self.gxg[:, ip]['v'], ip]['hd'].mean(),
                self.gxg[self.gxg[:, ip]['l'], ip]['hd'].mean(),
                self.gxg[self.gxg[:, ip]['h'], ip]['hd'].mean())


    def plot(self, ax=None, selection=[0, 1, 3, 4, 5], **kwargs):
        """Plot GXG.

        Parameters
        ----------
        selection : list
            list if indices to select the parcels for plotting
        nmax: int
            maximum number of graphs to plot
        """
        selection  = selection_check(selection)

        clrs = 'brgkmcy'

        for iclr, ip in enumerate(selection):
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

        lw = 0.5
        for iclr, ip in enumerate(selection):
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

    def __init__(self, dirs=None, IBOUND=None, parcel_data=None, tdata=None, gr=None):
        """Return Watbal object carrying the water budget for all cross sections in m/d.

        Parameters
        ----------
        dirs: Dir_struct object
            directory structure of GGOR
        IBOUND: numpy.ndarray (of ints(
            Modeflow's IBOUND array
        parcel_data: pd.DatFrame
            parcel property data and spacial data (table)
        tdata: pd.DataFrame
            time / meteo data. Only the tdata.index is required
        gr: fdm_tools.mfgrid.Grid

        Generates
        ---------
        self.W : np.recarray having labels ['RCH', 'EVT', 'GHB etc'] and values that
            are np.ndarrays with shape (1, nlay, nrow, nper)
            Each row has the cross sectional discharge in m/d.

        Note thta labels must be adapted if new CBC packages are to be included.

        @TO 20170823 as function
        @TO 20200907 turned into class
        """
        L = ['RCH', 'EVT', 'WEL', 'GHB', 'RIV', 'DRN', 'FLF', 'STO'] #Layer 0 labels

        self.case = os.path.basename(dirs.case)

        cbc_file = os.path.join(dirs.case_results, self.case  +'.cbc')
        print("\nReading binary cbc file '{}'.".format(cbc_file))
        print("Scanning this file may take a minute or so ...")

        self.CBC=bf.CellBudgetFile(cbc_file)

        dtype=[(lbl, float, (self.CBC.nlay, self.CBC.nrow, self.CBC.nper)) for lbl in L]

        self.W = np.zeros(1, dtype=dtype)
        self.gr = gr

        # Area of each cross section, made 3D (1, nlay, nrow, 1) compatible with W
        A_xsec = np.sum(gr.DxLay * gr.DyLay * IBOUND, axis=-1)[np.newaxis, :, :, np.newaxis]

        # Relative contribution of parcel tot total after reducing model parcel area to 1 m2
        #Arel = ((parcel_data['A_parcel'] / model_parcel_areas(gr, IBOUND)) / (parcel.data['A_parcel'].sum()))[:, np.newaxis]
        vals3D = np.zeros(gr.shape).ravel()
        print()
        for lbl in L:
            print(lbl, end='')
            if lbl in ['WEL', 'GHB', 'RIV', 'DRN']:
                vals = self.CBC.get_data(text=cbc_labels[lbl])
                for iper in range(self.CBC.nper):
                    vals3D[vals[iper]['node']-1] = vals[iper]['q']
                    self.W[lbl][0][:, :, iper] = np.sum(vals3D.reshape(gr.shape), axis=-1)
                    vals3D[:] = 0.
                    if iper % 100 == 0: print('.',end='')
                print(iper)
            elif lbl in ['RCH', 'EVT']:
                vals = self.CBC.get_data(text=cbc_labels[lbl])
                for iper in range(self.CBC.nper):
                    self.W[lbl][0][0, :, iper] = np.sum(vals[iper][1], axis=-1)
                    if iper % 100 == 0: print('.',end='')
                print(iper)
            else:
                vals = self.CBC.get_data(text=cbc_labels[lbl])
                for iper in range(self.CBC.nper):
                    self.W[lbl][0][:, :, iper] = np.sum(vals[iper], axis=-1)
                    if iper % 100 == 0: print('.', end='')
                print(iper)
            self.W[lbl] /= A_xsec # from m3/d to m/d

        # FLF lay 0 is the inflow of lay 1 and an outflow of lay 0
        self.W['FLF'][0][1, :, :] = +self.W['FLF'][0][0, :, :]
        self.W['FLF'][0][0, :, :] = -self.W['FLF'][0][1, :, :]

        print('Done. See self.CBC and self.W')


    def plot(self, parcel_data=None, tdata=None,
                    selection=None, sharey=False, ax=None):
        """Plot the running water balance of the GGOR entire area in mm/d.

        Parameters
        ----------
        parcel_data: pd.DataFrame
            parcel spatial and property data
        tdata: pd.DataFrame with time index
            time_index must correspond with MODFLOW's output
            if not specified, then CBB.times is used
        selection: list or sequence
            the water budget will be taken over the indices in selection.
            Use None for all.
        sharey: bool
            the y-axis of the two charts will be shared if True
        ax: plt.Axis object or None (axes will the be created)
            axes for figure.
        """
        m2mm = 1000. # from m to mm conversion

        selection = selection_check(selection)

        missing = set(watbal_label.keys()).difference(set(cbc_labels.keys()))
        if len(missing):
            raise ValueError("Missing labels = [{}]".format(', '.join(missing)))

        tindex = self.CBC.times if tdata is None else tdata.index

        # Sum over all Parcels. The watbal values are in mm/d. To sum over all
        # parcels multiply by their share of the regional area [-]
        Arel = (parcel_data['A_parcel'].values[selection] /
                parcel_data['A_parcel'].values[selection].sum()
                )[np.newaxis, :, np.newaxis] # (nlay, selected_parcels, nper)

        # summed over the entire selection!! (so nrow drops)
        dtype = [(lbl, float, (self.CBC.nlay, self.CBC.nper)) for lbl in cbc_labels]
        V = np.zeros(1, dtype=dtype)

        # From now in mm/d
        for lbl in cbc_labels: # note selection
            V[lbl][0] = np.sum(self.W[lbl][0][:, selection, :] * Arel * m2mm,
                            axis=-2) # also to mm/d

        clrs = [watbal_label[L]['clr'] for L in watbal_label]
        lbls = [watbal_label[L]['leg'] for L in watbal_label]

        if isinstance(selection, slice):
            ttl = ' Taken over parcels[{}:{}:{}]'.format(
                        selection.start, selection.stop, selection.step)
        else:
            ttl = ' Taken over parcels [{}]'.format(
                ', '.join(['{}'.format(s) for s in selection]))

        if ax is None:
            ax = newfig2(titles=(
                    'Water balance top layer.'   + ttl,
                    'Water balance botom layer.' + ttl), xlabel='time', ylabels=['mm/d', 'mm/d'],
                         size_inches=(14, 8), sharey=False, sharex=True)

        if not isinstance(ax, (list, tuple, np.ndarray)):
            raise ValueError("Given ax must be an sequence of 2 axes.")

        V0 = np.zeros((len(watbal_label), self.CBC.nper))
        V1 = np.zeros((len(watbal_label), self.CBC.nper))
        for i, lbl in enumerate(watbal_label):
            V0[i] = V[lbl][0, 0, :]
            V1[i] = V[lbl][0, 1, :]

        ax[0].stackplot(tindex, V0 * (V0>0), colors=clrs, labels=lbls)
        ax[0].stackplot(tindex, V0 * (V0<0), colors=clrs) # no labels
        ax[1].stackplot(tindex, V1 * (V1>0), colors=clrs, labels=lbls)
        ax[1].stackplot(tindex, V1 * (V1<0), colors=clrs) # no labels

        ax[0].legend(loc='best', fontsize='xx-small')
        ax[1].legend(loc='best', fontsize='xx-small')

        return ax


def run_modflow(dirs=None, parcel_data=None, tdata=None, laycbd=(1, 0), dx=1.,
                use_w_not_c=False):
        """Simulate GGOR using MODFLOW.

        Parameters
        ----------
        dirs: DirStruct object
            directory structure object, containing home and case information
        parcel_data: pd.DataFrame
            parcel data (spacial)
        tdata: pd.DataFrame
            meteo data used to generate stress periods
        use_w_not_c: bool
            use analytical ditch resistance w_ditch instead of real c_ditch
            c_ditch uses omega and effect of partial penetration of ditches
        """
        gr   = grid_from_parcel_data(
                    parcel_data=parcel_data, dx=dx, laycbd=laycbd)

        par = set_spatial_arrays(parcel_data=parcel_data, gr=gr)

        spd  = set_stress_period_data(tdata=tdata)


        bdd  = set_boundary_data(parcel_data=parcel_data,
                                     tdata=tdata,
                                     gr=gr,
                                     IBOUND=par['IBOUND'],
                                     use_w_not_c=use_w_not_c)

        # MODEL package parameter defaults for GGOR
        ipakcb = 53    # all cbc output to this unit
        chani = 1e-20  # no contact between adjacent rows in GGOR
        compact = True # means CBC output is of compact type
        nrchop = 3 # precip in first active layer
        nevtop = 3 # precip in first active layer
        layvka = 0 # vka valuea are vk
        mxiter = 200
        iter1 = 200
        hclose = 0.001
        rclose = 0.01

        fm = flopy.modflow

        case = os.path.basename(dirs.case)

        mf  = fm.Modflow(case, exe_name=dirs.exe_name, model_ws=dirs.case_results, verbose=True)

        #MODFLOW packages
        dis = fm.ModflowDis(mf, gr.nlay, gr.ny, gr.nx,
                            delr=gr.dx, delc=gr.dy, top=gr.Ztop[0], botm=gr.Zbot,
                            laycbd=gr.LAYCBD,
                            nper=spd['NPER'], perlen=spd['PERLEN'],
                            nstp=spd['NSTP'],
                            steady=spd['STEADY'])

        bas = fm.ModflowBas(mf, ibound=par['IBOUND'], strt=par['STRTHD'])

        lpf = fm.ModflowLpf(mf, hk=par['HK'], layvka=layvka, vka=par['VK'], chani=chani,
                            sy=par['SY'], ss=par['SS'],
                            laytyp=par['laytyp'], vkcb=par['VKCB'], ipakcb=ipakcb,
                            storagecoefficient=par['storagecoefficient'])

        ghb = fm.ModflowGhb(mf, stress_period_data=bdd['GHB'], ipakcb=ipakcb)
        riv = fm.ModflowRiv(mf, stress_period_data=bdd['RIV'], ipakcb=ipakcb)
        drn = fm.ModflowDrn(mf, stress_period_data=bdd['DRN'], ipakcb=ipakcb)
        wel = fm.ModflowWel(mf, stress_period_data=bdd['WEL'], ipakcb=ipakcb)

        rch = fm.ModflowRch(mf, nrchop=nrchop, rech=spd['RECH'], ipakcb=ipakcb)
        evt = fm.ModflowEvt(mf, nevtop=nevtop, evtr=spd['EVTR'],
                            surf=par['SURF'], exdp=par['EXDP'], ipakcb=ipakcb)

        oc  = fm.ModflowOc( mf, stress_period_data=bdd['OC'], compact=compact)

        #pcg = fm.ModflowPcg(mf, mxiter=200, iter1=200, hclose=0.001, rclose=0.001)
        sms = fm.ModflowSms(mf, mxiter=mxiter, iter1=iter1, hclose=hclose, rclosepcgu=rclose)

        packages = {'dis': dis, 'bas':bas, 'lpf':lpf, 'ghb':ghb, 'riv':riv, 'drn':drn,
                    'wel':wel, 'rch':rch, 'evt':evt, 'oc':oc, 'sms':sms}
        print('Pakages used:''[{}]'.format(', '.join(packages.keys())))

        # write data and run modflow
        mf.write_input()
        success, mfoutput = mf.run_model(silent=False, pause=False)

        print('Running success = {}'.format(success))
        if not success:
            raise Exception('MODFLOW did not terminate normally.')

        return par, spd, bdd, gr


def data_to_excel(dirs=None, data=None, fname=None):
    """Sace the current parcel_data DataFrame to an excel file.

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

    test=False

    # Parameters to generate the model. Well use this as **kwargs
    GGOR_home = os.path.expanduser('~/GRWMODELS/python/GGOR') # home directory
    case = 'AAN_GZK'

    #GGOR directory structure
    dirs = Dir_struct(GGOR_home, case=case)

    #Get the meteo data from an existing file or directly from KNMI
    meteo_data = knmi.get_weather(stn=240, start='20100101', end='20191231',
                                  folder=dirs.meteo)

    # Add columns "summer' and "hyear" to it"
    tdata = handle_meteo_data(meteo_data, summer_start=4, summer_end=10)

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
        # Bofek data, coverting from code to soil properties (kh, kv, sy)
        # The BOFEK column represents a dutch standardized soil type. It is used.
        # Teh corresponding values for 'kh', 'kv' and 'Sy' are currently read from
        # and excel worksheet into a pd.DataFrame (table)
        bofek = pd.read_excel(os.path.join(dirs.bofek, "BOFEK eenheden.xlsx"),
                              sheet_name = 'bofek', index_col=0, engine="openpyxl")

        # Create a GGOR_modflow object and get the upgraded parcel_data from it
        parcel_data = GGOR_data(defaults=defaults, bofek=bofek, BMINMAX=(5, 250),
                GGOR_home=GGOR_home, case=case).data

    # MODFLOW input arrays are int the returned dicts
    par, spd, bdd, gr =  run_modflow(dirs=dirs, parcel_data=parcel_data, tdata=tdata,
                                     use_w_not_c=True)

    #%% Get the modflow-computed heads and cell by cell flows

    heads = Heads_obj(dirs, tdata=tdata, IBOUND=par['IBOUND'], gr=gr)

    watbal = Watbal_obj(dirs,
                           tdata=tdata,
                           parcel_data=parcel_data,
                           IBOUND=par['IBOUND'],
                           gr=gr)

    #%% Open HDS file and plot heads (with GXG)
    if test:
        for tst in set(parcel_data['Test'].values):
            selection = list(parcel_data.index[parcel_data['Test'] == tst])
            test_vals_str = '{}'.format(', '.join(
                [str(tv) for tv in parcel_data[tst].iloc[selection]]))

            titles=['Parcel averaged heads, testvar {} in [{}]'.format(tst, test_vals_str),
                    'Parcel averaged heads, testvar {} in [{}]'.format(tst, test_vals_str)]
            ax = heads.plot(tdata=tdata,
                           parcel_data=parcel_data,
                           selection=selection,
                           titles=titles,
                           size_inches=(14, 8))
            if False:
                ax = watbal.plot(parcel_data=parcel_data,
                                 tdata=tdata,
                                 selection=selection[0],
                                 sharey=True)
    else:
        selection = [0, 1, 2, 3]
        titles=['Parcel averaged heads', 'Parcel averaged heads']
        ax = heads.plot(tdata=tdata,
               parcel_data=parcel_data,
               selection=selection,
               titles=titles,
               size_inches=(14, 8))
        ax = watbal.plot(parcel_data=parcel_data,
                         tdata=tdata,
                         selection=None,   # over all parcels
                         sharey=True)

#%%
    print('---- All done ! ----')
#%%
    ax = watbal.plot(parcel_data=parcel_data,
                         tdata=tdata,
                         selection=None,   # over all parcels
                         sharey=True)
