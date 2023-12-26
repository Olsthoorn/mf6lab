#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 22 10:58:55 2018

Analyical transient simulations of cross sections between parallel ditches.

Names of the different analytial solutions:

    L + [1|2] + [q]f] [+ W]
    This tells that the solution has 1 or 2 computed layers and the boundary
    condition in the underlying regional aquifer is eigher seepge (q) or head
    (f) and that it has or does not have entry/outflo resistance at the ditch.

    So we can name solutions as follows
    L1q, L2q, L1f L2f, L1qw, L1fw, L2f2, L2qw

@author: Theo
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def get_meteo(txtfile):

    data = pd.read_csv(txtfile, delim_whitespace=True, names=['date', 'P', 'E'],
                        index_col=0, parse_dates=True, dayfirst=True)
    data['P'] /= 1000. # to m/d
    data['E'] /= 1000. # to m/d
    return data


template = {'L': 100, # distance between the parallel ditches
           'bd': 1.0, # width of each ditch

           'z0': 0.0, # elevation of ground surfacee
           'zd': -0.6, # elevation of bottom of ditch
           'z1': -11, # elevation of bottom of cover layer
           'z2': -50, # elevation of bottom of second layer

           'Sy': 0.1, # specific yield of cover layer
           'S2': 0.001, # elastic storage coeff. for second layer

           'kh1': 2, # horizontal conductivity of cover layer
           'kv1': 0.2, # vertical conductivity of cover layer
           'kh2': 25, # horizontal conductivity of regional aquifer
           'kv2': 10, # vertical conductivity of regional aquifer

           'c': 100, # vertical resistance at bottom of cover layer
           'w': 0, # entry resistance of dich (as a vertical wall) [L]

           }

def simulate_single_Layer(name, props=None, time_data=None):
    '''Return results for solution 1 layer given upward seepage

    parameters
    ----------
        name: str
            solution name, one of ['lq1', 'lq1w']
            'l1'  : Single layer with no seepage.
            'l1q' : Single layer with upward seepage no ditch resistance
            'l1qw': Single layer with upward seepage with ditch resistance
            'l1f' : Single layer with given head in regional aquifer
            'l1fw': Same but with ditch resistance.
            'l2q  : Two aquifer system with given seepage from lower aquifer'
            'l2qw': Same but with ditch resistance
        ime_data: pd.DataFrame
            required fields: hLR, N, q or hLR

        The resistance between cover layer and regional aquifer is concentrated
        in the given cover-layer resistance 'c' at the bottom of the cover layer.
    '''

    L, c, mu = props['L'], props['c'], props['Sy']
    k, D, b  = props['kh1'], np.abs(props['z1'] - props['z0']), L / 2
    lamb     = np.sqrt(k * D * c)
    Lamb     = np.tanh(b / lamb) / (b / lamb)

    if name in ['l1q', 'l1f']:
        T = mu * c * (1 - Lamb) / Lamb
    elif name in ['l1qw', 'l1fw']:
        w = props['w']
        k = props['kh1']
        B = 1 / (k * b * w / lamb**2 + b/lamb / np.tanh(b/lamb) - 1) + 1
        T = mu * c / (B - 1)
    else:
        raise ValueError("Don't recognize solution name <{}>".format(name))

    hLR, N = time_data['hLR'].values, time_data['N'].values,

    if name in ['l1f', 'l1fw']: # head is given
        Q = - (hLR - time_data['f'].values) / c
    else: # seepage given
        Q = time_data['q'].values

    # Times in days since start of day given by the first
    times = time_data.index - time_data.index[0] + pd.Timedelta(1, 'D')

    # t includes start point one day before first time in days
    t = np.hstack((0, np.asarray(np.asarray(times,
                                        dtype='timedelta64[D]'), dtype=float)))

    DT = np.diff(t)
    h = np.zeros(len(t)) # initialize h
    h[0] = hLR[0] # intiial head at start of day 1.

    # Integration
    for i, (dt, hlr, n , q) in enumerate(zip(DT, hLR, N, Q)):
        e = np.exp(-dt/T)
        h[i + 1] = hlr + (h[i] - hlr) * e + (n + q) * T/mu * (1 - e)

    # Keep only h[1:] and put this in time_data as column 'h'
    time_data['h'] = h[1:]

    return time_data


class Solution:
    '''Analytic solution object. Allows simulation head betweeen parallel
    ditches in a one or two-layer aquifer system.
    The outcome is the average head in the cross section over time and its
    running water budget. Or the head as a function of x in the cross section
    for the input data averaged over time according to the input.

    template = {'b': 100, # distance between the parallel ditches
               'bd': 1.0, # width of each ditch

               'z0': 0.0, # elevation of ground surfacee
               'zd': -0.6, # elevation of bottom of ditch
               'z1': -11, # elevation of bottom of cover layer
               'z2': -50, # elevation of bottom of second layer

               'Sy': 0.1, # specific yield of cover layer
               'S2': 0.001, # elastic storage coeff. for second layer

               'kh1': 2, # horizontal conductivity of cover layer
               'kv1': 0.2, # vertical conductivity of cover layer
               'kh2': 25, # horizontal conductivity of regional aquifer
               'kv2': 10, # vertical conductivity of regional aquifer

               'w':   20, # ditch resistance as that of a vertical wall

               }
    '''

    def __init__(self, properties=template):
        '''Return an instance of an analytical solution only storing name and properties.

        parameters
        ----------
            template: dict
                a dict containing the properrties. The necessary properties
                are given in the example tamplate in this class. Not all
                properties are used by all solutions. Unsued properties
                may be omitted from the actual template.
        '''
        self.properties = template


    def solve(self, solution_name=None, time_data=None):
        '''Return average heads versus time and water budget versus time

        parameters
        ----------
            solution_name: str
                the name of the solutions, legal names are
                [L1q, L2q, L1f L2f, L1qw, L1fw, L2f2, L2qw]

            time_data: input pd.DataFrame containing a time series as defined.
                data's index must be pd.Timestamps
                columns must include hLR, P, E, I, q in which
                    hLR: water level elevations in the ditches [m]
                    P: Precipitation [m/d]
                    E: Evapotranspiration [m/d]
                    If recharge is given, one can omit E or set it to zeros.
                If interception is known include column
                    I: Interception [m/d]
                Depending on the solution, columns must include one of
                    q: Upward seepage from lower aquifer
                    f: Head in lower aquifer

        The time stamps mark dates at the end of which the head is to be computed.

        This is irrespective of the duration between two successive time stamps.

        Time steps are computed by the difference of successive time stamps.

        The flows and dictch levels are assumed constant druing the time steps.
        The values pertain to the end of the data of the time stamp given in
        the index of the data series.

        This scheme leaves the head at the start of the first time step
        undetermined. This initial head lies in principle before the first
        time stamp of the index. The siplest solution is to choose the
        ditch water elevation of the first tme step as its value.

        The unknown time stamp of this value is taken as the beginning of the
        data geiven by the first time stamp. Hence the first head is computed
        assuming the head at the beginning of the first date equals hLR and
        the flows given for the first time step apply.

        This way we'll compute a head at every time stamp, including the first.

        '''

        # We'll only use N (net recharge as meteo input series)
        # Verify presence of column 'N' or try to compute it
        if not 'N' in time_data.columns:
            time_data['N'] = time_data['P']
            try:
                time_data['N'] -= time_data['E']
            except:
                pass
            try:
                time_data['N'] -= time_data['I']
            except:
                pass
        else:
            raise ValueError("Can't compute recharge N from meteo data")

        # Veriry presence of 'q' in time_data
        if not 'q' in time_data.columns:
            raise ValueError('need seepage "q" as one of the data columns.')

        # Verify presence of 'hLR' in time_data
        if not 'hLR' in time_data.columns:
            raise ValueError('need hLR ditch water elevation in input data.')

        # Verify valid solution names
        if solution_name not in ['l1q', 'l1qw']:
            raise ValueError('Unknown solution name.')

        # simulate
        data = simulate_single_Layer(name=solution_name,
                                    props=self.properties, time_data=time_data)
        return data


if __name__ == '__main__':

    home = '/Users/Theo/GRWMODELS/python/GGOR/'

    metfile = os.path.join(home, 'meteo/PE-00-08.txt')
    data    = get_meteo(metfile)

    data['P'] = 0.005
    data['E'] = 0.

    data['hLR'] = 0. # add 'hLR' to data
    data[  'q'] = 0. # add '  q' to data
    data[  'f'] = 0. # add 'phi' to data


    # generate solution index pass its properties
    mySolution = Solution(properties=template)

    # solve using a specific solution and meteo data
    # copy() is necessary !
    data1 = mySolution.solve(solution_name='l1q', time_data=data.copy())
    data2 = mySolution.solve(solution_name='l1qw', time_data=data.copy())
    data3 = mySolution.solve(solution_name='l1f', time_data=data.copy())
    data4 = mySolution.solve(solution_name='l1fw', time_data=data.copy())


    # This shows that both solutions yield the same values of w == 0.
    fig, ax = plt.subplots()
    ax.set_title('Mean groundwater head in cross section, solution "l1q"')
    ax.set_xlabel('time')
    ax.set_ylabel('head, elevation [m]')
    ax.grid()
    ax.plot(data1.index, data1['h'], 'rx', label='no ditch resistance')
    ax.plot(data2.index, data2['h'], 'b.', label='with ditch resistance')
    ax.plot(data3.index, data2['h'], 'b.', label='with ditch resistance')
    ax.plot(data4.index, data2['h'], 'b.', label='with ditch resistance')
    ax.legend()
