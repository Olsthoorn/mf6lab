
# %%
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# %%
one_foot = 0.3
dt_start = np.datetime64('2018-01-01')
dt_end   = np.datetime64('2020-04-01')

#%%
def get_time_series(dirs, names=('Lakes', 'WVP1', 'WVP2')):
    """Return time series with given names from csv files in data_folder.
    
    Parameters
    ----------
    drs: instantiation of mf6tools.Dirs
        holds project directory Structure
    names: list of str
        names of the data series in the csv files.
    
    Returns
    -------
    dict with names as keys dicts, where the last dict holds the finromation
    necessary for mdflow to deal with the data as a time series. It has the keys
    defined by mf6 (filename, time_series_namerecord, timeseries, interpolation_methodrecord and sfacrecord)
    """
    
    time_series = dict()

    for name in names:
        series = pd.read_csv(os.path.join(dirs.Data, name + '.csv'), header=0, index_col=0)
        series[series.columns[0]] *= one_foot
        series.index = dt_start + np.round(
                            (dt_end - dt_start) / np.timedelta64(1, 'D') * series.index
                            ) * np.timedelta64(1, 'D')
        
        ts = np.vstack((
            np.asarray(series.index - series.index[0]) / np.timedelta64(1, 'D'),
            series.values.T)).T

        time_series[name] = {
            'filename': os.path.join(dirs.case, 'GWF', name + '.ts'),
            'time_series_namerecord': name,
            'timeseries': [(t, h) for t, h in ts],
            'interpolation_methodrecord': 'linear',
            'sfacrecord': 1.1
        }
            
    return time_series

if __name__ == '__main__':
    
    data_folder = '/Users/Theo/GRWMODELS/python/Pennink-Model/cases/VergroesenFremont/Data'
    os.chdir(data_folder)
    
    time_series = get_time_series(data_folder)
    
    fig, ax = plt.subplots()
    fig.set_size_inches(12, 5)
    ax.set_title("Time Series, Fremont (Ca), Quarry Lakes, Recr. Park")
    ax.set_xlabel("date time")
    ax.set_ylabel("elevation [m]")
    ax.grid(True)
    
    for name, series in time_series.items():
        ts = np.asarray(series['timeseries'])
        date = np.datetime64('2018-01-01') + ts.T[0] * np.timedelta64(1, 'D')
        ax.plot(date, ts[:, 1], '.-', label=name)
    ax.legend()
    
    plt.show()
