
import os
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

class Meteo():
    """Class to handle meteo data."""
    
    def __init__(self, met_csv_file=None):
        """Get the meteo time series data from csv file."""
        
        if met_csv_file is None:
            met_csv_file = get_meteo_path()
        
        assert os.path.isfile(met_csv_file), f"No such file {met_csv_file}"
        
        met = pd.read_csv(met_csv_file, skiprows=16, header=None)
        
        met.columns = ['STN','YYYYMMDD',   'DR',   'RH',  'RHX', 'RHXH', 'EV24']
        
        met.index = [__class__.to_timestamp(p) for p in met['YYYYMMDD']]
        
        met = met[['RH', 'EV24']]
        
        met['RH'] = __class__.met_col_to_m(met['RH'])
        met['EV24'] = __class__.met_col_to_m(met['EV24'])
        met['RCH'] = met['RH'] - met['EV24']
        
        self.data = met
        
    @property
    def recharge(self):
        return self.data['RH'] - self.data['EV24']
    
    @staticmethod
    def to_timestamp(p):
        """Return timestamp from int like 20251218.
        
        to_timestamp(20251218) --> pd.Timestap("2025-12-18")
        """
        p = str(p)
        return pd.Timestamp(f"{p[:4]}-{p[4:6]}-{p[6:8]}")

    @staticmethod
    def met_col_to_m(metcol):
        """Conver meteo column data from 10ths of mm to m"""
        values = metcol.values / 10000
        values[values < 0] = 0.00005
        return values

def get_meteo_path():
    """Return path to the current meteo file."""
    P = Path(os.getcwd())
    try:
        home = os.path.join(*P.parts[:P.parts.index('AAN_GZK') + 1])
    except ValueError:
        home =  os.path.join(*P.parts[:P.parts.index('GGOR') + 1])
        home = os.path.join(home, 'cases', 'AAN_GZK')
    
    meteo = os.path.join(home, 'meteo')

    assert os.path.isdir(meteo), 'Meteo dir not found.'

    return os.path.join(meteo, os.listdir(meteo)[0])


if __name__ == '__main__':
    #metfile = get_meteo_path()
    #meteo = Meteo(metfile)
    # meteo.data.plot(lw=0.5)
    #meteo.recharge.plot(lw=0.5)
    #plt.show()
    from scipy.special import erfc
        
    n = np.logspace(0, 2, 101)[1:]
    
    fig, ax = plt.subplots()
    fig.suptitle("Required terms in series")
    ax.set_title("Remaining exp series")
    ax.set(xlabel='n', ylabel='rest of series', xscale='log')
    for tau in [0.001, 0.01, 0.1, 1]: 
        y = 0.25 * np.sqrt(np.pi / tau) * erfc((2 * n + 1) * np.sqrt(tau))
        ax.plot(n, y, label=f"t/T={tau}")
    ax.grid(True)
    ax.legend()
    parts = list(Path(os.getcwd()).parts).index('GGOR')
    pth = os.path.join(os.getcwd(), 'doc', 'images')
    fig.savefig(os.path.join(pth, 'needed_terms.png'))
    
    plt.show()