# %% Question from Toine Vergroesen regarding groundwater behavior around Quarry Lakes Recreational Area, Fremont, California (2023/12/08)

"""
See for situation and images the folder images/*.xlsx

Tekst van Toine Vergroesen in email van 2023/12/08:

Beste Theo,
 
In een project in Californië kom ik een (nog) onverklaarbaar fenomeen tegen (zie bovenste schematische figuur onder deze email).
Het grondwater wordt vnl gevoed door infiltratievijvers, oude grindwinningen, die insnijden in wvp1, maar niet tot de bodem daarvan.
 
De gemeten variatie in stijghoogte in peilbuizen in het gebied rondom deze vijvers (zie onderste plaatje onder deze email) is in wvp2 overal beduidend groter dan in wvp1.
De stijghoogte zelf is wezenlijk lager in wvp2 (rode getallen) dan in wvp1 (zwarte getallen).Er zit dus een flinke weerstand tussen wvp1 en wvp2, een laag van ca 10 meter (30 ft) dikte.
 
De kD in wvp1 is 2 tot 5 keer zo hoog als in wvp2.
Het grondwater in beide pakketten is zoet.
Er wordt onttrokken uit beide pakketten, maar de windebieten variëren relatief weinig in de tijd.
Het variabele patroon van het waterniveau in de vijvers werkt sterk door in alle meetpunten, verschillen in neerslagoverschot zie je niet of nauwelijks terug in de metingen.
 
Heb jij een idee hoe zou zoiets kunnen?
En (hoe) zou je dat met een model kunnen simuleren?
 
Kan het te maken hebben met het gewicht van het water in de vijvers dat doorwerkt in het confined wvp2?
 
Ik ben benieuwd naar jouw ideeën hierover en hoop dat je tijd hebt om hierover na te denken.
 
Vriendelijke groeten,
Toine Vergroesen
 
"""

# TO 090312

import os
import sys
import numpy as np
from fdm import Grid
from mf6lab import mf6tools
import get_obs

HOME = '/Users/Theo/GRWMODELS/python/Pennink-Model/'

LENGTH_UNITS = 'meters'
TIME_UNITS = 'days'
FEET = 0.3 # [m]

dirs = mf6tools.Dirs(HOME)

sim_name = 'Fremont'
dirs = dirs.add_case(sim_name)

# Parameters workbook
params_wbk = os.path.join(dirs.case, sim_name + '.xlsx')

# chdir to case directory
os.chdir(dirs.case)
sys.path.insert(0, dirs.case)


# %% === tdis ==========  Period Data:
start_date_time = '2018-01-01' # Must be a string.

perDF = mf6tools.get_periodata_from_excel(params_wbk, sheet_name='PER')
period_data = [tuple(sp) for sp in perDF[['PERLEN', 'NSTP', 'TSMULT']].values]

# Stress period start times
sp_start_times = np.datetime64(start_date_time) + np.cumsum(
                  np.hstack((0., np.array([p[0] for p in period_data])[:-1]))
                  ) * np.timedelta64(1, 'D')

Simtdis = {'perioddata': period_data,
           'nper': len(period_data),
           'start_date_time': start_date_time,
           'time_units': TIME_UNITS,
           }

# %% === DIS ========== Grid
Rlakes = 725 # [m] Omtrek 5200 m, opp 1.63e m2.
idLake, idRand = 2, 3

xGr = np.hstack((0., np.logspace(1, 4, 41), 10001., Rlakes)) # [m] radiaal model
zGr = np.hstack((0., -np.cumsum(np.array([20, 40, 10, 30]))))

gr = Grid(xGr, [-0.5, 0.5], zGr, axial=True) # mesh houskeeping

IDOMAIN = gr.const(1, dtype=int)
IDOMAIN[gr.XM <= Rlakes] = idLake
IDOMAIN[:, :, -1]        = idRand

Gwfdis = {'gr': gr,
          'idomain': IDOMAIN,
          'length_units': LENGTH_UNITS}


# %% ==== STO ===============
Gwfsto = {'sy': gr.const(np.array([0.2,  0.2,  0.1,  0.1 ])[:, np.newaxis, np.newaxis]),
          'ss': gr.const(np.array([1e-3, 1e-3, 1e-6, 1e-6])[:, np.newaxis, np.newaxis])}

# %% === Gwfnpf =========== Horizontal and vertical conductivity
Gwfnpf = {  'k':   gr.const(np.array([0.1, 10, 1e-6, 10])[:, np.newaxis, np.newaxis]),
            'k33': gr.const(np.array([0.1, 10, 1e-6, 10])[:, np.newaxis, np.newaxis])}

# %% === Gwfic ======== Initial condictions (head)


strthd = gr.const(np.array([14., 14., 7., 0.])[:, np.newaxis, np.newaxis] * FEET)

Gwfic = {'strt': strthd}

# %% === CHD, Fixed head period data (Only specify the first period)

time_series = get_obs.get_time_series(dirs, names=['Lakes', 'WVP1', 'WVP2'])

lrc_lakes = gr.LRC(gr.NOD[gr.XM <= Rlakes], astuples=True)
lrc_wvp1  = gr.LRC(gr.NOD[1, :, gr.xm > 2000][:, 0].ravel(), astuples=True)
lrc_wvp2  = gr.LRC(gr.NOD[3, :, gr.xm > 4000][:, 0].ravel(), astuples=True)

period_data = []
for (name, ts), lrcs in zip(time_series.items(), [lrc_lakes, lrc_wvp1, lrc_wvp2]):
      period_data += [(lrc, name, name) for lrc in lrcs]
      
tseries  = []
for name, ts in time_series.items():
      tseries.append(ts)
      
Gwfchd = {'stress_period_data': {0: period_data},
          'timeseries': tseries,
          'boundnames': True,
          'maxbound': 30,
}

# %% === WEL, required to lower head to measured  value
# Location of wells
# Rwvp1, Rwvp2
# Extraction by wells

# %% === OC ====

Gwfoc = {'head_filerecord':   os.path.join(dirs.SIM, "{}Gwf.hds".format(sim_name)),
         'budget_filerecord': os.path.join(dirs.SIM, "{}Gwf.cbc".format(sim_name)),
         'saverecord': [("HEAD", "ALL"), ("BUDGET", "ALL")],
}

print('Done mf_adapt')

if __name__ == '__main__':
   
   print(dirs)     