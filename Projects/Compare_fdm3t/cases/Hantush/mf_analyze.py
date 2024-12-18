# %% Analyzing output of Pennik Series 2 tests
# TO 090223
# %%
import os
import numpy as np
import matplotlib.pyplot as plt
import flopy
import mf_adapt
import settings
from analytic.hantush_convolution import Wh
from etc import color_cycler
from src import mf6tools
from fdm.fdm3t import fdm3t, Grid, dtypeGHB, dtypeQ, dtypeH
from fdm.mf6_face_flows import get_struct_flows # (flowjas, grb_file=None, verbose=False)

use_models, use_packages  = mf6tools.get_models_and_packages_from_excel(
                                                mf_adapt.params_wbk, sheet_name='NAM')

SAVE_ANIMATION = True
cmap = 'binary_r'

colors = ['gray', 'black']

dirs     = mf_adapt.dirs
grb_file = os.path.join(dirs.GWF, mf_adapt.sim_name + 'Gwf.dis.grb')
gr       = mf_adapt.gr
pr       = mf_adapt.pr # properties from settings
               
# %% load the unformatted files with the heads, the concentration and the budget terms
sim = flopy.mf6.MFSimulation.load(
    sim_name=mf_adapt.sim_name, version='mf6', sim_ws=dirs.SIM)

time_units = sim.tdis.time_units.get_data()

sim_name = sim.name.lower()
gwf = sim.get_model(f'{sim_name}Gwf') # list(sim.model_names)[0])

# %% load the simulation results
headsObj = gwf.output.head()
budObj   = gwf.output.budget()

# %% the date times at the end of each time step
kstpkper = np.array(headsObj.get_kstpkper())
sim_time =  np.array(headsObj.get_times())

#  %% Get the flow right face for computing the stream function:
flowjas = budObj.get_data(text='FLOW-JA')

# fflows: dict with face flows (keys: 'fff', 'frf', 'flf') as in MF5 
fflows = get_struct_flows([flowjas], grb_file=grb_file, verbose=False)

# %% Get suitable levels for contouring the stream function psi [m2/d]
hmin, hmax, hInactive = np.unique(headsObj.get_alldata())[[0, -2, -1]]

# %% Set up initial plot of cross section with contouring
fig, ax = plt.subplots()
fig.set_size_inches(10,8)
ax.grid(True)
ax.set_xlabel('1/u')
ax.set_ylabel('s / (Q / (4 pi kD)')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim(1e-3, 1e1)

ttl =  ax.set_title(settings.title)

hds = headsObj.get_alldata()

ir = np.arange(gr.nx, dtype=int)[gr.xm <= settings.r_][-1]

fq = np.zeros(1, dtype=dtypeQ)
Q = 4 * np.pi * settings.kD
fq['I'], fq['q'] = gr.NOD[-1, 0, 0], Q
fq = {0: fq}

grAx = settings.grAx
use_ghb = grAx.nz < 2
kr = settings.props['kr'][:, np.newaxis, np.newaxis] * grAx.const(1.0)
kz = settings.props['kz'][:, np.newaxis, np.newaxis] * grAx.const(1.0)
ss = settings.props['ss'][:, np.newaxis, np.newaxis] * grAx.const(1.0)
hi = grAx.const(0.)
idomain = grAx.const(1, dtype=int)

clrs = color_cycler()
# Plot graph for each rho
for iy, rho in zip(range(gr.ny), settings.rhos):
    clr = next(clrs)
    Ir_ = gr.NOD[-1, iy, ir]
    
    # Modflow
    ax.plot(settings.um1[1:], hds[:, -1, iy, ir], '.', color=clr, label=f'MF6, rho = {rho:.4g}')
    
    # Analytic
    s = Wh(1/settings.um1, rho)[0]
    ax.plot(settings.um1, s, color=clr, label=f'Wh_analytic, rho={rho:.4g}')
    
    
    # Fdm3t
    t = grAx.xm[ir] ** 2 * settings.S / (4 * settings.kD) * settings.um1 # u=r^2S/(4kDt), um1=4kDt/(r^2S), t= um1 r^2 S/(4kD) = um1 (r/beta)^2
    L = grAx.xm[ir] / rho
    ctop = L ** 2  / settings.kD
        
    if use_ghb:
        ghb = np.zeros(grAx.nx * grAx.ny, dtype=dtypeGHB)
        ghb['I'], ghb['h'], ghb['C'] = grAx.NOD[-1].ravel(), 0., grAx.Area.ravel() / ctop
        ghb = {0: ghb}
        c = None
        fh = None
    else:
        fh = np.zeros(grAx.nx * grAx.ny, dtype=dtypeH)
        fh['I'], fh['h'] = grAx.NOD[0].ravel(), 0.
        fh = {0: fh}
        c = gr.const(0.)[:-1]
        c[0] = ctop 
        ghb = None       
        
    out = fdm3t(gr=grAx, t=t, k=(kr, kr, kz), ss=ss, c=c, fh=fh, fq=fq, ghb=ghb, hi=hi, idomain=idomain, epsilon=settings.props['epsilon'])
        
    # um1 = 4 * settings.kD * t  / (settings.S * grAx.xm[ir] ** 2)
    wh = out['Phi'][:, -1, 0, ir] # /(Q  / (4 *np.pi * settings.kD))
        
    ax.plot(settings.um1, wh, 'x',   color=clr, label=f"fdm3t, rho={rho:.4g}")
    
# logo:
ax.text(0.85, 0.05, str(np.datetime64('today')),
        fontsize=10, fontweight='normal', va='bottom', transform=fig.transFigure) # ax.transAxis
ax.legend()

fstr= f"_nt{len(t)-1}_gr{'-'.join([str(s) for s in gr.shape])}_eps{settings.props['epsilon']*100:.0f}perc"

plt.savefig(os.path.join(dirs.images, sim.name + f"{fstr}"+'.png'))

plt.show()
