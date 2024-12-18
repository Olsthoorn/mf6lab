# %% Analyzing output of Pennik Series 2 tests
# TO 090223
# %%
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import lfilter
from scipy.special import exp1
import flopy
import mf_adapt
import settings
from analytic.hantush_convolution import Wh
from etc import color_cycler
from src import mf6tools
from fdm.fdm3t import fdm3t, dtypeGHB, dtypeQ, dtypeH
from fdm.mf6_face_flows import get_struct_flows, get_structured_flows_as_dict # (flowjas, grb_file=None, verbose=False)

use_models, use_packages  = mf6tools.get_models_and_packages_from_excel(
                                                mf_adapt.params_wbk, sheet_name='NAM')

dirs     = mf_adapt.dirs
grb_file = os.path.join(dirs.GWF, mf_adapt.sim_name + 'Gwf.dis.grb')
gr       = mf_adapt.gr
               
# %% load the unformatted files with the heads, the concentration and the budget terms
sim = flopy.mf6.MFSimulation.load(
    sim_name=mf_adapt.sim_name, version='mf6', sim_ws=dirs.SIM)

time_units = sim.tdis.time_units.get_data()

sim_name = sim.name.lower()
gwf = sim.get_model(f'{sim_name}Gwf') # list(sim.model_names)[0])

# %% load the simulation results
headsObj = gwf.output.head()
budObj   = gwf.output.budget()

hds = headsObj.get_alldata()

# %% the date times at the end of each time step
kstpkper = headsObj.get_kstpkper()
sim_time =  np.array(headsObj.get_times())

#  %% Get the flow right face for computing the stream function:
# fflows: dict with face flows (keys: kstpkper with tuples 3 3D arrays: frf, fff and flf') as in MF5 
fflows = get_structured_flows_as_dict(budObj, grb_file=grb_file, verbose=False)

# %% Get suitable levels for contouring the stream function psi [m2/d]
hmin, hmax, hInactive = np.unique(headsObj.get_alldata())[[0, -2, -1]]

# %% Set up initial plot of cross section with contouring
# fdm3t
Icyl = gr.NOD[-1, :, 0]
cellid = gr.I2LRC(Icyl)
fh = np.zeros(len(Icyl), dtype=dtypeH)
fh['I'], fh['h'] =Icyl, settings.hb 
fh = {0: fh}

grAx = settings.grAx
kr = settings.props['kr'][:, np.newaxis, np.newaxis] * grAx.const(1.0)
kz = settings.props['kz'][:, np.newaxis, np.newaxis] * grAx.const(1.0)
ss = settings.props['ss'][:, np.newaxis, np.newaxis] * grAx.const(1.0)
hi = grAx.const(0.)
hi[:, : ,0] = settings.hb
idomain = grAx.const(1, dtype=int)

out = fdm3t(gr=grAx, t=settings.t, k=(kr, kr, kz), ss=ss, fh=fh, hi=hi, idomain=idomain, epsilon=settings.props['epsilon'])

FRF = np.zeros((len(settings.t) - 1, gr.nz, gr.ny, gr.nx))
for i, ksp in enumerate(kstpkper):
    FRF[i] = fflows[ksp][0]['frf']

HDS = headsObj.get_alldata()

# Plotting
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(14, 10))

# Flows Q(r)
ax1.grid(True)
ax1.set_xlabel('t [d]')
ax1.set_ylabel('Q m3/d')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylim(1e-1, 1e5)
ax1.set_ylim(1e2, 1e5)
ax2.grid(True)
ax2.set_xlabel('t [d]')
ax2.set_ylabel('head change [m] m')
#ax2.set_xscale('log')

ttl =  ax1.set_title(settings.title)

ax1.plot(settings.t[1:], FRF[:, 0, 0, 0], '--', color='k', label=f"MF6,   Q[{gr.xm[0]:.4g}] m")
ax2.plot(settings.t[1:], HDS[:, 0, 0, 0], '--', color='k', label=f"MF6,   Q[{gr.xm[0]:.4g}] m")

clrs = color_cycler()
for ir in range(1, gr.nx, 20):
    clr = next(clrs)
    Q = FRF[:, 0, 0, ir]
    ax1.plot(settings.t[1:], Q, '.-', color=clr, label=f"MF6,   Q[{gr.xm[ir]:.4g}] m")
    ax1.plot(settings.t[1:], out['Qx'][:, 0, 0, ir], 'x-', color=clr, label=f"fdm3t: Q[{gr.xm[ir]:.4G}] m")
    print(ir, Q[0], out['Qx'][0, 0, 0, 0])
    
    ax2.plot(settings.t[1:],    HDS[:, -1, 0, ir], '.-', color=clr, label=f"MF6:   r={gr.xm[ir]:.4g} m")
    ax2.plot(settings.t, out['Phi'][:, -1, 0, ir], 'x-', color=clr, label=f"fdm3t: r={gr.xm[ir]:.4g} m")

# logo:

fstr= f"_nt{len(settings.t)-1}_gr{'-'.join([str(s) for s in gr.shape])}_eps{settings.props['epsilon']*100:.0f}perc"

# Approximate by convolution
# Q that generates the desired head change at the radius of the cylinder at any time
# For convolution we must use equal time steps
tau = np.linspace(0, 1000., 1001)
u = gr.x[1] ** 2 * settings.S / (4 * settings.kD * tau)
Q = settings.hb * (4 * np.pi * settings.kD) / Wh(u[1:], 0.)[0]
ax1.plot(tau[1:], Q, 'o', color='k', mfc='none', label="appr. by convolution")

# Compute the head at the cylinder by well in center with Q computed above
SR = np.hstack((0., 1 / (4 * np.pi * settings.kD) * Wh(u[1:], 0.)[0]))
BR = np.hstack((0., SR[1:] - SR[:-1]))
s = lfilter(BR, 1, Q)
ax2.plot(tau[1:], s, 'o', color='k', mfc='none', label="appr. by convolution")

ax1.text(0.85, 0.05, str(np.datetime64('today')),
        fontsize=10, fontweight='normal', va='bottom', transform=fig.transFigure) # ax.transAxis

ax1.legend(fontsize=8)
ax2.legend(fontsize=8)

plt.savefig(os.path.join(dirs.images, sim.name + f"{fstr}"+'.png'))

plt.show()
