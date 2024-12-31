# %% Analyzing output of Burrgeman's solution 223_02
# 
# This workout of the solution not only contains the results of Modflow.
# It also constains the results of fdm3t, ttim, the direct integration in
# Bruggemans's solution and the Laplace transform.
# 
# Concolusions are:
#  MF6, ffdm3t and ttim yield essentially the same results, where fdm3t is a bit
# more accurate than MF6 when fdm3t applies epsilon < 1 (i.e. 0.6), where MF6
# only has epsilon = 1 embedded. Fdm3t is close to ttim.
# Numerical back-transformation of Laplace with Graver-Stehfest's algorithm only
# works for Q, not for heads. Other methods have not been tried.
# Direct integration as in the analytical solution given by Bruggeman, only
# works if integrated over 100 logcycles of u, or -100<z<4, where z = exp(u).
# The approximate method by convolution always works and has a reasonble results,
# that should always work in practice.
#
# TO 20241224
# %%
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import lfilter
import flopy
import ttim
import mf_adapt
import settings
from analytic.hantush_convolution import Wh
from etc import color_cycler
import Brug223_02 as b223
from src import mf6tools
from fdm.fdm3t import fdm3t, dtypeH
from fdm.mf6_face_flows import get_structured_flows_as_dict # (flowjas,grb_file=None, verbose=False)

use_models, use_packages  = mf6tools.get_models_and_packages_from_excel(
                                                mf_adapt.params_wbk, sheet_name='NAM')

dirs     = mf_adapt.dirs
grb_file = os.path.join(dirs.GWF, mf_adapt.sim_name + 'Gwf.dis.grb')
gr       = mf_adapt.gr

times = settings.t
        
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

out = fdm3t(gr=grAx, t=times, k=(kr, kr, kz), ss=ss, fh=fh, hi=hi, idomain=idomain, epsilon=settings.props['epsilon'])

FRF = np.zeros((len(times) - 1, gr.nz, gr.ny, gr.nx))
for i, ksp in enumerate(kstpkper):
    FRF[i] = fflows[ksp][0]['frf']

HDS = headsObj.get_alldata()

# TTIM
tmin, tmax = 1e-2, 1e3
rw = settings.props['r'][1]
ml = ttim.ModelMaq(kaq=settings.props['kr'],
                   z=-np.cumsum(np.hstack((0, settings.props['D']))),
                   Saq= settings.props['ss'], tmin=tmin, tmax=tmax)
hwel = ttim.HeadWell(ml, rw=rw, tsandh=[(tmin, settings.hb)], label='well')
ml.solve()

# For analytic stuff
kD = np.sum(settings.props['kr'] * settings.props['D'])
S  = np.sum(settings.props['ss'] * settings.props['D'])
hb = settings.hb
u = np.logspace(-8, 4, 100000)
z = np.log10(u)
R = gr.x[1]
a, b = z[[0, -1]]


# Plotting
fig, (ax2, ax1) = plt.subplots(2, 1, sharex=True, figsize=(14, 10))

# Flows Q(r)
ax1.grid(True)
ax1.set_xlabel('t [d]')
ax1.set_ylabel('Q m3/d')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylim(1e-1, 1e5)
ax2.set_ylim(-1, 2)
ax2.grid(True)
ax2.set_xlabel('t [d]')
ax2.set_ylabel('head change [m] m')
#ax2.set_xscale('log')

suptitle = settings.title + '\n' + f"hb={hb:.4g} m, R={R:.4g} m, S={S:.4g}, kD={kD:.4g} m"
fig.suptitle(suptitle)

ttl =  ax2.set_title("Verloop van de stijghoogten voor verschillende afstanden r")
ttl =  ax1.set_title("Verloop van het debiet voor verschillende afstanden r")
Ittim = np.arange(len(times))[times > 10 * tmin]
tttim = times[Ittim]
Qttim =  -ml.elementdict['well'].discharge(tttim)[0]

ax1.plot(tttim, Qttim, 's', color='k', mfc='none', label=f"TTIM,   Q[{gr.xm[0]:.4g}] m")
ax1.plot(times[1:], FRF[:, 0, 0, 0], '--', color='k', label=f"MF6,   Q[{gr.xm[0]:.4g}] m")
ax2.plot(times[1:], HDS[:, 0, 0, 0], '--', color='k', label=f"MF6,   Q[{gr.xm[0]:.4g}] m")

Nstehfest = 10

clrs = color_cycler()

for ir in range(1, gr.nx, 20):
    if gr.xm[ir] > 10000:
        continue
    
    r, R = gr.xm[ir], gr.x[1]
    
    Fiquad2 = np.zeros_like(times)
    Qquad2 = np.zeros_like(times)
    for it, t in enumerate(times):
        Fiquad2[it] = b223.Fint_z2(z[0], z[-1], hb=hb, r=gr.xm[ir], R=R, S=S, kD=kD, t=t)
        Qquad2[it]  = b223.Qint_z2(z[0], z[-1], hb=hb, r=gr.xm[ir], R=R, S=S, kD=kD, t=t)
    
    args = (hb, r, R, S, kD)
    
    clr = next(clrs)
    # Flows
    Qmf6 = FRF[:, 0, 0, ir]
    ax1.plot(times[1:], Qmf6, '.-', color=clr, label=f"MF6,   Q[{gr.xm[ir]:.4g}] m")
    ax1.plot(times[1:], out['Qx'][:, 0, 0, ir], 'x-', color=clr, label=f"fdm3t: Q[{gr.xm[ir]:.4G}] m")

    #Qttim =  -ml.elementdict['well'].discharge(times)[0]
    qrttim =  ml.disvec(gr.xm[ir], 0., tttim, layers=[0])[0][0]
    Qrttim = 2 * np.pi * gr.xm[ir] * qrttim
    ax1.plot(tttim, Qrttim, '*', color=clr, mfc='none', label=f"TTIM,  Q[{gr.xm[ir]:.4g}] m")
    
    ax1.plot(times, b223.Fback(b223.qhat, times, Nstehfest, args), '.-', color=clr, mfc='none', label=f"laplace, r={gr.xm[ir]:.4g} m")
    
    ax1.plot(times, Qquad2,'s-', color=clr, mfc='none', label=f"Q quad2, r={gr.xm[ir]:.4g}")

    # Heads
    ax2.plot(times[1:],    HDS[:, -1, 0, ir], '.-', color=clr) # , label=f"MF6:   r={gr.xm[ir]:.4g} m")
    ax2.plot(times, out['Phi'][:, -1, 0, ir], 'x-', color=clr) # , label=f"fdm3t: r={gr.xm[ir]:.4g} m")
    
    httimr = ml.head(gr.xm[ir], 0, tttim, layers=0)[0]
    ax2.plot(tttim, httimr, '*', color=clr, mfc='none') #, label=f"fdm3t: r={gr.xm[ir]:.4g} m")
    ax2.plot(times, b223.Fback(b223.fhat, times, Nstehfest, args), '.-', color=clr, mfc='none') #, label=f"laplace    r={gr.xm[ir]:.4g} m")
    
    ax2.plot(times, Fiquad2, 's-', color=clr, mfc='none') # , label=f"Fi quad2, r={gr.xm[ir]:.4g}")
    
# logo:

fstr= f"_nt{len(times)-1}_gr{'-'.join([str(s) for s in gr.shape])}_eps{settings.props['epsilon']*100:.0f}perc"

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

fig.legend(bbox_to_anchor=(0.15, 0.45, 0.2, 0.4), fontsize='small')


plt.savefig(os.path.join(dirs.images, sim.name + f"{fstr}"+'.png'))

plt.show()
