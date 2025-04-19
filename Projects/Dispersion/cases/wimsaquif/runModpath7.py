### === Run Modpath7 ===

# %% === Imports ===

import os
import numpy as np
import matplotlib.pyplot as plt
import flopy
import settings

pr = settings.props
gr = pr['gr']
sim_name = pr['sim_name']

# === 1. Settings and paths ===
dirs = settings.dirs

os.chdir(dirs.case)

# Workspaces voor GWF en MP7
sim_ws = dirs.SIM
mp7_ws = dirs.MP7

for path in [sim_ws, mp7_ws]:
    if not os.path.isdir(path):
        raise FileNotFoundError(f"Workspace does not exist: {path}")

mf6_exe = '/Users/Theo/GRWMODELS/USGS/mf6.4.2_mac/bin/mf6' # mf6 executable
if not os.path.isfile(mf6_exe):
    raise FileNotFoundError(f"Can't open file <{mf6_exe}>")

mp7_exe = '/Users/Theo/GRWMODELS/USGS/modpath_7_2_001/bin/mpath7gf.mac' # mp7 executable
if not os.path.isfile(mp7_exe):
    raise FileNotFoundError(f"Can't open file <{mp7_exe}>")

# %% === 2. Load existing MODFLOW 6 model ===

sim = flopy.mf6.MFSimulation.load(sim_name=sim_name,
                                  version="mf6",
                                  exe_name=mf6_exe,
                                  sim_ws=sim_ws)

# Copy de tdis file van SIM directory naar de GWF directory
# while removing the line starting with START_DATE_TIME

infile = os.path.join(dirs.SIM, sim_name + '.tdis')
outfile = os.path.join(dirs.GWF, sim_name + '.tdis')
with open(infile, 'r') as f:
    lines = f.readlines()
with open(outfile, 'w') as f:
    for line in lines:
        if not line.strip().upper().startswith("START_DATE_TIME"):
            f.write(line)

# Get the flow model object
gwf_model_name = sim.model_names[0]
flowmodel = sim.get_model(gwf_model_name)

assert isinstance(flowmodel, flopy.mf6.ModflowGwf), "flowmodel is not a GWF model"
print(f"Flow model name: {flowmodel.name}, workspace: {flowmodel.model_ws}")


# %% === 3. Create and configure MODPATH7 model ===

mp7_model_name = sim_name + 'mp7'

mp7 = flopy.modpath.Modpath7(
    modelname=mp7_model_name,
    flowmodel=flowmodel,
    exe_name=mp7_exe,
    model_ws=mp7_ws,
)

# %% === Create a MODPATH7 object ===

os.chdir(mp7_ws)

mpbas = flopy.modpath.Modpath7Bas(
    mp7,
    porosity=0.35,
)

# %% === 4. Define particles and simulation options ===

psi = np.load(os.path.join(dirs.GWF, "psi.npy"))
psi_max = psi[0, 0]
num_particles = 1000

# Psi values marking the regions of dPsi
psi_ = np.linspace(0, psi_max, num_particles + 1)

# Psi values for the release of the particles
psim_ = 0.5 * (psi_[:-1] + psi_[1:])

Zp = np.interp(psim_, psi[::-1, 0], gr.Z[::-1, 0, 0])[::-1]
Ip = np.interp(Zp, gr.Z[::-1, 0, 0], np.arange(gr.nlay + 1)[::-1])
Iz = np.asarray(Ip, dtype=int)
localz = 1 - (Ip - Iz)

particlelocs = [ (int(iz), gr.nrow //  2, 0) for iz in Iz]
num_particles = len(particlelocs)
particleids = np.arange(num_particles, dtype=int)

def check_particle_location(lrc_tuple, grid):
    k, i, j = lrc_tuple
    assert 0 <= k < grid.nlay, f"Invalid layer index: {k}"
    assert 0 <= i < grid.nrow, f"Invalid row index: {i}"
    assert 0 <= j < grid.ncol, f"Invalid col index: {j}"

for loc in particlelocs:
    check_particle_location(loc, gr)

p = flopy.modpath.ParticleData(partlocs=particlelocs,
                               structured=True,
                               particleids=particleids,
                               localx=0.5,
                               localy=0.5,
                               localz=localz,
                               timeoffset=0.,
                               drape=0.
                               )
                            
pg1 = flopy.modpath.ParticleGroup(particlegroupname='PG1',
                                  filename='PG1.ptcls',
                                  releasedata=0.0,
                                  particledata=p)

# %% === Generate the mpsim package ===

mpsim = flopy.modpath.Modpath7Sim(
    mp7,
    simulationtype="pathline", # valid ["endpoint", "pathline", "combined", "timeseries"]
    trackingdirection="forward", # "backward"
    weaksinkoption="pass_through",
    weaksourceoption="pass_through",
    budgetoutputoption="no",
    budgetcellnumbers=None,
    traceparticledata=None,
    referencetime=[0, 0, 0.0],
    stoptimeoption="extend",
    timepointdata=None,
    zonedataoption="off",
    zones=0,
    particlegroups=pg1
)

# %% === 5. Run Modpath7 ===

mp7.write_input()  # Writes the MODPATH input file

# --- Corrigeer het .mpnam-bestand naar relatieve paden (nodig voor gfortran) ---
mpnam_path = os.path.join(mp7.model_ws, f"{mp7_model_name}.mpnam")
with open(mpnam_path, "r") as f:
    lines = f.readlines()

corrected_lines = []
for line in lines:
    parts = line.strip().split()
    if len(parts) == 2:
        label, filepath = parts
        rel_path = os.path.relpath(os.path.abspath(filepath), start=mp7.model_ws)
        corrected_lines.append(f"{label:<10s} {rel_path}\n")
    else:
        corrected_lines.append(line)
        
with open(mpnam_path, "w") as f:
    f.writelines(corrected_lines)
# --------------------------------------------------------------------------------

# === run modpath ===
mp7.run_model()    # Runs MODPATH7 simulation

#%%  === 6. Extract results ===
fpth = os.path.join(mp7_ws, mp7_model_name + ".mppth")
fend = os.path.join(mp7_ws, mp7_model_name + ".mpend")

if not os.path.isfile(fpth):
    raise FileNotFoundError(f"Can't fine file <{fpth}>.")
if not os.path.isfile(fend):
    raise FileNotFoundError(f"Can't fine file <{fend}>.")

pl = flopy.utils.PathlineFile(fpth)
ep = flopy.utils.EndpointFile(fend)
# Get all pathline data as a list

data_pl = pl.get_alldata()
data_ep = ep.get_alldata()

# %% === Compute and plot breakthrough curves at given observation points x_obs ===
xR = 10 ** int(np.log10(gr.x[-1]))
xObs = np.arange(xR / 10, xR + 1., xR / 10)

dtype = np.dtype([('id', '<i4'), ('xObs', '<f8', len(xObs)), ('time', '<f8', len(xObs))])

pass_times = np.zeros(len(data_pl), dtype=dtype)

for i, arr in enumerate(data_pl):
    pass_times['id'][i] = arr['particleid'][0]
    pass_times['xObs'][i] = xObs
    pass_times['time'][i] = np.interp(xObs, arr['x'], arr['time'])
 
fpt ='pass_times_' + pr['k_field_pars']['name'] + '.npy'
print(f"Saving  pass_times array to <{fpt}>")
np.save(fpt, pass_times)
 
# Fig header line:
kfield_str = pr['k_field_pars']['k_field_str']

# === Plot the cumulative curves:
fig, ax = plt.subplots(figsize=(10, 6))
ax.set(title="Cumulative breakthrough curves.\n" + kfield_str, xlabel="time [d]", ylabel="Frac")
for j, xo in enumerate(xObs):
    idx = np.argsort(pass_times['time'][:, j])    
    t = pass_times['time'][:, j][idx]
    ax.plot(t, psim_ / psi_.max(), label=f"x = {xObs[j]} m")
ax.grid()
ax.legend()
# ax.set_xlim(ax.get_xlim()[0], 112000.)
fig.savefig(os.path.join(dirs.images, "cumulatives_" + pr['k_field_pars']['name'] +".png"))

# === Plot the cumulative curves:
fig, ax = plt.subplots(figsize=(10, 6))
ax.set(title="Histograms of breakthrough curves.\n" + kfield_str, xlabel="time [d]", ylabel="Frac")
for j, xo in enumerate(xObs):
    idx = np.argsort(pass_times['time'][:, j])    
    t = pass_times['time'][:, j][idx]
    ax.hist(t, bins=30, alpha=0.5, ec='k', label=f"x = {xObs[j]} m")
ax.grid()
ax.legend()
# ax.set_xlim(ax.get_xlim()[0], 112000.)
fig.savefig(os.path.join(dirs.images, "histograms_" + pr['k_field_pars']['name'] + '.png'))

# %% === Plotting the pathlines ===
fig, ax = plt.subplots(figsize=(10, 6))
ax.set(title="Cross section particle tracking, pathlines.\n" + kfield_str,
       xlabel="x [m]", ylabel="z [m]")

mm = flopy.plot.PlotCrossSection(ax=ax, model=flowmodel, line={"row": 0}, extent=None)

# mm.plot_grid(lw=0.5)

# Use a colormap to assign a color to each particle
colors = plt.cm.viridis(np.linspace(0, 1, len(data_pl)))

for i, xz in enumerate(data_pl):
    if np.mod(i, 5) == 0:
        ax.plot(xz['x'], xz['z'], label=f'particle {i}', color=colors[i])

ax.legend(loc='best', title="Particle ID")
fig.savefig(os.path.join(dirs.images, "pathlines_" + pr['k_field_pars']['name'] + '.png'))

# Plot the k_field (ln(k_field))
fig, ax = plt.subplots(figsize=(10, 6))
ax.set(title=f"ln(k)field, vertical cross section, shape={gr.shape}\n" +
           kfield_str, xlabel='x [m]', ylabel='z [m]')
kplot = ax.imshow(np.log(pr['k_field']), aspect='auto', extent=pr['extent'])
plt.colorbar(kplot, label="ln(k)", orientation='horizontal', ax=ax)
fig.savefig(os.path.join(dirs.images, "kfield_" + pr['k_field_pars']['name'] + '.png'))
plt.show()

# %%

