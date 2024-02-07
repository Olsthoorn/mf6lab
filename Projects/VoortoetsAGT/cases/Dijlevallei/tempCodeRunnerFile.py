# === load the model with simulation results. =====
sim = flopy.mf6.MFSimulation.load(sim_name=sim_name,
                                  version='mf6',
                                  sim_ws=dirs.SIM)

# === load simulation results =====
gwf = sim.get_model('{}Gwf'.format(sim.name).lower()) # list(sim.model_names)[0])
headsObj = gwf.output.head()
budObj   = gwf.output.budget()

# === Use head from last time step (steady) =====
h = headsObj.get_data(kstpkper=headsObj.get_kstpkper()[-1])
