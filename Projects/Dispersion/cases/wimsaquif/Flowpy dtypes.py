# Flowpy dtypes

# %%
import flopy
import numpy as np

# ========= MF5 ========================
model = flopy.modflow.Modflow("test")
chd = flopy.modflow.ModflowChd(model)
ghb = flopy.modflow.ModflowGhb(model)
drn = flopy.modflow.ModflowDrn(model)
riv = flopy.modflow.ModflowRiv(model)
wel = flopy.modflow.ModflowWel(model)

# "MF5"
dtypes = {
    'chd': np.dtype([('k', '<i8'), ('i', '<i8'), ('j', '<i8'), ('shead', '<f4'), ('ehead', '<f4')]),
    'ghb': np.dtype([('k', '<i8'), ('i', '<i8'), ('j', '<i8'), ('bhead', '<f4'), ('cond', '<f4')]),
    'drn': np.dtype([('k', '<i8'), ('i', '<i8'), ('j', '<i8'), ('elev', '<f4'), ('cond', '<f4')]),
    'riv': np.dtype([('k', '<i8'), ('i', '<i8'), ('j', '<i8'), ('stage', '<f4'), ('cond', '<f4'), ('rbot', '<f4')]),
    'wel': np.dtype([('k', '<i8'), ('i', '<i8'), ('j', '<i8'), ('flux', '<f4')]),
}


#%%
print("MF5")
for bnm, bnd in zip(['chd', 'ghb', 'drn', 'riv', 'wel'],
                     [chd,   ghb,   drn,   riv,   wel]):
    empty_array = bnd.get_empty(1)  # get empty array for 1 record
    print(f"{bnm}  : \n", empty_array.dtype)

# ======== MF6 =====================

# MF6 has other dtypes. However MF6 does not allow entering recarrays
# directly. Which I think is an inconsistence and a great pity.
# Instead it wants the data entered as lists like so:

# # Convert the structured numpy array into a list of tuples
# wel_data_list = [
#     (-100.0, 0),  # Well at cell 0 with a flux of -100.0
#     (-200.0, 4)   # Well at cell 4 with a flux of -200.0
# ]

print("MF6")
dtypes = {
'chd'  : np.dtype([('cellId', '<i8'), ('shead', '<f4'), ('ehead', '<f4')]),
'ghb'  : np.dtype([('cellId', '<i8'), ('bhead', '<f4'), ('cond', '<f4')]),
'drn'  : np.dtype([('cellId', '<i8'), ('elev', '<f4'), ('cond', '<f4')]),
'riv'  : np.dtype([('cellId', '<i8'), ('stage', '<f4'), ('cond', '<f4'), ('rbot', '<f4')]),
'wel'  : np.dtype([('cellId', '<i8'), ('flux', '<f4')]),
}
print(dtypes)



# %% Getting flopy to issue the MF6 dtypes (used chatGPT to do this iteratively)

# Model workspace
ws = "./modflow6_test"
name = "example_model"

# Create simulation and time discretization
sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name="mf6")
tdis = flopy.mf6.ModflowTdis(sim, time_units='DAYS', perioddata=[(1.0, 1, 1)])

# Create groundwater flow model
gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

# Define grid
nlay, nrow, ncol = 1, 3, 3
dis = flopy.mf6.ModflowGwfdis(
    gwf, 
    length_units="meters",
    nlay=nlay, nrow=nrow, ncol=ncol,
    delr=10.0, delc=10.0,
    top=10.0, botm=0.0
)

# Define initial conditions
ic = flopy.mf6.ModflowGwfic(gwf, strt=5.0)

# Node property flow (for hydraulic conductivity)
npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=1.0)

# --- STEP 1: Create WEL package without data to inspect dtype ---
wel = flopy.mf6.ModflowGwfwel(gwf)

# --- STEP 2: Get the dtype of WEL package stress_period_data ---
wel.stress_period_data.set_data({0: [(0, 100.0)]})  # Simple well at cell 0, with a flux of 100
wel_dtype = wel.stress_period_data.dtype
print("WEL dtype:", wel_dtype)

# River package
riv = flopy.mf6.ModflowGwfriv(gwf)
riv.stress_period_data.set_data({0: [(0, 1.0, 0.5, 0.25)]})  # River at cell 0 with stage and flux

# Drain package
drn = flopy.mf6.ModflowGwfdrn(gwf)
drn.stress_period_data.set_data({0: [(7, 1.0, 0.5)]})  # Drain at cell 7 with conductance and stage

# Constant Head package
chd = flopy.mf6.ModflowGwfchd(gwf)
chd.stress_period_data.set_data({0: [(0, 1.0)]})  # Constant head at cell 0 with a specific value

# Now that the stress_period_data is initialized for each package, you can access dtype
riv_dtype = riv.stress_period_data.dtype
drn_dtype = drn.stress_period_data.dtype
chd_dtype = chd.stress_period_data.dtype

print("RIV dtype:", riv_dtype)
print("DRN dtype:", drn_dtype)
print("CHD dtype:", chd_dtype)


# --- STEP 3: Build structured array for WEL using dtype ---
# Let's place wells in cellids 0 and 4 with fluxes -100 and -200
# Define a structured array for WEL
wel_dtype = np.dtype([('q', '<f8'), ('cellid', '<i8')])
wel_data_array = np.array([
    (-100.0, 0),
    (-200.0, 4)
], dtype=wel_dtype)

# Convert the structured array to a list
wel_data_list = wel_data_array.tolist()

print(wel_data_list)

# Internally flopy uses pandas but does not directly accept a DataFrame
# import pandas as pd
# 
# Create a pandas DataFrame for the stress_period_data
# wel_data_df = pd.DataFrame([
    # (-100.0, 0),  # Well at cell 0 with a flux of -100.0
    # (-200.0, 4)   # Well at cell 4 with a flux of -200.0
# ], columns=["q", "cellid"])

# Assign pandas DataFrame to WEL stress_period_data
# wel.stress_period_data.set_data({0: wel_data_df})

# Check if the data was set correctly
print(wel.stress_period_data)

# Assign list of tuples to WEL stress_period_data
wel.stress_period_data.set_data({0: wel_data_list})

# Assign structured array to WEL stress_period_data
# wel.stress_period_data.set_data({0: wel_data_array})

# Check if the data was set correctly
print(wel.stress_period_data)

# Optional: Add output control and solution
oc = flopy.mf6.ModflowGwfoc(
    gwf,
    budget_filerecord=f"{name}.cbc",
    head_filerecord=f"{name}.hds",
    saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]
)

ims = flopy.mf6.ModflowIms(sim, print_option='SUMMARY')
sim.write_simulation()


# %%

# dtype = [('k', int), ('i', int), ('j', int), ('flux', float)]
# data = np.array([
#     (0, 10, 20, -500.0),
#     (0, 11, 21, -600.0),
# ], dtype=dtype)

# wel.stress_period_data = {0: data}
