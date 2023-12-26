SimPkgs = { # first arg is sim
        'Simsim' : "MFSimulation",
        'Simgwf' : "ModflowGwf",
        'Simgwt' : "ModdlowGwt",
        'Simgnc' : "Sim gnc (ghost node correction) package",
        'Simims' : "Sim ims (solver package",
        'Simmvr' : "Sim mover package",
        'Simnam' : "Sim nam",
        'Simtdis': "Sim tdis package",
}


GwfPkgs = { # first arg is gwf
        'Gwtgwf' : "ModflowGwf",  # Same as 'Simgwt', first arg is sim not gwf
        'Gwfgwt': "ModflowGwf", # First arg is sim, not gwf
        'Gwfexc' : "Exchange, Modflowgwfgwt - GwfGwt exchange"
        'Gwfdis': "Rectilinear Discretization Input File",
        'Gwfdisu': "Unstructured Discretization Input File",
        'Gwfdisv': "Discretization by Vertices Input File",
        'Gwfbuy': "Buoyancy Package",
        'Gwfchd': "Time-Variant Specified Head Option *",
        'Gwfcsub': "Compaction and Subsidence Package",
        'Gwfdrn': "Drain Package *",
        'Gwfevt': "Evapotranspiration Package *",
        'Gwfghb': "General-Head Boundary Package *",
        'Gwfgnc': "Ghost-Node Correction Package",
        'Gwfhfb': "Horizontal Flow Barrier Package",
        'Gwfic': "Initial Conditions Package",
        'Gwflak': "Lake Package *",
        'Gwfmaw': "Multi-Aquifer Well Package *",
        'Gwfmvr': "Water Mover Package",
        'Gwfnpf': "Node Property Flow Package",
        'Gwfobs': "Observations Option",
        'Gwfoc': "Output Control Option",
        'Gwfrch': "Recharge Package *",
        'Gwfriv': "River Package *",
        'Gwfsfr': "Streamflow Routing Package *",
        'Gwfsto': "Storage Package",
        'Gwfuzf': "Unsaturated Zone Flow Package *",
        'Gwfvsc': "Viscosity Package",
        'Gwfwel': "Well Package *",  
        }

GwtPkgs = { # first arg is gwt
        'Gwtgwt' : "ModdlowGwt", # Same as 'Simgwt', first arg is sim not gwt
        'Gwtdis': "Spatial Discretization structured grid",
        'Gwtdisu': "Spatial Discretization unatrutured",
        'Gwtdisv': "Spatial Discretization with in layer unstructured",
        'Gwtadv': "Advection",
        'Gwtcnc': "Constant Concentration",
        'Gwtdsp': "Dispersion",
        'Gwtfmi': "Flow Model Interface",
        'Gwtic': "Initial Conditions",
        'Gwtist': "Immobile Storage and Transfer",
        'Gwtlkt': "Lake Transport",
        'Gwtmst': "Mobile Storage and Transfer",
        'Gwtmvt': "Mover Transport",
        'Gwtmwt': "Multi-Aquifer Well Transport",
        'Gwtobs': "Model Observations",
        'Gwtoc': "Output Control",
        'Gwtsft': "Streamfow Transport",
        'Gwtsrc': "Mass Source Loading",
        'Gwtssm': "Source-Sink Mixing",
        'Gwtuzt': "Unsaturated Zone Transport",
        }


Mf6Models ={
    'Simsim': "sim_name='sim', version='mf6', exe_name='mf6.exe', sim_ws='.', verbosity_level=1, continue_=None, nocheck=None, memory_print_option=None, write_headers=True, lazy_io=False",
    
    'Simgwf': "modelname='model', model_nam_file=None, version='mf6', exe_name='mf6', model_rel_path='.', list=None, print_input=None, print_flows=None, save_flows=None, newtonoptions=None, packages=None",
    
    'Simgwt': "modelname='model', model_nam_file=None, version='mf6', exe_name='mf6', model_rel_path='.', list=None, print_input=None, print_flows=None, save_flows=None, packages=None",
    
    'Simgnc': "loading_package=False, print_input=None, print_flows=None, explicit=None, numgnc=None, numalphaj=None, gncdata=None, filename=None, pname=None",
    
    'Simims': "loading_package=False, print_option=None, complexity=None, csv_output_filerecord=None, csv_outer_output_filerecord=None, csv_inner_output_filerecord=None, no_ptcrecord=None, ats_outer_maximum_fraction=None, outer_hclose=None, outer_dvclose=None, outer_rclosebnd=None, outer_maximum=None, under_relaxation=None, under_relaxation_gamma=None, under_relaxation_theta=None, under_relaxation_kappa=None, under_relaxation_momentum=None, backtracking_number=None, backtracking_tolerance=None, backtracking_reduction_factor=None, backtracking_residual_limit=None, inner_maximum=None, inner_hclose=None, inner_dvclose=None, rcloserecord=None, linear_acceleration=None, relaxation_factor=None, preconditioner_levels=None, preconditioner_drop_tolerance=None, number_orthogonalizations=None, scaling_method=None, reordering_method=None, filename=None, pname=None",
    
    'Simmvr': "loading_package=False, print_input=None, print_flows=None, modelnames=None, budget_filerecord=None, budgetcsv_filerecord=None, maxmvr=None, maxpackages=None, packages=None, perioddata=None, filename=None, pname=None",
    
    'Simnam': "loading_package=False, continue_=None, nocheck=None, memory_print_option=None, maxerrors=None, print_input=None, tdis6=None, models=None, exchanges=None, mxiter=None, solutiongroup=None, filename=None, pname=None",
    
    'Simtdis': "loading_package=False, time_units=None, start_date_time=None, ats_perioddata=None, nper=1, perioddata=((1.0,1,1.0),), filename=None, pname=None",

}


GwfSignature = {
'Gwfgwf': "modelname='model', model_nam_file=None, version='mf6', exe_name='mf6', model_rel_path='.', list=None, print_input=None, print_flows=None, save_flows=None, newtonoptions=None, packages=None",

# Dynamic exchange GWF-GWT
'GwfExc': "loading_package=False, exgtype='GWF6-GWT6', exgmnamea=None, exgmnameb=None, filename=None, pname=None",

'Gwfdis': "loading_package=False, length_units=None, nogrb=None, xorigin=None, yorigin=None, angrot=None, nlay=1, nrow=2, ncol=2, delr=1.0, delc=1.0, top=1.0, botm=0.0, idomain=None, filename=None, pname=None, parent_file=None",

'Gwfdisu': "loading_package=False, length_units=None, nogrb=None, xorigin=None, yorigin=None, angrot=None, nlay=None, ncpl=None, nvert=None, top=None, botm=None, idomain=None, vertices=None, cell2d=None, filename=None, pname=None",

'Gwfdisv': "loading_package=False, length_units=None, nogrb=None, xorigin=None, yorigin=None, angrot=None, nlay=None, ncpl=None, nvert=None, top=None, botm=None, idomain=None, vertices=None, cell2d=None, filename=None, pname=None, parent_file=None",

'Gwfbuy': "loading_package=False, hhformulation_rhs=None, denseref=1000.0, density_filerecord=None, dev_efh_formulation=None, nrhospecies=None, packagedata=None, filename=None, pname=None, parent_file=None",

'Gwfchd': "loading_package=False, auxiliary=None, auxmultname=None, boundnames=None, print_input=None, print_flows=None, save_flows=None, timeseries=None, observations=None, maxbound=None, stress_period_data=None, filename=None, pname=None, parent_file=None",

'Gwfdrn': "loading_package=False, auxiliary=None, auxmultname=None, auxdepthname=None, boundnames=None, print_input=None, print_flows=None, save_flows=None, timeseries=None, observations=None, mover=None, maxbound=None, stress_period_data=None, filename=None, pname=None, parent_file=None",

'Gwfevt': "loading_package=False, fixed_cell=None, auxiliary=None, auxmultname=None, boundnames=None, print_input=None, print_flows=None, save_flows=None, timeseries=None, observations=None, surf_rate_specified=None, maxbound=None, nseg=None, stress_period_data=None, filename=None, pname=None, parent_file=None",

'Gwfevta': "loading_package=False, readasarrays=True, fixed_cell=None, auxiliary=None, auxmultname=None, print_input=None, print_flows=None, save_flows=None, timearrayseries=None, observations=None, ievt=None, surface=0.0, rate=0.001, depth=1.0, aux=None, filename=None, pname=None, parent_file=None",

'Gwfghb': "loading_package=False, auxiliary=None, auxmultname=None, boundnames=None, print_input=None, print_flows=None, save_flows=None, timeseries=None, observations=None, mover=None, maxbound=None, stress_period_data=None, filename=None, pname=None, parent_file=None",

'Gwfgnc': "loading_package=False, auxiliary=None, auxmultname=None, boundnames=None, print_input=None, print_flows=None, save_flows=None, timeseries=None, observations=None, mover=None, maxbound=None, stress_period_data=None, filename=None, pname=None, parent_file=None",

'Gwfhfb': "loading_package=False, print_input=None, maxhfb=None, stress_period_data=None, filename=None, pname=None, parent_file=None",

'Gwfic': "loading_package=False, strt=1.0, filename=None, pname=None, parent_file=None",

'Gwflak': "loading_package=False, auxiliary=None, boundnames=None, print_input=None, print_stage=None, print_flows=None, save_flows=None, stage_filerecord=None, budget_filerecord=None, package_convergence_filerecord=None, timeseries=None, observations=None, mover=None, surfdep=None, time_conversion=None, length_conversion=None, nlakes=None, noutlets=None, ntables=None, packagedata=None, connectiondata=None, tables=None, outlets=None, perioddata=None, filename=None, pname=None, parent_file=None",

'Gwfmaw': "loading_package=False, auxiliary=None, boundnames=None, print_input=None, print_head=None, print_flows=None, save_flows=None, head_filerecord=None, budget_filerecord=None, no_well_storage=None, flow_correction=None, flowing_wells=None, shutdown_theta=None, shutdown_kappa=None, timeseries=None, observations=None, mover=None, nmawwells=None, packagedata=None, connectiondata=None, perioddata=None, filename=None, pname=None, parent_file=None",

'Gwfmvr': "loading_package=False, print_input=None, print_flows=None, modelnames=None, budget_filerecord=None, maxmvr=None, maxpackages=None, packages=None, perioddata=None, filename=None, pname=None, parent_file=None",

'Gwfnam': "loading_package=False, list=None, print_input=None, print_flows=None, save_flows=None, newtonoptions=None, packages=None, filename=None, pname=None, parent_file=None",

'Gwfnpf': "loading_package=False, save_flows=None, alternative_cell_averaging=None, thickstrt=None, cvoptions=None, perched=None, rewet_record=None, xt3doptions=None, save_specific_discharge=None, save_saturation=None, k22overk=None, k33overk=None, icelltype=0, k=1.0, k22=None, k33=None, angle1=None, angle2=None, angle3=None, wetdry=None, filename=None, pname=None, parent_file=None",

'Gwfoc': "loading_package=False, budget_filerecord=None, head_filerecord=None, headprintrecord=None, saverecord=None, printrecord=None, filename=None, pname=None, parent_file=None",

'Gwfrch': "loading_package=False, fixed_cell=None, auxiliary=None, auxmultname=None, boundnames=None, print_input=None, print_flows=None, save_flows=None, timeseries=None, observations=None, maxbound=None, stress_period_data=None, filename=None, pname=None, parent_file=None",

'Gwfrcha': "loading_package=False, readasarrays=True, fixed_cell=None, auxiliary=None, auxmultname=None, print_input=None, print_flows=None, save_flows=None, timearrayseries=None, observations=None, irch=None, recharge=0.001, aux=None, filename=None, pname=None, parent_file=None",

'Gwfriv': "loading_package=False, auxiliary=None, auxmultname=None, boundnames=None, print_input=None, print_flows=None, save_flows=None, timeseries=None, observations=None, mover=None, maxbound=None, stress_period_data=None, filename=None, pname=None, parent_file=None",

'Gwfsfr': "loading_package=False, auxiliary=None, boundnames=None, print_input=None, print_stage=None, print_flows=None, save_flows=None, stage_filerecord=None, budget_filerecord=None, package_convergence_filerecord=None, timeseries=None, observations=None, mover=None, maximum_picard_iterations=None, maximum_iterations=None, maximum_depth_change=None, unit_conversion=None, nreaches=None, packagedata=None, connectiondata=None, diversions=None, perioddata=None, filename=None, pname=None, parent_file=None",

'Gwfsto': "loading_package=False, save_flows=None, storagecoefficient=None, ss_confined_only=None, iconvert=0, ss=1e-05, sy=0.15, steady_state=None, transient=None, filename=None, pname=None, parent_file=None",

'Gwfuzf': "loading_package=False, auxiliary=None, auxmultname=None, boundnames=None, print_input=None, print_flows=None, save_flows=None, wc_filerecord=None, budget_filerecord=None, package_convergence_filerecord=None, timeseries=None, observations=None, mover=None, simulate_et=None, linear_gwet=None, square_gwet=None, simulate_gwseep=None, unsat_etwc=None, unsat_etae=None, nuzfcells=None, ntrailwaves=7, nwavesets=40, packagedata=None, perioddata=None, filename=None, pname=None, parent_file=None",

'Gwfwel': "loading_package=False, auxiliary=None, auxmultname=None, boundnames=None, print_input=None, print_flows=None, save_flows=None, auto_flow_reduce=None, timeseries=None, observations=None, mover=None, maxbound=None, stress_period_data=None, filename=None, pname=None, parent_file=None",
}

GwtSignature = {
'Gwtgwt' : "modelname='model', model_nam_file=None, version='mf6', exe_name='mf6', model_rel_path='.', list=None, print_input=None, print_flows=None, save_flows=None, packages=None",

'Gwtdis': "loading_package=False, length_units=None, nogrb=None, xorigin=None, yorigin=None, angrot=None, nlay=1, nrow=2, ncol=2, delr=1.0, delc=1.0, top=1.0, botm=0.0, idomain=None, filename=None, pname=None, parent_file=None",

'Gwtdisu': "loading_package=False, length_units=None, nogrb=None, xorigin=None, yorigin=None, angrot=None, vertical_offset_tolerance=0.0, nodes=None, nja=None, nvert=None, top=None, bot=None, area=None, idomain=None, iac=None, ja=None, ihc=None, cl12=None, hwva=None, angldegx=None, vertices=None, cell2d=None, filename=None, pname=None, parent_file=None",

'Gwtdisv': "loading_package=False, length_units=None, nogrb=None, xorigin=None, yorigin=None, angrot=None, nlay=None, ncpl=None, nvert=None, top=None, botm=None, idomain=None, vertices=None, cell2d=None, filename=None, pname=None",

'Gwtadv': "loading_package=False, scheme=None, filename=None, pname=None, parent_file=None",

'Gwtcnc': "loading_package=False, auxiliary=None, auxmultname=None, boundnames=None, print_input=None, print_flows=None, save_flows=None, timeseries=None, observations=None, maxbound=None, stress_period_data=None, filename=None, pname=None, parent_file=None",

'Gwtdsp': "loading_package=False, xt3d_off=None, xt3d_rhs=None, diffc=None, alh=None, alv=None, ath1=None, ath2=None, atv=None, filename=None, pname=None, parent_file=None",

'Gwtfmi': "loading_package=False, save_flows=None, flow_imbalance_correction=None, packagedata=None, filename=None, pname=None, parent_file=None",

'Gwtic': "loading_package=False, strt=0.0, filename=None, pname=None, parent_file=None",

'Gwtist': "loading_package=False, save_flows=None, sorption=None, first_order_decay=None, zero_order_decay=None, cim_filerecord=None, cimprintrecord=None, cim=None, thetaim=None, zetaim=None, decay=None, decay_sorbed=None, bulk_density=None, distcoef=None, filename=None, pname=None, parent_file=None",

'Gwtlkt': "loading_package=False, flow_package_name=None, auxiliary=None, flow_package_auxiliary_name=None, boundnames=None, print_input=None, print_concentration=None, print_flows=None, save_flows=None, concentration_filerecord=None, budget_filerecord=None, timeseries=None, observations=None, packagedata=None, lakeperioddata=None, filename=None, pname=None, parent_file=None",

'Gwtmst': "loading_package=False, save_flows=None, first_order_decay=None, zero_order_decay=None, sorption=None, porosity=None, decay=None, decay_sorbed=None, bulk_density=None, distcoef=None, sp2=None, filename=None, pname=None, parent_file=None",

'Gwtmvt': "loading_package=False, print_input=None, print_flows=None, save_flows=None, budget_filerecord=None, filename=None, pname=None, parent_file=None",

'Gwtmwt': "loading_package=False, flow_package_name=None, auxiliary=None, flow_package_auxiliary_name=None, boundnames=None, print_input=None, print_concentration=None, print_flows=None, save_flows=None, concentration_filerecord=None, budget_filerecord=None, timeseries=None, observations=None, packagedata=None, mwtperioddata=None, filename=None, pname=None, parent_file=None",

'Gwtnam': "loading_package=False, list=None, print_input=None, print_flows=None, save_flows=None, packages=None, filename=None, pname=None, parent_file=None",

'Gwtoc': "loading_package=False, budget_filerecord=None, concentration_filerecord=None, concentrationprintrecord=None, saverecord=None, printrecord=None, filename=None, pname=None, parent_file=None",

'Gwtsft': "loading_package=False, flow_package_name=None, auxiliary=None, flow_package_auxiliary_name=None, boundnames=None, print_input=None, print_concentration=None, print_flows=None, save_flows=None, concentration_filerecord=None, budget_filerecord=None, timeseries=None, observations=None, packagedata=None, reachperioddata=None, filename=None, pname=None, parent_file=None",

'Gwtsrc': "loading_package=False, auxiliary=None, auxmultname=None, boundnames=None, print_input=None, print_flows=None, save_flows=None, timeseries=None, observations=None, maxbound=None, stress_period_data=None, filename=None, pname=None, parent_file=None",

'Gwtssm': "loading_package=False, print_flows=None, save_flows=None, sources=None, filename=None, pname=None, parent_file=None",

'Gwtuzt': "loading_package=False, flow_package_name=None, auxiliary=None, flow_package_auxiliary_name=None, boundnames=None, print_input=None, print_concentration=None, print_flows=None, save_flows=None, concentration_filerecord=None, budget_filerecord=None, timeseries=None, observations=None, packagedata=None, uztperioddata=None, filename=None, pname=None, parent_file=None",
}

def verify_keys(set1, set2):
    set1 = set(set1)
    set2 = set(set2)
    return set1.difference(set2)


def sort_keys(oldObj):    
    newObj = {}
    for k in sorted(oldObj.keys()):
        newObj[k] = oldObj[k]
    return newObj

def print_sorted(dictObj):
    for k, datum in dictObj.items():
        print('\'{}\': "{}",'.format(k, datum))
        #print()
     
def print_for_excel(dictObj):   
    for key, keys_values in dictObj.items():
        key_val_list = keys_values.split(', ')
        for item in key_val_list:
            k, v = item.split('=')
            if k == 'save_flows':
                v = 'True'
            print('{}\t{}\t{}'.format(key, k, v))
        print()
        
print_for_excel(Mf6Models)

def gwf_template(Pkgs):
    for p in Pkgs.keys():
        print(
 '''### {0} ==================
    if '{0}' in use_packages:
        logging.info('{0}')
        model_dict['{0}'].update()
        fp_packages['{0}'] = flopy.mf6.Modflow{0}({1}, **model_dict['{0}'])
 '''.format(p, p[:3].lower())
        )

def print_pkgs(Pkgs):
    for k, pkg in Pkgs.items():
        print('{}    {}    "{}"'.format(k, k[3:], pkg))

#print_for_excel(Mf6Models)
print_for_excel(GwfSignature)
#print_for_excel(GwtSignature)


#gwf_template(GwfPkgs)        
#gwf_template(GwtPkgs)

#for k in Mf6Models:
#    print(k)

#print_pkgs(GwfPkgs)
#print_pkgs(GwtPkgs)
