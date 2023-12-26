""" MFLAB script or function that is to become the backbone of Modflow simulation in python

The original mflab was built on Matlab. Because I have no access to Matlab since my emeritat,
I switched to Python after 2015. Python is free and, therefore, generally accessible to students
as well as anyone who's interested in programming.

This mflab.py file is to be the backbone of the Pyton versoin of mflab. It may start as a script,
but will become a module with the essential functionality in it.

Modelfow modelling will make use of Modflow 6 as well as related models, such as the USGS version
of the old mt3dms for transport and seawat for adding density flow to the Modflow.

Use will be made of flopy to generate the input files for Modflow etc. and to
reada out the results. mflab will facilitate model creation by automating much
of the use of flopy, with direct focus on the purpose of the models.

The old mflab in Matlab started in 2008. The current project in November 2023.

Copyright 2009-2023 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
under free software foundation GNU license version 3 or later

Each modflow model together with its transport model shall be contained in a single directory.

Like the old mflab, to generate a model and readout and show the results of it, three files
are requried:

mf_parameters.xlsx --> f"{name}.xlx" where name is the name of the current model.

mf_adapt.py
mf_setup.py
mf_analyze.py

mf_aparameters.xlsx is an Excel workbook in which parameters for the current model are specified.
It has several worksheets in it. MFLOW contains the parameters used by Modflow with their values.
This sheet hardly needs adaptation as most of these parameters are defauls and will be ok for the
current model. Other worksheets are PER, LAY, ... that require adpation to the current model. The
sheet PER defines the stress periods to be used. The sheet LAY defines the layer-wide properties
of the current model.
# TODO: to be extended.

mf_adapt.py should be adapted by the user for every model, as it specifies the current model.
When ready, mf_setup is run, whic takes mf_adapt and generates the input for Modflow while
writing it to disk to the files that will be grabbed by Modflow etc.

When ready, Modflow etc. is automatically run.

After Modflow etc. have (successfully finished), mf_analyze.py is run to show the results.
Like mf_adapt.py, mf_analyze.py is adapted by the user to analyze and visualize the model results
in the way required by the user.

The user should not have to change mf_setup.py.
"""
import numpy as np
import flopy
import logging
import mf_adapt
from mf6lab import mf6tools
logging.basicConfig(level=logging.WARNING, format=' %(asctime)s - %(levelname)s - %(message)s')

NOT = np.logical_not
AND = np.logical_and
OR  = np.logical_or

# Ssame for PER and LAY

gr = mf_adapt.gr

def mf_setup():

    #params_wbk = os.path.join(mf_adapt.HOME, mf_adapt.sim_name + '.xlsx')

    use_models, use_packages  = mf6tools.get_models_and_packages_from_excel(
                                                mf_adapt.params_wbk, sheet_name='NAM')
    
    model_dict = mf6tools.get_mf6_params_from_excel(
                                    mf_adapt.params_wbk, sheet_name='SIM6')

    if 'Gwf' in use_models:
        model_dict.update(mf6tools.get_mf6_params_from_excel(
                                    mf_adapt.params_wbk, sheet_name='GWF6'))
    if 'Gwt' in use_models:
        model_dict.update(mf6tools.get_mf6_params_from_excel(
                                    mf_adapt.params_wbk, sheet_name='GWT6'))

    fp_packages=dict()

### === S I M U L A T I O N =========================================

### sim =====================
    logging.info("sim")
    model_dict['Simsim'].update(sim_name=mf_adapt.sim_name, sim_ws=mf_adapt.dirs.SIM)
    sim = flopy.mf6.MFSimulation(**model_dict['Simsim'])
    fp_packages['Simsim'] = sim

### tdis ====================
    logging.info("Simtdis")
    model_dict['Simtdis'].update(**mf_adapt.Simtdis)
    fp_packages['Simtdis'] = flopy.mf6.ModflowTdis(sim, **model_dict['Simtdis'])
            
### ======== F L O W ================================================

    if 'Gwf' in use_models:
    ### Gwf =====================
        if True: # Groundwater flow model, use name from 'sim').
            logging.info("Gwf")
            Gwf_model_name = sim.name + 'Gwf'
            model_dict['Gwfgwf'].update(modelname=Gwf_model_name,
                                    model_rel_path=mf_adapt.dirs.GWF,
                                    exe_name=sim.exe_name,
                                    )
            gwf = flopy.mf6.ModflowGwf(sim, **model_dict['Gwfgwf'])
            fp_packages['Gwfgwf'] = gwf
            
            ### Ims =====================
            model_dict['Gwfims'].update()
            logging.info("Gwfims")
            Gwfims = flopy.mf6.ModflowIms(sim, **model_dict['Gwfims'])
            fp_packages['Gwfims'] = Gwfims
            sim.register_ims_package(Gwfims, [Gwf_model_name])    
    

    ### Gwfdis ==================
        if 'Gwfdis' in use_packages:
            logging.info('Gwfdis')
            model_dict['Gwfdis'].update(**mf_adapt.Gwfdis)
            gr = model_dict['Gwfdis'].pop('gr')
            model_dict['Gwfdis'].update(nlay=gr.nlay, nrow=gr.nrow, ncol=gr.ncol,
                                        delr=gr.dx, delc=gr.dy,
                                        top=gr.Z[0], botm=gr.Z[1:]                              
            )
            fp_packages['Gwfdis'] = flopy.mf6.ModflowGwfdis(gwf, **model_dict['Gwfdis'])

    ### Gwfdisu ==================
        if 'Gwfdisu' in use_packages:
            logging.info('Gwfdisu')
            model_dict['Gwfdisu'].update(**mf_adapt.Gwfdisu)
            fp_packages['Gwfdisu'] = flopy.mf6.ModflowGwfdisu(gwf, **model_dict['Gwfdisu'])
    
    ### Gwfdisv ==================
        if 'Gwfdisv' in use_packages:
            logging.info('Gwfdisv')
            model_dict['Gwfdisv'].update(**mf_adapt.Gwfdisv)
            fp_packages['Gwfdisv'] = flopy.mf6.ModflowGwfdisv(gwf, **model_dict['Gwfdisv'])
    
    
    ### Gwfbuy ==================
        if 'Gwfbuy' in use_packages:
            logging.info('Gwfbuy')
            model_dict['Gwfbuy'].update(**mf_adapt.Gwfbuy)
            fp_packages['Gwfbuy'] = flopy.mf6.ModflowGwfbuy(gwf, **model_dict['Gwfbuy'])
    
    ### Gwfchd ==================
        if 'Gwfchd' in use_packages:
            logging.info('Gwfchd')
            timeseries = mf_adapt.Gwfchd.pop('timeseries', None)
            model_dict['Gwfchd'].update(**mf_adapt.Gwfchd)
            fp_packages['Gwfchd'] = flopy.mf6.ModflowGwfchd(gwf, **model_dict['Gwfchd'])
            Gwfchd = fp_packages['Gwfchd']
            if timeseries:
                Gwfchd.ts.initialize(**timeseries.pop(0))
                while timeseries:
                    Gwfchd.ts.append_package(**timeseries.pop(0))
    
    ### Gwfcsub ==================
        if 'Gwfcsub' in use_packages:
            logging.info('Gwfcsub')
            model_dict['Gwfcsub'].update(**mf_adapt.Gwfcsub)
            fp_packages['Gwfcsub'] = flopy.mf6.ModflowGwfcsub(gwf, **model_dict['Gwfcsub'])
    
    ### Gwfdrn ==================
        if 'Gwfdrn' in use_packages:
            logging.info('Gwfdrn')
            model_dict['Gwfdrn'].update(**mf_adapt.Gwfdrn)
            fp_packages['Gwfdrn'] = flopy.mf6.ModflowGwfdrn(gwf, **model_dict['Gwfdrn'])
    
    ### Gwfevt ==================
        if 'Gwfevt' in use_packages:
            logging.info('Gwfevt')
            model_dict['Gwfevt'].update(**mf_adapt.Gwfevt)
            fp_packages['Gwfevt'] = flopy.mf6.ModflowGwfevt(gwf, **model_dict['Gwfevt'])
    
    ### Gwfghb ==================
        if 'Gwfghb' in use_packages:
            logging.info('Gwfghb')
            model_dict['Gwfghb'].update(**mf_adapt.Gwfghb)
            fp_packages['Gwfghb'] = flopy.mf6.ModflowGwfghb(gwf, **model_dict['Gwfghb'])
    
    ### Gwfgnc ==================
        if 'Gwfgnc' in use_packages:
            logging.info('Gwfgnc')
            model_dict['Gwfgnc'].update(**mf_adapt.Gwfgnc)
            fp_packages['Gwfgnc'] = flopy.mf6.ModflowGwfgnc(gwf, **model_dict['Gwfgnc'])
    
    ### Gwfhfb ==================
        if 'Gwfhfb' in use_packages:
            logging.info('Gwfhfb')
            model_dict['Gwfhfb'].update(**mf_adapt.Gwfhfb)
            fp_packages['Gwfhfb'] = flopy.mf6.ModflowGwfhfb(gwf, **model_dict['Gwfhfb'])
    
    ### Gwfic ==================
        if 'Gwfic' in use_packages:
            logging.info('Gwfic')
            model_dict['Gwfic'].update(**mf_adapt.Gwfic)
            fp_packages['Gwfic'] = flopy.mf6.ModflowGwfic(gwf, **model_dict['Gwfic'])
    
    ### Gwflak ==================
        if 'Gwflak' in use_packages:
            logging.info('Gwflak')
            model_dict['Gwflak'].update(**mf_adapt.Gwflak)
            fp_packages['Gwflak'] = flopy.mf6.ModflowGwflak(gwf, **model_dict['Gwflak'])
    
    ### Gwfmaw ==================
        if 'Gwfmaw' in use_packages:
            logging.info('Gwfmaw')
            model_dict['Gwfmaw'].update(**mf_adapt.Gwfmaw)
            fp_packages['Gwfmaw'] = flopy.mf6.ModflowGwfmaw(gwf, **model_dict['Gwfmaw'])
    
    ### Gwfmvr ==================
        if 'Gwfmvr' in use_packages:
            logging.info('Gwfmvr')
            model_dict['Gwfmvr'].update(mf_adapt.Gwfmvr)
            fp_packages['Gwfmvr'] = flopy.mf6.ModflowGwfmvr(gwf, **model_dict['Gwfmvr'])
    
    ### Gwfnpf ==================
        if 'Gwfnpf' in use_packages:
            logging.info('Gwfnpf')
            model_dict['Gwfnpf'].update(**mf_adapt.Gwfnpf)
            fp_packages['Gwfnpf'] = flopy.mf6.ModflowGwfnpf(gwf, **model_dict['Gwfnpf'])
    
    ### Gwfobs ==================
        if 'Gwfobs' in use_packages:
            logging.info('Gwfobs')
            model_dict['Gwfobs'].update(**mf_adapt.Gwfobs)
            fp_packages['Gwfobs'] = flopy.mf6.ModflowGwfobs(gwf, **model_dict['Gwfobs'])
    
    ### Gwfoc ==================
        if 'Gwfoc' in use_packages:
            logging.info('Gwfoc')
            model_dict['Gwfoc'].update(**mf_adapt.Gwfoc)
            fp_packages['Gwfoc'] = flopy.mf6.ModflowGwfoc(gwf, **model_dict['Gwfoc'])
    
    ### Gwfrch ==================
        if 'Gwfrch' in use_packages:
            logging.info('Gwfrch')
            model_dict['Gwfrch'].update(**mf_adapt.Gwfrch)
            fp_packages['Gwfrch'] = flopy.mf6.ModflowGwfrch(gwf, **model_dict['Gwfrch'])
    
    ### Gwfriv ==================
        if 'Gwfriv' in use_packages:
            logging.info('Gwfriv')
            model_dict['Gwfriv'].update(**mf_adapt.Gwfriv)
            fp_packages['Gwfriv'] = flopy.mf6.ModflowGwfriv(gwf, **model_dict['Gwfriv'])
    
    ### Gwfsfr ==================
        if 'Gwfsfr' in use_packages:
            logging.info('Gwfsfr')
            model_dict['Gwfsfr'].update(mf_adapt.Gwfsfr)
            fp_packages['Gwfsfr'] = flopy.mf6.ModflowGwfsfr(gwf, **model_dict['Gwfsfr'])
    
    ### Gwfsto ==================
        if 'Gwfsto' in use_packages:
            logging.info('Gwfsto')
            model_dict['Gwfsto'].update(**mf_adapt.Gwfsto)
            fp_packages['Gwfsto'] = flopy.mf6.ModflowGwfsto(gwf, **model_dict['Gwfsto'])
    
    ### Gwfuzf ==================
        if 'Gwfuzf' in use_packages:
            logging.info('Gwfuzf')
            model_dict['Gwfuzf'].update(**mf_adapt.Gwfuzf)
            fp_packages['Gwfuzf'] = flopy.mf6.ModflowGwfuzf(gwf, **model_dict['Gwfuzf'])
    
    ### Gwfvsc ==================
        if 'Gwfvsc' in use_packages:
            logging.info('Gwfvsc')
            model_dict['Gwfvsc'].update(**mf_adapt.Gwfvsc)
            fp_packages['Gwfvsc'] = flopy.mf6.ModflowGwfvsc(gwf, **model_dict['Gwfvsc'])
    
    ### Gwfwel ==================
        if 'Gwfwel' in use_packages:
            logging.info('Gwfwel')
            model_dict['Gwfwel'].update(**mf_adapt.Gwfwel)
            fp_packages['Gwfwel'] = flopy.mf6.ModflowGwfwel(gwf, **model_dict['Gwfwel'])


    ### ===== T R A N S P O R T =================================================

    if 'Gwt' in use_models:
    ### Gwt ===============================
        if True: # Groundwater transport model, use name from 'sim').
            logging.info("Gwtgwt")
            Gwt_model_name = sim.name + 'Gwt'
            model_dict['Gwtgwt'].update(modelname=Gwt_model_name,
                                    model_rel_path=mf_adapt.dirs.GWT,
                                    exe_name=sim.exe_name,
                                    )
            gwt = flopy.mf6.ModflowGwt(sim, **model_dict['Gwtgwt'])
            fp_packages['Gwtgwt'] = gwt
            
            ### ims =====================
            model_dict['Gwtims'].update(complexity='MODERATE')  # SIMPLE | MODERATE | COMPLEX
            logging.info("Gwtims")
            Gwtims = flopy.mf6.ModflowIms(sim, **model_dict['Gwtims'])
            fp_packages['Gwtims'] = Gwtims
            sim.register_ims_package(Gwtims, [Gwt_model_name])
            
    ### Gwtdis ==================
        if 'Gwtdis' in use_packages:
            logging.info('Gwtdis')
            model_dict['Gwtdis'].update(nlay=gr.nlay, nrow=gr.nrow, ncol=gr.ncol,
                        delr=gr.dx, delc=gr.dy,
                        top=gr.Z[0], botm=gr.Z[1:],
                        length_units=mf_adapt.LENGTH_UNITS,
                        idomain=mf_adapt.IDOMAIN,                                 
            )
            fp_packages['Gwtdis'] = flopy.mf6.ModflowGwtdis(gwt, **model_dict['Gwtdis'])
    
    ### Gwtdisu ==================
        if 'Gwtdisu' in use_packages:
            logging.info('Gwtdisu')
            model_dict['Gwtdisu'].update()
            fp_packages['Gwtdisu'] = flopy.mf6.ModflowGwtdisu(gwt, **model_dict['Gwtdisu'])
    
    ### Gwtdisv ==================
        if 'Gwtdisv' in use_packages:
            logging.info('Gwtdisv')
            model_dict['Gwtdisv'].update()
            fp_packages['Gwtdisv'] = flopy.mf6.ModflowGwtdisv(gwt, **model_dict['Gwtdisv'])
            
    ### Gwtadv ==================
        if 'Gwtadv' in use_packages:
            logging.info('Gwtadv')
            model_dict['Gwtadv'].update(scheme=mf_adapt.SCHEME)
            fp_packages['Gwtadv'] = flopy.mf6.ModflowGwtadv(gwt, **model_dict['Gwtadv'])
    
    ### Gwtcnc ==================
        if 'Gwtcnc' in use_packages:
            logging.info('Gwtcnc')
            model_dict['Gwtcnc'].update(stress_period_data=mf_adapt.CONSTCONC,
                                        maxbound = len(mf_adapt.CONSTCONC))
            fp_packages['Gwtcnc'] = flopy.mf6.ModflowGwtcnc(gwt, **model_dict['Gwtcnc'])
    
    ### Gwtdsp ==================
        if 'Gwtdsp' in use_packages:
            logging.info('Gwtdsp')
            model_dict['Gwtdsp'].update(**mf_adapt.DISPERSIVITIES)
            fp_packages['Gwtdsp'] = flopy.mf6.ModflowGwtdsp(gwt, **model_dict['Gwtdsp'])
    
    ### Gwtfmi ==================
        if 'Gwtfmi' in use_packages:
            logging.info('Gwtfmi')
            model_dict['Gwtfmi'].update(**mf_adapt.Gwtfmi)
            fp_packages['Gwtfmi'] = flopy.mf6.ModflowGwtfmi(gwt, **model_dict['Gwtfmi'])
    
    ### Gwtic ==================
        if 'Gwtic' in use_packages:
            logging.info('Gwtic')
            model_dict['Gwtic'].update(strt=mf_adapt.STRTC)
            fp_packages['Gwtic'] = flopy.mf6.ModflowGwtic(gwt, **model_dict['Gwtic'])
    
    ### Gwtist ==================
        if 'Gwtist' in use_packages:
            logging.info('Gwtist')
            model_dict['Gwtist'].update()
            fp_packages['Gwtist'] = flopy.mf6.ModflowGwtist(gwt, **model_dict['Gwtist'])
    
    ### Gwtlkt ==================
        if 'Gwtlkt' in use_packages:
            logging.info('Gwtlkt')
            model_dict['Gwtlkt'].update()
            fp_packages['Gwtlkt'] = flopy.mf6.ModflowGwtlkt(gwt, **model_dict['Gwtlkt'])
    
    ### Gwtmst ==================
        if 'Gwtmst' in use_packages:
            logging.info('Gwtmst')
            model_dict['Gwtmst'].update(**mf_adapt.Gwtmst)
            fp_packages['Gwtmst'] = flopy.mf6.ModflowGwtmst(gwt, **model_dict['Gwtmst'])
    
    ### Gwtmvt ==================
        if 'Gwtmvt' in use_packages:
            logging.info('Gwtmvt')
            model_dict['Gwtmvt'].update()
            fp_packages['Gwtmvt'] = flopy.mf6.ModflowGwtmvt(gwt, **model_dict['Gwtmvt'])
    
    ### Gwtmwt ==================
        if 'Gwtmwt' in use_packages:
            logging.info('Gwtmwt')
            model_dict['Gwtmwt'].update()
            fp_packages['Gwtmwt'] = flopy.mf6.ModflowGwtmwt(gwt, **model_dict['Gwtmwt'])
    
    ### Gwtobs ==================
        if 'Gwtobs' in use_packages:
            logging.info('Gwtobs')
            model_dict['Gwtobs'].update()
            fp_packages['Gwtobs'] = flopy.mf6.ModflowGwtobs(gwt, **model_dict['Gwtobs'])
    
    ### Gwtoc ==================
        if 'Gwtoc' in use_packages:
            logging.info('Gwtoc')
            model_dict['Gwtoc'].update(**mf_adapt.GWTOC)
            fp_packages['Gwtoc'] = flopy.mf6.ModflowGwtoc(gwt, **model_dict['Gwtoc'])
    
    ### Gwtsft ==================
        if 'Gwtsft' in use_packages:
            logging.info('Gwtsft')
            model_dict['Gwtsft'].update()
            fp_packages['Gwtsft'] = flopy.mf6.ModflowGwtsft(gwt, **model_dict['Gwtsft'])
    
    ### Gwtsrc ==================
        if 'Gwtsrc' in use_packages:
            logging.info('Gwtsrc')
            model_dict['Gwtsrc'].update()
            fp_packages['Gwtsrc'] = flopy.mf6.ModflowGwtsrc(gwt, **model_dict['Gwtsrc'])
    
    ### Gwtssm ==================
        if 'Gwtssm' in use_packages:
            logging.info('Gwtssm')
            model_dict['Gwtssm'].update()
            fp_packages['Gwtssm'] = flopy.mf6.ModflowGwtssm(gwt, **model_dict['Gwtssm'])
    
    ### Gwtuzt ==================
        if 'Gwtuzt' in use_packages:
            logging.info('Gwtuzt')
            model_dict['Gwtuzt'].update()
            fp_packages['Gwtuzt'] = flopy.mf6.ModflowGwtuzt(gwt, **model_dict['Gwtuzt'])



### ==== Dynamic exchange ===========================================
        
    if 'Gwf' in use_models and 'Gwt' in use_models:
        logging.info("Dynamic exchange GWF-GWT active")
        model_dict['Gwfexc'].update(exgtype='GWF6-GWT6',
                                    exgmnamea=model_dict['Gwfgwf']['modelname'],
                                    exgmnameb=model_dict['Gwtgwt']['modelname'])
        exchange = flopy.mf6.ModflowGwfgwt(sim, **model_dict['Gwfexc'])
        fp_packages['Gwfexc'] = exchange


    return fp_packages, model_dict, use_models, use_packages


def run_modflow(sim=None):
        """Simulate GGOR using MODFLOW.

        Parameters
        ----------
        sim: flopy.mf6.MFSimulaition object
            completely defined model.
        """
        # Write simulation
        sim.write_simulation(silent=mf_adapt.user_settings['silent'])
        
        # Run simulation
        success, buff = sim.run_simulation(silent=mf_adapt.user_settings['silent'])

        print('Running success = {}'.format(success))
        if not success:
            print(buff)
            logging.critical("Buffer printed because MODFLOW did not terminate normally.")
            raise Exception('MODFLOW did not terminate normally.')
        return success



if __name__ == '__main__':
    
    fp_packages, model_dict, use_models, use_packages = mf_setup()
    
    success = run_modflow(sim=fp_packages['sim'])
        
    print(mf_adapt.sim_name)
    