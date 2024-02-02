# mf6lab (follow-up of mflab)

Groundwater modeling using **flopy** and **Modflow 6** (**MF6**)

`Modflow 6` changed the world for groundwater modelers. It promises a way into the future by being completely modular and extensible with new packages, and the ability to combine with more and diverse models with which it can exchange data.

Modflow 6 integrates most of the Modflow developments made in the past, which often led to separate versions of Modflow. One of the advances of **MF6** is to allow for unstructured grids. But also internally **MF6** has been modernized and improved. For instance, the use of Newton solutions was used to solve instabilities with non-linear problems. The new **mf6** also soles the problem with rewetting of dry cells that have so long plagued geneartions of Modflow users because it often prevented to reach convergence.

Furthermore, **md6** contains both the flow model and the transport model (previously MT3DMS was used for this) and can handle variable density and viscosity (replacing the previous *seawat* model).

**Modflow** can now be conveniently interacted with by **python cripting** through **flopy**, which was also newly developed by a team of people at the USGS together with scholars from over the world.

**Mf6lab** builds on **floy** and **MF6** by standardizing setup and analysis of groundwater models with **flopy** and **m6**, thus streamlining its use in an overlying framework. **mf6lab** is applied in many projects. More will follow as I will use **mf6lab** for all new groundwater modelling projects in the future.

Examples can be found under **mf6lab/projects/<project>/cases/<case>**

The idea central to **mf6lab** is to have a backbone file mf_setup.py that never changes which creates a dictionary with the data for all packages to be used in the simulation a standard file ***run.py*** which invokes it and the files `mf_adapt.py`, `mf_analyze.py` and `<case_name>.xlsx` which are specfiic for each model. Note however, that for most models the visualization of the simulation results, which is wat `mf_analyze.py` does is essentially the same for most models. `mf_adapt` is used to set the parameters for the packages to be used. It can conveniently import a dedicated file `settings.py` in which certain parameters are stored together and read from the Excel workbook file `<case_name>.xlsx` which has some sheets to guide the model. The first sheet is named **NAM**. This is where you tell wich flow and transport packages to use for the current simulation. You do this by setting the  value of appropriate cells to 1 or zero. The transport model will be automatically invoked when one or more transport packages are switched to on. The second sheet is **SIM6** it holds the parameters and their default values of the simulation packages. The third sheet is named **GWF6**. It holds the parameters of all availabe packages with their default values for the flow model. As long as your are satisfied with them, you don't have to change them, but you can do it here. Alternatively you can overwrite any parameter in `mf_adapt.py`, whatever is more natural to you. In any case, the sheet gives a complete overview of all parameters and their default values. mf_adapt.py will swallow this entire sheet to generate a dictionary with the parameters for each of the available packages. The next fourth sheet is named **GWT6**. It has the same structure as the previous sheets and holding all parametes and their default values pertaining to the tranport model. This too will seldom be changed. The fifth sheet is called **PER**, it specifies the stress periods. It is a table with at least the following header **IPER** **PERLEN** **NSTEPS** **TSMULT** needed to tell **mf6** about the stress periods and subdivision for this simulatioin. One can add as many columns as desired to this table, which may then be utilized in `mf_adapt.py'. For instnace, in which stress period some action takes place. The sheet will be read in by mf_adapt.py and the resulting pandas DataFrame is yours to use at your liking. Only the mentioned columns are requied for **mf6**. A further sheet ofthen useful is **LAY** to specify layers and their properties. This sheet may be use if deemed useful. For instance, one may specify for each layer its conductivities, storativities, dispersivities as well as some code and geologic description. The sheet can be read in by `mf_adapt.py` and use to provide data to prepare the input, but data in the resulting DataFrame can also be used to set parameters for visualization such as itmes like formation names to be used in the legend. This is completely up to you. 

Often it is convenient to create an exta file like `settings.py` to store divers simulation parametes on one place, and to set up, test and visualize the network. This file can be imported by 'mf_adapt.py' and by `mf_analyze.py` to gain access to the parameters in it from a central place. Wether to use such an extra file `settings.py` is up to the user, but it turns out to be convenient in many circumstances as a central place to find the parameters to fine-tune a given siulation-run.

Thus, having prepared `mf_adapt.py` run `run.py` to let flopy generate the **mf6** files for this simulation. The smiulation is automatically started. It's progress can be followed on the screen. When successfully finished run `mf_analyze.py` to read, analyze, visualize, print and store the results. That's it.

For a new run change a few parameters which may be convenientlu stored in `settings.py`, launch run.py again and afterwards launch `mf_analyze.py` again to get the new results.

Some actions of done is:
switch from steady state to transient --> just swich on the Gwfsto pacakge in the NAM sheet of the Excel file and run again.
to switch on the `Gwfchd` package ---> just switch this package on in the `NAM` sheet of the Excel file.
to use the transport model ---> switch on the approprate Gwt-packages in the `NAM` sheet. Of course, in `mf_adapt.py` you set the appropriate parameters for these packages.
to change the stress periods, time steps and tsmult--> change them in the `PER` sheet of the Excdel file.
to change layer properties ---> change them in the `LAY` sheet of the Excel file.
to use density flow --> switch on the Gwfbuy package in the NAM shet of the Excel file.

The actual input tha **MF6** is genearted in a number of files by **flopy**. They can be inspected in the folder ***GWF** of the current case. The files for the input generated by **flopy** for the transport simulation can be found in the folder **GWT** of the current case. The files for the input of the current simulaiton under which the flow and the tranport model resied can be found in the filder **SIM** of the current case. This folder also receives the output files that **MF6** generated i.e. the heads (`<case>.hds`) and budget (`<case>.cbc`) files and possibly the concentration (`<case>.unc`) files.

To know which packages are needed in a specific model and whether or not the groundwar transport model is to be run, **mf6lab** reads the sheet NAM in the Excel file `<case>.xlsx`. Of course one can add as many special purpose sheets to this Excel file as one desires for specific actions or data in this simulation.

See the **mf6lab/projects** fiolders for the different projects. Under each project see the cases folder for the different individual cases. Note that every case has the same folders (more can be added of desired or removed if not needed). Gernally the visualization outcome is stored in the folder **Images** (i.e. **fm6lab/Projects/<project>/cases/<case>/Images**).

@TOlsthoorn 20231226, 20240202

Theo Olsthoorn (prof. dr.ir. TN), emeritus prof TUDelft, eritus hydrologist at Waternet Amsterdam.
Heemstede 2024/02/02


TO 20231225
