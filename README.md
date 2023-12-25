# mf6lab
Groundwater modeling using FloPy and Modflow 6 (MF6)

Modflow 6 changed the world of groundwater modeling. It promises a way into the future by being completely modular and extensible with new packages, but also with combining with more and diverse models which can exchange data with the grounwater models.

It integrates most of the developments made in the past and which often led to separate version of Modflow. One of the advances is to allow for unstructured grids. But also internally, the use of Newton solutions and solving the everlasting problem with rewetting of cells that so often prevented convergence.

Furthermore md6 contains both the flow model, the transport model and can handle variable density and viscosity.

Modflow itself can now be conveniently addressed by python cripting through flopy.

Mf6lab builds on floy and modflow 6 in making setup and analysis of the modeling with flopy/m6 more standardized and, therefore, more streamlined and convenient. It is applied in many projects. More will follow as I will use mf6lab for all new groundwater modelling projects in the future.

The idea central to mf6lab is for each model to create two files named mf_adapt and mf_analyze.
mf_adapt provides the data specific to the model (soil properties, recharge etc).
then mf_run is run, which dalls mf_setup. mf_setup imports mf_adapt thorugh which it gets the data pertaining to the current model and then instantiates the flopy.mf6 modules, after which mf6 is called and run.
After finishing run mf_analyze to analyze and visualize the results.
