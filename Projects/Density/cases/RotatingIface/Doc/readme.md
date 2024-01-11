# The Rotating Inteface Problem (Modflow Development Team (2023) Example 53)

@TO 20240111

The problem consists of a cross-sectional box filled with three fluids of different densities (fig. 53–1). The initial boundaries between the fluids are not horizontal, and thus, the fluids rotate. Bakker and others (2004) and Langevin and others (2003) compared simulated velocities at the onset of rotation with velocities obtained using an analytical solution. This same comparison with simulated velocities is also shown here for MODFLOW 6.

Regarding updates see Preface of Modflow Development Team (2023):

The document presents example problems for MODFLOW 6. The MODFLOW 6 program can be downloaded from the USGS for free. The performance of MODFLOW 6 has been tested with these examples and other tests. Future testing, however, might reveal errors that were not detected in the test simulations. Updates are routinely made to the MODFLOW 6 program and to these examples. Users can check for updates on the MODFLOW Web page (https://doi.org/10.5066/F76Q1VQV).

Experience
Small time steps are required or else the solution will not be stable and will not converge. Trouble often occurs when increasing time step sizes are used. Then the NO_PCT option in `ims` package may help as does switching the `sto` package on to dampen the solution, but this is not the correct way as the problem lies in too big time steps. The general (`Peclet` number) rule should be applied, i.e. the water must in no cell pass more than one cell in a sing time step. An optimal way is to use a single period with many time steps, i.e. 10000 time steps, and then set the save frequency in the `saverecord` in the `oc` package. This reduces storage and post-processing time and still yields a smooth solution, suitable for animation. Example of the

For the flow model:
Gwfoc = {'head_filerecord':   os.path.join(dirs.SIM, "{}Gwf.hds".format(sim_name)),
         'budget_filerecord': os.path.join(dirs.SIM, "{}Gwf.cbc".format(sim_name)),
         'saverecord': [("HEAD", "FREQUENCY", 100), ("BUDGET", "FREQUENCY", 100)],
}

For the transport model:
Gwtoc = {
      'concentration_filerecord' : os.path.join(dirs.SIM, '{}Gwt.ucn'.format(sim_name)),
      'budget_filerecord':         os.path.join(dirs.SIM, '{}Gwt.cbc'.format(sim_name)),
      'saverecord' : [("CONCENTRATION", "FREQUENCY", 100), ("BUDGET", "FREQUENCY", 100)],
}

combined with

'IPER', 'PERLEN', 'NSTP', 'TSMULT' = [[0, 10000, 10000, 1.0]]

i.e. one stress period with PERLEN 100000 and 100000 time steps with `tsmult=1.0`.

Then setting the save frequency to 100 results in 100 records (every 100 days) in the output, which is about optimal for this problem.

Further, `sto` package is off so quasi steady state solutions are computed that only change due to the change of the density distribution.

Boundary conditions are set using chd.

Either
 1) left and right of the model both head 0 and cR and cL as concentrations.

Or
 2) lower left cell with head 0 and cL as concentration


***References**

Modflow Developement Team (2023) Modflow 6 - Example Problems. Modflow 6 Development Team, with contributions form Chieh Ying Chen and Mike Toews (2023/06/26)

Bakker, M., Oude Essink, G.H., and Langevin, C.D., 2004, The rotating movement of three immiscible fluids—a benchmark problem: Journal of Hydrology, v. 287, no. 1, p. 270 – 278, https://doi.org/10.1016/j.jhydrol.2003.10.007.

Langevin, C.D., Shoemaker, W.B., and Guo, W., 2003, MODFLOW-2000 the U.S. Geological Survey Modular Ground-Water Model–Documentation of the SEAWAT-2000 Version with the Variable-Density Flow Process (VDF) and the Integrated MT3DMS Transport Process (IMT): U.S. Geological Survey Open-File Report 03-426, 43 p., accessed July 25, 2019, at https://pubs.er.usgs.gov/publication/ofr03426.


