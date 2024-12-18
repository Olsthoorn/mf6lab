# The purpose of this case Hantush

The idea is to compare the accuracy of fdm3t with Modflow and an analytic solution for which the Hantush well function is used.

The Hantush well function is computed for a number of rho values.

In Modflow a 1 layer 7 row model is made to compute the results for 7 rho values at once (in each row one rho value). This is done by setting k22 to 1e-20 and using the correct conductance for each row and column.
Axial symmetry is achieved by multiplying kr, kz and ss by 2 pi r and using C = pi (r[1:]^2 - r[:-1]^2) / c als concuctance, wich c = (r / rho)^2 / kD.

Furthermore, "alternative_k_averaging" is set to "logarithmic" in Modflow to get the highest accuracy for the Axial case. fdm3t does this automatically and more sophisticately.

For fdm3t an axial grid is used.

The results of MF, fdm3t and the Hantush well function are similar and practially equal when sufficient r-values are used. I.e. >= 20 per logcycle seems ok. With 10 r values per log cycle, differences occur for large rho values.

fdm3t is more accurate than MF6 when less than 20 r-values per log cycle are used.


@TO 20241218
