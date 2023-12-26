**GGOR** is acronym for *Desired Groundwater and Surface Water Regime* (*Gewenst Grond- en Oppervlaktewater Regime*)

The **GGOR tool** simulates the transient development of *shallow groundwater simultaneously in a large number of parcels*, which are influenced by precipitation surplus, vertical seepage exchange with a regional aquifer and verying level of the surface water.

In the GGOR each parcel is treated as as a **cross section of a shallow aquifer bounded by two parallel ditches** with equal but possibly varying surface water level. Each parcel is assumed to be *symmetric about the center line between its bounding ditches*.The parcels exchange water with the underlying regional aquifer, which is goverend by their head difference. The parcels receive a time-varying rechare (precipitaion minus evapotranspiration) and the regional aquifer may be subject to injection or extraction, which will cause seepage either upward (positive) or downward (negative).

The properties of the parcels, i.e. their conductivity, vertical resistance, specific yield or elastic storativity and resistance between aquifer and ditches differ.

The GGOR tool accepts a time series of precipitation and evapotranspiration, varying ditch levels, injection/extraction into/from the regional aquifer and the data pertaining to all individual parcels, extracted from a large GIS file.

The simulations are done using the numerical model *ModflowUSG*, which allows to obtain simulated time series for any point in any of the parcels, and by *analytical solutions* but only for the *averaged head in the parcels between the ditches*, yielding time-variing parcel-averaged head and flows. Simulating the heads and flows both numerically and analytically allows comparing the two methods, of which the analytical should be faster with many simultaneous parcels, and the numerical more accurate because it requires less assumptions, but with similar outcomes. The outcomes can't be exactly the same, because the time varying analytic solution, in fact, uses successive steady-state solutions and is, therefore, an approximation. Of course, the numeric solution is an approximation as well, but it can be made more accurate by refining the space and time resolution.

The GGOR analytical solutions implemented and tested in the notebook in this directory provides tools to analyze the transient development of groundwater heads in the cover layer of the first layer of an aquifer system inwhich the shallow groundwater is affacted by precipiation surplus and by seepage (upward possitve) by which it interacts with a more regional aquifer.

@T.N.Olsthoorn 209-01-21