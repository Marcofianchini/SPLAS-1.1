# SPLAS 1.1 #

SPLAS 1.0 (Sea PLAnkton Simulator) is a size-structure resolving N-P-Z (nutrient-phytoplankton-zooplankton) model, based on two simple assumptions: 
fundamental plankton physiological traits scale with organism size and grazing by zooplankton is structured by an optimal predator-prey length ratio
(along with a complex log-gaussian preference function in this case).
The model is based on empirical allometric relationships used to parameterize all major plankton characteristics (taken from available literature).
70 size-classes are defined for phytoplankton, in the range 1-200 μm, and there are 70 matching size-classes for zooplankton,
in the range 2-1000 μm. Here phytoplankton are considered all autotrophic, uptaking Nitrogen and growing according to an internal dynamics (cellular quota).
The cellular quota is very important as it allows to link phytoplankton biomass with the density of cell, which is also an important derived variable of the model,
allowing for a conversion from Nitrogen to Carbon biomass and fluxes.
Zooplankton are all heterotrophic, feeding on phytoplankton through a size-based preference. No mixotrophy is represented inside SPLAS.
Despite the idealized approach, the 0-D version of the model is able to catch the overall pattern of plankton size distribution and diversity, and reproduces correctly the bottom-up and top-down interactions,
showing an emergent pattern of plankton size-structure consistent with theoretical and model prediction, e.g. Ward et al. (2013).
SPLAS was also extended to a 1-D framework (adding the vertical dimension) with further complications as detritus compartments, light effects and turbulent vertical diffusion.
With 1-D SPLAS we obtained sensible vertical profiles of phytoplankton total biomass along with depth-inhomogeneous size-structures 
(which are not possible to reproduce without light inclusion).
We also obtained a reasonable interval of estimates for Carbon export (central value about 68 mgC m-2 day-1).
For further details, see the PhD Thesis by Marco De Pasquale (2018) "The role of plankton size in community structure, biodiversity and biogeochemical fluxes: a modeling approach". It is present in this repository.
If you would like to use this code for further development or modifications/adaptations to your scope, please give credits to this thesis.
Marco put all of his effort in developing this stuff, so please give the right credit to this work.
Thank you! :) 

SPLAS 1.1 - inclusion of light input (izero),the possibility to increase/decrease light input linearly in time (delta) and the possibility to make it normally distributed(choosing the standard deviation,sd, parameter)
For further details ask to fianchini.marco@libero.it 
### Set-up of the model ###

Technical notes:

* Fortran 2003 mathematical libraries + Python 3.5 I/O interface
* Makefile compiles Fortran libraries and then python3.5 calls the main program (main.py)
* 2 possible configurations or sets of dependencies: 
	* intel environment (compile with ifort)
		- intel/pe-xe-2016--binary
 	    - intelmpi/5.1.3--binary
 		- python/3.5.0
 		- mkl/11.3.3--binary

	* gnu environment (compile with gfortran)
		- mpi4py/1.3.1--openmpi--1.8.2--gnu--4.8.3
		- lapack/3.5.0--gnu--4.8.3
		- blas/3.5.0--gnu--4.8.3
		- python/3.5.0
		- openmpi/1.8.2--gnu--4.8.3
		
* Some profiling has been performed using vtune/16.3 and code has been optimized to a certain extent

Few further notes:

* SPLAS is already tested and validated for the 0-D part
* 1-D part could be further developed and improved (link with biogeochemistry and C export should be enforced in particular!),
  but it works well and gives rise to interesting patterns (although still far from reality)

### Contribution guidelines ###

Please refer to the documentation cited above (and here recalled):

* For further details see the PhD Thesis by Marco De Pasquale, "The role of plankton size in community structure, biodiversity and biogeochemical fluxes: a modeling approach" (2018).
  If you would like to use this code for further development or modifications/adaptations to your scope, please give credits to this work.
* For the technical part, try to use the libraries reported in the section above ("Set-up of the model") which make up the two possible environments
  that were used to conduct the experiments reported in the thesis. Since such experiments were finished at the end of 2017, if using those libraries is not possible for some reason,
  or if updates with breaking changes have been released meantime...well, good luck with that!
  UPDATE 05/2021 it still works but the warning for the future remains. 

### Contacts ###

* For any issue or question, admin of this project are: Fianchini Marco, Marco De Pasquale
* Contact by E-mail: fianchini.marco@libero.it , depasquale.marco@gmail.com
