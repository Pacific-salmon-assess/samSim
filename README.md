# samSim
## Closed-loop simulation package for salmon CUs

-----  

**Authors: Cameron Freshwater**  
**Date: 2019-02-01 (ONGOING)**

-----

### Summary
This repository contains the necessary files to run stochastic closed-loop simulations parameterized with Pacific salmon stock-recruitment data. The principal function, `recoverySim()`, is intended to test the performance of different management procedures (broadly a mix of harvest control rules and assessment methods) across operating models representing distinct ecological hypotheses. A suite of performance metrics are generated that allow analysts to evaluate different management procedures ability to achieve multiple, interacting conservation- and catch-based objectives. In short, the model is intended to provide a framework for the quantitative component of a management strategy evaluation.

The focal unit of the simulated dynamics are conservation units (CUs) - genetically distinct aggregates of one or more spawning populations that are unlikely to recolonize in the event of extirpation. Under Canada's Wild Salmon Policy these will be the target of future rebuilding strategies. Thus the `samSim` package is well suited to evaluating the rebuilding potential for depleted CUs in a mixed-stock context, but should not be used to evaluate the dynamics of subpopulations within CUs and should only be used to evaluate multiple management units (distinct aggregates of CUs managed quasi-independently) with care.

A summary of relevant files and how to run a simulation are provided below. Most functions contain relatively detailed documentation (and sometimes functioning examples). Details on the operating model (biological dynamics and fishery interactions) and the management procedures (harvest control rule and assessment process) will be provided in a vignette to come.

-----

### Files
All files are stored in the following directories:

#### data
Includes a collection of loose .rda files that are necessary to run the examples associated with certain `samSim` functions. Also includes four subdirectories:

  - *fraserDat* - includes .csv files necessary to run Fraser River sockeye salmon simulations. Currently also includes an "unused" directory that holds raw or deprecated files. 
  - *area3Chum* - includes .csv files necessary to run Area 3 chum salmon simulations. *Note*: migrated directly from `salmon-sim` repo and has not been recently used. May be removed from package unless necessary for integrated examples.
  - *manProcScenarios* - includes example .csv files to run simulations focused on evaluating the performance of different management procedures across a small range of operating models.
  - *opModelScenarios* - includes example .csv files to run simulations focused on evaluating the performance of different oparing models across a small range of management procedures.
  
#### man
Includes the .rd files used to populate help files for each function. Created automatically via `roxygen`.

#### outputs
Directory generated automatically by running `recoverySim()`. Contains a *diagnostics* directory that includes diagnostic plots and *simData* that includes output data files summarizing performance.

#### R
Includes all functions including `recoverySim()` as well as necessary helper and post-processing functions.

#### reports
Includes supplementary documentation describing model development and input data. **NOTE Should probably be removed.**

#### Rmd
Includes Rmarkdown files that include an example simulation run, as well as descriptions model structure and parameterization.

#### src
Includes scripts necessary for several helper C++ functions.

------

### Running a simulation

Simulations are run by installing the samSim package and using the `recoverySim()` function. Parameter values are passed to the function using a series of .csv files with the `simPar` and `cuPar` arguments being most critical. `simPar` contains parameter values that are shared among CUs and define a given scenario (e.g. species, simulation length, OM and MP characteristics). `cuPar` contains parameters that are CU-specific including at a minimum names and SR model type, but typically stock-recruit parameters as well. Details of how to pass a suite of scenarios to the simulation model are provided in the Rmd directory.


