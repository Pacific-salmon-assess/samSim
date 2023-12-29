# samSim
## Closed-loop simulation package for salmon CUs - Southern BC chinook populations



### Summary
This repository contains the necessary files to run stochastic closed-loop simulations parameterized with Pacific salmon stock-recruitment data. The principal function, `recoverySim()`, is intended to test the performance of different management procedures (broadly a mix of harvest control rules and assessment methods) across operating models representing distinct ecological hypotheses. A suite of performance metrics are generated that allow analysts to evaluate different management procedures ability to achieve multiple, interacting conservation- and catch-based objectives. In short, the model is intended to provide a framework for the quantitative component of a management strategy evaluation.

The focal unit of the simulated dynamics are conservation units (CUs) - genetically distinct aggregates of one or more spawning populations that are unlikely to recolonize in the event of extirpation. Under Canada's Wild Salmon Policy these will be the target of future rebuilding strategies. Thus the `samSim` package is well suited to evaluating the rebuilding potential for depleted CUs in a mixed-stock context, but should not be used to evaluate the dynamics of subpopulations within CUs and should only be used to evaluate multiple management units (distinct aggregates of CUs managed quasi-independently) with care.

A summary of relevant files and how to run a simulation are provided below. Most functions contain relatively detailed documentation (and sometimes functioning examples). Details on the operating model (biological dynamics and fishery interactions) and the management procedures (harvest control rule and assessment process) will be provided in a vignette to come.

[After installing the package development prerequisites](https://support.rstudio.com/hc/en-us/articles/200486498) install `samSim` with:

```{r}
devtools::install_github("Pacific-salmon-assess/samSim")
```

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
Directory generated automatically by running `recoverySim()`. Contains a *diagnostics* directory that includes diagnostic plots and *simData* that includes output data files summarizing performance. Note that in practice this specific output directory will not be used by the analyst because its contents are only populated when the source code is ran or when the function is run within the `samSim` package, which is not recommended. **However** an equivalent directory is automatically generated by `samSim` when it is first run in a new working directory (see *Running a simulation* below for clarification).

#### R
Contains `recoverySim()` function as well as necessary helper and post-processing functions.

#### reports
Includes supplementary documentation describing model development and input data.

#### Rmd
Includes Rmarkdown files that include an example simulation run, as well as descriptions model structure and parameterization.

#### src
Includes scripts necessary for several helper C++ functions.

------

### Running a simulation

Simulations are run by installing the samSim package and using the `recoverySim()` function. Generally this should occur in a fresh working directory (e.g. a new `.Rproj`), which will be automatically generate an `outputs` directory and necessary subdirectories. Parameter values are passed to the function using a series of .csv files with the `simPar` and `cuPar` arguments being most critical. `simPar` contains parameter values that are shared among CUs and define a given scenario (e.g. species, simulation length, OM and MP characteristics). `cuPar` contains parameters that are CU-specific including at a minimum names and SR model type, but typically stock-recruit parameters as well. See *Input file details* below. Details of how to pass a suite of scenarios to the simulation model are provided in `Rmd/exampleSimRun.Rmd`.

------

### Input file details

#### `simPar`
 A detailed descrption of the contents of the `simPar` file can be found by accessing ?simParexample 
`simPar` is a .csv file that contains the input parameters that characterize a specific simulation run, but which are *shared* among CUs. Each row represents a unique scenario (i.e. combination of operating model and management procedure). Generally it is easiest to create multiple `simPar` input files, each of which contain a coherent analysis (e.g. one input focusing on the effects of different harvest control rules across changing productivity regimes, a second input examining the effects of survey effort), but this is not strictly necessary. Contents include:
  
  - `scenario` - scenario name
  - `nameOM` - operating model name
  - `nameMP` - management procedure name
  - `keyVar` - focal variable of the analysis; subjective since typically multiple variables will differ among scenarios, but should be a focal point of main figures. Currently can be one of the following arguments: `prodRegime`, `synch`, `expRate`, `ppnMix`, `sigma`, `endYear`, `adjustAge`, `mixOUSig`, `adjustForecast`, `adjustEnRoute`, `obsSig`, `obsMixCatch` (**NOTE these should eventually be defined explicitly**)
  - `plotOrder` - order in which grouped scenarios will be plotted (useful when keyVar is not an ordinal or numeric variable)
  - `species` - lower case species name (chum and sockeye have been tested robustly; pink and coho have not; chinook should be used with extreme caution since most stocks do not meet assumptions of the model)
  - `simYears` - length of the simulation period (excluding priming period)
  - `harvContRule` - harvest control rule (`TAM`, `fixedER`, `genPA`)
  - `benchmark` - biological benchmark used to assess conservation status (`stockRecruit`, `percentile`)
  - `canER` - total Canadian exploitation rate
  - `usER` - American exploitation rate (note can also be supplied as CU-specific value in `cuPars`)
  - `propMixHigh` - proportion of Canadian catch allocated to mixed-stock fisheries (can range from 0 to 1)
  - `enRouteMortality` - on/off switch for en route mortality
  - `constrain` - if `TRUE` and harvest control rule is TAM then mixed stock fisheries are constrained
  - `singleHCR` - single stock harvest control rule (`FALSE`, `retro`, `forecast`)
  - `moveTAC` - if `TRUE` and single stock quota from low-abundance CUs is re-allocated to other CUs
  - `prodRegime` - productivity regime (`low`, `lowStudT`, `med`, `studT`, `skew`, `skewT`, `decline`, `divergent`, `oneUp`,  `oneDown`, `high`)
  - `startYear` - indicates when a productivity decline (if specified by `prodRegime == "decline"`) should start
  - `endYear` - indicates when a productivity decline (if specified by `prodRegime == "decline"`) should end
  - `rho` - temporal autocorrelation coefficient in recruitment deviations
  - `arSigTransform` - if `TRUE` estimates of sigma from input are transformed so that they account for temporal autocorrelation
  - `correlCU` - the correlation among CUs in recruitment deviations
  - `corrMat` - if `TRUE` a custom correlation matrix is required to be passed as an input and is used to specify the covariance matrix for recruitment deviations
  - `mu_logCovar1` - mean of lognormal distribution for annual SR covariate (e.g., marine survival)
  - `sig_logCovar1` - log-normal variation for annural SR covariate
  - `sampCU_coef1` - TRUE/FALSE indicating whether CU-specific coefficients for the SR covariate should be sampled
  - `sigCU_coef1` - among-CU variation in coefficients for the SR covariate (only applied if sampCU_coef1 == TRUE)
  - `tauCatch` - logistic variation in CU-specific catches
  - `obsSig` - log-normal variation in spawner observation error 
  - `mixOUSig` - beta-distributed variation in mixed-stock fishery outcome uncertainty; input parameter represents the standard deviation used to calculate the location parameter
  - `singOUSig` - beta-distributed variation in single-stock fishery outcome uncertainty; input parameter represents the standard deviation used to calculate the location parameter
  - `obsMixCatch` - log-normal variation in mixed-stock catch observation error
  - `obsSingCatch` - log-normal variation in single-stock catch observation error
  - `obsAgeErr` - logistic variation in observed age error
  - `lowCatchThresh` - lower aggregate catch target (used as a performance metric)
  - `highCatchThresh` - upper aggregate catch target (used as a performance metric)
  - `extinctThresh` - quasi-extinction threshold
  - `adjustSig` - scalar on CU-specific sigmas (recruitment deviations)
  - `adjustAge` - scalar on `tauCatch`
  - `adjustEnRouteSig` - scalar on en route mortality rates


#### `CUPars`
`CUPars` are .csv files that contain CU-specific input parameters. Note that these parameters should *not* vary among simulation runs. Differences in operating models that involve CU-specific traits (e.g. population dynamics) can typically be introduced via options in the `simPar` file. Each row represents a specific CU. 

Mandatory contents include:

  - `manUnit` - management unit
  - `stkName` - CU name
  - `stk` - CU identification number (can be assigned arbitrarily or based on previous modeling exercises)
  - `model` - stock-recruit model used to forward simulate dynamics (`ricker`, `larkin`)
  - `minER` - minimum Canadian exploitation rate
  - `alpha` - productivity parameter for Ricker models
  - `beta0` - density-dependence parameter for Ricker models
  - `sigma` - recruitment variation for Ricker models
  - `meanRec2` - mean proportion of age-2 recruits
  - `meanRec3` - mean proportion of age-3 recruits
  - `meanRec4` - mean proportion of age-4 recruits
  - `meanRec5` - mean proportion of age-5 recruits
  - `meanRec6` - mean proportion of age-6 recruits
  - `medianRec` - median historical recruitment
  - `lowQRec` - 25th percentile historical recruitment
  - `highQRec` - 75th percentile historical recruitment

Optional contents include:
  
  - Necessary if modeling cyclic stocks
    - `domCycle` - integer to identify the dominant cycle line in Larkin stocks (`1, 2, 3, 4` or `NA`)
    - `tauCycAge` - logistinc variation in age structure
    - `larkAlpha` - productivity parameter for Larkin models
    - `larkBeta0` - density-dependence parameter for Larkin models
    - `larkBeta1` - lag-1 density-dependence parameter for Larkin models
    - `larkBeta2` - lag-2 density-dependence parameter for Larkin models
    - `larkBeta3` - lag-3 density-dependence parameter for Larkin models
    - `larkSigma` - recruitment variation for Larkin models
  - Necessary if American exploitation differs among CUs
    - `usER` - American exploitation rate
  - Necessary if modeling en route mortality
    - `meanDBE` - mean difference between estimates (a proxy for en route mortality)
    - `sdDBE` - interannual standard deviation of difference between estimates
  - Necessary if modeling TAM harvest control rule
    - `medMA` - median mortality adjustment used in TAM harvest control rule
  - Necessary if modeling forecast process
    - `meanForecast` - mean forecast relative to observed
    - `sdForecast` - interannual standard deviation of forecast
