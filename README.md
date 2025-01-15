# samSim <img src="man/figures/logo.png" align="right" height="138" alt="" />
## Closed-loop simulation package for salmon conservation units and stock management units



### Summary
This repository contains the necessary files to run stochastic closed-loop simulations parameterized with Pacific salmon stock-recruit data. The principal function, `genericRecoverySim()`, is intended to test the performance of different management procedures (broadly a mix of harvest control rules and assessment methods) across operating models representing distinct ecological hypotheses. A suite of performance metrics are generated that allow analysts to evaluate different management procedures ability to achieve multiple, interacting conservation- and catch-based objectives. In short, the model is intended to provide a framework for the quantitative component of a management strategy evaluation.

The focal unit of the simulated dynamics are conservation units (CUs) - genetically distinct aggregates of one or more spawning populations that are unlikely to recolonize in the event of extirpation. Under Canada's Wild Salmon Policy these will be the target of future rebuilding strategies. The CUs are generally grouped into Stock management units (SMUs) , which are groups of one or more CUs that are managed together to achieve an joint objective. The SMUs are also equivalent to the "Major Stocks" in the Fish Stocks Provisions under the Fisheries Act. The `samSim` package includes feature to evaluare management startegies applied at both the CU and SMU scales. 


A summary of relevant files and how to run a simulation are provided below. Most functions contain relatively detailed documentation (and sometimes functioning examples). 

Install `samSim` with:

```{r}
devtools::install_github("Pacific-salmon-assess/samSim")
```

-----

------

### Running a simulation

Simulations are run by installing the samSim package and using the `genericRecoverySim()` function. Generally this should occur in a fresh working directory (e.g. a new `.Rproj`), which will be automatically generate an `outputs` directory and necessary subdirectories. Parameter values are passed to the function using a series of .csv files with the `simPar` and `cuPar` arguments being most critical. `simPar` contains parameter values that are shared among CUs and define a given scenario (e.g. species, simulation length, OM and MP characteristics). `cuPar` contains parameters that are CU-specific including at a minimum names and SR model type, but typically stock-recruit parameters as well. See *Input file details* below. Details of how to pass a suite of scenarios to the simulation model are provided in `?genericRecoverySim`. 

------

### Input file details

#### `simPar`
 A detailed description of the contents of the `simPar` file can be found by accessing ?simParexample 


#### `CUPars`
`CUPars` are .csv files that contain CU-specific input parameters. A detailed description of the contents of the `simPar` file can be found by accessing ?CUParexample
