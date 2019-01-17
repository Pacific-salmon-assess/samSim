# Summary of modifications to C. Holt's chum model

### I) Addition of random number vectors within salmon-sim code (Dec. 11); **replaced Jan. 17**
-	Originally used complex system of random number calls in vectors to ensure 	appropriate number called regardless of extinction events
-	Removed `if(extinctAg == 1)` statement and instead modified model so that 	entire cycle run regardless of spawner abundance


### II) Modifications to make model usable across species (Dec. 13/14)
#### Input parameters
-   Include species argument within function
	-   Generation length (gen) and age class distribution used to represent 	differences
	-   Could potentially also specify stray rate and species-specific extinction       	thresholds
	-   Adding coho should be straight forward, CK more difficult
-	Note that nInit and total simulation length will change as a result of 	adjusting gen
-	Treat odd and even year cycle lines as independent

#### Options for making code flexible
1.	nAges argument represents species and interacts w/ external functions
2.	nAges assumed to be static and relative ppn of each group changes to reflect  	species diversity

**Option 2 selected for now**

#### Code adjustments
-	Add in species call separating pink from chum/sockeye/coho
-	Allow for age-2 returns; necessary for pink obviously but may also be useful 	for stocks such as Harrison sox
-	Estimate age class proportion, recBY, recRY accordingly


### III) CU-specific estimates within model (Jan. 3/4)
-	Calculate age structure, benchmarks, catches, forecasts and exploitation rates 	at level of individual CUs
-	Will require modifications for inputs because some (e.g. `ageStruc`) are now 	matrices
-	**Update Jan 10**: edited var-covar matrix to allow for CU-specific sigma 	estimates
-	**Update Jan 24**: added CU-specific estimates of `tauAge`

#### Code adjustments
-	Largely consisted of changing outputs from vectors/matrices to arrays
-	E.g. `ppnAges, Smsy, estRicB, ppnAge2Ret4` now have additional dimension
	-	All data stored as years X CUs X trials


### IV) Change model outputs (Jan. 8/9; edited: Jan. 17, Jan. 24, Mar. 23)
-	Diagnostic figures currently include
	-	1) 
		-	a. Randomly draw one CU from one trial and plot: SR relationships and 	spawner abundance relative to median and quantile values
		-	b. For all CUs plot (true and observed): alpha, S, obsS, recBY, recRY, 	Canadian catch, US catch, single catch, exploitation rate, sMSY, sGen, 	(saved as `plotOneTrial` function)
	-	2) Plot median and 20-80% CI limites for aggregate: spawners, observed 	spawners, recruits (BY/RY), obs recruits (BY/RY), catch, obs catch, 	synchrony, obs syncrhony, benchmarks, estimated benchmarks (saved as 	`plotMedian` function) 
-	Data exported as one .csv of aggregate values (nTrials X nVariables) and one 	list of CU-specific matrices (nTrials X nCUs X nVariables)
	-	1) `aggDat.csv` includes: 
		-	General info: model run, harvest control rule, trial number and 	targetER
		-	Aggregate spawner abundance (mean and SD across years per trial)
		-	Aggregate spawner abundance in first and last two generations in TS
		-	Aggregate catch (mean and SD across yrs per trial)
		-	Observed aggreagte catch (mean and SD)
		-	Mean aggregate catch in first and last two gens
		-	Mean realized and observed realized exploitation rate
		-	Proportion of years aggregate is above upper/lower BMs
		-	Ppn of years aggregate catch is above aggregate threshold
		-	Ppn of CUs above upper/lower BMs for last generation of TS
		-	Ppn of CUs extinct at end of sim
	-	2) `cuDat.Rdata` includes CU-specific values for:
		-	General info: harvest control rule, CU name, MU name, and targetER
		-	Spawner abundance (mean and SD)	
		-	Alpha (mean and SD)
		-	Beta or B0 for Larkin (mean and SD)
		-	Total true catch (mean and SD)
		-	Total observed catch (mean and SD)
		-	Realized exploitation rate and observed ER
		-	Ppn years above upper/lower BMs
		-	Ppn years **observed** abundance is above **estimated** BMs
		-	Ppn years with COSEWIC or WSP decline
		-	Ppn years above Canadian and single-CU catch thresholds


### V) Model inputs (Jan. 15; **updated Feb. 12**, **updated Feb. 16**, **updated March 6**)
*Required files*; **Additional details in inputParSummary.md**

1.	`...simInputs.csv` contains simulation specific parameters that are shared 	among CUs
	-	Can be modified manually by entering different values for specific 	simulations
1. `...CUPars.csv` contains CU-specific SR parameters and age estimates, as 	well as entries for different initial abundances and exploitation rates
	-	Age proportions are estimated using all available data; ageTau excluded 	age classes that are never observed (i.e. age 2 for non-Harrison CUs)
 
*Optional files*

1. To account for uncertainty in underlying SR relationship, parameters can now 	be input as medians (in cuPars.csv) or sampled from two additional input 	csvs - containing either 1000 samples per stock for Ricker or Larkin pars 		(e.g. `fraserDat/rickerMCMCPars.csv`)
 	-	To account for correlations among parameters, only one sample is drawn per 	stock, which contains each par necessary
1. `...RecDatTrim.csv` contains stock and age-specific recruitment data used to 	set initial age proportions and prime model with abundance TS
1. `...CatchDatTrim.csv` contains sector specific catch data 
1. **Fraser only**  Added file containing cycle line specific fishery reference 	points drawn from FRSSI (`tamRefPts.csv`)
	

### VI) Divided catches across multiple fisheries (Jan. 19 and 23; updated Feb. 7)
-	Assumed to be captured in either US, mixed Canadian or single-CU Canadian
-	True catches estimated using three different ERs
	-	`usER` and `canER` are in `sims.csv`; vector `singER` is in 	`...CUPars.csv`
	-	Also incorporate two levels of outcome uncertainty: lower for US/mixed Can, 	higher for single CU (reflects differences between commercial and FN/rec)
-	Observed catches for US and mixed Canadian fisheries
	-	Incorporate two levels of catch observation errors: lower for US/mixed Can 	fisheries, higher for single CU fisheries (same rationale as above)
	-	Incorporate uncertainty when splitting observed aggregate catches (e.g. 	`obsCanCatch`) into CU-specific catch estimates; `singObsCatch` does not 	because these fisheries would assume all captured individuals originate in 	that CU
-	Proportional contribution of each CU is estimated using `ppnCatchErr` 	function, which is closely related to `ppnAgeErr` function, and based on each year's true 	proportional contribution plus error associated with mis-assignment of 	aggregate catches to specific CUs (see catchTau below for details)
-	Added an en-route migration mortality term `migMort` that represents fish that 	escape the marine American/Canadian fisheries, but die prior to spawning and 	terminal fisheries (could be moved if this jars with allocation issues)
	-	Mortality rate `enRouteMR` is CU-specific and varies across years based on 	`enRouteSig`
	-	In the case of Fraser River mean M is the mean proportional difference 	between estimates (pDBE) since 2000 and sigma is the standard deviation in 	pDBE over the same period
		-	Full dataset in `data/fraserDat/postSeasonMortality_2017.csv`
	-	En route mortality is observed with an error equivalent to observing 	spawner abundance `obsSig`
	

### VII) Adjusted MSY calculation (Feb. 2)
-	Following Scheuerell 2016 PeerJ, replaced linear approximation for sMSY 	suggested by Hilborn (1985) `(ricA/ricB)*(0.5-0.07ricA)` with explicit solution 	generated with Lambert W function: `(1 - lambert_W0(exp(1-ricA)))/ricB`


### VIII) ageTau parameter selection (Jan. 19 and Feb. 5)
-	Evaluated performance of four ageTau structures using `exploreTau.R` and output 	files in `simData/exploreTau/exploreAgeTau`; input files in `data/tauDummySim`
	-	1) High tau value (2.0) that was uniform across all CUs
	-	2) Low uniform tau value (0.5)
	-	3) Tau estimated across full TS using `tauEstimator.R`; mean=1.17
	-	4) Tau estimated for each cycle line; these values were 0.1-0.5 lower with 	mean=0.98
	-	Full and cycle line specific estimates for Fraser River CUs saved in `data/	fraserTauEstimates.csv`
-	Compared observed to generated age proportion TS `ageProportion.pdf` and 	compared SD of each age class's proportional abundance
	-	Tau based on cycle lines are somewhat higher than observed, but closer than 	low, high or full TS tau
	-	**Will use cycle line specific estimates of tauAge for now**
-	Examined correlations between cycle line specific abundances, age structure 	(i.e. proportion of age-4 individuals) and SD in age structure
	-	Tau values only strongly correlated with interannual variability in age 	structure; suggests that there is substantial variation among cycle lines 	in proportion of age-4 individuals (seems consistent with DD hypotheses)
-	**Update May 8**: average, year-specific deviations in true vs. observed 	recruit abundance; however total observed abundance (i.e. sum of recruits 	observed across years doesn't change)
	-	0.1 = 1.02 (+/- 0.3)
	-	0.3 = 1.04 (+/- 0.4)
 	-	0.6 = 1.09 (+/- 0.5)
	-	0.9 = 1.16 (+/- 0.7)
	-	1.5 = 1.39 (+/- 1.5)


### IX) catchTau parameter selection (Jan. 19; Feb. 5; Feb. 19)
-	catchTau estimated using CU-specific return data and script 	`tauEstimator.R`, which is based on `get.mv.logistic.tau.R` from C Holt
-	Since misassignment is influenced by many factors (quality of genetic baseline, 	sample size in GSI runs, representativeness of test fisheries), accurate 	estimates of this uncertainty do not exist. As a proxy used two alternative 	methods
	-	1a) Estimate tauCatch based on interannual variation in observed CU-specific 	proportions within each cycle year; assumption is that observed deviations 	between years could represent misassignment rates if management is assuming 	set proportion will return
		-	Values between 0.3 and 0.6 (mean=0.5) depending on cycle year
	-	1b) As above but instead of using spawner abundance used total catch
		-	Values between 0.3 and 0.7 (mean=0.475)	
	-	2) Estimate tauCatch based on intra-annual variation in observed CU-	specific proportions across test-fishery surveys; assume that variation at 	smaller temporal scales reflects variation at larger temporal scales
		-	Limited dataset (2010-2014) and for now only included data from 	Johnstone Strait TF, but estimated values are similar (0.3-0.7, 	mean=0.42)
-	Additionally conducted simulation study with range of values (0.1, 0.5, 1.0 and 	2.0) for tauCatch, with all other sources of variability set to 0
	-	Example output figures in `diagnostics/tauExplore` show true/observed CU 	specific contributions to aggregate catch and percent error
	-	Values from 0.5-1 produced spikes in error of ~10% and seemed reasonable
-	**Will use value of 0.5 for tauCatch for now**	


### X) Replace first priming loop with true recruit/spawner abundance (Feb. 6,7)
-	Instead of using constant spawner abundances in years 1 through 5 allow SR 	data to be used as input; used to generate initial values for S, recBY, and 	ppnAge
-	Intended to reflect actual status of populations in recovery context and also 	necessary for cyclic stocks
-	Input data should be .csv with following columns: stk, yr, rec2, rec3, rec4, 	rec5, and ets
-	**EDIT**
	-	To allow BMs to be calculated shifted to using obs spawner abundances for 	entire time series available (SR datasets are currently trimmed and 	cleaned within model)
	-	Also made the assumption that observed and true abundances/benchmarks are 	the same during initial priming period. This allows BMs to be integrated 	seamlessly between priming and simulation component 
	-	Saved original script as `SimpleSimWCatch.R` and renamed `recoverySim.R`
-	**EDIT 2**
	-	Added marine and First Nations catch data to represent Canadian catches 	and single CU catches prior to sim period
-	**Note**
	-	May need additional modifications to account for differences in available 	input data
	-	Add argument so that if a CU specific benchmark cannot be estimated in a 	given year it takes the last value non-NA value (e.g. estSGen[nPrime,k])


### XI) Add `mortCalc` function (Feb. 8)
-	Wrapper function that sequentially calculates usCatch, canCatch, en route 	mortality, and singleCatch providing output as a list


### XII) Add preliminary HCRs (Feb. 8, 13, 15;)
-	Add simple HCR as alternatives to fixed ER
-	Intended to reflect FRSSI TAM rules 
	-	When forecasted recruitment below lower BM, canadian ER and singER each = 	0.05 for all MUs except late which gets 0.1 for each
	-	When between upper and lower BMs, Canadian and single CU ERs are calculated 	to allow for escapement equal to the lower BM after accounting for **MU-	specific average** en route mortality and US exploitation
	-	When above upper BM Canadian ER is 0.3 and single CU ER is 0.3 (normally 0.6 	for TAM)
-	Currently using cycle line and MU-specific fishery reference points provided by 	AMH as BMs (`tamRefPts.csv`)
-	**NOTE LARGELY DEPRECATED VIEW MAY 29 NOTE**


### XIII) Account for gappy input data (Mar. 9)
-	If stock-recruit TS is incomplete in recent years initial error term for Ricker 	AR operating model cannot be estimated; similarly recBY are not available to 	prime recRY
-	Modify B. Davis infill function to estimate recruits and spawners in last 5 	years where data are missing using proportional contribution (geometric mean) 	from observations over past 25 years (or length of TS, whichever longer)
	-	Added as an intermediate for loop between observed historical data and 	closed-loop forward sim


### XIV) Incorporate realistic estimates of forecast error (Mar. 16, Mar. 31, Apr. 27)
-	Use pre-season, in-season (3 days after 50% migration date) and post-season 	(meeting after TAC date) of run size from 2006-2017 for each MU as potential 	estimates of forecast error in the model 
	-	See `camGenDat/forecastEstimates.csv` and `prepFRParInputs.R` for details
-	Post-season run size estimates come from `fraserDat/camGenDat/frTotalCatch.csv`
-	Currently forecast error is normally distributed with mean equal to mean in-	season estimate relative to post-season run size estimate (e.g. mu = 1.1 if 	mean in-season estimate 10% larger than observed) and sigma equal to sd among 	years
-   **Update**: replaced MU-specific distributions with:
	-   mu = 1.20 for all MUs except summer (~mean value among these MUs)
	-   mu = 0.85 for summer because underestimates seem common for this MU)
	-   sd = 0.15 slightly larger than observed for 3 grouped MUs
	-   Modified `mortCalcTAM` to constrain target Canadian FR to 0.01 when 	realized ER is over 100%
	-   Tried to apply TAC by applying ER to forecasted recruits rather than true 	(see `mortCalcTam` in `simFunctions.R`); however high level uncertainty, 	biased towards overfishing seemed to result in unrealistically large levels 	of extirpation
		-   Moderating OU seems like a more tractable, equivalent process
-	**NOTE**: currently there is no equivalent process for Nass chum model because 	a) we lack historical estimates of forecast error (likely because no forecast 
-	takes place) and b) TAM rules are unlikely to be applied to these CUs, 	therefore not necessary


### XV) Generate multi-simulation run outputs (Mar. 21, 22; Apr 23)
-	`pmFunctions.R` contain functions to produce output plots using summary lists 	of CU-specific and aggregate performance metrics
-	Figures
	1.	Trade-off dot plots that contrast catch PMs w/ various conservation PMs, but 	can easily incorporate greater range
		-	Option to include inter-trial variability as confidences intervals is 	functional, but produces very "busy" figures for cu-specific plots and 	should 	likely be modified
	2.	Continuous trade off plot (e.g. Walters Skeena) that show escapement and 	catch on the left axis, extirpation risk and benchmark status on the right
	3.	Dot plots outlining aggregate and CU-specific multi-trial status across 	proportional performance metrics (e.g. ppn of CUs above upper BM at end of 	trial)
	4.	Multi-operating model tradeoff plots are produced separately 	(`compareOMsPlot`) and plot agg PMs vs. MPs for multiple OMs simultaneously


### XVI) Switch from internal MVT AR to simple MVT Ricker model (Mar. 30, edited Apr. 17 w/ new function)
-	Using `Ricker.MVT` from `simFunctions.R`

		Ricker.MVT<-function(S,a,b,cv,error){	
 			N<-length(a)
  			R<-NA
  			if(a>=0){
    			if(b!=0) R <- S*exp(a-b*S)*exp(err)
    			if(b==0) R <- S*exp(error)
  			}
  			if(a<0) R<-S*exp(a)*exp(error)*exp(err)
  			return(R)	
		}

-	Justification: input SR parameters currently being estimated with 	hierarchical, non-AR models
-	Currently structure still contains AR-model specific objects (e.g. `err` and 	`arRicSig`), but could be removed depending on future model goals


### XVII) Pooled observation error (Mar. 31, Apr. 10)
-	Add `obsErrDat` dataframe which contains observation error estimates
-	Observations of mixed CU catch and forecast errors are MU-specific, while 	spawners and single CU catches are CU-specific (drawn individually but from a 	shared distribution)
	-	In case of Fraser en route mortality shares the *same* error as spawner 	observations
	-	Recruit observations do not have their own observation error because they 	are simply sum of observed spawners and catch
-	Re-estimated each year and not stored


### XVIII) Ran profiler to ID bottlenecks (Apr. 10)
-	Used package `profvis`
-	Little room for obvious improvement; largest processing times associated with 	1) constructing `ppnAges` matrix, 2) estimating y-intercept, slope and rSq with 	internal SR relationship, and 3) estimating slope of decline


### IXX) Streamline outputs generated by simRun.R (Apr. 14, Apr. 17)
-	Add parallel processing
	-	Allows different scenarios to be run in parallel using `parallel`, 	`doParallel`, and `foreach` packages
	-	Note that bugs are most likely to appear when passing user defined 	functions that are nested within `parlapply` call; be sure to specify 	**all** necessary libraries, objects, and external functions using 	`clusterExport` or 	`cluserEvalQ`
-	Change output functions to access combination of directories (based on	 	scenarios run simultaneously) and subdirectories (based on specific OMs) where 	diagnostics, data, and output figures are now stored


### XX) Changed external hierarchical model used to estimate SR pars for chum (Apr. 17)
-	Instead of being passed observed recruits, passed log(R/S)
-	Intended to make it consistent w/ B. Connor's PSF model, but parameter 	estimates 	still don't converge
-	Compared to previous model structure, parameter estimates are basically 	identical except TauA is an order of magnitude larger


### XXI) Separated outcome uncertainty (Apr. 18)
-	Previously OU in fisheries was uniform across CUs within a year, i.e. if CUs 	had same target ERs in a fishery they would have same realized ERs
-	Now draw `mixOutErr`, `migMortErr` and `singOutErr` from `qnorm(runif(nCU, 	0.0001, 0.9999), 0, fisheryOUSigma)`


### XXII) Add weighted CV and agg CV estimates (Apr. 20)
-	Following Thibaut and Connolly 2013, sqrt(synch) * wtdCV provide an index of 	aggregate variability
-	Added 10 year moving window estimates to explore how different OMs/MPs 	contribute to component variability/synchrony


### XXIII) Identify optimal number of trials (Apr. 24; Apr. 29)
-	Sampled 25-1500 trials worth of PM outputs for Fraser River frSoxVaryProd 	simulation runs
-	Median values in aggregate and CU-specific PMs appear to stabilize at 250-500 	trials, but variance can fluctuate; suggest running 500 for first checks, 	1000-1500 final
-	Outputs in `outputs/miscSensitivity/minTrials`


### XXIV) Misc. productivity adjustments (May 7 - 8; May 15; May 21)
-	Add varyA == FALSE term to function inputs
	-	Used to constrain Ricker/Larkin alpha term to be the same across trials
	-	Currently set to mean value of all CUs in the simulation
-	Add realProd matrix to function outputs
	-	Contains realized productivity estimates (i.e. recBY/S) - used to provide a 	more meaningful representation of synchrony (covary more strongly with 	correlCU; see simulationRunNotes.R for additional details)
	-	Added to aggregate output plots as well
-	Change productivity operating models and do not sample from posterior but 	calculate point values
	-	High productivity: 90th percentile of alpha values + associated beta and 	sigma
	-	Moderate productivity: 50th percentile
	-	Low productivity: 10th percentile 
	-	Decline from medium to low w/ start and end date 
-	**Edit May 21**
	-	Adjust above so that only alphas are sampled at different prod regimes, 	betas/sigmas remain at median	


### XXV) Add allocation variables (May 15; changed July 25)
-	Added three input variables 
	-	`ppnMixHigh` - proportion of Canadian catch allocated to mixed stock 	fisheries by default (currently set at 0.8)
	-	`ppnMixLow` - as above but lower value intended to test whether outcomes 	improve when singleCU allocations increase; only used when following true
	-	`varyAllocation` - if TRUE then ppnMixHigh is replaced with ppnMixLow when 	threshold values reached
-	Currently when `varyAllocation` is true the HCR control rule switches to 	ppnMixLow if agS < agSGen in 3 of previous 5 years
	-	Eventually this should at the very least be replaced by observed and 	estimated
	-	With TAM rule could be directly incorporated so that allocation shifts 	based on ref pts, but currently just passed same ppnMix value as fixedER 	would be
-	**NEW RULE** Replaced previous version since agSGen is too easily met; now 	ppnLow is used whenever the average CU within an MU is not above its lower BM in 	at least 50% of previous generation
		  
    	meanStatus <- ifelse(length(CUs) > 1, 
      						 mean(apply(lowerBM[(y - 1):(y - gen), CUs], 2, sum)),
       						 sum(lowerBM[(y - 1):(y - gen), CUs]))
		poorStatusMU[y, CUs] <- ifelse(meanStatus > 0.5 * gen, 1, 0)


### XXVI) Divide catch observations by MU (May 18, 21)
-	Switch observations of US and Canadian mixed stock fisheries to occur at MU 	level so that mis-allocations are realistic
-	Prevents MUs that may contribute to large portion of run, but rarely intercepted 	(e.g. Early Stuart) from making up a disproportionate portion of observed catch 	due to catchTau
-	Calculated with `calcObsCatch`
-	As a result also downweight catchTau to 0.1 from 0.5


### XXVII) Cap recruitment (May 25)
-	Limit recruitment to be no more than 3x maximum observed historic recruitment
-	If SR time series not available cap at 5x 75% quantile (passed in CU pars file)


### XXVIII) Adjust TAM rule and exploitation rate calculations (May 29; edit Aug 16)
-	Forecasts
	-	To ensure consistency between HCRs, forecasts are now made in all MPs
	-	In data limited scenarios, median and SD estimates for forecast accuracy are 	set to 1 and 0 respectively, i.e. forecasts perfectly represent reality and 	outcome uncertainty should be adjusted upwards to compensate
-	Instead of adjusting forecasts downward based on median en route mortality, use 	actual management adjustment which is calibrated from same data (pDBE dataset)
	-	pMA * esc goal = numerical management adjustment to add to esc.goal
	-	pDBE * (run size – total catch) = number of fish that will be counted at 	Mission but not counted on spawning grounds (i.e., projected DBE)
	-	pMA = 1 / (1 + pDBE) - 1
	-	Currently using median pMA with interannual variance since 2000, pulled 	from `postSeasonMortality_2017.csv`
-	**NOTE** In fixed exploitation rate HCRs, forecasted recruitment cannot be 	moderated via pMA because this requires an escapement goal which is increased by 	the pMA	(adjusted escapement goal * harvest rate = TAC under TAM); this may lead 	to issues where TAM rule is more conservative than fixed independent of changes 	in harvest rate with abundance, but I don't know how to work around it
-	Pseudo-code
	1.	Forecast of recruitment is made (at MU level because this is scale at which 	PSC data available)
	2.	If the species is not sockeye, amTAC and canTAC are calculated simultaneously 	from forecasted recruitment and prespecified exploitation rates
	3.	If the species is sockeye and HCR = fixedER, a totalTAC is calculated based 	on Canadian ER and Americans get 16.5% of this total, Canadians get remainder
	4.	If the species is sockeye and HCR = TAM, totalTAC depends on forecasted 	abundance
		- If forecast below lower reference point (RP), totalTAC = minimumER * 	forecast	
		- If forecast above LRP and below URP, escapement target (by default LRP) is 	adjusted upwards using MU-specific median management adjustment (MA); if 	forecasted recruitment is above the adjusted target, an ER is calculated 	based on the difference between forecast and target; this ER is only used 	if it is a) larger than minimum ER and b) forecast is greater than 	adjusted target
		- If forecast above URP, escapement target adjusted upwards using 				(1 - maxER)*foreRec; ER selected as above 
	5.	Canadian TAC (totalTAC - amTAC; still at MU level) is distributed 	proportionally between mixed and single CU fisheries creating a mixTAC and 	singleTAC
	6.	mixTAC and singleTACs are used to generate target (mixHR and singleHR) and 	realized harvest rates (mixHR + outcome uncertainty and singleHR + outcome 	uncertainty) by dividing by true recruitment at the MU level. 
		-	This step is necessary because outcome uncertainty is parameterized 	for exploitation rates, but not catches. 
		-	MixHR and singleHR are at the MU level (because TAC is), but because 	outcome uncertainty varies among CUs, the realized HRs vary among CUs
		-	TACs are divided by true recruitment, rather than forecasted 	recruitment, because outcome uncertainty incorporates deviations between 	target and realized harvest rates
	7.	Finally harvest rates are applied to each CU to produce catches.
-	Relevant functions are `calcTAC` and `calcHarvRate`, replacing `mortCalc` and 	`tamERCalc` in `simUtilityFunc.R`
-	**NOTE AUG 16** Change forecast component so that forecasting error differs 	among MUs


### XXIX) Add custom correlation matrices (June 5)
-	Add argument to simPars input so that if `corrMat == "custom"` model looks for 	an input correlation matrix to use instead of applying static value passed from 	`simPars$correlCU`
-	For Fraser River this correlation matrix generated from stock-specific Ricker/	Larkin models (chosen as above based on PSC report and fit using EFF not ETS) 	fit to all available data; however correlation matrix looks only at residuals 	since 1990
	-	Details in `fraserCovariance.R`


### XXX) Parameterize OU and obsCatch for Nass (June 7)
-	Use deviations between two catch estimates (NBSSR model vs. hail + otolith data; 	2000-2014) to get estimate of SD in observation error 
-	Use deviations between ERs of two methods and target ER (10%) to estimate outcome 	uncertainty
-	Conducted simulation tests w/ different levels of obsSig and ouSig to identify 	values that approximated historical patterns
-	Details in `ncChumDoc.md` and `simulationRunNotes.md`, as well as referenced R 	scripts and output files


### XXXI) Add single stock ER HCR (June 11; edited July 30)
**OBSOLETE SEE NEW NOTE FOR SEP 25**
-	Single CU ERs coded to remove fixed proportion of Canadian TAC (currently fixed 	at 20%)
-	Added simple management lever so that this TAC is only exploited when 	recruitment has greater than lower benchmark in at least 50% of years in 	previous generation
-	Currently uses true recruitment and true sGen, but would be more realistic to 	switch to observed and estimated respectively

        singTAC[y, k] <- ifelse(sum(lowerBM[(y-gen):(y-1), k]) > (0.5 * gen), tacs[[3]], 0)
-	Note that second `ifelse()` used to define singCatch is separate and helps 	constrain realized ERs < 1
	

### XXXII) Identified bottlenecks part II (June 20)
-	After discussion with S Anderson decide additional performance required
	1.	Remove aggregate diagnostic plot - very inefficient and rarely useful since 	most information in single trial plot
	2.	Remove moving window estimates of synchrony/CV due to high overhead with 	rollapplyr; instead export a list of 2 arrays (nYears x nCUs x nTrials) for 	recBY and spawners that can be used to generate these posthoc
	3.	Replace standard `lm()` with `quickLM()` which is a wrapper for a C++ 	equivalent that is trimmed down and much faster
-	Overall runtime decreases by ~40%


### XXXIII) Re-introduce temporal autocorrelation (June 21)
-	Concerned that forward simulated populations don't have sufficiently strong 	recruitment deviations so reintroduce autocorrelation
-	Use `modRickerAR1.MVT`, which is modified version of C. Holts function but 	modified to account for CU-specific recruitment deviations being estimated 	outside model
	-	When rho = 0, results are equivalent to `Ricker.MVT`
-	Details of simulation tests w/ different ranges of rho are in 	simulationRunNotes.md
-	Reducing sigma helps, but makes dynamics overly deterministic


### XXXIV) Add proportion of fisheries open PM (July 19)
-	Add a new PM to track how frequently the fishery is closed
-	Calculated for each MU and follows different rules depending on species 		and statistical area
	-	In Fraser, an MU's fishery is open if the true recruit abundance 		`recRY` in a given year is above the cycle-specific fishery reference 		point ID'd by TAM rule
	-	in Nass, the MU's fishery is open if the **previous** years Canadian 		catch (i.e. mix catch + single catch) did not exceed 0.1 (target 		maximum exploitation rate for SA)
-	Final PM is the proportion of fisheries that are open in a given year; this 	is also used to calculate the proportion of years that **all** fisheries 		are open in `aggDat.csv`
-	**EDIT Aug 14** Additional fishery PMs and adjustments 
	-	AMH	suggested adding Fraser specific fishery PMs representing TAC values 	corresponding to stressed (500k TAC across all MUs) and healthy fisheries 	(1 million across all MUs)
		-	Added as inputs (`lowCatchThresh` and `highCatchThresh` to simInput 	files)
		-	Note that these values incorporate forecast error when TAM HCR is		being used	
	-	However given current framework these thresholds are met in vast majority 	of years and are unlikely to be sensitive to different MPs/OMs
	-	As alternative modified ppnFisheryOpen stat so that it no longer flips on/	off, but instead a mean proportion of MUs that are open each year is 	estimated for full sampling period


### XXXV) Add skewness operating model option (July 24)
-	**NOTE UPDATED OCTOBER SEE SECTION 41**
-	Some evidence that declines in sockeye productivity might be triggered by 	synchronous negative recruitment deviations rather than declines in mean 	productivity
	-	Plots of histograms indicate some CUs have marginally non-normal residuals
-	To simulate specify a time period (equivalent to method used for `decline` 	scenario); during these years recruitment deviations are drawn from a skewed 	multivariate t-distribution
-	Currently generated using `sn` package with `rmst(n = 1, xi = rep(0, nCU), 	alpha = rep(-0.5, nCU), nu = 7, Omega = corMat)`
-	Generally results in higher median performance but also greater variance among	trials


### XXXVI) Adjust BMs for cyclic benchmarks (Aug 17)
-	Change cyclic BMs so that status is only assessed on dominant cycle lines 	otherwise status assumed to be the same as that of previous year
-	Requires estimating status in priming loop as well


### XXXVII) Add harvest constraints to pseudo-TAM rule (Aug 30)
-	Because there is considerable overlap in MU migration phenology (and hence 	exploitation rates), FRSSI will constrain harvest by reducing TAC 
-	In reality, this is explicitly modeled based on MU-specific overlaps and how 	close each MU is to its upper fishery reference point 
-	We simplified the process, reducing TAC by 25% unless one of the following 
	-	1) All adjacent MUs are above their upper FRP (after incorporating the management 	adjustment for en route mortality) 
	-	2) A MU's TAC is calculated based on minimum exploitation rate 
-	See `constrain` and `calcTAC` functions in simUtilityFuc.R for details	


### XXXVIII) Added catch PMs (Aug 30)
-	`catchChange` calculates the relative difference in total catch between 	years, then inverts the value to provide an index of catch stability (i.e. 	lower is better)
-	`lowCatchThresh` and `highCatchThresh` are Fraser specific PMs that represent 	minimum TACs corresponding to higher levels of management strain (500k and 1 	million fish)


### XXXIX) Adjust calcTAC and calcHarvRate functions (Sep 25)
Goal was to ensure that TACs that are calculated based on TAM rule and 	relative allocation to mixed vs. single-CU fisheries result in harvest rates 	that **do not** decline as fishery progresses

1.	Calculate TAC based on forecasted abundance at MU level relative to 	fishery reference points (cycle/MU-specific), after accounting for 	management 	adjustment for en route mortality `calcTAC` function
2.	Multiply `calcTAC` output by true proportional contribution of each CU 	(for 	mixed stock fisheries) or forecasted proportion (for single CU fisheries)
	-	Accounts for fact that this is a managed process only in the latter 	case
3.	Convert adjusted TACs into harvest rates so that outcome uncertainty can be applied, then back convert to TACs; otherwise en route mortality results in single CU fisheries underharvesting
	-	To ensure that harvest rates result in catches equivalent to TACs, the 	denominator of the HR calculation (i.e. recruit abundance) needs to 	decrease as we move from mixed to single stock fisheries
	-	See `calcHarvRate` for details
4.	Subtract realized TAC from recruits with the constraint that the remaining 

Potential issues
-	Because outcome uncertainty applied at multiple steps (both mixed stock 	fisheries, en route mortality) it is very possible that realized catch rates 	will not reflect proportions that are pre-specified
-	Management adjustment in TAM rule may not be sufficient to prevent overharvest 	in single stock CUs


### XXXIX) Remove bias corrections (Sep 27)
-	Observed strong effect of bias corrections on recruitment deviations (i.e. -	sigma^2/2 to center mean on true value)
	-	Detailed comparison of OMs with and without bias correction are in 	`outputs/miscSensitivity/biasCorrection`
-	Ultimately removed from both recruitment deviations and other log-normal based 	distributions (e.g. observation error) because the original parameter values 	were not estimated with bias estimators (and medians not means are being used 	to seed anyways) 
-	C. Holt's opinion on issue from email on Sep25
    -	Sean, yes, the Ricker parameters are estimated with normal errors using the linearized form of he Ricker. However, recruitment is forward simulated with log-normal errors. That step of the forward simulation is similar to the one used in Anderson et al. 2015 (portfolio paper). Likewise, we initially included a –sig2/2 correction when simulating recruitment because of the log-normal errors, to ensure that the arithmetic mean of the error term =1 when exponentiated .  However, we noticed that very few salmon simulation papers actually do this correction. Papers estimating Ricker parameters and using them subsequently to predict/forecast recruitment generally include a +sig2/2 correction to the recruitment predictions or the Ricker a value itself to provide the arithmetic mean prediction (expected value) instead of the median prediction (e.g., Dorner et al. 2009, Haeseker et al. 2005, Parken et al. 2006). My understanding was that if  variable x is log-normally distributed (e.g., recruitment), then in the distribution of log(x) is normal with mean=mu and SD=sig  (like the linearized Ricker model). In this case, then exp(mu) = median (x)  AND the expected value or arithmetic mean of (x) = exp(mu+sig2/2).  Others have justified omitting the bias corrections because they have simulated, estimated, and reported medians exclusively and not worried about providing expected values. I guess I’m trying to assess if this is consistent with FRSSI and SBC Chinook MSE thinking, and/or any other simulation/estimation projects any of you have been involved with.  `


### XXXX) Add correlated en route mortality (Oct 3)
-	Attempt to increase synchrony in recruit abundance by add correlations in en 	route mortality deviations (multivariate normal)
-	Use both a uniform strong correlation (0.75) and a custom covariance matrix 	intended to represent pairwise relationships among CUs
-	Although both options resulted in correlated migration mortality *rates* 	neither produced strong correlations in en route losses, spawner abundances, or 	recruit abundances
	-	E.g. mean pairwise correlation in mortality rates increased from 0.18 to 	0.21
-	Unsure why, although likely due to dilution effects moving from deviations to 	mortality rates to en route losses
-	Accept minimal correlations among recruits for now


### XXXXI) Update skewed and student-t operating model (Oct 3)
-	Use skewed normal and skewed student-t versions of Ricker/Larkin models to 	generate parameters for alternative productivity regime operatin models
-	Parameters estimated in synchSalmon repo using skewedModels.R script; stan 	models written by S. Anderson
-	Median skewness parameter estimate = 0.85; to represent a more severe scenario 	use 75th percentile = ~0.65
-	Weak evidence for heavy tails except for Weaver however given relatively short 	time series not unexpected; use nu = 2.5 or 3 to represent heavy tail scenario 	(i.e. catastrophic OM)


### XXXXII) Update single stock HCR (Oct 15)
-	Addition 1)
	-	Secondary harvest control rule - single stock TAC is only removed if a new 	CU-specific OCP is passed
	-	For Ricker models this is median spawner abundance in previous generation 	must be greater than lower BM
	-	For Larkin models single stock harvest only occurs on dominant cycle lines 	and abundance of previous dominant cycle line has to exceed upper BM
-	Addition 2)
	-	Secondary HCR can also be used to re-assign TAC from low status to high 	status stocks
	-	CU X, Y and Z are in the same MU; CU X is below its lower RP, CU Y is 	above its upper RP and CU Z is between the two
	-	CU Z receives its own single stock TAC
	-	The single stock TAC from CU X is shifted to CU Y, which also receives its 	own TAC; since CU Y may not be able to sustain the full TAC from X, it 	receives a proportion based on its abundance in the run
		-	For simplicity's sake, if CU Y makes up 25% of the run, it will receive 	25% of any other CU's TAC
-	Addition 3)
	-	Add toggle so that single stock fisheries can occur before or after en route 	mortality
	-	Allows sensitivity analyses to identify relative efficacy of fisheries in 	different locations without explicitly endorsing a specific allocation 	scheme


### XXXXIII) Add AFE to TAM rule (Oct 16)
-	Add Aboriginal Fisheries Exclusion to Fraser TAM rule
-	This guarantees 400k fish to Canadian TAC prior to US TAC being calculated
	-	I.e. AM TAC is calculated with US ER assuming that aggregate abundance is 	400k smaller and that those 400k fish are distributed proportionally across 	MUs


### XXXXIV) Add forecast based secondary HCR (Nov 9)
-	As alternative to retrospective HCR (i.e. based on median obs S over prev 	generation), add one based on forecasted spawner abundance
-	Calculated using CU-specific forecasts of recruitment, adjusted downward by 	mean mortality rate if ER occurs before the fishery does, minus amTAC minus 	mixTAC
-	Reference points are still the same, i.e. TAC is only taken if forecast is 	above lower BM


### XXXXV) Switched outcome uncertainty to additive (Jan 4)
- Following C. Holt's suggestion OU is now additive (i.e. 1 + rnorm) rather than multiplicative
- Also considered switching forecast error to additive but it appears to be commonly applied multiplicatively (e.g. Catalano and Jones 2014)


### XXXXVI) Switched outcome uncertainty to beta distribution (Jan 16)
-	Parameterized harvest rate outcome uncertainty for Fraser River sockeye using deviations between mid-season TAC and run size estimates vs. post-season catch and run size estimates
	-	Could not parameterize TAC because target TAC is often 0 creating a bimodal distribution that is complicated to replicate
	-	Details in `salmon-sim/estimateOU.R`
-	As a result modified `samSim/recoverySim.R` to apply OU on ER rather than TAC by adding a `calcRealCatch` function that receives target TAC as inputs
	-	Also switched from additive normal to beta distribution following S. Anderson's advice
	-	Parameterized w/ data in `estimateOU.R` script
	-	Requires back-calculation; see `synchSalmon/ouSimulationTests.Rmd` for details


### XXXXVII) Remove forecast uncertainty from TAM rule calculations (Jan 16)
-	K Holt noted that using forecasted recruitment to set TAC and then applying outcome uncertainty from observed data would result in inflated error rates
-	Note that forecasts of recruitment are still generated for use in single stock harvest control rules