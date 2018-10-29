# Summary of Input Parameters #

*Describes current primary input parameters for recoverySim.R and, where applicable, how they were derived*

### 1. simPar - defines simulation processes
-	**NOTE: justification for specific sim values (largely sigma values) are 	summarized in `reports/refParameterValues.xlsx` and 	`parameterJustification.docx`**
-	nameOM - operating model name
-	nameMP - management procedure name 
-	species - currently tested with `sockeye` and `chum`; `pink` and `coho` will be 	supported eventually
-	timeVaryA - set to `FALSE` (i.e. stable), `decrease`, `increase` or  `cycle`
	-	Currently specified at simulation run level, but could be swapped for CU-	specific values
-	simYears
-	harvContRule - `fixedER` or `TAM`
-	canER - Canadian mixed CU fishery (i.e. non-terminal marine)
-	corrrelCU - correlation among CUs in recruitment deviations
-	tauCatch - interannual variation in age structure; estimated multiple ways 	using proxy data (details in `fraserRiverDoc.md` and `exploreTau.R`)
	-	Not possible to parameterize with Nass data
-	obsSig - estimated spawner abundance error
-	mixOUSig - mixed (i.e. marine) fisheries outcome uncertainty (assumed to be the 	same for both US and Canadian mixed)
-	singOUSig - single fisheries outcome uncertainty
	-	Currently specified at simulation run level, but could be swapped for CU-	specific values
-	obsMixCatchSig - estimated catch error in mixed-CU fisheries
-	obsSingCatchSig - estimated catch error in single-CU fisheries
-	ageErr - error assigning recruits to a brood year
-	extinctThresh - pseudo-extinction threshold
-	minA - minimum productivity level at which declines stabilize

### 2. cuPar - defines CU-specific parameters
-	manUnit
-	stkName
-	stk
-	cuNameOM/cuNameMP - vector that is only used when multiple simulations are run 	with varying CU-specific traits (e.g. differences in single CU exploitation 	rate)
-	model - `ricker` or `larkin` (latter Fraser only)
-	usER - US exploitation rate; occurs before Canadian or single-CU exploitation
-	canER - Canadian mixed-stock exploitation rate
	-	**Note** currently overridden by simPar input, but could be adjusted
-	tauCycAge - error term used to generate interannual deviations from mean age
	-	For the Fraser these values are the mean of cycle line specific values 	estimated in `exploreTau.R`
	-	For the Nass these values could not be directly estimated (age structure is 	not regularly estimated), so **arbitrary** moderate value used
-	alpha - median Ricker's alpha
	-	For Fraser provided by AMH/FRSSI group (mean values from rickerLL - summary 	stats file)
	-	For Nass estimated externally using `ncStockRecModels.R` and then median 	calculated from 1000 samples
-	beta0 - median Ricker's beta; values as above for alpha with exception of Lower 	Nass CU (both median and sampled values are not estimated in SR model, but 	based on median spawner abundance (see script for details))
-	sigma - median Ricker's sigma; values estimated as above
-	meanRecX - mean number of recruits in each age class (estimated from time 	series for Fraser, static values used for Nass)
-	medianRec/QRec - median, lower, and upper quartiles for recruit abundance; 	added to some output plots
-	meanDBE/sdDBE - mean and interannual variation in "difference between 	estimates", a proxy for en route mortality 
	-	Fraser values estimated since 2000 from `postSeasonMortality_2017.csv`
	-	Nass fixed at NA
-	meanForecast - mean forecast error 
	-	Fraser estimated from data provided by PSC `ManagementTables2005-2017.xls`
		-	Given relatively small sample sizes, pooled most data, but kept summer 	estimate separate because persistently conservate, while other MUs' 	forecasts are overestimated
	-	Nass fixed at NA
-	sdForecast - variation in in-season forecast error; standard deviation among 	years	