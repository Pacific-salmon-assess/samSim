# Future modules to add to closed-loop simulation model

### 1) Increase ecological complexity
-	Add indicator streams that are surveyed then expanded to provide CU-specific 	estimates of abundance; goal is to better represent fine-scale changes in 	survey effort
-	Moved to a high priority due to relevance to north coast chum systems

### 2) Temporal changes in age structure
-	Alter age at maturity (presumably towards older individuals) through time; 	would be modeled in ways equivalent to productivity changes

### 3) Mechanistic drivers of variation in productivity
-	Instead of directly introducing changes in productivity, link changes in alpha 	to environmental drivers (e.g. SST) that can then be manipulated based on IPCC 	forecasts
-	Excluded for now due to likely complexity

### 4) Changes to carrying capacity and depensation
-	Incorporate different SR relationships/operating models to reflect changes in 	habitat capacity (i.e. beta term changes through time) or depensation
-	Relatively straightforward but postponed to minimize complexity of outputs, as 	well as because current focus on ERs as management levers and depensation can 	be partially accounted for with quasi-extinction thresholds

### 5) Additional sources of uncertainty
-	Currently need to limit sources of variance to make working with model 	tractable, but assumptions should be noted and potentially addressed at later 	date
-	Likely sources
	-	1) CU-specific deviations in recruitment forecasts: `recErr` term could be 	drawn for each CU in every year instead of once per year
	-	2) CU-specific deviations in outcome uncertainty; `singOUSig` moved to 	parsCU file
	-	3) CU-specific deviations in catch/spawner observation error; `obsCatchSig`

### 6) Changes to internal management model
-	Change internal function used to estimate stock status from simple `lm` to 	hierarchical (i.e. mixed effects model) accounting for correlations among 	stocks
-	For sockeye case study estimate BMs and SR parameters using Larkin model - will 	likely have to be Bayesian

### 7) Additional adjustments to HCRs
-	Add fine-scale adjustments to HCRs based on factors external to species model 	(e.g. reduce exploitation of Fraser sockeye when tail will overlap with coho)

### 8) Account for unsampled CUs
-	Given evidence that MUs such as Nass frequently have gappy TS would be worth 	adding missed surveys (note that this is currently cooked into chum CSAS model, 	but has not been adjusted for new model's structure and will be incompatible 	with many current outputs)
