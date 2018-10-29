## Simplified Closed Loop Simulation Model##
**Diverged from primary model June 25, 2018**

###Major differences from recoverySim.R###
-	Simulation input parameters
	-	Pre-specified initial abundances (time steps defined as 	obsLag 	and set equal to generation length + 2) instead of 	time series of observed spawners
	-	Alpha drawn from a shared distribution
	-	Single beta/sigma values shared among subpopulations (Eric's 	request)
-	Removed
	-	Larkin model options
	-	Single CU fishery (now just American and overall Canadian)
	-	En route mortality
	-	Observations of component (i.e. subpop) catches which had 	considerable knock on effects (more detail below)
	-	Certain benchmarks (e.g. catch thresholds, synchrony) and 	objects that are used in higher level plotting functions  	(e.g. HCR names, exploitation rate calcs)
-	Catch observations are estimated directly, not with calcObsCatch, 	which partitions catches at MU level, then splits to specific CUs 	based on their proportional abundance
-	Changes to observation submodel
	-	In real world spawner abundances can be estimated at 	subpopulation level (e.g. indicator streams), but catches 	cannot unless baseline is comprehensive and funding extensive
	-	Therefore estimating recruitment at subpop level requires 	subdividing catch w/ arbitrary expansions (e.g. same 	proportions observed in streams); this could be added to the 	model, but it seems unlikely that management decisions will 	be made based on SR relationships at the subpop level, 	instead subpop-specific observations removed entirely
	-	Knock on effects mean that following are also not 	**observed** at subpopulation level (i.e. still exist in 	process model):
		-	American and Canadian catches
		-	Recruitment
		-	Estimated SR benchmarks
		-	Exploitation rate
-	Added subsampling of spawner abundance and spawner bias back into 	model (**NOTE appears to be functional but have not been 	rigorously tested**)
	-	Sample a set proportion of population
		-	E.g. if ppnSampled = 0.6, first 3 of 5 subpops always the 	ones observed
	-	Expansion factor calculated with `calcExpFactor()` after a 	prespecified number of observations of **full** aggregate 	(defaults to three equivalent to south coast chum model); more 	details around following code chunk
	
		    if (y == gen * 2) {
    			expFactor[n] <- calcExpFactor(obsS, nSamplePops, gen)
    		}
-	Replaced `plotOneTrial()` function used to generate diagnostic 	plots w/ simplified version
	-	Note that aggregate "true" SR curve is not plotted because 	these parameters are unknown (recruits generated at subpop 	level, observed at aggregate)


### Summary of inputs
-	simParFull: .csv of simulation wide (i.e. common to all subpops)
	-	Scenario: grouping variable used to define multiple OMs to be 	run simultaneously
	-	nameOM/nameMP: name of operating model and management 	procedure respectively
	-	species
	-	nPops: number of subpopulations to simulate
	-	simYears: number of years to forward simulate **AFTER** priming period, `nPrime`, which is a function of species-specific generation length
	-	initialAbundance: number of fish per subpopulation used to generate starting values
	-	usER/canER: American and Canadian exploitation rates
	-	model: stock-recruit model used to forward sim dynamics; **currently only Ricker is functional**
	-	alpha: mean of alpha distribution from which subpop specific 	values are drawn
	-	alphaSig: sd of alpha distribution from which subpop specific 	values are drawn
	-	sigma: sd of recruitment deviations; used to generate 	**diagonal** of multivariate normal distribution's error 	matrix
	-	correlPop: correlation among subpops in rec deviations; used 	to generate **off-diagonal** of MVN error matrix
	-	corrMat: if TRUE, then the simulation looks for a custom 	correlation matrix to pass instead of calculating one using 	correlPop and sigma
	-	tauCycAge: parameter used to generate variation in age-at-	maturity
	-	meanRecX: mean proportion of recruits return at a given age
	-	obsSig: sd of deviations between true and observed spawner 	abundances
	-	outcomeSig: sd of deviations between target and realized 	exploitation rates
	-	obsCatchSigma: sd of deviations between true and observed 	catches
	-	obsSBias: used to bias S observations up/down
	-	ppnPopSample: proportion of subpopulations sampled each year
	-	sampleProb: probability that sampling occurs in a given year
	-	ageErr: sd of deviation between true and observed age 	proportions
	-	extinctThreshold: quasi-extinction threshold; below this 	number recruitment and spawners go to zero within the model
-	EXAMPLE_prodCorrMatrix.csv
	-	This is an example of a custom input correlation matrix for an 	aggregate w/ 19 subpopulations (the Fraser)
	-	Off-diagonal element represent mean pairwise correlations 	among subpops (i.e. covariance)
	-	Diagonal represents variance elements; **these can be set to 	anything because they are replaced by sigma within the model**


### Summary of outputs
-	Single trial diagnostic figure (code in diagnosticFunc.R)
	-	Randomly pulls one trial and plots 
	-	Time series of spawner/recruits (one subpop)
	-	SR relationship for the aggregate w/ true and observed (plus 	estimated SR
	-	Time series of subpop specific variables, followed by 	aggregate variables
-	Aggregate diagnostic figure (code in diagnosticFunc.R)
	-	Shows similar time series of aggregate attributes, but with 	medians and 50% PI
	-	Finishes w/ spaghetti plot showing subset of trials colored 	based on final status relative to true benchmarks
-	List of subPop-specific summary statistics (`cuList`; saved as 	`paste (nameOM, nameMP, "cuDat.RData", sep = "_")`)
-	List of aggregate attributes over time (`agTSList`; saved as 	`paste(nameOM, nameMP, "aggTimeSeries.RData", sep = "_")`)
-	.csv of aggregate summary statistics (`aggDat`; saved as `paste	(nameOM, nameMP, "aggDat.csv", sep = "_")`)


