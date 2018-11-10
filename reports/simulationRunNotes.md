# Notes on Simulation Runs #
## Organized reverse-chronologically ##

#### Run single stock scenarios w/ secondary HCRs (Nov. 9)
-	Run retrospective w/ no secondary HCR, as well as forecasting and 	retrospective options
-	See singleStock/mdReports/secondaryHCRs.Rmd for details

#### Run single stock scenarios w/ productivity (Nov. 9)
-	Run retrospective secondary HCR w/ skewed normal or low alpha
-	See singleStock/mdReports/productivity.Rmd for details

#### Run single stock scenarios w/ uncertainty (Nov. 8)
-	Manipulate tauCatch, obsS, obsCatch, OU, and ageErr independently 	to evaluate relative impact on proportion allocation effects
-	See singleStock/mdReports/uncertainty.Rmd for details

#### Rerun full spectrum of single stock scenarios - (Oct 25)
-	Easiest to look at early catch and early spawner abundance
-	Reference case (no mixed stock constraints, no secondary harvest 	control rule, no moving TAC and no en route mortality) mix and 	single stock TAC will differ because mixTAC is split by true 	proportions, singleTAC by forecast proportions
	-	Furthermore, realized catches will differ due to uncorrelated 	outcome uncertainty
-	Aggregate
	-	Early basically as expected
		-	Constraint results in pos. correlation between mix and cons
		-	Single HCR results in opposite, pos. correlation between 	single and mix
		-	Single + constraint dampens differences between mixed and 	single resulting in increased S and decreased C overall
		-	Moving TAC further dampens differences, reducing aggregate 	spawner abundance considerably
		-	Adding ER mortality before harvest completely removes 	spawner abundances, but catches obviously still decline
	-	Differences between late and early results
		-	Relative benefits of mixed ppn (constraint only) weakens 	suggesting that catches benefit over long term
		-	Other trends typically weeaken suggesting that there is a 	high probability of less catch, but probability of greater 	spawner abundance is more dubious
		-	However impacts on status are generally stronger
-	CU-specific
	-	Constrain only: generally single bigger than mixed (as expected), 	but due to patterns described above possible for individual years 	to have mix catches bigger than single
	-	Single HCR only: generally mixed bigger than single as expected 	(note that over entire simulation period some weird patterns emerge 	due to abundance being driven down)
	-	Single HCR + constraint: primarily driven by single HCR so that 	catches maximized at high mixed proportion
	-	Single + cons + moveTAC: similar to above but magnifies differences 	among CUs (e.g. healthy CUs like Chilko differ from depleted)
	-	Single + cons + ER before fishery: stocks with high ER and 	exploitation rates benefit from mixed stock fisheries, others 	resemble above (e.g. EStud); others are just noisy (e.g. Chilko)


#### Compare MPs that do or do not reallocate TAC - (Oct 19)
-	MPs incorporate same rules as below, but TAC that is not taken is re-	allocated to CUs within the same MU that are above their upper biological 	benchmark
-	**Note in all of these simulation runs 100% of en route mortality occurs 	before single CU fishery and HCR is based on TRUE status**
-	Figures in `outputs/varyAllocation/mortalityOrder/`
-	CU-specific trends
	-	For several CUs moving TAC increases the contrast between different 	levels of allocation 
	-	I.e. single-stock fisheries have lower catches and spawner abundance when 	TAC is moved around
	-	Supplementary analyses suggest that this is driven by overfishing  	relatively early in the time series, resulting in permanent declines in 	recruit abundance; typically driven by scenario where high TAC coincides 	with poor return and high en route mortality
	-	I.e. TAC based on en route mortality is still risky
	-	Specific patterns vary among CUs based on how frequently they are the ones 	receiving other populations TAC


#### Compare MPs that vary allocation - Part 2 (Oct 16)
-	Includes more strict secondary HCR and toggle for timing of en route 	mortality, but single stock TAC is not re-assigned
-	Mechanisms influencing catch/conservation trade-offs vary across fisheries - 	constraint in mixed stock fishery, status/cycle line based closures in single
-	Overall catches in single stock fisheries remain lower than expected due to 	relatively frequent closures
-	Figures in `outputs/varyAllocation/mortalityOrder/`
-	CU-specific trends
	-	Mortality before fishery results in:
		-	a) flat trend (e.g. EStu) of stable spawner abundance across 	allocations, but decreasing catches; suggests minimal mixed stock 	harvest does not impact recovery (may produce modest improvements in 	BM based PMs)
		-	b) linear trend (e.g Bowron) suggests catches decline overall (due to 	frequent closures), but combination of mortality + TAC results in 	relatively 	fewer spawners reaching the grounds when they do occur
		-	c) clump (e.g. Chilko) suggests that for healthy stocks with low 	levels of en route mortality, it shouldn't matter
	-	Mortality after fishery results in:
		-	a) moderate increases in spawner abundance to dramatic declines in 	catch in stocks where closures are common (EStu)
		-	b) very tightly clumped for stocks where closures ~= constraints 	(Harrison, Chilko, Bowron)
		-	c) linear trend where closures are common, but don't result in 	improved performance (Late)
-	Aggregate trends
	-	Mortality before fishery results in reductions in median agg catch and 	small increase in spawner abundance, but decrease in ppn CUs above BM
		-	I.e. moderate gains in conservation performance are uneven and likely 	outweighed by reductions in catches
	-	Mortality after fishery results in linear trend with particularly bad 	performance in ppn above BM and ppn extinct
		-	I.e. overharvest still common even w/ rules?
	-	Notably flex assignment performs particularly poorly

#### Compare MPs that vary allocation (static) (Sep. 7)
-	Allocation to mixed stock fisheries varied from 0 to 1
	-	TAM MP w/ standard overlap constraints + secondary closures based on CU 	specific abundance; specifically no single CU TAC unless above true 	lower BM in 2 of previous 4 years
	-	Three operating models: reference, skewed productivity, and high 	synchrony
-	Outputs in `varyAllocation`
-	CU-specific trends were variable
	-	E.g. As expected catches and spawner abundance declined w/ mixed stock 	for EStu, but for Chilko spawner abundance positively correlated w/ mix 	and catch negatively
	-	Driven by complex interactions between number of CUs in an MU and a 	CU's status; e.g. can catch more fish w/ single stock fisheries in 	Chilko because it has good status and no constraints, EStu the opposite
	-	Generally weak effects on cons PMs except for Cultus, which benefited 	considerably from lower mixed stock
-	At aggregate level catch and spawner abundance, ppnCUs were positively 	correlated 	with mixed stock proportion 
-	**NOTE Sep. 10** Reran without CU-specific HCR which resulted in decreases in 	catch, but not in spawner abundance
	-	Likely due to the same ERs being applied before en route mortality in the 	case of mixed stock fisheries, and afterwards in single stock


#### Compare changes in productivity via direct effects on alpha (standard regimes) via deviations using skewed multivariate t (July 30; update Aug 16)
-	Scenarios listed in `data/opModelScenarios/fraserOMInputs_varySkew.csv`
-	Effects on PMs were as expected, i.e. lower productivity, long regimes w/ 	skewed t (30 yrs) or high frequency of drawing skewed t (0.4)
-	However time series of productivity have deviations that are too severe w/ 	skewed t relative to observed 
	-	Saved in `miscSensitivity/skewness`
-	Should either tweak student-t parameters or try skewed normal
	-	Weakened student-t a little but still strong; try skewed normal
- **Update Aug 16** ran models w/ skewed normal and low productivity
	- PM outputs are in `miscSensitivity/skewness`
	- Interactions complex but generally ref > skewN > skewT > lowProd
	- Skewed relationships also typically have stronger synchrony effects
	- Skewed N shows nearly identical patterns to normal just less steep
	- Also includes histograms comparing observed (grey), simulated w/ normal (red), and simulated w/ skewed (green)

#### Compare covariance and sigma effects on performance metrics - June 29 and July 3
-	Continuation of June 8 entry focusing on Fraser
-	Trial 1 - effects of uncertainty
	-	Ran models w/ 0 uncertainty (except for sigma which is necessary to create 	treatments) and had similar results
-	Trial 2 - effects of HCR, ER  and productivity
	-	Prod and exploitation independently had no strong effect, but now combined 	(0.65 exploitation and low productivity regime
	-	Effects largely the same across synchrony treatments, but PMs uniformly low 	so perhaps to severe
	-	Increasing ER alone resulted in generally similar trends, but as noted some 	PMs appeared dampened due to uniformly low response
-	Trial 3 - effects on observed rather than true PMs are identical
-	Trial 4 - increase synchrony and CVc simultaneously; add custom correlation 	matrix
	-	Custom correlation matrix, which mimics recently observed correlations among 	Fraser CU had similar effect to mod low OM (i.e. rho/sig = 0.3)
	-	Trends were stronger particularly at higher values rho/sig > 0.7
	-	**Note** sigma value used in synchrony treatments may be too low, which makes 	sense given the mean estimated across CUs is ~1; adjust accordingly by 	setting default sigma to 0.9 and allowing it to vary (i.e. for synchrony 	treatments)
-	Trial 5 - sigma adjusted upwards to realistic levels - ref value is 0.9 w/ range 	from 0.75 - 1.5
	-	Trends strengthen a little bit, but still seem low
-	Trial 6 - extend simulation length from 55 to 125 years and run w/ both synch and 	sigma increased
	-	Increasing simulation length only strongly influenced uncertainty - slightly 	tightened bounds
-	**NOTE CU-specific PMs are not uniform!**
	

#### Compare strength of autocorrelation - June 22
-	Mike Hawkshaw and Sean Anderson suggest that lack of autocorrelation may result in 	insufficiently strong recruitment deviations and reduced effects of synchrony
-	Additionally meta-analysis by J. Thorson (CJFAS 2014) indicates moderate 	autocorrelation in salmonids (~0.4)
-	Run simulations across suite of OMs w/ rho ranging from 0-0.9 and other 		parameters at moderate levels or reference levels; MP is fixedER at 0.4
	-	`data/opModelScenarios/X_varyAR.csv`
	-	`scripts/scenarioSpecific/runModelRobustOM_A.R`
-	Generated plots of time series of spawners, recruitsBY and logRS, as well as 	autocorrelation plots (`outputs/miscSensitivity/autocorrelation`)
-	Chum
	-	Overall having no autocorrelation seemed to strongly dampen variability, 	while high AR resulted in spikes that were far to large
	-	Time series of raw data (CU-specific) suggests medianAR (rho = 0.5) 
	-	ACF plots (summed values) suggest lowAR (rho = 0.1)
	-	Baseline (little uncertainty) vs. ref had negligible effect on the patterns
	-	Thorson's estimate may be most appropriate for this species
	-	**Update June 29**
		- 	Quick analysis of residuals from simple SR models indicates rho varies 	between 0.4 and 0.6 with median of 0.5; will use that as default rho
-	Sockeye	
	-	Overall less consistent results across different metrics 
	-	logRS: best supported by lowAR
	-	recBY: best supported by low/medAR
	-	S: best supported by low/medAR
	-	Difficulties in replicating cyclic trends produced in Larkin stocks
	-	**Update June 26**
		- 	Quick analysis of residuals from simple SR models indicates rho varies 	between ~0 and ~0.4 w/ mean of about 0.2; will use that as default rho
- Reruns conducted **June 25**
	-	Reran with subset of Ricker stocks only; low-med AR seemed to be most 	favored
	-	Reran Larkin stocks 
		-	W/out exploitation (i.e. natural dynamics only) both with and without 	AR; cycles still weak in both cases
		-	W/out en route mortality; cycles still weak
		-	W/ high exploitation (0.7); cycles still weak
	-	Reran Nass 
		-	W/ moderate and no AR; moderate increases in contrast among 	treatments when AR included
			-	Note that variability in spawners greater than variability		recruits (potentially due to more sources uncertainty)
		-	W/ high exploitation (0.7) and AR; still no strong link between 	synchrony and PMs



#### Compare covariance and sigma effects on synchrony and aggregate variability - June 8
-	Look at spawners, recruits (RY) and logRS
-	Parameters (`nassOMInputs_varycorr`) were 0, 0.1, 0.3, 0.7 and 1; 
-	Figures in `miscSensitivity/synchTests`
-	First round results w/ chum and minimal uncertainty from other variables (0.1 	across the board); saved as `chumBaseline`
-	Changes in metrics
	-	Synchrony
		-	Covariance: generally greatest variation in logRS (0.35-0.9 median) then 	rec and spawners were similar (0.5-0.75)
		-	Sigma: all treatments had similar effect resulting in lower bound of 	covariance treatment
	-	Wtd CV
		-	Covariance: all treatments had similar effects on spawners and rec 	(~0.5), as well as logRS (~1.1)
		-	Sigma: strong treatment effect on spawners and recruits (0.35 - 0.8), and 	even stronger effect on logRS (0.25-3.5)
	-	Agg CV
		-	Covariance: patterns reflected differences in synchrony (as expected), 	with logRS still showing largest difference between treatments
		-	Sigma: like covariance, patterns mirrored effects of wtd CV but muted 
-	Performance metrics
	-	Covariance: no effect on PMs
	-	Sigma: moderate effects on ppn of CUs above BMs (higher sigmas reduce 	proportion) and strong positive effects on median (spawner/catch) abundance
		-	Latter seems counterinuitive at first, but makes sense because CUs are 	no longer strongly regulated, which allows for extremely large positive 	recruitment deviations
-	Repeat above first w/ chum, but w/ reference parameter value set; then with 	sockeye with both limited (i.e. pars at 0.1) and then with reference set
	1.	Chum reference set: higher levels of stochasticity meant that there was an 	increase in the ppn of CUs above BMs when synchronous and with higher 	variance; other PMs remained the same
	2.	Sox base set
		-	Generally similar patterns to chum in TS, although mean synchrony shifted 	down when correlations low (likely due to more CUs), resulting in greater 	"spread" across treatments as well as greater influence on aggregate CV
		-	Much higher weighted CV in spawners/recruits, moderately greater in 	logRS; this was amplified when CU-specific sigma was varied, once again 	resulting in greater spread
		-	Counterintuitively, absolute value of weighted CV in logRS, but not 	spawner or recruitment **Potentially due to fewer CUs?**
		-	Compared to chum, high correlations result in moderate increase in 	variation in spawners
		-	High sigma results in similar trends in ppn CUs above BMs and in spawner 	abundance
	3.	Sox reference set: similar patterns in time series and pms from reference to 	baseline
-	Repeat sox analysis w/ low productivity (5th percentile of alpha distribution)
	-	Effect of synchrony still nil and effect of sigma unchanged


#### Revised sensitivity analysis using updated model - June 7
-	Ran simulations with Nass dataset varying obsSig, obsMixCatch, mixOUSig and 	singOUSig across 4 levels low (0.1), medLow (0.3), medHigh (0.7) and upper 	(1)
	-	All other par values held at reference levels
	-	Ran simulations w/ fixed ER = 0.3, then repeated w/ fixed ER = 0.6
-	Generated time series plots (`outputs/miscSensitivity/uncertaintyTS`) and 	spreadsheet of SD of variation in relevant error outputs, e.g. when varying 	obsSig, compared observed spawners to true (`outputs/miscSensitivity/	uncertaintyTS`)
-	Increasing ER introduced trends at higher outcome uncertainties, presumably 	due to lower ERs as abundance declined; **however SDs were similar**
-	Difficult to determine what is a realistic level of uncertainty however 	based on SD in observation error and outcome uncertainty 2000-2014 (see 	`ncCHumDoc.md` for details), suggest using following sigmas
	-	obsSig = 0.6: ~30% deviations in observed spawners
	-	obsCatchSig = 1: 75% deviations in observed catch
	-	obsMixOU/obsSingOU = 0.7: ~25% deviations in target exploitation rate


#### Compare age error - May 7
-	ageErr = 0.0, 0.1, 0.3, 0.6, 0.9
-	Basically no effect on management procedures because:
	-	a) it only influences obsRecBY
	-	b) the absolute abundance of recruits is unaltered, they're just misallocated 	to years. As a result, metrics that are integrated over "many" years (i.e. 	most PMs) are unaffected


#### Compare spawner obs error - Apr. 25
-	obsSig = 0.0, 0.1, 0.4, 0.6 (reference is 0.25)
-	Strongest impacts
	-	varObsSpawners - overestimate relative to true, most severe above 0.4
	-	ppnYrsEstUpper - fairly severe underestimates (i.e. less conservative) of 	status 
	-	ppnCUEstUpper - very severe underestimates of status
	-	ppnCUEstLower - fairly severe underestimates of status


#### Compare covariance - Apr. 19; May 7
-	correlCU = 0.0, 0.3, 0.6, 0.9
-	Note run with preliminary catch thresholds, if these become focus then 	simulations should be repeated
-	Strongest impacts
	-	sdS - moderate increases in variance
	-	ppnYrsUpp - moderate decreases in true status
	-	meanSynch - increases with covariance but surprisingly moderate impact (phi 	increases by 0.1 with strongest treatment) 
	-	endSynch - increases with covariance and stronger effect than on meanSynch 
-	Moderate increases in variability in catch and decreases in CU-specific S
-	**Update**
-	Observed increases smaller than expected, i.e. correlCU = 0.9 results in 	moderately stronger covariance in spawner abundance
-	Compare base run (plots in `varySimPar/prodCov`) with alternative OMs
	-	H1: due to exploitation introducing additional coherence; **no evidence** 	dropping ERs to 0.0 reduced synchrony, but not range of values 
		-	`prodCov_lowER`
	-	H2: due to stochasticity from varying productivity; 
		-	a) **no evidence** when CUs differed in prod, but prod stable among 
			trials (though variance decreased) 
			-	`prodCov_lowERmedianPars`
		-	b ) **no evidence** when CUs had same prod across trials and other 	parameters static
			-	`prodCov_lowERoneAlpha`
	-	H3: due to stochasticity from unique sigmas; **no evidence** although 	correlation matrix has much stronger covariance terms, doesn't translate to 	noticeable increase in synch
	-	H4: due to examining synchrony in spawner abundance, rather than realized 	productivity (true recBY/S)
		-	a) when productivity and sigma are shared among CUs (CU-specific 	recruitment deviations are still drawn each year), this results in large 	increase in synchrony range
		-	b) when median CU-specific productivity and sigma values are used, range 	narrows, but is still considerably larger than with just spawner 	abundance


#### Compare estimated and true parameters - Feb. 12
-	Independent of any observation error (except tauCatch) sGen and sMSY are 	consistently underestimated because ric A and ric B are overestimated, i.e. 	intrinsic productivity and DD effects estimated too high
-	Interestingly if there is any substantial error in assigning captured fish to 	CUs (e.g. tauCatch=0.5), it's effectively impossible to identify CU-specific 	ERs due to noise


#### Explore different observation error effects - Feb. 1
-	Input .csv's in data/observationSim and generated figures/data in outputs/	diagnostics/observationSim
-	Ran with subset (first 6) of CUs
-	Ran trials with low-moderate exploitation rates (0.05, 0.1, 0.1 respectively) 	and either low (0.05) or moderate (0.25) values for: `obsSig, 			mixOUSig/singOUSig, obsMixCatch/obsSingCatch, ageErr`
	-	`CorrelCU, tauCatch, arRicRho, obsRSig` set to 0; all bias parameters also 	set to 0
	-	CU-specific dynamics assumed based on in press Fraser CSAS report
-	Notable patterns based on differences between test (no error), low, and 	moderate values
	-	ageErr: minor deviations in estProd, obsRecBY, and sGen
	-	obsSig: 
		-	Minor deviations in obsS, obsRecBY, obsRecRY, benchmarks (increased 	contrast among CUs)
		-	Introduces considerable noise in obsER
		-	Little effect on aggregate indicators
	-	OU (both mixed and single manipulated simultaneously):
		-	Introduces variation in spawner abundance rather than between true and 	observed per se
		-	Introduces noise in true exploitation rate
		-	In aggregate, increasing OU moderately positively correlated with 	inter-trial variation in observed recruitment and aggregate catch, but 	most effects moderate
	-	Catch error (both mixed and single manipulated simultaneously):
		-	Deviations to obsRecBY and obsRecRY are proportional to magnitude of ER
		-	Introduces noise in obs exploitation rate
		-	*When all ERs are equal across CUs these are coherent deviations, but 	will become less synchronized as singleERs drift from mixed and from 	one another*
		-	Some differences in aggregate catch, but fairly minor


#### Explore different parameterizations for tauAge and tauCatch - Jan. 31
-	Based on `exploreTau.R`; input .csv's in data/tauDummySim directory and 	generated figures/data in outputs/diagnostics/tauExplore or outputs/data/	tauExplore
-	Runs included no source of error except tauAge and tauCatch and moderate ERs 	(0.1 for all fisheries)
-	tauAge simulations run with:
	-	tauEst - time series wide estimates of tau (calculated with 	`tauEstimator.R`)
	-	tauCycle - cycle line specific estimates of tau (calc as above)
	-	tauLow - uniformly low value (0.5)
	-	tauHigh - uniformly high value (2.0)
-	Compared generated time series for each to observed and estimated overall 	variance across TS for each age class; cycle line specific values were closer 	for most stocks, even those where Larkin model not top ranked
-	**Conclusion**: use tauCycle to generate age variation
-	tauCatch simulations bounded by 0.1 and 2.0 with annual percent error ranging 	from less than 1% in former to 10-15% in latter
-	**Conclusion**: more subjective since real world estimates don't seem to exist, 	but a value of 0.5-1.0 seems to produce reasonable assignment rates for well 	studied regions with many CUs (e.g. Fraser) or poorly studied regions with 	fewer (e.g. Nass chum)








