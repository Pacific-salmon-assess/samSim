# North Coast Chum Case Study Documentation #

### Data available (data/northCoastDat)
Unless otherwise noted all North Coast chum data provided by Shaun Davies and Katie Beach (DFO)

#### area3OtolithResults.xlsx
1. *Worksheet 1: Otolith results*
	-	Provides marked vs. unmarked catch breakdown 2011-2016
	-	Data:
		-	Year, gear, area, week, estimates of catch (unsure what this value 	represents precisely), proportion of readable and marked otoliths, 
		-	Depending on year ~69-91% of catch
2.	*Worksheet 2: weekly var*
	-	Week to week variability in proportion of catch allocated to hatcheries
	-	Typically stable for 3-4 weeks before declining 
3.	*Worksheet 6: Hatchery ppn (3)*
	-	Age specific catches from 2011-2016
	-	Almost all age 4 except for 2016 when catches were dominated by age-5 	(totals were low but not abnormally so)
4.	Other worksheets have additional data that is likely not relevant – e.g. 	hatchery-specific estimates, data on individual fish
	-	**What is AlexData worksheet?**


#### northCoastChumAgeStruc.xls
-	Summary of results from otolith aging studies conducted on SA 2e, 2w, 3, 4, 5, 	6, 7, 8 and 10 chum
-	Data provided:
	-	Management subarea, sample location, wild/enhanced, year, frequency at age, 	proportion at age, brood year at age, sample size
-	For area 3 samples collected from two subareas and nine sample locations, 	collected from 2008-2012, with 2-151 fish per sampling location
-	Generally age 4s are at least 75% but some years/locations ~50%


#### northCoastCommercialCatch_2011to2015.xls
-	Provides daily, sub-area specific catch estimates for north coast statistical 	areas from 2005-2015
-	Includes seine, gillnet and troll fisheries
-	Data:
	-	Estimate type (e.g. in-season), area opening descriptor, target species, 	fishing date, area/subarea, hours open, number of vessels, species specific 	catches retained and released
-	**Haven’t examined closely yet**


#### area3ChumEnglishWB.xlsx
-	Provides estimates of stream-specific survey effort (# years surveyed within 	decade) and returns (not age-specific) from 1950-2015
-	Data:
	-	CU, area, reviewer, survey assessment type, return estimates per year, and 	decadal averages
-	61 streams have data across all 4 CUs
-	**NOTE LARGE GAPS IN COVERAGE AMONG CUs**
	-	Lower Nass: several years with no indicator streams surveyed prior to 	1992; 1993-2008 **one** stream surveyed in **one** year
	-	Portland Canal-Observatory: some coverage in all years
	-	Portland Inlet: some coverage in all years


#### summariesTRTC.xlsx
-	Appears to provide data used to generate K. English reports on north coast 	salmon populations (i.e. all species totals as well as specifics for 	statistical areas)
-	**Basically all data is a function of common ERs at the SA level, applied to 	CU-specific estimates of escapements generated from indicator streams w/ 	expansions**
-	**Used to generate catch data for our purposes** (no estimates of FN catch 	repavailable)
-	Worksheets include *Pivot* that combines species totals and subsets of SA 	data, *TotalPlots* which provides species specific totals, *StatAreaPlots*, 	and	*ConservationUnitPlots*; last two most relevant
-	Data for SAPlots and CUPlots worksheets runs from 1954-2014
	-	ID variables
	-	Q1-Q3 columns represent data quality
		-	Q1 is categorical index of survey quality (e.g. visual surveys vs. 	mark-recapture vs. counting fence)
		-	Q2 is survey execution reflecting EF1 (see below) and the proportion 	of indicators nsampled
		-	Q3 is survey representation reflecting EF2 and strength of correlation 	between index streams and CU as a whole
	-	Total index escapement + 2 expansion factors = observed escapement
		-	EF1 corrects for missing indicators, EF2 corrects for non-indicator 	streams
	-	Observed escapement + 1 expansion factor = total escapement
		-	EF3 corrects for observer efficiency (assumed to always be 1.5)
	-	Canadian harvest and exploitation rates + total escapement = TRTC (total 	return to Canada)
	-	Total harvest and total ER included US harvest/ERs (which aren't shown); 	total harvest + total escapement = total run
	-	**Note** all ERs are assumed to be static across CUs
-	Some escapement estimates available for all years in Portland Inlet and 	Portland Canal-Observatory, but not in Lower Nass

#### outputAge.xlsx
-	Provides age-specific recruitment estimates by brood year based on total 	returns to each statistical area and **mean** estimates of age composition
	-	In SA3 8.9% age-3, 76% age-4, 14.4% age-5 and 0.06% age-6
-	Not clear how these average values are generated, but presumably from recent 	otolith studies

#### chumSRDatPSF.csv
-	Brood table for north/central coast chum stocks provided by PSF (E. Hertz) 	used to generate SR parameter estimates
	-	Original data in `NBBR/Age Tables/CM OUTPUT AGE StatArea.xlsx`
-	Assumes age structure stable through time
-	Produced as output from NBBR model by K. English
	
#### psfMCMCPars.out
-	Stock recruitment parameter estimates for 3 Skeena and 2 Nass CUs provided by 	E. Hertz
-	Fit using jags model `MultiStockHierPSF.txt` which was replicated in our 	notation as `MultiStockHierPSFRep.txt`
-	Priors for these CUs were also provided by E. Hertz; modified in script to 	represent beta (log(1/prSmax)) and tau (1/prCV)
	
		CU				prSmax	prCV		
		Skeena_Estuary	1118	1		
		Lower_Skeena	20378	10		
		Middle_Skeena	2544	1		
		Portland_C_Obs	42230	10		
		Portland_Inlet	24524	10	
-	Analysis comparing estimates using their priors/model and our own is in 	`ncStockRecModels.R`; some details below in Apr. 30 update

#### northCoastDFOData_KBeach.xlsx
-	Catch, ER, and escapement estimates that appear to account for otolith (i.e. 	estimate of AK hatchery chum contribution)
-	Data from 2000-2017
-	Estimates differ substantially from summariesTRTC.xlsx


----

### Estimating stock-recruit parameters
**First pass on fitting SR models - Feb. 13** 

-	Received useable SR dataset 
-	Began working to estimate SR parameters using standard and hierarchcial Ricker 	models
	-	Based on JAGS scripts and associated R functions written by Brooke Davis for retro SC chum project
	-	Minor edits by CF with SR fitting script `ncStockRecModels.R` and `stockRecFunctions.R`
-	Modifications due to small sample sizes for Lower Nass CU
	-	Included 3 Skeena CUs in hierarchical model
	-	Smax priors are a function of maximum observed abundance over entire time 	series, not just years with full recruit data
		-	I.e. in `stockRecFunctions.R` escapement time series used to calc sMax, while paired escapement/recruitment used to estimate SR
	-	Removing Lower Nass and running model with 5 CUs instead of 6 does not strongly influence parameter estimates; however estimates of alpha for Lower Nass do not show as much shrinkage as expected
-	Stock Key:	
	`if(stk==1){ CU<-"Portland Inlet" };  if(stk==2){CU<-"Lower Nass" };
	if(stk==3){ CU<-"Portland Canal-Observatory" }`

**Validate age structure (RESOLVED) - Mar. 8** 

-	Double checked age-specific returns as follows; seem to be correct:
	
		for(i in 1:nrow(newChumDat)){
  		age3Esc <- newChumDat$Age3[i+3]*newChumDat$Escape[i+3]
  		age3Catch <- newChumDat$Age3[i+3]*(newChumDat$Total.ER[i+3]*(newChumDat$Escape[i+3]/(1-newChumDat$Total.ER[i+3])))
  		newChumDat$rec3[i] <- age3Catch+age3Esc
  		age4Esc <- newChumDat$Age4[i+4]*newChumDat$Escape[i+4]
  		age4Catch <- newChumDat$Age4[i+4]*(newChumDat$Total.ER[i+4]*(newChumDat$Escape[i+4]/(1-newChumDat$Total.ER[i+4])))
 		 newChumDat$rec4[i] <- age4Catch+age4Esc
  		age5Esc <- newChumDat$Age5[i+5]*newChumDat$Escape[i+5]
  		age5Catch <- newChumDat$Age5[i+5]*(newChumDat$Total.ER[i+5]*(newChumDat$Escape[i+5]/(1-newChumDat$Total.ER[i+5])))
  		newChumDat$rec5[i] <- age5Catch+age5Esc
  		age6Esc <- newChumDat$Age6[i+6]*newChumDat$Escape[i+6]
  		age6Catch <- newChumDat$Age6[i+6]*(newChumDat$Total.ER[i+6]*(newChumDat$Escape[i+6]/(1-newChumDat$Total.ER[i+6])))
  		newChumDat$rec6[i] <- age6Catch+age6Esc
		}

**Adjustments to stock recruit models - Mar. 16**

-	To make chum SR estimates consistent with Fraser River modeling strategy had 	to combine process and observation variance estimates into single value (i.e. 	sigma)
-	Commented out tauv and trimmed code moderately; saved as 	`MultiStockHierCF.txt`
-	These models fit more efficiently (neff much larger), gave reasonable 	precision estimates and behaved as predicted, i.e. data limited CU shrank more 	than others
-	Given small sample size of Lower Nass CU (i.e. stk 2), didn't use fitted 	values
	-	Instead that stock in `nassChumMCMCPars.csv` contains alpha values from the posterior hyperdistribution for alpha, beta0 value is equal to median from *escapement* time series and sigma is shared with Portland Canal-Observatory (larger sigma Portland Inlet)

**Discrepancies with PSF estimates - Apr. 30 (RESOLVED MAY 31)**

-	Attempted to get PSF and our own stock recruit estimates to converge
-	Ran four model combinations to try to identify the discrepancy
	1.	Normal model - our standard priors, all six CUs and 	`MultiStockHierSimple.txt` model
	2.	Trimmed model - our standard priors, five CUs and 	`MultiStockHierSimple.txt` model
	3.	PSF replicate model - PSF priors, five CUs and `MultiStockHierPSFRep.txt` 	model
	4.	PSF model - PSF priors, five CUs and `MultiStockHierPSF.txt` model 	(basically a copy of B. Connors model with very minor changes
-	We couldn't easily compare model 4 to our outputs because of differences in 	variable names, however 3 and 4 provided nearly identical summary statistics
-	Model 3 and 4 outputs were **NOT** consistent w/ parameter distributions (i.e. 	raw output data in `data/northCoastDat/psfMCMCPars.out`
-	They were close to model 2 and slightly less similar to model 1, suggesting 	our model structure is acceptable
-	For now using `MultiStockHierSimple.txt` model with our own priors
-	**Resolved May 31** divergence due to error by CF in calculating total number 	of recruits; etimates now similar

**Issues w/ data limited productivity - June 1**

-	When running low productivity scenarios based on 10th percentile of alpha, 	noticed that L Nass CU (data limited) had negative productivity values, i.e. 	well below replacement
-	Pattern due to using muAlpha as proxy for this CU's alpha; former has a highly 	uncertain distribution with a substantial amount of distribution below 0
-	Since there is no evidence to suggest that the Lower Nass has 	disproportionately lower productivity than other CUs suggest
	-	a) using informative priors to constrain hyper-distribution to more 	realistic values (e.g. change prior on muAlpha's sigma from dunif(0, 100) 	to dunif(0, 10) 
	-	b) sample other two CUs distributions with replacement to produce a mean 	area 3 value (**currently selected this option for both alpha and sigma**)
		
**Concerns around forward simulated variability - June 5**

-	During meeting with S. Davies, K. Beach, AM Huang and K Holt several commented 	that forward simulated variability appeared reduced relative to historic
-	Appears to largely be due to taking medians (individual trials still quite 	variable)
-	Although sigma values are low relative to Fraser, this seems unlikely to be 	driving the pattern as they are consistent with values in other chum SR studies 	(e.g. Collie et al. 2012 also varied between 0.3 and 0.7)
-	May also be influenced by:

	1. How data are generated, i.e. as model outputs (indeed individual samples of 	simulations show expected levels of variability)
	2. Less interannual variability in exploitation rates
	3. Unrealistically low levels of outcome uncertainty

**Catch comparison and outcome uncertainty - June 7**

-	Used two sources of catch data to help parameterize catch obs error (difference 	between method 1 and method 2) and outcome uncertainty (deviations between 	estimated ER and target ER)
-	Catch methods
	-	1) NBSSR model outputs (mod) - primary source of catch data for north coast 	chum
		-	Distinct American and Canadian harvest rates produced using expanded 	escapement and adjusted exploitation rates from other fisheries (either 	pink or sockeye depending on time period)
	-	2) Hail data (otoHigh and otoLow) - available for subset of years (2000 - 	2017)
		-	Consists of commercial catch data + est. release mortality + Nisgaa 	marine catch 
		-	Available at area level and majority assumed to be Alaskan hatchery 	origin
		-	To calculate potential SA3 origin harvest multiply by the proportion of 	thermally marked fish observed in catch as well as a dummy variable 	representing the ppn that may be migrating into SA3 (rather than 4 or 	6); assume high scenario (0.9; OtoHigh) and low (0.1; OtoLow)
-	otoHigh and mod catches track each other fairly closely, but begin to diverge 	~2010; mean relative difference is 1.76 +/- 3.14 (i.e. otolith catches are 2x 	higher on average)
	-	otoLow exploitation is much lower than either
	-	**Suggest using 1.5 as sigma for catch**
-	Outcome uncertainty is scaled relative to target ER of 10%
	-	modER typically ~10% lower, otoHighER ~25% high; both sigmas ~0.75
	-	**Suggest using 0.75 as sigma for mixedOU**
-	Summary figures are in `miscSensistivity/catchComparison/catchTrends.pdf` and 	analysis in `misc/ncCatchComparison.R`
-	Notes on simulations conducted examining different levels of obsSig and OUSig 	are in `simulationRunNotes.md`

**Update stock recruitment models - July 23**
-	Given evidence of autocorrelation in residuals, decided to re-estimate SR 	parameters using an AR Ricker model (MultiStockHierAR.txt)
-	Prior for rho ~ dnorm(0, 1)
	-	Resulted in realistic estimates for all parameters
	-	Increasing sd in prior to 10 resulted AR effects completely disappearing 
	-	Limited examples to compare to, but generally analysts seem to use 	uninformative prior w/ broad uniform distribution
-	At same time changed prior on hyperdistribution for alpha to be set on 	TauA_global ~ dgamma(1, 1)
	-	I.e. sigA_global no longer defined
-	Data limited CU's parameters estimated by drawing 1000 samples from dataset 	that contains the other 5; beta based on 1/median observed escapement