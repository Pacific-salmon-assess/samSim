# Fraser River Case Study Documentation #

### Data available (data/fraserDat)
Unless otherwise noted all Fraser River sockeye data provided by Ann-Marie Huang (DFO); data quality summaries available online (http://open.canada.ca/data/en/dataset/3d659575-4125-44b4-8d8f-c050d6624758)

#### `sockDat.xls`
-	Contains stock-specific estimates of age-specific recruitment, effective 	female spawners (EFF), total spawners (ETS), and estimates of juvenile 	abundance
-	For now used to estimate CU-specific age structure and provide historic  	spawner abundance in priming period; SR parameters taken from `rickerLL - 	mcmc_pars.csv` files
-	Currently only complete (i.e. total recruits for all age classes) through 2011
-	**NOTE** whether stock forward simulated w/ Ricker or Larkin based on 2017 WSP 	Reassessment document

Stock key

	if(stk==1){ stk.name<-"E. Stuart" };  if(stk==2){stk.name<-"L. Stuart" };
	if(stk==3){ stk.name<-"Stellako" }; if(stk==4){ stk.name<-"Bowron" };  
	if(stk==5){ stk.name<-"Raft" }; if(stk==6){ stk.name<-"Quesnel" }; 
	if(stk==7){ stk.name<-"Chilko" };  if(stk==8){ stk.name<-"Seymour" }; 
	if(stk==9){ stk.name<-"L. Shuswap"}; if(stk==10){ stk.name<-"Birkenhead" }; 
	if(stk==11){ stk.name<-"Cultus" }; if(stk==12){ stk.name<-"Portage" };  
	if(stk==13){ stk.name<-"Weaver" }; if(stk==14){ stk.name<-"Fennell" }; 
	if(stk==15){ stk.name<-"Scotch" };  if(stk==16){ stk.name<-"Gates" };  
	if(stk==17){ stk.name<-"Nadina" };  if(stk==18){ stk.name<-"Upper Pitt" };  
	if(stk==19){ stk.name<-"Harrison" }

Model key

    if (stk == 1 | stk == 2 | stk == 6 | stk == 8 | stk == 9) {
	  fraserCUPars$model[i] <- "larkin"
	}	else {
	  fraserCUPars$model[i] <- "ricker"
	}

#### `rickerLL - summary stats.csv` (or larkin equivalent)
-	Contains median + quartile values for stock-recruitment parameters (alpha, beta 	and sigma; plus b1, b2 and b3 for Larkin models)
-	Use medians for initial model runs
-	Based on relationship between recruits and total spawners (ETS), not effective females (EFF); PSM accounted for in estimates of ETS --> EFF --> rec

#### `rickerLL - mcmc_pars.csv` (or larkin equivalent)
-	1000 SR parameter sets for each stocks
-	From AMH: "For something that is “close enough”, I’d say take the median alpha 	from the summary stats and find a par set with that alpha in it and use that. 	For something closer, I’ve included some code that will sort through and winnow 	down the 1000 par sets to ~10-20 sets based on how close each par set is to the	median alpha, b0, etc."
-	Based on relationship between recruits and total spawners (ETS), not effective 	females (EFF); PSM accounted for in estimates of ETS --> EFF --> rec 
-	**Model samples directly from these files to select a coherent set of larkin or 	ricker pars each sim run**

#### `postSeasonMortality_2017.csv`
-	Includes time series of en route (pDBE or proportional difference between 	Mission and spawning grounds abundance estimates) and prespawn mortality that 	are used to generate meanDBE and sdDBE (representing mean and variation in en 	route mortality since 2000) vectors in `fraserCUPars.csv` input
-	pMA values used to inform management adjustments for forecasts
-	Note 1/0 column that is used to exclude wonky years

#### `tamRefPts.csv`
-	Fishery reference points courtesy of AMH that are current values for 2010-2014, 	used to generate HCR within model
-	Units are 1000s of fish
-	medDBE value represents median difference between estimate post-2000 (calculated 	in `prepFRParInputs.R` and described in detail below); used to adjust forecasts 	and account for expected en route mortality

#### `camGenDat\frTotalCatch.csv`
-	CU-specific estimates of marine, first nations, and total catch (as well as 	exploitation rates)
-	Used to estimate catchTau, calculate forecast error rates, and also passed to 	simulation model for inclusion in diagnostic plots
-	**Note**: marine catches are technically a mix of US and Canadian commercial 	fisheries, therefore their inclusion in figures as "Can Catch" is not strictly 	correct

#### `camGenDat\americanERs.csv`
-	MU-specific estimates of American catches and total run size from Fraser River 	panel reports from 1996 to 2016
-	MU-specific estimates of ERs were <0.0005 for early Stuart and ~0.06 for the 	other three MUs
-	These are lower than management targets (16.5% for all MUs but E Stu) because 	management targets are applied to fish after accounting for en route mortality 	losses; the model now reflects this and 16.5% ERs are default 

#### `camGenDat/forecastEstimates.csv`
-	MU-specific estimates of pre-season, in-season (3 days after 50% migration 	date) and post-season (meeting after TAC date) of run size from 2006-2017
-	Error calculated as estimate/total run, with total run size coming from 	`frTotalCatch.csv`
-	Use mean and sd of in-season values to parameterize forecast error in closed-	loop model; added to fraserCUPars.csv
-	Raw data provided as .xls file `ManagementTables2005-2017` by PSC; summaries 	available in FRP annual reports

#### `camGenDat/cultusSpawnerEsc.csv`
-	Estimates of adult returns to Cultus since EFF and ETS data are missing from 	sockDat file due to captive breeding program
-	Data were copied from *Recovery Objectives* tab of `Cultus Matrix 2017 _psm 	assumptions_ (22-Feb-2017).xlsx` (not saved on github)
-	NOTE that SR parameters for this CU were not estimated post-2000 due to captive 	breeding program; as a result ETS is likely uncoupled from recruit abundance

#### Used to generate fraserCUPars.csv, which contains CU-specific SR parameters, age estimates, initial abundances and exploitation rates, and is one of two input files for simulation