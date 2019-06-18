**Intro**
- Fisheries management is typically complicated by pervasive uncertainty and competing objectives. 
	- Paritcularly true of Pacific salmon stocks 
		- Extreme variability in abundance
		- Complicated, highly migratory life history 
		- Evidence of regime-like shifts in productivity
		- High genetic diversity, but predominantly mixed-stock fisheries
		- Diverse array of stakeholders
		- Mix of healthy and highly depleted populations
- Closed-loop simulations can be used to identify management strategies that are robust to stochastic recruitment variation, future changes in productivity regimes, and outcome/implementation uncertainty.
	- Closed-loop simulations is also well suited to evaluating how well a given strategy balances conservation and socio-economic objectives. 
	- When the quantitative simulation process is paired with repeated dialogue between stakeholders, managers and analysts --> MSE
- Closed-loop simulations characterized by two main components
	1) Operating model - represents a hypothesis about the ecological dynamics of the harvested stocks, as well as how the fishery interacts with them
	2) Management procedure - represents the assessment process, as well as the harvest control rule that sets exploitation rates
	- In each year a new generation is generated based on spawner abundance using a stock-recruit relationship (e.g. Ricker); this abundance is observed with error within the model, providing an estimate of abundance, which is fed into the HCR; when individuals recruit into the fishery some portion is harvested based on the HCR, the remainder are allowed to spawn
	- Process is iterated through time to span the management period of interest and across numerous Monte Carlo trials to capture the uncertainty across multiple parameters
	- Running simulations across different scenarios (i.e. a unique combination of one operating model and one management procedure)
- Each scenario assessed against a suite of performance metrics that represent different conservation- or catch-based objectives
	- Can vary in ecological scale and time frame; e.g. median catch of one stock over 5 years, proportion of stocks above their biological benchmark in 20 years
- Closed-loop simulation models can be used to evaluate a wide range of management questions
	- At the simplest scale can be used to evaluate how abundance/recovery varies across different exploitation rates
	- Alternatively can be used to evaluate more detailed HCRs (single-stock fisheries), the impacts of regime shifts, or changes in assessment strategy

**Generic model details**
- The `samSim` package is used to conduct closed-loop simulations focused on aggregates of Pacific salmon CUs
	- Basic requirements are estimates of CU-specific stock-recruit parameters and age at maturity, to drive the operating model, and a fully terminal fishery, because the model cannot account for harvest of immature fish
	- Model initialized with observed time series of spawner abundance so that status at the beginning of the simulation reflects most recent assessment
- Model outline
	- Recruitment simulated with Ricker (or Larkin) stock recruit models, which incorporate interannual variability in rec devs
		- May also covary among populations and introduce temporal autocorr
	- American fisheries (if relevant) are applied based on mean exploitation rates 
	- Three different harvest control rules can be used to regulate Canadian catch: fixed exploitation rates, an abundance based exploitation rate, or a TAM rule (a Fraser sockeye salmon-specific abundance-based rule)
	- Canadian fisheries are applied, following the harvest control rule, with catches deviating from target due to outcome uncertainty
	- Remaining recruits allowed to spawn, producing subsequent generation
- Key model outputs
	- Diagnostic plots
	- Files for both raw data and summarized performance metrics
	- Performance metric figures: continuous line plots, dot plots, and tradeoff plots

**Additional components**
- The basic model is intended to provide an indication of how changes in exploitation rate will influence the long term status of a stock aggregate.
	- However many additional components of the simulation model can be adjusted to explore how these patterns change under different ecological assumptions or, alternatively, how more nuanced harvest control rules perform
- Changes in productivity (i.e. the average number of recruits per spawner) can be modified by:
	1) By passing arguments in the simPar input file; options include a scalar multiplier, decline, divergent trends, or one CU going up/down
	2) Manipulating alpha directly in external stock recruit analyses (e.g. kalman filter)
- Different probability distributions impacting recruitment deviations (skewed, student-t) and different levels of covariance
- Single-stock fisheries based on proportional allocation
- En route mortality including mean, interannual variation and timing relative 
- Observation subcomponent (spawner, catch, age-at-maturity, forecast) that informs assessment; decoupled from true dynamics unless linked by harvest control rule 