### Fate Mapping Quantifies the Dynamics of B Cell Development and Activation throughout Life

https://www.sciencedirect.com/science/article/pii/S2211124720313656

#### coding
The models directory has stan files for different models.
The scripts directory has Rscripts that call stan models to fit on the data.

sh files are designed to run the Rscripts in a batch mode in order to submit them as parallel jobs across multiple nodes on the local cluster at Yate's lab at CUMC, New York.
The Rscript named "stan_plots.R" in the /scripts dircetory was used to complie and combine parallel jobs that ran on different nodes within the local cluster. The combined fitted object was then used to make inferences about parameter estimates, model suitability and plots for post-predictive analyses.

#### Modelling 
We aimed to quantify the dynamics of various subsets within mature B cell population and to understand the rules of replacement of old cells by that of new ones within each subset.

FM cells were assumed to circulate freely in the lymphatic system, and so we pooled the numbers of these cells recovered from spleen and lymph nodes when modelling their dynamics. We model Spleen and LN GC cells separately assuming that circulating FM cells feed into both these populations with a constant rate over time.

The timecourses of total cell counts, Normalised chimerism and the porportions of Ki67hi cells in host and donor compartments were used to confront and compare various models.

Models tested here:

    Simple Homoegenous model (SHM): constant rates of loss and division.
    Time-dependent model (TDM): Either loss or division varies with host-age.
    Kinetic heterogeneity model (KHM): Fast and Slow (with respect to their division and turnover) subsets within host and donor compartments.
    Incumbent model (INC): Fast and Slow (pre-existing) subsets only within host.

We also tried variations withing each model by comparing different sources, by varying the rate of change of source influx with time and by keeping the rate of loss of Ki67 (beta) either free or fixed (1/3.5 days). 

