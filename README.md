This repo is associated with the manuscript, "Weighted Bayesian adaptive randomization for faster optimal treatment identification in N-of-1 trials".

It contains the following files:

```
.
├── BAR.R
├── WBAR.R
├── helpers.R
├── pipeline.R
├── compareRAR.R
└── weightingfunctions.R
```

- `pipeline.R`: main function for users. Used for running simulation studies for the paper. This sets up the simulation parameters, prepares the Stan model, and runs the replications. 
- `BAR.R`: contains an function of an implementation of standard Bayesian adaptive randomization (Thompson Sampling) with an MCMC sample
- `WBAR.R`: contains an function of our implementation of our weighted Bayesian adaptive randomization algorithm with an MCMC sample
- `helpers.R`: contains some helper functions used by the algorithm function
- `compareRAR.R`: a function for simulating an N-of-1 trial with different adaptive randomization algorithms. To be used with `pipeline.R`.
- `weightingfunctions.R`: contains functions for implementing different weighting functions to be used with WBAR
