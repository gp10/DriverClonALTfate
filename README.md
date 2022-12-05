# DriverClonALTfate :: Dynamics of Driver mutant Clones with ALTered progenitor cell fate
Computational methods to infer the altered progenitor behavior of a particular driver mutant in mouse epithelium. This set of methods allows to determine which mechanism drives the competitive advantage of driver mutant clones over WT counterparts based on modelling and fitting of experimental lineage tracing data in transgenic mice.

This repository contains code that was originally used in the following manuscript:
  > Herms A, Colom B, Piedrafita G, et al. (2020) Competitive fitness of Pik3ca mutant clones in normal esophagus is determined by the metabolic state of wild type cells. _submitted_

### Graphical abstract
![GraphicalAbstract](https://github.com/gp10/DriverClonALTfate/blob/master/Graphical_abstract_DriverClonALTfate.png)

### Overview
The computational methods included in this repository to infer WT and mutant cell behavior are based on the **Single Progenitor (SP) model** paradigm of cell renewal in murine squamous epithelium. It is assumed that both WT and mutant clones rely on a single population of progenitor cells undergoing stochastic fates but with different dynamical attributes. While the division rate (_Lambda_) is assumed invariant between WT and mutant progenitors based on experimental data, unknown parameters are the probability of a symmetric division outcome leading to two dividing or two differentiating cells (_r_), the stratification rate (_Gamma_) (which is related with the proportion of proliferating cells in the basal layer (_rho_), and in the case of the mutant, a possible fate imbalance between the two possible symmetric division outcomes (_Delta_).

In order to determine the likeliest mechanism of mutant clonal advantage, the code includes Monte-Carlo methods to simulate the SP model dynamics in the basal compartment under given parameter conditions. A script called **Analysis-BasalCloneDynamics-logLike.m** allows to define a set of parameter conditions for both WT and mutant progenitors and fit the corresponding simulated clonal data to the experimental clone sizes at the different time points, computing log-likelihood values. This operation would need to be repeated iteratively for different combinations of parameter values to determine the maximum likelihood estimate for WT and mutant populations. In practice, the `Datasets` folder contains files that allow loading the log-likelihood values obtained for each genotype after a large parameter grid search in a supercomputing facility, so as to avoid computation time to user. The **Analysis_MLE_WT** and **Analysis_MLE_Mut** scripts then allow to retrieve the maximum likelihood estimate (MLE) on the SP model parameters for WT and mutant, respectively.

The Monte-Carlo simulators are extended to explore total clone sizes including dynamics in the suprabasal compartment. In **Analysis_TotalCloneDynamics** one can fit simulated data onto experimental total clone sizes and study the dynamics of floating clones having no basal attachment.

### Main scripts
- **Analysis-BasalCloneDynamics-logLike.m** : main script to run simulations of basal clone sizes under specific parameter conditions, plot fits on experimental basal clone data (average basal clone size and distributions of the no. of basal cells / clone) and compute corresponding log-likelihood values.
- **Analysis-MLE-WT.m** : script used to compute and plot best-fit SP-model parameter values (on WT basal clone dynamics) from a matrix of all log-likelihood values obtained by parameter grid search changing _rho_ and _r_. A maximum likelihood estimate (MLE) is obtained with 95% CI on parameter values based on Likelihood-Ratio (LR)-test. Results from the LR-test are plotted as a heatmap on the 2D parameter space.
- **Analysis-MLE-Mut.m** : script used to compute and plot best-fit SP-model parameter values (on Mutant basal clone dynamics) from a matrix of all log-likelihood values obtained by parameter grid search changing _Delta_, _rho_ and _r_. A maximum likelihood estimate (MLE) is obtained with 95% CI on parameter values based on Likelihood-Ratio (LR)-test. Results from the LR-test are plotted as a heatmap on the 3D parameter space.
- **Analysis-TheoreticalSuprabTotalRatioScaling.m** : script plotting how the proportion of suprabasal cells per clone is inferred to scale (decrease) when introducing a progenitor fate imbalance (_Delta_).
- **Analysis-TotalCloneDynamics.m** : script to run simulations of total clone sizes under specific parameter conditions (typically bounded to MLE values computed above for basal layer), plot fits on experimental total clone data (average total clone size and distributions of the no. of total cells / clone) and validate the inferred clonal dynamics in the suprabasal compartment. It includes an analysis of the proportion of floating clones.

### Dependencies
- Monte Carlo simulators :
  - MonteCarloSimulator-SP-basal.m : Monte Carlo simulator of basal clone dynamics under the SP model with balanced fates.
  - MonteCarloSimulator-SP-total-varMu.m : Monte Carlo simulator of basal & total clone dynamics under the SP model with balanced fates (a time-dependent shedding rate is assumed).
  - MonteCarloSimulator-SP-total-varMu-evenStart.m : Monte Carlo simulator of basal & total clone dynamics under the SP model with balanced fates (a time-dependent shedding rate is assumed) - version considering an even initial induction condition not restricted to progenitor cells.
  - MonteCarloSimulator-SP-basal-unbalance.m : Monte Carlo simulator of basal clone dynamics under the SP model with unbalanced fates (_Delta_).
  - MonteCarloSimulator-SP-total-unbalance-varMu.m : Monte Carlo simulator of basal clone dynamics under the SP model with unbalanced fates (_Delta_) (a time-dependent shedding rate is assumed).
  - MonteCarloSimulator-SP-total-unbalance-varMu-evenStart.m : Monte Carlo simulator of basal & total clone dynamics under the SP model with unbalanced fates (_Delta_) (a time-dependent shedding rate is assumed) - version considering an even initial induction condition not restricted to progenitor cells.

- Functions for specific calculations:
  - calculate-AvgCloneSize.m : function called to calculate the average no. of (basal or total) cells / surviving clone over time, from simulation data.
  - calculate-CloneSizeDist.m : function called to calculate the relative frequecies of clones of different (basal or total) sizes at the different time points, either for simulation or experimental data.
  - calculate-FloatingClones.m : function called to calculate the proportion of floating clones over time, either for simulation or experimental data.
  - calculate-logLike-bootstrapping.m : function called to calculate log-likelihood value estimates for subsets (built by random permutation) of experimental data (bootstrapping) using experimental and simulated basal clone sizes as input.
  - calculate-MLE-2D.m : function called to compute the maximum likelihood estimate (MLE; with 95% CI) on the SP-model parameters _rho_ and _r_ for the WT population.
  - calculate-MLE-3D.m : function called to compute the maximum likelihood estimate (MLE; with 95% CI) on the SP-model parameters _Delta_, _rho_ and _r_ for the Mutant population.

- Plotting scripts:
  - plot-AvgCloneSize-experimental.m : computes and plots the average (basal or total) experimental clone sizes over time.
  - plot-AvgCloneSize-simulated.m : plots the average (basal or total) clone sizes over time from the SP-model simulations.
  - plot-CloneSizeDist.m : plots the (basal or total) clone size distributions at the different time points overlying simulation fits on experimental frequencies.
  - plot-FloatingClones-experimental.m : plots the percentage of floating clones over time from the experimental data.
  - plot-FloatingClones-simulated.m : plots the percentage of floating clones over time from the SP-model simulations.
  - plot-MLE-2D.m : plots a 2D histogram of the likelihood ratio (LR)-test value for the different parameter sets {_rho_, _r_} within the domain of explored parameter space (used for WT).
  - plot-MLE-3D.m : plots a 3D histogram of the likelihood ratio (LR)-test value for the different parameter sets {_Delta_, _rho_, _r_} within the domain of explored parameter space (used for Mutant).

- Others:
  - SelectModelParamVal.m : allows loading specific parameter values of interest corresponding to MLE values found in the manuscript for WT or Mutant progenitors.
  - size2freq.m : calculates the frequency histogram (distribution) of experimental or simulated clone sizes from their individual sizes.
  - size2freqbinned.m : calculates the frequency histogram (distribution) for experimental or simulated clone sizes binned in categories increasing in size in powers of 2.
  - freq2squeeze.m : allows reassigning the few outlier clones into the maximum ordinary binned size category instead of discarding them.
  - logLike-calc.m : computes the log-Likelihood match of simulated vs. experimental clone size distributions.
  - subsampling-experim-clones.m : partitions experimental clones into random subsets (sampling by random permutation from the pool of experimental data), for the purpose of ensuring robustness in MLE calculations.
  - subsampling-simulated-clones.m : partitions simulated clones into random subsets (sampling by random permutation from the pool of simulated data), for the purpose of confidence interval calculations.
  - subsampling-clones-overtime.m : partitions simulated clones into random subsets (sampling by random permutation from the pool of simulated data), for the purpose of confidence interval calculations (alternative version).
  - rotateXLabels.m : allows rotating x-axis tick labels. Imported from MathWorks. See license details in: Ben Tordoff (2020). rotateXLabels( ax, angle, varargin ) (https://www.mathworks.com/matlabcentral/fileexchange/27812-rotatexlabels-ax-angle-varargin), MATLAB Central File Exchange. Retrieved July 6, 2020.

- `Datasets` folder : contains experimental lineage tracing data and previously calculated log-likelihood values for a large collection of parameter values both for WT and mutant populations.

### Requirements
- Matlab (built on version R2016b and tested on R2016b & R2020b)
- No dependencies needed
- Runs on all OS


