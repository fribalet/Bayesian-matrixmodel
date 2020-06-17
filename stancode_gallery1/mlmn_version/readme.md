# multi-level multinomial (mlmn) version of the gallery1 code

**Note**: Code has not yet been fully tested. Some changes were made to the model parameterization and priors; these are detailed below. The goal of these changes was to parameterize the Dirichlet-Multinomial model in a way that consistently converges when applied to SeaFlow data. Most of the changes were aimed at reducing the number of MCMC samples that were rejected due to violating parameter constraints, which was somewhat prevalent during the initial warm-up phase. Initial size distribution estimation often returned nan values, so the parameterization for this was completely redefined in a more flexible manner. The monotone division samples often caused deltas to exceed their upper bound in some samples, hence an explicit upper bound was added to the deltas and the monotone division was reparameterized to be fully constrained within these bounds. Priors for gamma and rho were loosened to cover the full possible range of these parameters, though with most of the probability mass concentrated at low values. Other parameter bounds were tweaked as well to avoid violating any model constraints (details below). The result was convergence on both SeaFlow and Zinser data.

**These changes are contained in the file** `matrixmodel_mlmultinom_estinilnorm2_monodelta2_respv2_normparam_trackgrowth_xval2.stan`.

Python code can be found in `experimental_zinser_seaflow_20200616_gallery1_mlmn_truecount.ipynb` (uncompiled) and `experimental_zinser_seaflow_20200616_gallery1_mlmn_truecount.ipynb`.

### INITIAL SIZE DISTRIBUTION:
Changed estimation of initial size distribution:
- Removed parameter "w_ini_mu"
- Removed parameter "w_ini_sigma"
- Moved parameter "w_ini" from transformed parameter block to parameter block, parameterized as simplex (sum-to-one vector)
- No prior specified here, hence STAN uses a diffuse prior (could make more informative)

### MONOTONE DIVISION:
- Renamed parameters "delta_max" to "delta"
- Renamed "delta_max_incr" to "delta_incr" and reparameterized as simplex
- Added scaling parameter "delta_max," define "delta" as cumulative sum of "delta_incr" times "delta_max"
- Again, didn't specify priors for either of these; could definitely come up with more informative priors especially for "delta_max"

### PARAMETER CONSTRAINTS:
- Added upper bound of 1.0/dt_days to "delta" and upper bound of 1.0/dt_norm to "rho_max"
- Increased upper bound of "gamma_max" to 1.0/dt_norm

### PRIORS:
- Changed priors for "gamma_max" and "rho_max" to truncated normal distributions to cover the full range
- Added truncation to "E_star" to explicitly constrain the prior to be nonnegative
- Changed "sigma" prior from exponential to lognormal to get a somewhat broader prior, added truncation with minimum of 1 (sigma of 1 corresponds to observing only a single particle, so it is a safe lower bound)

### MODEL LIKELIHOOD:
- Added 1 to all entries of "alpha" used as mean of Dirichlet distribution to avoid zero-valued entries (Dirichlet parameters must be strictly positive)
