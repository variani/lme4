blme/lme4 notes
===============

# What current blme implements

Both linear and generalized linear mixed models are available.

## Currently available covariance priors in blme
- none
  - flat prior
- correlation
  - $\Sigma = DRD$ where $D$ is diagonal and $R$ is a correlation matrix
  - place priors directly on the elements of $D$ and $R$
- spectral
  - $\Sigma = UDU^\top$ where $U$ is orthonormal $D$ is diagonal and positive
  - no priors on $U$ and independent and identical priors over eigenvalues
- direct
  - gamma, inverse gamma, wishart, or inverse wishart

## Currently available fixed effects priors in blme
- multivariate normal
  - only zero mean priors allowed
  - but full covariance matrix can be specified

## Currently available common scale parameter priors in blme
- seems like none are available, i.e. flat

# Where to look in new lme4 for required changes

1.  $\mathbf{R}_{X}^\top \mathbf{R}_{X} = \mathbf{V^\top V} - \mathbf{R}_{ZX}^\top \mathbf{R}_{ZX}$ gets changed to $\mathbf{R}_{X}^\top \mathbf{R}_{X} = \mathbf{V^\top V} + \mathbf{\Sigma_b}^{-1} - \mathbf{R}_{ZX}^\top \mathbf{R}_{ZX}$.  This happens in `predModule.cpp` in the `updateDecomp` function. Might want to think carefully about how to use the `rankUpdate` method here.  Can we spherize the fixed effects to make this simpler?  Maybe in the common scale case, we could really just pass fixed effects like we would random effects, but without updating the covariance parameters?

2.  In general I think that `updateDecomp` will be the main place we'll need to provide hooks for when adding a prior to fixed effects.  Where else will we need them?  In any case, this should be fairly dangerous (for me) given that I'm not very comfortable with `Eigen`.

3.  In terms of machinery for optimizing the common scale, when we choose to parameterize independently of it, I'm not sure where we need hooks.

