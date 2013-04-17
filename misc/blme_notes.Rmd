blme/lme4 notes
===============

Both linear and generalized linear mixed models are available.

# Currently available covariance priors in blme
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

# Currently available fixed effects priors in blme
- multivariate normal
  - only zero mean priors allowed
  - but full covariance matrix can be specified

# Currently available common scale parameter priors in blme
- seems like none are available, i.e. flat
