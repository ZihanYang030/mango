# **mango**
MANGO (MAppiNG multi-cell-type multi-scale spatial relationships in spatially resolved Omics) \
is a permutation-based computational framework for detecting higher-order spatial co-localization across multiple spatial scales. \
MANGO defines the local microenvironment of each reference cell using a k-nearest-neighbor neighborhood with Gaussian distance-decay weighting. Within each neighborhood, MANGO computes a multivariate local co-localization quotient (mLCLQ) that quantifies whether a predefined set of query cell types co-occurs more frequently than expected given their global spatial distributions. Significance is assessed through permutation testing. MANGO is implemented as an open-source R package, available. 

<img align="top" src="man/figures/mango_overview.png" alt="drawing" width="600"/>

## Install the Package
Install mango from GitHub:
```r
remotes::install_github("ZihanYang030/MANGO")
library(MANGO)
```

## License

mango is licensed under the GNU General Public License v3.0 or any later version.
