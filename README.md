
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastGFPCA

<!-- badges: start -->
<!-- badges: end -->

fast generalized functional principal components analysis

-   Authors: [Julia Wrobel](http://juliawrobel.com), Ciprian
    Crainiceanu, Andrew Leroux
-   License: [MIT](https://opensource.org/licenses/MIT). See the
    [LICENSE](LICENSE) file for details
-   Version: 0.9

## Installation

You can install the development version of fastGFPCA from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("julia-wrobel/fastGFPCA")


# build vignette with installation
devtools::install_github("julia-wrobel/fastGFPCA", build_vignettes = TRUE)
```

## Example

This example performs GFPCA on simulated binary data. More details on
the use of the package can be found in the vignette.

``` r
library(fastGFPCA)

# simulate data
set.seed(10001)
df_gfpca <- sim_gfpca(N = 50, J = 500, case = 1)$df_gfpca

# run gfpca
bin_mod <- fast_gfpca(Y = df_gfpca, overlap = TRUE, binwidth = 20, npc = 4, family = "binomial")
```
