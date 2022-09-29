
<!-- README.md is generated from README.Rmd. Please edit that file -->

# syntenet <img src="man/figures/logo.png" align="right" height="138" />

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/almeidasilvaf/syntenet)](https://github.com/almeidasilvaf/syntenet/issues)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check-bioc](https://github.com/almeidasilvaf/syntenet/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/almeidasilvaf/syntenet/actions)
[![Codecov test
coverage](https://codecov.io/gh/almeidasilvaf/syntenet/branch/master/graph/badge.svg)](https://codecov.io/gh/almeidasilvaf/syntenet?branch=master)
<!-- badges: end -->

The goal of `syntenet` is to infer synteny networks from whole-genome
protein sequence data and analyze them. Anchor pairs from synteny
analyses are treated as an undirected unweighted graph (i.e., a synteny
network), and users can perform:

-   **Synteny detection** using a native implementation of the [MCScanX
    algorithm](https://doi.org/10.1093/nar/gkr1293), a C++ program that
    has been modified and ported to R with Rcpp. This way, users do not
    need to install MCScanX beforehand, because `syntenet` has its own
    implementation of the same algorithm.
-   **Synteny network inference** by treating anchor pairs as edges of a
    graph;
-   **Network clustering** using the Infomap algorithm;
-   **Phylogenomic profiling**, which consists in identifying which
    species contain which clusters. This analysis can reveal highly
    conserved synteny clusters and taxon-specific ones (e.g., family-
    and order-specific clusters);
-   **Microsynteny-based phylogeny reconstruction** with maximum
    likelihood, which can be achieved by inferring a phylogeny from a
    binary matrix of phylogenomic profiles with IQTREE.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `syntenet` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("syntenet")
```

And the development version from
[GitHub](https://github.com/almeidasilvaf/syntenet) with:

``` r
BiocManager::install("almeidasilvaf/syntenet")
```

## Citation

Below is the citation output from using `citation('syntenet')` in R.
Please run this yourself to check for any updates on how to cite
**syntenet**.

``` r
print(citation('syntenet'), bibtex = TRUE)
#> 
#> To cite package 'syntenet' in publications use:
#> 
#>   Almeida-Silva F, Zhao T, Ullrich K, Van de Peer Y (2022). _syntenet:
#>   Inference And Analysis Of Synteny Networks_. R package version
#>   0.99.4, <https://github.com/almeidasilvaf/syntenet>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {syntenet: Inference And Analysis Of Synteny Networks},
#>     author = {Fabrício Almeida-Silva and Tao Zhao and Kristian K Ullrich and Yves {Van de Peer}},
#>     year = {2022},
#>     note = {R package version 0.99.4},
#>     url = {https://github.com/almeidasilvaf/syntenet},
#>   }
```

Please note that `syntenet` was only made possible thanks to many other
R and bioinformatics software authors, which are cited either in the
vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `syntenet` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

-   Continuous code testing is possible thanks to [GitHub
    actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
    through *[usethis](https://CRAN.R-project.org/package=usethis)*,
    *[remotes](https://CRAN.R-project.org/package=remotes)*, and
    *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)*
    customized to use [Bioconductor’s docker
    containers](https://www.bioconductor.org/help/docker/) and
    *[BiocCheck](https://bioconductor.org/packages/3.15/BiocCheck)*.
-   Code coverage assessment is possible thanks to
    [codecov](https://codecov.io/gh) and
    *[covr](https://CRAN.R-project.org/package=covr)*.
-   The [documentation website](http://almeidasilvaf.github.io/syntenet)
    is automatically updated thanks to
    *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
-   The code is styled automatically thanks to
    *[styler](https://CRAN.R-project.org/package=styler)*.
-   The documentation is formatted thanks to
    *[devtools](https://CRAN.R-project.org/package=devtools)* and
    *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.15/biocthis)*.
