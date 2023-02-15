
# MDShelper

<!-- badges: start -->
<!-- badges: end -->

The goal of `MDShelper` is to provide some convenience functions when conducting multidimensional-scaling (MDS) analysis. So far the package contains the following functions (use `?function_name` for more information):

- `cross_validation_MDS()`:  Runs a cross-validation procedure to find the best fitting number of dimensions of the MDS space.

- `BIC_MDS()`: Computes the BIC for a specific number of dimensions of the MDS space according to Lee (2001)


- `gen_data_MDS()`: Simulates data by generating a matrix with the true underlying dimensions and item coordinates, as well as the corresponding true pairwise distances

- `fill_mat()`: Transforms a distance matrix with only lower triangle entries (rest `NA`) to full matrix with 0 on the diagonals

- `compute_dists()`: Wrapper around a `Rcpp`-function (`pair_dist_cpp`)  to compute pairwise distances between two vectors.



## Installation

You can install the development version of MDShelper from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dizyd/MDShelper")
```


## Example

This is a basic example which shows you how generate pairwise distance data based on a defined number of underlying dimensions and running a cross-validation procedure to try to recover this number of dimensions

``` r
library(MDShelper)

## Generate pairwise distance matrix
sim_data <- gen_data_MDS(ndims = 4, n = 16)

## Run Cross-Validation
cross_validation_MDS(sim_data$dist_mat, max_dim = 5)


```

