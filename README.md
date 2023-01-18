
# MDShelper

<!-- badges: start -->
<!-- badges: end -->

The goal of MDShelper is to ...

## Installation

You can install the development version of MDShelper from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dizyd/MDShelper")
```

## Functions 

So far the package contains the following functions (use `?function_name` for more information):

## Example

This is a basic example which shows you how generate pairwise distance data based on a defined number of underlying dimensions and running a cross-validation procedure to try to recover this number of dimensions

``` r
library(MDShelper)

## Generate pairwise distance matrix
sim_data <- gen_data_MDS(ndims=4, n=16)

## Run Cross-Validation
cross_validation_MDS(sim_data$dist_mat, max_dim = 5)


```

