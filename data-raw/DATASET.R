## code to prepare `DATASET` dataset goes here

set.seed(1234)
df <- gen_lat_dim(ndims=4,n=16)

usethis::use_data(df, overwrite = TRUE, internal = TRUE)
