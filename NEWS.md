---
title: "NEWS"
output: html_document
---

# MDShelper 0.0-7

* All functions that fit MDS solutions (`loo_CV_MDS()`, `holdout_CV_MDS()`,
  `cross_validation_MDS()`, `compute_BIC_MDS()`) have new arguments `eps`
  and `itmax`, passed to `smacof::mds()`, with the defaults `eps = 1e-8` and
  `itmax = 1e5` (smacof: 1e-6 and 1000).
* Why: with missing cells, smacof's Torgerson start fills them with a
  constant. In as many or more dimensions than the data have, this start
  already fits the observed pairs almost perfectly, so with `eps = 1e-6` the
  algorithm stops after a few iterations and the held-out distance stays near
  the filled-in value. The cross-validation functions then favour too few
  dimensions. Results can therefore differ from 0.0-6, mainly for the
  cross-validation functions; on complete data (BIC) choices change rarely.
* Fitting takes longer (about 5-10 times for the cross-validation functions).

# MDShelper 0.0-6

* New `holdout_CV_MDS()` (hold-out cross-validation after Richie & Verheyen, 2020)
