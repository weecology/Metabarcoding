Plants by Treatment
================
Ellen Bledsoe,
March 2021

## Annuals by Sampling Period and Treatment

Running CCA (from Supp et al., 2012 and Portal-LTREB) to see if there
are significant differences in plant communities between controls and KR
exclosures for each season of diet
    sampling.

### Summer Annuals 2016

#### Bray-Curtis Distances & Heatmap

    ##   year season plot quads treatment control.4 control.11 control.14 control.17
    ## 1 2016 summer    4    16   control 0.0000000  0.4835996  0.3536797  0.2378190
    ## 2 2016 summer   11    16   control 0.4835996  0.0000000  0.4506154  0.3317974
    ## 3 2016 summer   14    16   control 0.3536797  0.4506154  0.0000000  0.2740507
    ## 4 2016 summer   17    16   control 0.2378190  0.3317974  0.2740507  0.0000000
    ## 5 2016 summer    3    16 exclosure 0.2332703  0.4505304  0.4168591  0.3600123
    ## 6 2016 summer   15    16 exclosure 0.4170132  0.1966651  0.3624720  0.2551381
    ## 7 2016 summer   19    16 exclosure 0.3025542  0.4418471  0.4240323  0.3564089
    ## 8 2016 summer   20    16 exclosure 0.3658103  0.4582654  0.2796710  0.3171195
    ## 9 2016 summer   21    16 exclosure 0.3086403  0.5235435  0.3734916  0.3580675
    ##        ex.3     ex.15     ex.19     ex.20     ex.21
    ## 1 0.2332703 0.4170132 0.3025542 0.3658103 0.3086403
    ## 2 0.4505304 0.1966651 0.4418471 0.4582654 0.5235435
    ## 3 0.4168591 0.3624720 0.4240323 0.2796710 0.3734916
    ## 4 0.3600123 0.2551381 0.3564089 0.3171195 0.3580675
    ## 5 0.0000000 0.4843415 0.2681849 0.3255163 0.2552856
    ## 6 0.4843415 0.0000000 0.4602162 0.4049219 0.4612882
    ## 7 0.2681849 0.4602162 0.0000000 0.2184779 0.2412041
    ## 8 0.3255163 0.4049219 0.2184779 0.0000000 0.2895958
    ## 9 0.2552856 0.4612882 0.2412041 0.2895958 0.0000000

![](Plants_by_Treatment_files/figure-gfm/annual_summer2016_heatmap-1.png)<!-- -->

#### CCA by Treatment

![](Plants_by_Treatment_files/figure-gfm/annual_summer2016_cca-1.png)<!-- -->

    ## [1] "Variance Inflation Factor (<10 is fine)"

    ## summerannuals_2016_trtexclosure 
    ##                               1

    ## [1] "Variance explained:"

    ## [1] 0.1117901

    ## Permutation test for cca under NA model
    ## Marginal effects of terms
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: cca(formula = summerannuals_2016_sp ~ summerannuals_2016_trt)
    ##                        Df ChiSquare     F Pr(>F)
    ## summerannuals_2016_trt  1   0.05687 0.881  0.617
    ## Residual                7   0.45185

### Winter Annuals 2017

#### Bray-Curtis Distances & Heatmap

    ##   year season plot quads treatment control.4 control.11 control.14 control.17
    ## 1 2017 winter    4    16   control 0.0000000  0.4287844  0.5515273  0.3958126
    ## 2 2017 winter   11    16   control 0.4287844  0.0000000  0.3470930  0.2995290
    ## 3 2017 winter   14    16   control 0.5515273  0.3470930  0.0000000  0.4049866
    ## 4 2017 winter   17    16   control 0.3958126  0.2995290  0.4049866  0.0000000
    ## 5 2017 winter    3    16 exclosure 0.3268529  0.4888661  0.5193372  0.5416008
    ## 6 2017 winter   15    16 exclosure 0.4925699  0.3684930  0.4495784  0.3261664
    ## 7 2017 winter   19    16 exclosure 0.4290136  0.3654876  0.5450190  0.4302616
    ## 8 2017 winter   20    16 exclosure 0.4873352  0.4075804  0.3912079  0.3939218
    ## 9 2017 winter   21    16 exclosure 0.3170895  0.4774594  0.4732551  0.4452816
    ##        ex.3     ex.15     ex.19     ex.20     ex.21
    ## 1 0.3268529 0.4925699 0.4290136 0.4873352 0.3170895
    ## 2 0.4888661 0.3684930 0.3654876 0.4075804 0.4774594
    ## 3 0.5193372 0.4495784 0.5450190 0.3912079 0.4732551
    ## 4 0.5416008 0.3261664 0.4302616 0.3939218 0.4452816
    ## 5 0.0000000 0.5686514 0.4734896 0.4524324 0.3295173
    ## 6 0.5686514 0.0000000 0.4720061 0.3889509 0.4356436
    ## 7 0.4734896 0.4720061 0.0000000 0.4424799 0.4020361
    ## 8 0.4524324 0.3889509 0.4424799 0.0000000 0.3708186
    ## 9 0.3295173 0.4356436 0.4020361 0.3708186 0.0000000

![](Plants_by_Treatment_files/figure-gfm/annual_winter2017_heatmap-1.png)<!-- -->

#### CCA by Treatment

![](Plants_by_Treatment_files/figure-gfm/annual_winter2017-1.png)<!-- -->

    ## [1] "Variance Inflation Factor (<10 is fine)"

    ## winterannuals_2017_trtexclosure 
    ##                               1

    ## [1] "Variance explained:"

    ## [1] 0.1674206

    ## Permutation test for cca under NA model
    ## Marginal effects of terms
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: cca(formula = winterannuals_2017_sp ~ winterannuals_2017_trt)
    ##                        Df ChiSquare      F Pr(>F)  
    ## winterannuals_2017_trt  1   0.11540 1.4076  0.081 .
    ## Residual                7   0.57389                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Summer Annuals 2017

#### Bray-Curtis Distances & Heatmap

![](Plants_by_Treatment_files/figure-gfm/annual_summer2017_heatmap-1.png)<!-- -->

#### CCA by Treatment

![](Plants_by_Treatment_files/figure-gfm/annual_summer2017-1.png)<!-- -->

    ## [1] "Variance Inflation Factor (<10 is fine)"

    ## summerannuals_2017_trtexclosure 
    ##                               1

    ## [1] "Variance explained:"

    ## [1] 0.1090373

    ## Permutation test for cca under NA model
    ## Marginal effects of terms
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: cca(formula = summerannuals_2017_sp ~ summerannuals_2017_trt)
    ##                        Df ChiSquare      F Pr(>F)
    ## summerannuals_2017_trt  1   0.05155 0.8567  0.701
    ## Residual                7   0.42123

## All Plants by Sampling Period and Treatment

### Summer All Plants 2016

![](Plants_by_Treatment_files/figure-gfm/all_summer2016-1.png)<!-- -->

    ## [1] "Variance Inflation Factor (<10 is fine)"

    ## summer_2016_trtexclosure 
    ##                        1

    ## [1] "Variance explained:"

    ## [1] 0.1326964

    ## Permutation test for cca under NA model
    ## Marginal effects of terms
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: cca(formula = summer_2016_sp ~ summer_2016_trt)
    ##                 Df ChiSquare     F Pr(>F)
    ## summer_2016_trt  1   0.10162 1.071  0.328
    ## Residual         7   0.66418

### Winter All Plants 2017

![](Plants_by_Treatment_files/figure-gfm/all_winter2017-1.png)<!-- -->

    ## [1] "Variance Inflation Factor (<10 is fine)"

    ## winter_2017_trtexclosure 
    ##                        1

    ## [1] "Variance explained:"

    ## [1] 0.1477049

    ## Permutation test for cca under NA model
    ## Marginal effects of terms
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: cca(formula = winter_2017_sp ~ winter_2017_trt)
    ##                 Df ChiSquare      F Pr(>F)
    ## winter_2017_trt  1   0.14903 1.2131  0.146
    ## Residual         7   0.85992

### Summer All Plants 2017

![](Plants_by_Treatment_files/figure-gfm/all_summer2017-1.png)<!-- -->

    ## [1] "Variance Inflation Factor (<10 is fine)"

    ## summer_2017_trtexclosure 
    ##                        1

    ## [1] "Variance explained:"

    ## [1] 0.1229296

    ## Permutation test for cca under NA model
    ## Marginal effects of terms
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: cca(formula = summer_2017_sp ~ summer_2017_trt)
    ##                 Df ChiSquare      F Pr(>F)
    ## summer_2017_trt  1   0.07275 0.9811  0.536
    ## Residual         7   0.51907
