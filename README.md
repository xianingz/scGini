# scGini
scGini implements a Negative-Binomial model for recovering the true Gini index for individual single cell based on scRNA-seq data, which is an estimation of the stemness level of the single cell.

## Installation
The latest version of scGini R package can be installed from Github.  
```r
devtools::install_github("xianingz/scGini")
library(scGini)
```

## Goodness of fit  
`NBFitPlot` function is used to visualiza the goodness of fit of the Negative-Binomial model for individual single cell.
```r
NBFitPlot(scdat, seed=1)
```

## Calculation of corrected Gini index  
`scGini` function calculates the corrected Gini index by using the raw counts scRNA-seq data as input (without any processings like normalization). Parallel computing can be enabled by setting the `ncore` parameter.
```r
scGini(dat, ncore=1)
```


## Reference

