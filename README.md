# BiModel <img src='man/logo/BiModel_logo.png' align="right" height="140" />
In presented repository the binomial distribution mixture model is applied.

## Installation
You can install the package from [GitHub](https://github.com/) with:
``` r
# install.packages("devtools")
devtools::install_github("karowid617/BiModel")
```

## Manual
Loading of library and examplary data.
``` r
library(BiModel)
data(example)
```
Using the examplary data the following code run the BMM decomposition where 2 clusters will be searched based on K-means initalization.

``` r
BMM_res<-BernoulliEM(example$nouli_data, K = 2,start_ini = 20,ini = "kmeans",m_iter = 3000,eps = 1e-40)
```

## References
Polanski, A., Marczyk, M., Pietrowska, M., Widlak, P., Polanska, J., 2018. Initializing the em algorithm for univariate gaussian, multi-component, heteroscedastic mixture models by dynamic programming partitions. International Journal of Computational Methods 15, 1850012.

