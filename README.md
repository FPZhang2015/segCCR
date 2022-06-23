## segCCR
Implementation of estimation and testing for the segmented correspondence curve regression. It fits the segmented correspondence curve regression by grid search method. A test procedure is proposed to test the existence of the change point in the segmented correspondence curve regression.

## Installing
The package can be installed from github. R package **devtools** is required.

```
library(devtools)
install_github("qunhualilab/segCCR")
library(segCCR)
```

## Example

```
data(ChIPseq)

m = 100
tm <- seq(0.01, 0.999, length.out = m)
nx = nlevels(factor(ChIPseq$x))
par.ini = c(0.5, 2, 1, rep(0.1, 2*(nx-1)))  # initial value
## estimate
fit = segCCR(data = ChIPseq,
           par.ini = par.ini,
           tm=tm,
           NB = 5)      

## test the existence of a change point
fit.test = segCCR_test(ChIPseq,
           tm=tm,
           NB = 10)
```

## Citation

Please cite the following paper:Zhang, F. and Li, Q. (2022).Segmented correspondence curve regression for quantifying
covariate effects on the reproducibility of high-throughput experiments. Biometrics. 
