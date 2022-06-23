## segCCR
Implementation of estimation and testing for the segmented correspondence curve regression. It fits the segmented correspondence curve regression by grid search method. 
A test procedure is proposed to test the existence of the change point in the segmented correspondence curve regression.

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
## estimate
m = 100
tm <- seq(0.01, 0.999, length.out = m)
nx = nlevels(factor(ChIPseq$x))
par.ini = c(0.5, 2, 1, rep(0.1, 2*(nx-1)))  # initial value
nx = nlevels(factor(ChIPseq$x))
par.ini = c(0.5, 2, 1, rep(0.1, 2*(nx-1)))  # initial value
fit = segCCR(data = ChIPseq,
           par.ini = par.ini,
           tm=tm,
           NB = 5)           
```


## Citation

Please cite the following paper: A semi-parametric statistical model for integrating gene expression profiles across different platforms
