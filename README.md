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
x = expr
plot(x[,1],x[,2],xlim=c(-8,8),ylim=c(-8,8),cex = .4);
idr.out = IDR.3component(x = x)
```
