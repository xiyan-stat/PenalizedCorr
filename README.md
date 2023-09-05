# PenalizedCorr

The PenalizedCorr package is used to estimate the autocorrelation and partial autocorrelation function using penalized methods. This is a tutorial example for the R package from the paper **Penalised M-estimation for Autocorrelation Colin Gallagher, Xiyan Tan**

## Quick start guide
### load package

```{r}
# Load the package
library(PenalizedCorr)

# Read or create time series data
## Example 1: univariate time series generated from AR(1) with parameter 0.5
data1 = arima.sim(n=100, model = list(ar=0.5))
# visualize data
plot(data1)
```
<img src="fig/example1.png" width="864" /> 

