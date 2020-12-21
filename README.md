# BTSM: Bayesian Time Series Models

<!-- badges: start -->
[//]: # "[![CRAN](http://www.r-pkg.org/badges/version/BGVAR)](https://cran.r-project.org/package=BGVAR)"
[//]: # "[![month](http://cranlogs.r-pkg.org/badges/BGVAR)](https://www.r-pkg.org/pkg/BGVAR)"
[//]: # "[![total](http://cranlogs.r-pkg.org/badges/grand-total/BGVAR)](https://www.r-pkg.org/pkg/BGVAR)"
<!-- badges: end -->

This repository should serve as a compendium of various Bayesian time series models. It consists of estimation procedures and various functions for analysing time series data efficiently. This toolbox does not claim to be complete and will be expanded by demand of its users.

## Installation

The latest development version can be installed from GitHub.

``` r
devtools::install_github("mboeck11/BTSM")
```

Note that Mac OS needs gfortran binary packages to be installed. See also: https://gcc.gnu.org/wiki/GFortranBinaries.

## Usage

This package focuses on multivariate time series models estimated in a Bayesian fashion. The table below summarises available estimation function for different models and which tool functions are available for each model.

| Model                                                 | Estimation  | IRFs | Predictions | FEVD | HD  |
|-------------------------------------------------------|-------------|------|-------------|------|-----|
| Bayesian Vector Autoregression                        | `bvar()`    | yes  | yes         | yes  | yes |
| Bayesian Vector Error Correction Model                | `bvec()`    | yes  | no          | no   | no  |
| Bayesian Interacted Vector Autoregression             | `bivar()`   | yes  | no          | no   | no  |
| Bayesian Panel Vector Autoregressions                 | `bpvar()`   | no   | no          | no   | no  |
| Bayesian Interacted Panel Vector Autoregression       | no          | no   | no          | no   | no  |
| Time-varying Parameter Bayesian Vector Autoregression | `tvpbvar()` | no   | no          | no   | no  |

We also provide S3methods for each function as summarized in the next table

| Model                                                 | Function    | `print()` | `summary()` | `logLik()` | `coef()` | `vcov()` | `resid()` |
|-------------------------------------------------------|-------------|-----------|-------------|------------|----------|----------|-----------|
| Bayesian Vector Autoregression                        | `bvar()`    | yes       | yes         | yes        | no       | no       | no        |
| Bayesian Vector Error Correction Model                | `bvec()`    | no        | no          | no         | no       | no       | no        |
| Bayesian Interacted Vector Autoregression             | `bivar()`   | no        | no          | no         | no       | no       | no        |
| Bayesian Panel Vector Autoregressions                 | `bpvar()`   | no        | no          | no         | no       | no       | no        |
| Bayesian Interacted Panel Vector Autoregression       | no          | no        | no          | no         | no       | no       | no        |
| Time-varying Parameter Bayesian Vector Autoregression | `tvpbvar()` | no        | no          | no         | no       | no       | no        |

## References

Arias, J.E., Rubio-Ramírez, J.F., and D.F. Waggoner (2018) Inference based on structural vector autoregressions identified with sign and zero restrictions: Theory and applications *Econometrica: Journal of the Econometric Society*, Vol. 86(2), pp. 685–720.

Bitto, A. and S. Frühwirth-Schnatter (2019) Achieving shrinkage in a time-varying parameter model framework. *Journal of Econometrics*, Vol. 210(1), pp. 75-97.

Doan, T., Litterman, R., and C. Sims (1984) Forecasting and conditional projection using realistic prior distributions *Econometric reviews*, Vol. 3(1), pp. 1–100.

George E. I., Sun D., and S. Ni (2008) Bayesian stochastic search for VAR model restrictions. *Journal of Econometrics*, Vol. 142(1), pp. 553–580. 

Giannone D., Lenza M., and G. E. Primiceri (2015) Prior selection for vector autoregressions. *Review of Economics and Statistics*, Vol. 97(2), pp. 436–451.

Huber, F. and M. Feldkircher (2016) Adaptive Shrinkage in Bayesian Vector Autoregressive Models. *Journal of Business and Economic Statistics*, Vol. 37(1), pp. 27-39.

Huber F., Kastner, G. and M. Feldkircher (2019) Should I stay or should I go? A latent threshold approach to large‐scale mixture innovation models. *Journal of Applied Econometrics*, Vol. 34, pp. 621-640.

Jarocinski, M. (2010) Responses to monetary policy shocks in the east and west of Europe: a comparison. *Journal of Applied Econometrics*, Vol. 25, pp. 833-868.

Kilian L., and H. Lütkepohl (2017) Structural Vector Autoregressive Analysis. *Cambridge University Press*.

Mertens K., and M. O. Ravn (2013) The dynamic effects of personal and corporate income tax changes in the United States *The American Economic Review*, Vol. 103(4), pp. 1212–47.

Rubio-Ramirez, J.F., Waggoner, D.F., and T. Zha (2010) Structural vector autoregressions: Theory of identification and
algorithms for inference *The Review of Economic Studies*, Vol. 77(2), pp. 665–696.

Sa, F., Towbin, P. and T. Wieladek (2014) Capital inflows, financial structure and housing booms. *Journal of the European Economic Association*, Vol. 12(2), pp. 522-546.

Sims, C. A. (1980) Macroeconomics and reality. *Econometrica: Journal of the Econometric Society*, Vol 48(1), pp. 1–48.

Towbin, P. and S. Weber (2013) Limits of floating exchange rates: The role of foreign currency debt and import structure. *Journal of Development Economics*, Vol. 101, pp. 179-194.

