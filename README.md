# Finite Sample Smeariness of Fréchet Means on the Circle

This is a complementary package to the work "Finite Sample Smeariness of Frechet Means with Application to Climate" where we develop a framework for more reliable statistical inference with Frechet means on the circle or the torus. In the paper we formalize the concept of finite sample smeariness, which negatively affects the reliability of asymptotic quantile based tests. Under the presence of finite sample smeariness - which can be tested - we recommend the usage of a bootstrap based approach, otherwise the asymptotic quantile based approach should be employed. 


# Package
We provide an R package `FSS` which implements the following functions:
+ `clt.hotelling.one.sample.test` and `clt.hotelling.two.sample.test`: Tests for equality of Fréchet means via an asymptotic quantile based approach. 
+ `bootstrap.hotelling.one.sample.test` and `bootstrap.hotelling.two.sample.test`: Tests for equality of Fréchet means via a bootstrap resampling approach
+ `FSS.test`: Test for the presence of finite sample smeariness. 
+ `guideline.one.sample.test` and `guideline.two.sample.test`: A  procedure which automatically chooses the appropriate test depending on the presence of finite sample smeariness. 

This package continues to be under development, and has only been tested on Mac OS 14.2.1 with R 4.0.3. 

# Installation
This package may be installed as follows, using the `devtools` R package. If you do not have the `devtools`
package installed, you may install it using the command `install.package("devtools")`.
```r
library(devtools)
devtools::install_github("hundrieser/FSS")
```
Alternatively, the file FSS_1.0.0.tar.gz’ can be downloaded and the package can be installed with the command.
```r
install.packages("/PathToFile/FSS_1.0.0.tar.gz", repos = NULL, type = "source")
```

# Code for Paper
In addition to the package we also provide code for every simulation performed in the manuscript (see Folder: `CodeForPaper`). This includes:
+ `VarianceModulationCurves`: The simulations on variance modulation curves in case of different parametric distributions (Figures 1-3 in manuscript)
+ `TestingForPresenceFSS`: The simulations on the test for the presence of finite sample smeariness (Figure 5 in manuscript)
+ `TestingForEqualityOfMeans`: The simulation on the performance of the asymptotic quantile and bootstrap based tests (Table 1 and Figures 6-8 in manuscript). Further we provide the code for the variance variance modulation at difference sample sizes for the considered distributions (Tables 3 and 4). 
+ `WindAnalysis`: The code for the wind data analysis (Figures 9,10, and in Appendix Figures 2,3). The underlying wind data for the cities Goettingen and Basel is not provided in the folder but can be accessed from meteoBlue. 