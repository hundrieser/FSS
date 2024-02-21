\name{FSS.test}
\alias{FSS.test}
\title{Test for presence of finite sample smeariness
}
\description{
Performs a test for the presence of finite sample smeariness for circular data. 
}

\usage{
FSS.test(x.sample, significance.level = 0.05, bootstrap.rep = 1000)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x.sample}{A numeric vector of angles in radian consisting of $n$ entries in $[-\pi,\pi)$. 
}
  \item{significance.level}{
    A numeric value between 0 and 1. The significance level for the test.
}
\item{bootstrap.rep}{
    An integer value. The number of bootstrap replicates used to estimate the variance of the empirical intrinsic circular mean.}
}

\details{
 The test is based on comparing the variance of bootstrapped empirical intrinsic means with the variance of the data set divided by n. In case this ratio, called the (empirical) variance modulation is close to one, the dataset is not affected by the presence of finite sample smeariness. If the variance modulation is however much larger than one, then the test will likely reject. 

 If the dataset is not univariate but multivariate, then for each component a test for the presence of finite sample smeariness is to be conducted in conjunction with a multiple testing procedure according to Benjamini-Hochberg. This way, it becomes clear which entries are affected by the presence of finite sample smeariness. 
}


\value{
  \item{test.statistic}{The value of the test statistic. This is the (estimated) variance modulation ratio of the empirical intrinsic mean. }
  \item{test.result}{If equal to 0, then the null hypothesis is not rejected. If equal to 1, then the null hypothesis is rejected. }
  \item{p.value}{The p-value of the test. If the p-value is smaller than the significance level, then the null hypothesis is rejected.}
}

\references{
  "Finite Sample Smeariness of Frechet Means with Application to Climate" (2020) by S.Hundrieser, B. Eltzner, and S.F. Huckemann.
}
\author{
S. Hundrieser
}

\seealso{See also \code{\link{clt.hotelling.one.sample.test.Rd}} and \code{\link{clt.hotelling.two.sample.test.Rd}} for tests on equality of intrinsic means which still perform reliable in the presence of finite sample smeariness.}
\examples{
set.seed(3)
  sample1 = runif(30,-pi/2,pi/2)        # Not affected by finite sample smeariness
  sample2 = runif(30,-3*pi/4,3*pi/4)    # Possibly affected by finite sample smeariness
  sample3 = runif(30,-pi+0.1,pi - 0.1)  # Likely affected by finite sample smeariness


  FSS.test(sample1, significance.level = 0.05, bootstrap.rep = 1000)
  # test.sttatistic = 0.888444
  # test.result = 0
  # p.value = 0.9979765
  # This sample is not affected by the presence of finite sample smeariness.

  FSS.test(sample2, significance.level = 0.05, bootstrap.rep = 1000)
  # test.sttatistic = 1.601373
  # test.result = 1
  # p.value = 3.655793e-05
  # This sample is affected by the presence of finite sample smeariness.

  FSS.test(sample3, significance.level = 0.05, bootstrap.rep = 1000)
  # test.sttatistic = 7.442801
  # test.result = 0
  # p.value = 0
  # This sample is affected by the presence of finite sample smeariness.
}
