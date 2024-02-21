\name{guideline.one.sample.test}
\alias{guideline.one.sample.test}

\title{Guideline for One-Sample Test of Circular data
}
\description{
  Assesses the presence of finite sample smeariness (FSS) for circular data and performs a suitable one-sample Hotelling tests for circular data - under presence of FSS the bootstrap based test is employed, otherwise the CLT based test is used.
}

\usage{
guideline.one.sample.test(x.sample, x.mu, significance.level = 0.05, bootstrap.rep = 1000, comments = TRUE)
}

\arguments{
  \item{x.sample}{A numeric vector of angles in radian consisting of $n$ entries in $[-\pi,\pi)$. 
}
  \item{x.mu}{ A numeric value in $[-\pi,\pi)$.  This is the element of the circle for which it is tested whether the underlying population mean of the sample is equal to it.
}
  \item{significance.level}{
    A numeric value between 0 and 1. The significance level for the test.
}
  \item{bootstrap.rep}{
    An integer value. The number of bootstrap replications to approximate the null distribution of the test statistic.
}
\item{comments}{
    A Boolean value. If TRUE, the function prints comments on which test is conducted.
}
}

\value{
  \item{test.statistic}{The value of the test statistic considered.}
  \item{testresult}{If equal to 0, then the null hypothesis is not rejected for the prescribed significance level. If equal to 1, then the null hypothesis is rejected.}
  \item{p.value}{The p-value of the test considered. If the p-value is smaller than the significance level, then the null hypothesis is rejected.}
}

\references{
  "Finite Sample Smeariness of Frechet Means with Application to Climate" (2020) by S.Hundrieser, B. Eltzner, and S.F. Huckemann.
}

\author{
S. Hundrieser
}

\seealso{See also \code{\link{clt.hotelling.one.sample.test.Rd}} for the one-sample Hotelling based on the central limit theorem and \code{\link{bootstrap.hotelling.one.sample.test.Rd}} for the one-sample Hotelling based on bootstrapping.}

\examples{
  set.seed(3)
  sample1 = runif(30,-pi/2,pi/2)        # Not affected by finite sample smeariness
  sample2 = runif(30,-3*pi/4,3*pi/4)    # Possibly affected by finite sample smeariness
  sample3 = runif(30,-pi+0.1,pi - 0.1)  # Likely affected by finite sample smeariness

  guideline.one.sample.test(x.sample = sample1, x.mu = 0, significance.level = 0.05, comments = TRUE)
  # Sample is not affected by the presence of finite sample smeariness.
  # CLT based one-sample Hotelling test is performed
  # test.statistic = 0.05851197
  # test.result = 0
  # p.value = 0.8088635

  guideline.one.sample.test(x.sample = sample2, x.mu = 0, significance.level = 0.05, comments = TRUE)
  # Sample is affected by the presence of finite sample smeariness.
  # Bootstrap based two-sample Hotelling test is performed
  # test.statistic = 1.352163
  # test.result = 0
  # p.value = 0.186

  guideline.one.sample.test(x.sample = sample3, x.mu = 0, significance.level = 0.01, comments = TRUE)
  # Sample is affected by the presence of finite sample smeariness.
  # Bootstrap based two-sample Hotelling test is performed
  # test.statistic = 6.304613
  # test.result = 0
  # p.value = 0.031
}