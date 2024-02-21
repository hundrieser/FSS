\name{guideline.two.sample.test}
\alias{guideline.two.sample.test}

\title{Guideline for Two-Sample Test of Circular/Toroidal data
}
\description{
  Assesses the presence of finite sample smeariness (FSS) for circular data and performs a suitable two-sample Hotelling tests for circular data - under presence of FSS in either data set the bootstrap based test is employed, otherwise the CLT based test is used.
}
\usage{
guideline.two.sample.test(x.sample, y.sample, significance.level = 0.05, bootstrap.rep = 1000, comments = TRUE)
}

\arguments{
  \item{x.sample}{A numeric vector of angles in radian consisting of $n$ entries in $[-\pi,\pi)$.
}
  \item{y.sample}{A numeric vector of angles in radian consisting of $n$ entries in $[-\pi,\pi)$.
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

\seealso{See also \code{\link{clt.hotelling.two.sample.test.Rd}} for the two-sample Hotelling based on the central limit theorem and \code{\link{bootstrap.hotelling.two.sample.test.Rd}} for the two-sample Hotelling based on bootstrapping.}

\examples{
  set.seed(3)
  x.sample1 = runif(30,-pi/2,pi/2)        # Not affected by finite sample smeariness
  y.sample1 = runif(30,-pi/2,pi/2)        # Not affected by finite sample smeariness
  x.sample2 = runif(30,-3*pi/4,3*pi/4)    # Possibly affected by finite sample smeariness
  y.sample2 = runif(30,-3*pi/4,3*pi/4)    # Possibly affected by finite sample smeariness
  x.sample3 = runif(30,-pi+0.1,pi - 0.1)  # Likely affected by finite sample smeariness
  y.sample3 = runif(30,-pi+0.1,pi - 0.1)  # Likely affected by finite sample smeariness

  guideline.two.sample.test(x.sample = x.sample1, y.sample = y.sample1, significance.level = 0.05, comments = TRUE)
  # Both samples are not affected by the presence of finite sample smeariness.
  # CLT based two-sample Hotelling test is performed.
  # test.statistic = 0.7561354
  # test.result = 0
  # p.value = 0.3845406

  guideline.two.sample.test(x.sample = x.sample2, y.sample = y.sample2, significance.level = 0.05, comments = TRUE)
  # At least one sample is affected by the presence of finite sample smeariness.
  # Bootstrap based two-sample Hotelling test is performed.
  # test.statistic = 1.694562
  # test.result = 0
  # p.value = 0.107

  guideline.two.sample.test(x.sample = x.sample3, y.sample = y.sample3, significance.level = 0.01, comments = TRUE)
  # At least one sample is affected by the presence of finite sample smeariness.
  # Bootstrap based two-sample Hotelling test is performed
  # test.statistic = 1.663219
  # test.result = 0
  # p.value = 0.276
}