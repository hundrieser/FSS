\name{clt.hotelling.one.sample.test}
\alias{clt.hotelling.one.sample.test}
\title{Asymptotic Quantile (CLT) based One-Sample Hotelling Test
}

\description{
  Performs a one-sample Hotelling test for circular or toroidal based on the Central Limit Theorem (CLT) for the empirical intrinsic mean. 
}

\usage{
clt.hotelling.one.sample.test(x.sample, x.mu, significance.level = 0.05)
}

\arguments{
  \item{x.sample}{ A numeric $m \times n$-matrix of angles in radians with entries in $[-\pi,\pi)$. The number of columns $m$ is the dimension of the torus, and the number of rows $n$ is the number of points on the torus. 
}
  \item{x.mu}{ A numeric vector with $m$ entries of angles in radians with entries in $[-\pi,\pi)$. This is the element of the torus for which it is tested whether the underlying population mean of the sample is equal to it.
}
  \item{significance.level}{
    A numeric value between 0 and 1. The significance level for the test.
}
}

\details{
    The test is based on the asymptotic distribution of the empirical intrinsic mean. The test is valid for large sample sizes and provided no finite sample smeariness is present in the data. The latter can be tested using the function \code{FSS.test}.
}

\value{
  \item{test.statistic}{The value of the test statistic.}
  \item{testresult}{If equal to 0, then the null hypothesis is not rejected. If equal to 1, then the null hypothesis is rejected. }
  \item{p.value}{The p-value of the test. If the p-value is smaller than the significance level, then the null hypothesis is rejected.}
}

\references{
  "Finite Sample Smeariness of Frechet Means with Application to Climate" (2020) by S.Hundrieser, B. Eltzner, and S.F. Huckemann.
  "Nonparametric Inference on Manifolds" (2012) by A. Bhattacharya and R. Bhattacharya
  "Large sample theory of intrinsic and extrinsic sample means on manifolds II." (2005) by  A. Bhattacharya and V. Patrangenaru.
}

\author{
S. Hundrieser
}

\seealso{See also \code{\link{clt.hotelling.one.sample.test.Rd}} for a similar test based on bootstrap resampling and \seealso{See also \code{\link{clt.hotelling.two.sample.test.Rd}} for the two-sample test CLT based test.}}

\examples{
  set.seed(3)
  sample1 = runif(30,-pi/2,pi/2)        # Not affected by finite sample smeariness
  sample2 = runif(30,-3*pi/4,3*pi/4)    # Possibly affected by finite sample smeariness
  sample3 = runif(30,-pi+0.1,pi - 0.1)  # Likely affected by finite sample smeariness

  clt.hotelling.one.sample.test(x.sample = sample1, x.mu = 0, significance.level = 0.05)
  # test.statistic = 0.05851197
  # test.result = 0
  # p.value = 0.8088635

  clt.hotelling.one.sample.test(x.sample = sample2, x.mu = 0, significance.level = 0.05)
  # test.statistic = 2.085583
  # test.result = 0
  # p.value = 0.1486955

  clt.hotelling.one.sample.test(x.sample = sample3, x.mu = 0, significance.level = 0.01)
  # test.statistic = 41.03194
  # test.result = 1
  # p.value = 1.497616e-10

  # The p-value of the CLT based test is for sample2 and sample3 
  # smaller compared to that of the bootstrap-based test. This is a 
  # consequence of the presence of finite sample smeariness in the 
  # data. The first dataset sample1 is not affected by finite sample 
  # smeariness and the p-value is similar to that of the 
  # bootstrap-based test.
}
