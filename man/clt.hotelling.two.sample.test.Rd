\name{clt.hotelling.two.sample.test}
\alias{clt.hotelling.two.sample.test}

\title{Asympototic Quantile (CLT) based Two-Sample Hotelling Test
}
\description{
  Performs a two-sample Hotelling test for circular or toroidal based on the Central Limit Theorem (CLT) for the empirical intrinsic mean. 
}

\usage{
clt.hotelling.two.sample.test(x.sample, y.sample, significance.level = 0.05)
}

\arguments{
  \item{x.sample}{ A numeric $m \times n$-matrix of angles in radians with entries in $[-\pi,\pi)$. The number of columns $m$ is the dimension of the torus, and the number of rows $n$ is the number of points on the torus. 
}
  \item{y.sample}{ A numeric $m \times n$-matrix of angles in radians with entries in $[-\pi,\pi)$. The number of columns $m$ is the dimension of the torus, and the number of rows $n$ is the number of points on the torus. 
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

\seealso{See also \code{\link{clt.hotelling.two.sample.test.Rd}} for a similar test based on bootstrap resampling and \seealso{See also \code{\link{clt.hotelling.one.sample.test.Rd}} for the one-sample test CLT based test.}}

\examples{
  set.seed(3)
  x.sample1 = runif(30,-pi/2,pi/2)        # Not affected by finite sample smeariness
  y.sample1 = runif(30,-pi/2,pi/2)        # Not affected by finite sample smeariness
  x.sample2 = runif(30,-3*pi/4,3*pi/4)    # Possibly affected by finite sample smeariness
  y.sample2 = runif(30,-3*pi/4,3*pi/4)    # Possibly affected by finite sample smeariness
  x.sample3 = runif(30,-pi+0.1,pi - 0.1)  # Likely affected by finite sample smeariness
  y.sample3 = runif(30,-pi+0.1,pi - 0.1)  # Likely affected by finite sample smeariness

  clt.hotelling.two.sample.test(x.sample = x.sample1, y.sample = y.sample1, significance.level = 0.05)
  # test.statistic = 0.7561354
  # test.result = 0
  # p.value = 0.3845406

  clt.hotelling.two.sample.test(x.sample = x.sample2, y.sample = y.sample2, significance.level = 0.05)
  # test.statistic = 3.212809
  # test.result = 0
  # p.value = 0.07306396

  clt.hotelling.two.sample.test(x.sample = x.sample3, y.sample = y.sample3, significance.level = 0.01)
  # test.statistic = 24.16841
  # test.result = 1
  # p.value = 8.826825e-07

  # The p-value of the CLT based test is for x.sample2 vs. y.sample2  
  # and x.sample3 vs y.sample3 smaller compared to that of the 
  # bootstrap-based test. This is a consequence of the presence of 
  # finite sample smeariness in the data. The first datasets x.
  # sample1 and y.sample1 are not affected by finite sample 
  # smeariness and the p-value is similar to that of the 
  # bootstrap-based test.
}
