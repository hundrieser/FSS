\name{bootstrap.hotelling.two.sample.test}
\alias{bootstrap.hotelling.two.sample.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bootstrap based Two-Sample Hotelling Test
%%  ~~function to do ... ~~
}
\description{
 Performs a two-sample Hotelling test for circular or toroidal based on bootstrapping empirical intrinsic mean. 
}
\usage{
bootstrap.hotelling.two.sample.test(x.sample, y.sample, significance.level = 0.05, bootstrap.rep = 1000)
}

\arguments{
   \item{x.sample}{ A numeric $m \times n$-matrix of angles in radians with entries in $[-\pi,\pi)$. The number of columns $m$ is the dimension of the torus, and the number of rows $n$ is the number of points on the torus. 
}
  \item{y.sample}{ A numeric $m \times n$-matrix of angles in radians with entries in $[-\pi,\pi)$. The number of columns $m$ is the dimension of the torus, and the number of rows $n$ is the number of points on the torus. 
}
  \item{significance.level}{
    A numeric value between 0 and 1. The significance level for the test.
}
 \item{bootstrap.rep}{
    An integer value. The number of bootstrap replications to approximate the null distribution of the test statistic.
}
}

\details{
  The test is based on the bootstrap resampling techniues of the empirical intrinsic mean. The test is valid for large sample sizes and even remains valid if finite sample smeariness is present in either dataset. It should be used when the test for the presence of finite sample smeariness (\code{FSS.test}) rejects.
}

\value{
  \item{test.statistic}{The value of the test statistic.}
  \item{testresult}{If equal to 0, then the null hypothesis is not rejected for the prescribed significance level. If equal to 1, then the null hypothesis is rejected.}
  \item{p.value}{The p-value of the test. If the p-value is smaller than the significance level, then the null hypothesis is rejected.}
}

\references{
  "Finite Sample Smeariness of Frechet Means with Application to Climate" (2020) by S.Hundrieser, B. Eltzner, and S.F. Huckemann.
}

\author{
S. Hundrieser
}

\seealso{See also \code{\link{clt.hotelling.two.sample.test.Rd}} for a similar test based on the central limit theorem, which is computationally more efficient but is not as reliable in the presenco of finite sample smeariness.}

\examples{
set.seed(3)
  x.sample1 = runif(30,-pi/2,pi/2)        # Not affected by finite sample smeariness
  y.sample1 = runif(30,-pi/2,pi/2)        # Not affected by finite sample smeariness
  x.sample2 = runif(30,-3*pi/4,3*pi/4)    # Possibly affected by finite sample smeariness
  y.sample2 = runif(30,-3*pi/4,3*pi/4)    # Possibly affected by finite sample smeariness
  x.sample3 = runif(30,-pi+0.1,pi - 0.1)  # Likely affected by finite sample smeariness
  y.sample3 = runif(30,-pi+0.1,pi - 0.1)  # Likely affected by finite sample smeariness

  bootstrap.hotelling.two.sample.test(x.sample = x.sample1, y.sample = y.sample1, significance.level = 0.05)
  # test.statistic = 0.7933288
  # test.result = 0
  # p.value = 0.358

  bootstrap.hotelling.two.sample.test(x.sample = x.sample2, y.sample = y.sample2, significance.level = 0.05)
  # test.statistic = 1.691997
  # test.result = 0
  # p.value = 0.12

  bootstrap.hotelling.two.sample.test(x.sample = x.sample3, y.sample = y.sample3, significance.level = 0.01)
  # test.statistic = 1.582631
  # test.result = 0
  # p.value = 0.278

  # The p-value of the bootstrap-based test is for x.sample2 vs. 
  # y.sample2  x.sample3 vs y.sample3 larger compared to that of the 
  # CLT based test. This is a consequence of the presence of finite 
  # sample smeariness in the data. The first datasets x.sample1 and y.
  # sample3 are not affected by finite sample smeariness and the 
  # p-value is similar to that of the CLT based test.
}

