\name{fre.Mean}
\alias{fre.Mean}
\title{Frechet Mean on the Circle
}

\description{
  Compute the Frechet mean of a set of points on the circle. All points are assumed to be in radians with entries in $[-\pi,\pi)$. The function uses the algorithm described in "Intrinsic means on the circle: uniqueness, locus and asymptotics" by T. Hotz and S.F. Huckemann.
}

\usage{
fre.Mean(sample)
}

\arguments{
  \item{sample}{
    A numeric vector of angles in radians with entries in $[-\pi,\pi)$.
  }
}

\value{
  Outputs the intrinsic circular sample mean of the vector parametrized in $[-\pi,\pi)$.
}

\references{
"Intrinsic means on the circle: uniqueness, locus and asymptotics" (2015) by T. Hotz and S.F. Huckemann
}

\author{
  S. Hundrieser
}

\note{
  Computational complexity scales with order O(n log(n)) where n is the number of entries in the input vector.
}

\examples{
  fre.Mean(c(-pi/2, 0, pi/2))       # = 0
  fre.Mean(c(0, pi/2, pi))         # 1.570796 = pi/2
  fre.Mean(c(-pi/2 - 1, pi/2 + 1)) # -3.141593 = -pi
}
