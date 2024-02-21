\name{fre.Mean.Multivariate}
\alias{fre.Mean.Multivariate}
\title{Multivariate Frechet Mean on the Torus
}

\description{
 Compute the Frechet mean of a set of points on the m-dimensional Torus. All points are assumed to be in radians with entries in $[-\pi,\pi)^m$. The function relies on reducing the problem to the circle and then using the fre.mean function, which relies on the  algorithm by described in "Intrinsic means on the circle: uniqueness, locus and asymptotics" by T. Hotz and S.F. Huckemann.
}

\usage{
fre.Mean.Multivariate(sample)
}

\arguments{
  \item{sample}{
    A numeric $m \times n$-matrix of angles in radians with entries in $[-\pi,\pi)$. The number of columns $m$ is the dimension of the torus, and the number of rows $n$ is the number of points on the torus. 
  }
}

\value{
  Outputs the intrinsic toroidal sample mean of the vector parametrized in $[-\pi,\pi)^m$.
}

\references{
"Finite Sample Smeariness of Frechet Means with Application to Climate" (2020) by S.Hundrieser, B. Eltzner, and S.F. Huckemann.
"Intrinsic means on the circle: uniqueness, locus and asymptotics" (2015) by T. Hotz and S.F. Huckemann.
}

\author{
S. Hundrieser
}

\note{
Computational complexity scales with order O(m * n log(n)) where m is the number of columgs (dimension of the data) and n is the number of rows (number of data points) of the input matrix.
}

\examples{
  sample1 = matrix(c(-pi/2, 0, pi/2, 0, 0, 0, pi/2, 0, -pi/2), ncol = 3)
  sample2 = matrix(c(-pi/2, -pi, pi/2, -pi, -pi, -pi, pi/2, -pi, -pi/2), ncol = 3)
  sample3 = matrix(c(-pi/4, pi/4,0, 0, pi/4,-pi/4), ncol = 3)

  fre.Mean.Multivariate(sample1)   # = c(0,0,0)
  fre.Mean.Multivariate(sample2)   # = c(-pi,-pi,-pi)
  fre.Mean.Multivariate(sample3)   # = c(0,0,0)
}
