##### Summary: Simulates the test for equality of Frechet means under various distirbutions

# Load the 'circular' package, which provides methods for circular (directional) statistics
library(circular)

# Load the 'ggplot2' package, which provides a system for 'declaratively' creating graphics,
# based on "The Grammar of Graphics"
library(ggplot2)

# Load the 'latex2exp' package, which parses and converts LaTeX math formulas to R's
# plotmath expressions
library(latex2exp)

# Load the 'gridExtra' package, which provides miscellaneous functions for "grid" graphics,
# including the ability to arrange multiple grid-based plots on a page
library(gridExtra)

# Load the 'RColorBrewer' package, which provides color schemes for maps and other graphics
# designed by Cynthia Brewer as described at http://colorbrewer2.org
library(RColorBrewer)

# Set the working directory to the directory of the currently active RStudio document
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Select Simulation setting to recreate data from manuscript
SETTING = 1
# Settings 1, 2, 3, 4 generate the entries of Table 1 in manuscript (with R = 100000 and N = 100, 1000, 1000)
# Settings 11, 12, 13 generate the data for  simulation setting from Figure 6 (one-sample test) in manuscript with r = 0.0, while
# Settings 14, 15, 16 generate the data for  simulation setting from Figure 6 (one-sample test) in manuscript with r = 0.1.
# Settings 21, 22, 23 generate the data for  simulation setting from Figure 7 (two-sample test) in manuscript with r = 0.0, while
# Settings 24, 25, 26 generate the data for  simulation setting from Figure 7 (two-sample test) in manuscript with r = 0.1.
# Settings 31, 32, 33 generate the data for  simulation setting from Figure 8 (two-sample test) in manuscript with r = 0.0, while
# Settings 34, 35, 36 generate the data for  simulation setting from Figure 8 (two-sample test) in manuscript with r = 0.1.

# Sample sizes for samples X and Y
N <- 30

# Number of Bootstrap repetitions
B <- 1000

# Number of repetitions for the simulation of the test
R <- 1000

# Significance level for statistical tests
alpha <- 0.05

# If chosen, save the output of the simulations
SAVE.OUTPUT <- FALSE

# precision with which the x-axis is discretized for the simulation on the power of the test
I = 31


### Preliminary functions

# Custom modulo function to handle 2*pi as period
modulo <- function(x) {
    # Calculate the modulo of x with 2*pi as the divisor
    output <- x %% (2 * pi)

    # If any values in the output are greater than or equal to pi,
    # subtract 2*pi from them to bring them into the range -pi to pi
    output[output >= pi] <- output[output >= pi] - (2 * pi)

    # If any values in the output are less than -pi,
    # add 2*pi to them to bring them into the range -pi to pi
    output[output < -pi] <- output[output < -pi] + (2 * pi)

    # Return the adjusted output
    return(output)
}

# Correction term to calculate a correction term for the intrinsic distance on a circle or torus
corrector <- function(theta, xi) {
    # If theta is positive and xi is between -pi and theta-pi (inclusive),
    if ((theta > 0) && (-pi <= xi) && (xi <= theta - pi)) {
        # Return 1 as the correction term
        return(1)
    }

    # If theta is negative and xi is between pi+theta and pi (inclusive),
    if ((theta < 0) && (pi + theta <= xi) && (xi <= pi)) {
        # Return -1 as the correction term
        return(-1)
    }

    # If neither of the above conditions is met, return 0 as the correction term
    return(0)
}

# Compute distance vector between two vectors on a circle or torus
distance.vec <- function(theta.vec, omega.vec) {
    # Initialize a vector to hold the distances
    distance <- rep(0, length(theta.vec))

    # For each element in the vectors,
    for (i in 1:length(theta.vec)) {
        # Calculate the distance between the corresponding elements in theta.vec and omega.vec
        # Subtract 2*pi times the result of the error function applied to the same elements
        # Store the result in the distance vector
        distance[i] <- theta.vec[i] - omega.vec[i] - 2 * pi * corrector(theta.vec[i], omega.vec[i])
    }

    # Return the distance vector
    return(distance)
}

### Compute Frechet sample mean

# Calculates the local minimum of the Frechet function for a given set of candidate values.
# It takes two arguments: a vector of candidate values and a sample vector.
frechet.Local.Min <- function(candidates, sample) {
    # Get the number of candidates
    n <- length(candidates)
    # Calculate the mean of the sample
    sample_mean <- mean(sample)

    # Initialize a flag to check if it's the first iteration
    FirstTime <- TRUE

    # Loop over each candidate, except the last one
    for (i in 1:(n - 1)) {
        # Check if the candidate is a potential local minimizer in the positive direction
        if (candidates[i] >= 0 && sample[i] < candidates[i] - pi && candidates[i] - pi < sample[i + 1]) {
            # Calculate the Frechet value for this candidate
            frechetValue <- mean((sample - sample_mean)^2) - (2 * pi * i / n)^2 + 4 * pi * (i / n) * (pi + mean(sample[1:i]) - sample_mean)

            # If it's the first iteration, initialize the vectors of Frechet values and local minimizers
            # Otherwise, append the new Frechet value and local minimizer to the existing vectors
            if (FirstTime) {
                frechetValues <- c(frechetValue)
                localMinimisers <- c(candidates[i])
                FirstTime <- FALSE
            } else {
                frechetValues <- c(frechetValues, frechetValue)
                localMinimisers <- c(localMinimisers, candidates[i])
            }
        }

        # Check if the candidate is a potential local minimizer in the negative direction
        if (candidates[i] < 0 && sample[i] < candidates[i] + pi && candidates[i] + pi < sample[i + 1]) {
            # Calculate the Frechet value for this candidate
            frechetValue <- mean((sample - sample_mean)^2) - (2 * pi * (n - i) / n)^2 + 4 * pi * ((n - i) / n) * (pi - mean(sample[seq(i + 1, n)]) + sample_mean)

            # If it's the first iteration, initialize the vectors of Frechet values and local minimizers
            # Otherwise, append the new Frechet value and local minimizer to the existing vectors
            if (FirstTime) {
                frechetValues <- c(frechetValue)
                localMinimisers <- c(candidates[i])
                FirstTime <- FALSE
            } else {
                frechetValues <- c(frechetValues, frechetValue)
                localMinimisers <- c(localMinimisers, candidates[i])
            }
        }
    }

    # Check if the last candidate is a potential local minimizer
    if (candidates[n] >= 0 && candidates[n] - pi < sample[1] || candidates[n] < 0 && candidates[n] + pi > sample[n]) {
        # Calculate the Frechet value for this candidate
        frechetValue <- mean((sample - sample_mean)^2)

        # If it's the first iteration, initialize the vectors of Frechet values and local minimizers
        # Otherwise, append the new Frechet value and local minimizer to the existing vectors
        if (FirstTime) {
            frechetValues <- c(frechetValue)
            localMinimisers <- c(candidates[n])
        } else {
            frechetValues <- c(frechetValues, frechetValue)
            localMinimisers <- c(localMinimisers, candidates[n])
        }
    }

    # Combine the vectors of local minimizers and Frechet values into a matrix
    localMinimisersAndfrechetValues <- cbind(localMinimisers, frechetValues, deparse.level = 0)

    # Return the matrix of local minimizers and Frechet values
    return(localMinimisersAndfrechetValues)
}

# Calculates the Frechet mean of a sample.
# It takes a sample vector as input.
fre.Mean <- function(sample) {
    # Sort the sample in ascending order
    sample <- sample[order(sample)]
    # Get the number of elements in the sample
    n <- length(sample)
    # Calculate the Euclidean mean of the sample
    euclideanMean <- mean(sample)

    # If the sample has only one element or all elements are the same,
    # return the first element as the Frechet mean
    if ((n == 1) || (min(sample) == max(sample))) {
        return(sample[1])
    }

    # Generate a vector of candidate values for the Frechet mean
    # by adding multiples of 2*pi/n to the Euclidean mean and taking modulo 2*pi
    candidates <- modulo(euclideanMean + 2 * pi * (c(1:n)) / n)
    # Calculate the local minimum of the Frechet function for each candidate
    localMinimisersAndFrechetValues <- frechet.Local.Min(candidates, sample)

    # Find the candidate with the smallest Frechet value
    # and return it as the Frechet mean
    intrinsicSampleMeanAndFrechetValue <- as.vector(localMinimisersAndFrechetValues[which.min(localMinimisersAndFrechetValues[, 2]), ])

    return(intrinsicSampleMeanAndFrechetValue[1])
}

# Calculate intrinsic sample means for given data set
fre.Mean.Multivariate <- function(sample) {
    if (is.vector(sample)) {
        intrinsicSampleMean.vec <- fre.Mean(sample)
    } else {
        m <- length(sample[1, ])
        intrinsicSampleMean.vec <- rep(0, m)
        for (i in 1:m) {
            intrinsicSampleMean.vec[i] <- fre.Mean(sample[, i])
        }
    }

    return(intrinsicSampleMean.vec)
}

# Rotate data (that was gathered by bootstrap) so sample mean is at zero
rotateData <- function(sample.vec) {
    if (is.vector(sample.vec)) {
        sample.mean <- fre.Mean(sample.vec)
        rotated.sample.vec <- modulo(sample.vec - sample.mean)
    } else {
        n <- length(sample.vec[, 1])
        m <- length(sample.vec[1, ])
        sample.mean <- t(matrix(rep(fre.Mean.Multivariate(sample.vec), n), nrow = m, ncol = n))
        rotated.sample.vec <- modulo(sample.vec - sample.mean)
    }
    return(rotated.sample.vec)
}

# Calculate Variance of a sample
calculateVariance <- function(sample) {
    return(var(rotateData(sample)))
}

# Generate a bootstrap sample of size m from a given sample x
generate.Bootstrap.Sample <- function(x, m) {
    # If x is a vector,
    if (is.vector(x)) {
        # Generate a bootstrap sample by sampling with replacement from x, with size m
        bootstrap.sample <- sample(x, size = m, replace = TRUE)

        # Return the bootstrap sample
        return(bootstrap.sample)
    }

    # If x is not a vector (i.e., it's a matrix or dataframe),
    else {
        # Generate a bootstrap sample by sampling with replacement from the rows of x, with size m
        bootstrap.sample <- x[sample(seq(1, length(x[, 1])), m, replace = TRUE), ]

        # Return the bootstrap sample
        return(bootstrap.sample)
    }
}

# Compute the variance of the mean of bootstrap samples via k repetitions
get.Variance.Bootstrap.Sample.Mean <- function(sample, b = B) {
    # Initialize a vector to store the means of the bootstrap samples
    bootstrapped.sample.means <- rep(0, b)

    # Loop over the number of bootstrap samples
    for (i in 1:b) {
        # Generate a bootstrap sample from the input sample and compute its intrinsic sample mean
        bootstrapped.sample.means[i] <- fre.Mean(generate.Bootstrap.Sample(x = sample, m = length(sample)))
    }

    # Compute and return the variance of the means of the bootstrap samples
    return(calculateVariance(bootstrapped.sample.means))
}



### Testing methods

# Circular/Toroidal one-sample Hotelling-test using asymptotic CLT with p.value as output
clt.hotelling.one.sample.test <- function(x.sample, x.mu, significance.level = 0.95) {
    x.sample <- as.matrix(x.sample)
    x.mu <- as.matrix(x.mu)

    x.mean <- fre.Mean.Multivariate(x.sample)

    nx <- length(x.sample[, 1])
    m <- length(x.sample[1, ])

    sigma.x <- calculateVariance(x.sample)

    distance.means.vec <- distance.vec(x.mean, x.mu)

    t1 <- nx * t(distance.means.vec) %*% solve(sigma.x) %*% (distance.means.vec)

    p.value <- (1 - pchisq(t1, m))

    if (t1 > qchisq(1 - significance.level, m)) {
        return(list(test.result = 1, p.value = p.value))
    } else {
        return(list(test.result = 0, p.value = p.value))
    }
}

# Circular/Toroidal two sample Hotelling-test using asymptotic CLT with p.value as output
clt.hotelling.two.sample.test <- function(x.sample, y.sample, significance.level = 0.05) {
    x.sample <- as.matrix(x.sample)
    y.sample <- as.matrix(y.sample)

    x.mean <- fre.Mean.Multivariate(x.sample)
    y.mean <- fre.Mean.Multivariate(y.sample)

    nx <- length(x.sample[, 1])
    ny <- length(y.sample[, 1])
    m <- length(y.sample[1, ])

    sigma.x <- calculateVariance(x.sample)
    sigma.y <- calculateVariance(y.sample)

    sigma <- (nx + ny) * (1 / (nx) * sigma.x + 1 / (ny) * sigma.y)
    # sigma = ((nx - 1)*sigma.x + (ny - 1)*sigma.y)/(nx + ny - 2)

    distance.means.vec <- distance.vec(x.mean, y.mean)

    t2 <- (nx + ny) * t(distance.means.vec) %*% solve(sigma) %*% (distance.means.vec)

    p.value <- (1 - pchisq(t2, m))

    if (t2 > qchisq(1 - significance.level, m)) {
        return(list(test.result = 1, p.value = p.value))
    } else {
        return(list(test.result = 0, p.value = p.value))
    }
}


# Circular/Toroidal one-sample Hotelling-test using Bootstrap methods with p.value as output
bootstrap.hotelling.one.sample.test <- function(x.sample, x.mu, significance.level = 0.05, bootstrap.rep = B) {
    x.sample <- as.matrix(x.sample)

    nx <- length(x.sample[, 1])
    m <- length(x.sample[1, ])

    x.mean <- fre.Mean.Multivariate(x.sample)

    d.vec.x <- matrix(0, nrow = length(x.sample[1, ]), ncol = bootstrap.rep)

    sigma.x <- matrix(0, nrow = length(x.sample[1, ]), ncol = length(x.sample[1, ]))

    for (i in 1:bootstrap.rep) {
        d.vec.x[, i] <- distance.vec(
            fre.Mean.Multivariate(generate.Bootstrap.Sample(x.sample, length(x.sample[, 1]))),
            x.mean
        )
        sigma.x <- sigma.x + d.vec.x[, i] %*% t(d.vec.x[, i])
    }

    #### How do you generate the statistic?

    sigma.x <- sigma.x / bootstrap.rep

    sigma <- sigma.x
    sigma.inv <- solve(sigma)

    t1.bootstrap <- rep(0, bootstrap.rep)
    for (i in 1:bootstrap.rep) {
        t1.bootstrap[i] <- t(d.vec.x[, i]) %*% sigma.inv %*% (d.vec.x[, i])
    }

    c_boundary <- quantile(t1.bootstrap, probs = 1 - significance.level, type = 2)

    t1 <- t(distance.vec(x.mean, x.mu)) %*% sigma.inv %*% (distance.vec(x.mean, x.mu))

    # Compute the p-value
    p.value <- (1 - ecdf(t1.bootstrap)(t1))

    # If the T^2 statistic is greater than the critical value, return 1 and the p-value; otherwise, return 0 and the p-value
    if (t1 > c_boundary) {
        return(list(test.result = 1, p.value = p.value))
    } else {
        return(list(test.result = 0, p.value = p.value))
    }
}


# Circular/Toroidal two-sample Hotelling-test using Bootstrap methods with p.value as output
bootstrap.hotelling.two.sample.test <- function(x.sample, y.sample, significance.level = 0.05, bootstrap.rep = B) {
    # Convert the samples to matrices
    x.sample <- as.matrix(x.sample)
    y.sample <- as.matrix(y.sample)

    # Compute the number of observations in each sample and the number of variables
    nx <- length(x.sample[, 1])
    ny <- length(y.sample[, 1])
    m <- length(x.sample[1, ])

    # Compute the multivariate intrinsic sample mean for each sample
    x.mean <- fre.Mean.Multivariate(x.sample)
    y.mean <- fre.Mean.Multivariate(y.sample)

    # Initialize matrices to store the distances and variances for each sample
    d.vec.x <- matrix(0, nrow = length(x.sample[1, ]), ncol = bootstrap.rep)
    d.vec.y <- matrix(0, nrow = length(x.sample[1, ]), ncol = bootstrap.rep)
    sigma.x <- matrix(0, nrow = length(x.sample[1, ]), ncol = length(x.sample[1, ]))
    sigma.y <- matrix(0, nrow = length(y.sample[1, ]), ncol = length(y.sample[1, ]))

    # Loop over the number of bootstrap replications
    for (i in 1:bootstrap.rep) {
        # Compute the distance and variance for the x sample
        d.vec.x[, i] <- distance.vec(
            fre.Mean.Multivariate(generate.Bootstrap.Sample(x.sample, length(x.sample[, 1]))),
            x.mean
        )
        sigma.x <- sigma.x + d.vec.x[, i] %*% t(d.vec.x[, i])

        # Compute the distance and variance for the y sample
        d.vec.y[, i] <- distance.vec(
            fre.Mean.Multivariate(generate.Bootstrap.Sample(y.sample, length(y.sample[, 1]))),
            y.mean
        )
        sigma.y <- sigma.y + d.vec.y[, i] %*% t(d.vec.y[, i])
    }

    # Compute the average variance for each sample
    sigma.x <- sigma.x / bootstrap.rep
    sigma.y <- sigma.y / bootstrap.rep

    # Compute the pooled variance and its inverse
    sigma <- sigma.x + sigma.y
    sigma.inv <- solve(sigma)

    # Initialize a vector to store the bootstrap T^2 statistics
    t2.bootstrap <- rep(0, bootstrap.rep)
    for (i in 1:bootstrap.rep) {
        t2.bootstrap[i] <- t(d.vec.x[, i] - d.vec.y[, i]) %*% sigma.inv %*% (d.vec.x[, i] - d.vec.y[, i])
    }

    # Compute the critical value for the T^2 statistic
    c_boundary <- quantile(t2.bootstrap, probs = 1 - significance.level, type = 2)

    # Compute the T^2 statistic for the original samples
    t2 <- t(distance.vec(x.mean, y.mean)) %*% sigma.inv %*% (distance.vec(x.mean, y.mean))

    # Compute the p-value
    p.value <- (1 - ecdf(t2.bootstrap)(t2))

    # If the T^2 statistic is greater than the critical value, return 1 and the p-value; otherwise, return 0 and the p-value
    if (t2 > c_boundary) {
        return(list(test.result = 1, p.value = p.value))
    } else {
        return(list(test.result = 0, p.value = p.value))
    }
}



# Function to apply a specific one- or two-sample test
performTests <- function(x.sample, y.sample, x.mu, test.method = c("clt1", "clt2", "bootstrap1", "bootstrap2"), significance.level) {
    # If the test method is "clt1", perform a one-sample Hotelling's T-squared test using the Central Limit Theorem (CLT)
    if (test.method == "clt1") {
        return(clt.hotelling.one.sample.test(x.sample, x.mu, significance.level)$test.result)
    }
    # If the test method is "clt2", perform a two-sample Hotelling's T-squared test using the Central Limit Theorem (CLT)
    else if (test.method == "clt2") {
        return(clt.hotelling.two.sample.test(x.sample, y.sample, significance.level)$test.result)
    }
    # If the test method is "bootstrap1", perform a one-sample Hotelling's T-squared test using the bootstrap method
    else if (test.method == "bootstrap1") {
        return(bootstrap.hotelling.one.sample.test(x.sample, x.mu, significance.level)$test.result)
    }
    # If the test method is "bootstrap2", perform a two-sample Hotelling's T-squared test using the bootstrap method
    else if (test.method == "bootstrap2") {
        return(bootstrap.hotelling.two.sample.test(x.sample, y.sample, significance.level)$test.result)
    }
}

# Test for presence finite sample smeariness
FSS.Test <- function(x, b = B, alpha = 0.05) {
    # Initialize a vector to hold bootstrap means
    mu_Bootstrap <- rep(0, b)

    # Get the length of the input vector x
    n <- length(x)

    # For each iteration up to b,
    for (i in 1:b) {
        # Generate a bootstrap sample from x and calculate its Frechet mean
        # Store the mean in mu_Bootstrap
        mu_Bootstrap[i] <- fre.Mean(generate.Bootstrap.Sample(x, n))
    }

    # Calculate the Frechet mean of collection mu_Bootstrap
    mu_mu_Bootstrap <- fre.Mean(mu_Bootstrap)

    # Calculate the variance of sample x
    var_n <- calculateVariance(x)

    # Calculate the variance of mu_Bootstrap
    var_n_bootstrap <- calculateVariance(mu_Bootstrap)

    # Calculate the mean of the fourth power of the distance between mu_mu_Bootstrap and each value in mu_Bootstrap
    W <- mean(distance.vec(rep(mu_mu_Bootstrap, b), mu_Bootstrap)^4)

    # Calculate the Bootstrap Variance Modulation (BVM)
    BVM <- n * var_n_bootstrap / var_n

    # Calculate the quantile of the standard normal distribution at 1-alpha
    phi_alpha <- qnorm(1 - alpha, 0, 1)

    # Calculate the threshold value h
    h <- n * phi_alpha / sqrt(b) * sqrt(W - var_n_bootstrap^2) / var_n

    # Calculate the p-value for the test
    p.value <- 1 - pnorm((BVM - 1) * var_n / sqrt(W - var_n_bootstrap^2) * sqrt(b) / n, 0, 1)

    # If BVM - 1 is greater than h, return 1; otherwise, return 0
    # This is the result of the Finite Sample Smeariness (FSS) test
    if (BVM - 1 > h) {
        return(list(testresult = 1, p.value = p.value))
    } else {
        return(list(testresult = 0, p.value = p.value))
    }
}


# Multple testing correction as in Benjamini-Hochberg
mult.test.correction <- function(p.values, alpha = 0.05) {
    # Define the number of tests
    N <- 19
    # Define a correction factor
    c <- 1

    # Initialize vectors to store the p-values and the row and column indices
    upperDiag <- rep(0, N)
    indexRow <- rep(0, N)
    indexCol <- rep(0, N)

    # Loop over the number of tests
    for (i in 1:N) {
        # Store the p-value and the row index
        upperDiag[i] <- p.values[i]
        indexRow[i] <- i
    }

    # Order the p-values in ascending order
    newOrdering <- order(upperDiag)
    # Create a data frame with the ordered p-values, the corresponding row indices, the order indices, and a sequence from 1 to N
    dataframe <- cbind(upperDiag[newOrdering], indexRow[newOrdering], newOrdering, seq(1, N))
    # Compute the new p-value threshold
    NewPValue <- dataframe[max(dataframe[, 4] * (dataframe[, 1] <= alpha * dataframe[, 3] / (c * N))), 1]

    # Initialize a vector to store the output
    output <- rep(0, N)
    # Loop over the number of tests
    for (k in 1:N) {
        # If the p-value is less than the new p-value threshold, set the corresponding element of the output to 1
        if (dataframe[k, 1] < NewPValue) {
            output[dataframe[k, 2]] <- 1
        }
    }

    # Return the output
    return(output)
}



### Sample generating mechanisms

# Calculates the density of a von Mises mixture consisting of two components.
# The components are located at 0 (weighted with alpha) and at -pi (weighted with (1-alpha)).
# The respective concentration parameters are given by upperParameter and lowerParameter.
vmm.Antipodal.Density <- function(upperParameter, lowerParameter, alpha, level) {
    # Calculate the density of the first component of the mixture
    # This is done by multiplying the weight (alpha) by the exponential of the negative concentration parameter (upperParameter)
    # and dividing by the modified Bessel function of the first kind of order 0 evaluated at the concentration parameter
    upperDensity <- alpha * exp(-upperParameter) / besselI(upperParameter, 0)

    # Calculate the density of the second component of the mixture
    # This is done in a similar way to the first component, but using the weight (1 - alpha) and the concentration parameter (lowerParameter)
    lowerDensity <- (1 - alpha) * exp(lowerParameter) / besselI(lowerParameter, 0)

    # Subtract the level from the sum of the densities of the two components
    output <- upperDensity + lowerDensity - level

    # Normalize the output by dividing by 2*pi and return it
    return(output / (2 * pi))
}

# Finds the lower parameter for a von Mises mixture (VMM) distribution to be smeary.
# A smeary distribution is one where the empirical Frechet mean admits a slower than parametric convergence rate to the true mean.
# The function takes two arguments: the upper parameter and alpha (the weight of the first component of the mixture).
find.Lower.Parameter.Smeary.VMM <- function(upperParameter, alpha) {
    # Define a function that calculates the density of the VMM distribution for a given lower parameter
    # The density is calculated using the vmm.Antipodal.Density function, with the level set to 1
    f <- function(x) vmm.Antipodal.Density(upperParameter, x, alpha, 1)

    # Use the uniroot function to find the root of the function f in the interval [0, 100]
    # The root is the value of the lower parameter that makes the density of the VMM distribution equal to 1
    # The tolerance for the root finding algorithm is set to 0.000001
    return(uniroot(f, interval = c(0, 100), tol = 0.000001)$root)
}

# Finds the lower parameter for a von Mises mixture (VMM) distribution that admits a density with a prescribed level at the antipode.
# The function takes three arguments: the upper parameter, alpha (the weight of the first component of the mixture), and the prescribed level.
find.Lower.Parameter.FSS.VMM <- function(upperParameter, alpha, level) {
    # If the density of the VMM distribution at the antipode (calculated using the vmm.Antipodal.Density function with the lower parameter set to 0)
    # multiplied by 2*pi is greater than the prescribed level, or if alpha is 1 (meaning the first component has all the weight),
    # return 0 as the lower parameter
    if ((vmm.Antipodal.Density(upperParameter, 0, alpha, 0) * 2 * pi > level) || (alpha == 1)) {
        return(0)
    }

    # Define a function that calculates the density of the VMM distribution for a given lower parameter
    # The density is calculated using the vmm.Antipodal.Density function, with the level set to the prescribed level
    f <- function(x) vmm.Antipodal.Density(upperParameter, x, alpha, level)

    # Use the uniroot function to find the root of the function f in the interval [0, 100]
    # The root is the value of the lower parameter that makes the density of the VMM distribution equal to the prescribed level
    # The tolerance for the root finding algorithm is set to 0.000001
    return(uniroot(f, interval = c(0, 100), tol = 0.000001)$root)
}

# Generate a sample from a smeary von Mises mixture (VMM) distribution.
# The function takes three arguments: the sample size (n), the upper parameter, and alpha (the weight of the first component of the mixture).
generate.Sample.Smeary.VMM <- function(n = N, upperParameter, alpha) {
    # Find the lower parameter for the VMM distribution to be smeary
    lowerParameter <- find.Lower.Parameter.Smeary.VMM(upperParameter, alpha)

    # Initialize a vector to hold the sample
    sample <- rep(0, n)

    # Draw a Bernoulli random variable to determine the number of data points from the first component of the mixture
    # This is done by generating a random number from a binomial distribution with parameters n and alpha
    numberOfUpperVonMises <- rbinom(n = 1, size = n, prob = alpha)

    # Generate the data points from the first component of the mixture
    # This is done by generating random numbers from a von Mises distribution with mean direction 0 and concentration parameter upperParameter
    # The modulo function is used to ensure the data points are in the range [0, 2*pi]
    sample[1:numberOfUpperVonMises] <- modulo(as.numeric(rvonmises(n = numberOfUpperVonMises, mu = circular(0), kappa = upperParameter)))

    # Generate the data points from the second component of the mixture
    # This is done in a similar way to the first component, but using the lower parameter and shifting the mean direction by pi
    sample[(numberOfUpperVonMises + 1):n] <- modulo(as.numeric(rvonmises(n = (n - numberOfUpperVonMises), mu = circular(0), kappa = lowerParameter) + pi))

    # Return the sample
    return(sample)
}

# Generate a sample from a von Mises mixture (VMM) distribution with a FSS.
# The function takes four arguments: the sample size (n), the upper parameter, alpha (the weight of the first component of the mixture), and the lower parameter.
generate.Sample.FSS.VMM <- function(n = N, upperParameter, alpha, lowerParameter) {
    # Draw a Bernoulli random variable to determine the number of data points from the first component of the mixture
    # This is done by generating a random number from a binomial distribution with parameters n and alpha
    numberOfUpperVonMises <- rbinom(n = 1, size = n, prob = alpha)

    # Initialize a vector to hold the sample
    sample <- rep(0, n)

    # If all data points are from the first component of the mixture,
    # generate the data points by generating random numbers from a von Mises distribution with mean direction 0 and concentration parameter upperParameter
    # The modulo function is used to ensure the data points are in the range [0, 2*pi]
    if (numberOfUpperVonMises == n) {
        sample <- modulo(as.numeric(rvonmises(n = n, mu = circular(0), kappa = upperParameter)))
    }
    # If all data points are from the second component of the mixture,
    # generate the data points in a similar way to the first component, but using the lower parameter and shifting the mean direction by pi
    else if (numberOfUpperVonMises == 0) {
        sample <- modulo(as.numeric(rvonmises(n = n, mu = circular(0), kappa = lowerParameter)) + pi)
    }
    # If there are data points from both components of the mixture,
    # generate the data points from each component separately and combine them into the sample
    else {
        sample[1:numberOfUpperVonMises] <- modulo(as.numeric(rvonmises(n = numberOfUpperVonMises, mu = circular(0), kappa = upperParameter)))
        sample[(numberOfUpperVonMises + 1):n] <- modulo(as.numeric(rvonmises(n = (n - numberOfUpperVonMises), mu = circular(0), kappa = lowerParameter)) + pi)
    }

    # Return the sample
    return(sample)
}

# Generate a sample from a specific distribution
generateSample <- function(distr = c("Uniform", "BimodalNormal", "VonMises", "BvMm"),
                           parameter,
                           n,
                           hole =0) {
    # Initialize a vector to hold the sample
    sample <- rep(0, n)

    # If the distribution is Uniform,
    if (distr == "Uniform") {
        # Generate the sample by generating random numbers from a uniform distribution with minimum and maximum given by the parameters
        # The modulo function is used to ensure the data points are in the range [0, 2*pi]
        sample <- modulo(runif(n, parameter[1], parameter[2]))
        return(sample)
    }
    # If the distribution is BimodalNormal,
    else if (distr == "BimodalNormal") {
        # Draw a Bernoulli random variable to determine the number of data points from each mode
        # This is done by generating a random number from a binomial distribution with parameters n and the ratio given by the parameters
        updown <- rbinom(n, 1, parameter[5])
        for (i in 1:n) {
            # If the Bernoulli random variable is 1, generate the data point from the first mode
            # This is done by generating a random number from a normal distribution with mean and standard deviation given by the parameters
            if (updown[i] == 1) {
                sample[i] <- modulo(rnorm(1, parameter[1], parameter[2]))
            }
            # If the Bernoulli random variable is 0, generate the data point from the second mode in a similar way
            else {
                sample[i] <- modulo(rnorm(1, parameter[3], parameter[4]))
            }
        }
    }
    # If the distribution is VonMises,
    else if (distr == "VonMises") {
        # If the third parameter is -1, generate the sample from a smeary VMM distribution
        if (parameter[3] == -1) {
            sample <-
                generate.Sample.Smeary.VMM(
                    n = n,
                    upperParameter = parameter[1],
                    alpha = parameter[2]
                )
        }
        # If the third parameter is not -1, generate the sample from a VMM distribution with FSS
        else {
            sample <- generate.Sample.FSS.VMM(
                n = n,
                upperParameter = parameter[1],
                alpha = parameter[2],
                lowerParameter = parameter[3]
            )
        }
    }
    # If the distribution is VonMises,
    else if (distr == "BvMm") {
        # Generate the sample from a VMM distribution with FSS.
        sample <- generate.Sample.FSS.VMM(
            n = n,
            upperParameter = parameter[1],
            alpha = 1,
            lowerParameter = 0
        )
        # randomly rotate the sample to the left or to the right
        sample <- modulo(sample + (2 * rbinom(n, 1, 0.5) - 1) * parameter[2])
    }

    # For each data point in the sample,
    for (i in 1:n) {
        # If the absolute value of the data point is greater than pi minus the hole size,
        # shift the data point by pi
        # This is done to ensure the data points are in the range [-pi, pi]
        if (abs(sample[i]) > pi - hole) {
            sample[i] <- modulo(sample[i] + pi)
        }
    }

    # Return the sample
    return(sample)
}

# Find the rejection probability based on r repetitions for a given test method applied on a sample from a specific distribution, where one distribution is rotated by rotation
determineRejectionProbability <- function(x.distr = c("Uniform", "BimodalNormal", "VonMises", "BvMm"), y.distr,
                                          x.parameter, y.parameter,
                                          x.mu, hole = 0,
                                          n.x, n.y, r = R, rotation, significance.level,
                                          test.method) {
    # Initialize the rejection frequency
    rejection.freq <- 0

    # Initialize the x and y samples
    x.sample <- rep(NA, n.x)
    y.sample <- rep(NA, n.y)

    # Loop over the number of repetitions
    for (i in 1:r) {
        # Generate the x sample
        x.sample <- generateSample(x.distr, x.parameter, n.x, hole)
        # If the test method is "clt1" or "bootstrap1", perform a one-sample test
        if (test.method == "clt1" || test.method == "bootstrap1") {
            rejection.freq <- rejection.freq + performTests(x.sample, y.sample, x.mu + rotation, test.method, significance.level)
        }
        # Otherwise, generate the y sample and perform a two-sample test
        else {
            y.sample <- modulo(generateSample(y.distr, y.parameter, n.y, hole) + rotation)
            rejection.freq <- rejection.freq + performTests(x.sample, y.sample, x.mu + rotation, test.method, significance.level)
        }
    }

    # Return the rejection probability
    return(rejection.freq / r)
}

# Generate a rejection probability curve for a given test method by applying the test on one or two samples of some specific distribution, where the first sample is rotated by multiple rotations
generateProbabilityCurve <- function(x.distr = c("Uniform", "BimodalNormal", "VonMises", "BvMm"), y.distr,
                                     x.parameter, y.parameter, x.mu = 0, hole, n.x, n.y, rotation, significance.level,
                                     test.method) {
    output <- rep(NA, length(rotation))
    for (i in 1:length(rotation)) {
        cat("(", i, "/", length(rotation), ")  - Simulation started ")
        output[i] <- determineRejectionProbability(x.distr, y.distr, x.parameter, y.parameter, x.mu, hole, n.x, n.y, r = R, rotation[i], significance.level, test.method)
        cat("- Simulation for finished for rotation = ", round(rotation[i], 2), ", Probability of rejection: ", output[i], "\n")
    }
    return(output)
}


# Generate a rejection probability curve by rotating the one distribtion for multiple test methods
generateProbabilityCurveForAllTests <- function(x.distr = c("Uniform", "BimodalNormal", "VonMises", "BvMm"), y.distr,
                                                x.parameter, y.parameter, hole = 0, n.x = 30, n.y = 30, intervals, significance.level = 0.05, performTests = c(1, 0, 0, 0)) {
    # The
    if (intervals == 0) {
        rotation <- c(0)
    } else if (intervals >= 1) {
        rotation <- pi * seq(-intervals , intervals , 1) / intervals
    } else {
        stop("intervals must be a positive integer")
    }

    clt1.output <- rep(NA, length(rotation))
    clt2.output <- rep(NA, length(rotation))
    bootstrap1.output <- rep(NA, length(rotation))
    bootstrap2.output <- rep(NA, length(rotation))

    if (sum(abs(performTests)) == 0) {
        stop("At least one test method must be selected")
    }

    if (performTests[1] == 1) {
        cat("Calculating rejection probabilities for one-sample CLT based Hotelling test with data x ~", x.distr, "\n")
        clt1.output <- generateProbabilityCurve(x.distr, y.distr, x.parameter, y.parameter,
            x.mu = 0, hole, n.x, n.y, rotation, significance.level,
            test.method = "clt1"
        )
        if(SAVE.OUTPUT){ write.table(clt1.output, paste("TESTSRESULTS_clt1_htest_x", x.distr,
            "_x.par", paste(x.parameter, collapse = "_"),
            "_n.x=", n.x, "_hole",round(10*hole), ".txt",
            sep = "", collapse = NULL
        ), row.names = FALSE, col.names = FALSE)}
    }

    if (performTests[2] == 1) {
        cat("Calculating rejection probabilities for two-sample CLT based Hotelling test with data x ~", x.distr, ", y ~", y.distr, "\n")
        clt2.output <- generateProbabilityCurve(x.distr, y.distr, x.parameter, y.parameter,
            x.mu = 0, hole, n.x, n.y, rotation, significance.level,
            test.method = "clt2"
        )
        if(SAVE.OUTPUT){write.table(clt2.output, paste("TESTSRESULTS_clt2_htest_x", x.distr,
            "_x.par", paste(x.parameter, collapse = "_"),
            "_n.x=", n.x, "_y", y.distr,
            "_y.par", paste(y.parameter, collapse = "_"),
            "_n.y=", n.y, "_hole",round(10*hole), ".txt",
            sep = "", collapse = NULL
        ), row.names = FALSE, col.names = FALSE)}
    }

    if (performTests[3] == 1) {
        cat("Calculating rejection probabilities for one-sample bootstrap based Hotelling test  with data x ~", x.distr, "\n")
        bootstrap1.output <- generateProbabilityCurve(x.distr, y.distr, x.parameter, y.parameter,
            x.mu = 0, hole, n.x, n.y, rotation, significance.level,
            test.method = "bootstrap1"
        )
        if(SAVE.OUTPUT){write.table(bootstrap1.output, paste("TESTSRESULTS_bootstrap1_htest_x", x.distr,
            "_x.par", paste(x.parameter, collapse = "_"),
            "_n.x=", n.x, "_hole",round(10*hole), ".txt",
            sep = "", collapse = NULL
        ), row.names = FALSE, col.names = FALSE)}
    }

    if (performTests[4] == 1) {
        cat("Calculating rejection probabilities for two-sample bootstrap based Hotelling test  with data x ~", x.distr, ",  y ~", y.distr, "\n")
        bootstrap2.output <- generateProbabilityCurve(x.distr, y.distr, x.parameter, y.parameter,
            x.mu = 0, hole, n.x, n.y, rotation, significance.level,
            test.method = "bootstrap2"
        )
        if(SAVE.OUTPUT){write.table(bootstrap2.output, paste("TESTSRESULTS_bootstrap2_htest_x", x.distr,
            "_x.par", paste(x.parameter, collapse = "_"),
            "_n.x=", n.x, "_y", y.distr,
            "_y.par", paste(y.parameter, collapse = "_"),
            "_n.y=", n.y, "_hole",round(10*hole), ".txt",
            sep = "", collapse = NULL
        ), row.names = FALSE, col.names = FALSE)}
    }

    return(list(clt1.output = clt1.output, clt2.output = clt2.output, bootstrap1.output = bootstrap1.output, bootstrap2.output = bootstrap2.output))
}

##### Simulation functions
simulationOnTesting <- function(setting, n = N) {
set.seed(1)

# Setting of Table 1
parameter.list= list(
    list(distr = "VonMises", parameter = c(3, 1/2, 0),   hole = 0),
    list(distr = "VonMises", parameter = c(3, 1/2, 1/4),   hole = 0),
    list(distr = "VonMises", parameter = c(3, 1/2, 1/2),   hole = 0),
    list(distr = "VonMises", parameter = c(3, 1/2, 3/4), hole = 0)
)
# Generate test results for Table 1
# Here R  was chosen as 100 000
if( floor(setting/10) == 0){
    s = setting %% 10
    # interval is chosen as 0 to ensure that we get the rejection probability under the null, i.e., if probability measures are identical
    return(generateProbabilityCurveForAllTests(x.distr = parameter.list[[s]]$distr, y.distr = parameter.list[[s]]$distr, x.parameter = parameter.list[[s]]$parameter, y.parameter = parameter.list[[s]]$parameter, hole = parameter.list[[s]]$hole, n.x = n, n.y = n, intervals = 0, significance.level = alpha, performTests = c(0,1,0,1)))
    print("finished!")
}

# Settings of Figure 6 and 7
parameter.list= list(
    list(distr = "VonMises", parameter = c(3, 1, 0),     hole = 0),
    list(distr = "VonMises", parameter = c(3, 1/2, 0),   hole = 0),
    list(distr = "VonMises", parameter = c(3, 1/2, 1/2), hole = 0),
    list(distr = "VonMises", parameter = c(3, 1, 0),     hole = 0.1),
    list(distr = "VonMises", parameter = c(3, 1/2, 0),   hole = 0.1),
    list(distr = "VonMises", parameter = c(3, 1/2, 1/2), hole = 0.1)
)

# Generate test results for Figure 6
if( floor(setting/10) == 1){
    s = setting %% 10
    return(generateProbabilityCurveForAllTests(x.distr = parameter.list[[s]]$distr, y.distr = parameter.list[[s]]$distr, x.parameter = parameter.list[[s]]$parameter, y.parameter = parameter.list[[s]]$parameter, hole = parameter.list[[s]]$hole, n.x = n, n.y = n, intervals = I, significance.level = alpha, performTests = c(1,0,1,0)))
}

# Generate test results for Figure 7
if( floor(setting/10) == 2){
    s = setting %% 10
    return(generateProbabilityCurveForAllTests(x.distr = parameter.list[[s]]$distr, y.distr = parameter.list[[s]]$distr, x.parameter = parameter.list[[s]]$parameter, y.parameter = parameter.list[[s]]$parameter, hole = parameter.list[[s]]$hole, n.x = n, n.y = n, intervals = I, significance.level = alpha, performTests = c(0,1,0,1)))
}

# Settings of Figure 8
parameter.list= list(
    list(distr = "BvMm", parameter = c(3, 1),   hole = 0),
    list(distr = "BvMm", parameter = c(3, 5/4), hole = 0),
    list(distr = "BvMm", parameter = c(3, 3/2), hole = 0),
    list(distr = "BvMm", parameter = c(3, 1),   hole = 0.1),
    list(distr = "BvMm", parameter = c(3, 5/4), hole = 0.1),
    list(distr = "BvMm", parameter = c(3, 3/2), hole = 0.1)
)

# Generate test results for Figure 8
if( floor(setting/10) == 3){
    s = setting %% 10
    return(generateProbabilityCurveForAllTests(x.distr = parameter.list[[s]]$distr, y.distr = parameter.list[[s]]$distr, x.parameter = parameter.list[[s]]$parameter, y.parameter = parameter.list[[s]]$parameter, hole = parameter.list[[s]]$hole, n.x = n, n.y = n, intervals = I, significance.level = alpha, performTests = c(0,1,0,1)))
}

 stop("Setting not found")
}

# Run the simulation with the specified setting and store the result in 'output'
output = simulationOnTesting(SETTING)

# # Print the output of the simulation
# print(output)

# If the setting corresponds to Table 1 (as determined by the first digit of the setting being 0), stop the execution and print an error message
if( floor(SETTING/10) == 0){
    stop("Choose a different setting to obtain plots")
}

# Define the labels for the tests
test.labels = c(" CLT based\n one-sample\n Hotelling Test", " CLT based\ntwo-sample\n Hotelling Test"," Bootstrap based\n one-sample\n Hotelling Test", " Bootstrap based\n two-sample\n Hotelling Test")

# Define the rotation values as a sequence from -pi to pi with 63 equally spaced intervals
rotation <- pi * seq(-I , I , 1) / I


# Combine the results of the four tests into a single data frame
plot.dataframe = rbind(
    # For each test, create a data frame with the rotation values, the rejection probabilities, and the test label
    data.frame( x = rotation, rejProb = output$clt1.output, test = test.labels[1]),
    data.frame( x = rotation, rejProb = output$clt2.output, test = test.labels[2]),
    data.frame( x = rotation, rejProb = output$bootstrap1.output, test = test.labels[3]),
    data.frame( x = rotation, rejProb = output$bootstrap2.output, test = test.labels[4])
)

# Remove the rows of the data frame where the rejection probability is NA
plot.dataframe <- plot.dataframe[!is.na(plot.dataframe$rejProb), ]

# Start a new PDF device for plotting, with specified width and height
pdf("ComparisonTests.pdf", width = 4, height = 2.5)

# Create a ggplot
ggplot(
    # Specify the data frame to use
    data = plot.dataframe,
    # Map the aesthetics: x to 'x', y to 'rejProb', and group to 'test'
    aes(x = x, y = rejProb, group = test)
) +
    # Set the limits of the y-axis
    ylim(0, 1) +
    # Set the limits of the x-axis
    xlim(-pi, pi) +
    # Add a title to the plot
    ggtitle(paste("Rejection probabilities for Setting ", SETTING, sep="")) +
    # Label the y-axis
    labs(y = "Rejection prob.", x = "x") +
    # Add lines to the plot, with color determined by 'test'
    geom_line(aes(col = test, group = test), size = 0.3) +
    # Remove minor grid lines
    theme(panel.grid.minor = element_blank()) +
    # Manually set the colors for the lines
    scale_colour_manual(values = c("red","blue"), labs(fill = "Test\nmethod"))+
    #   scale_colour_grey(start = 0.6, end = 0, labs(fill = "Test method")) +
    # Manually set the line types
    scale_linetype_manual(values = c("solid", "83")) +
    # Add a horizontal line at 'alpha', with dotted line type and black color
    geom_hline(yintercept = alpha, linetype = "dotted", col = "black") +
    # Use a white background
    theme_bw() +
    # Position the legend at the bottom and align the title to the left
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0)) +
    # Adjust the position of the x-axis text
    theme(axis.text.x = element_text(vjust = 0)) +
    # Set the x-axis ticks and labels
    scale_x_discrete(
        limit = c(-pi, -pi / 2, 0, pi / 2, pi),
        labels = c(
            expression(-pi),
            expression(paste(-pi, "/", 2, sep = "")),
            expression(paste(0, sep = "")),
            expression(paste(pi, "/", 2, sep = "")),
            expression(pi)
        )
    ) +
    # Adjust the position of the x-axis title
    theme(axis.title.x = element_text(vjust = -1))

# Close the PDF device
dev.off()
