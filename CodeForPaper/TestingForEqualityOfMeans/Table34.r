##### Summary: Generates Tables 3 and 4 of manuscript, depicting the variance modulation for  von Mises mixtures with disks of multiple radii cut out at around the antipode of the mean.


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

# Estimate order of smearyness for finite sample smearyness


### Parameters for initialisation of simulation
# Vector of sample sizes
N.vector <- c(30, 100, 300)

# Number of empirical Frechet means generated for approximation of the numerator in the variance modulation
R <- 100000


# Create a vector of hole values
hole.vec <- c(
    0,
    0.1
    # , 0.2,
    # pi / 2
)



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

# Generate a sample from a specified distribution.
# The function takes four arguments: the distribution type (distr), the parameters of the distribution (parameter), the sample size (n), and the hole size (hole).
generate.Sample <-
    function(distr = c("Uniform", "BimodalNormal", "VonMises", "BvMm"),
             parameter,
             n,
             hole) {
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


# Computes the variance modulation for given sample size 'n' from a specified distribution, where 'distr' is the distribution type, 'parameters' are the parameters for the distribution, and 'hole' is the hole size that is cut out near the antipode of the mean.
computeVarianceModulation <- function(n, distr, parameters, hole) {
    # Set the seed for reproducibility
    set.seed(1)

    # Generate a sample of size 'R' from the specified distribution
    sample <- generate.Sample(distr = distr, parameter = parameters, n = R, hole = hole)

    # Approximate the numerator of the variance modulation by calculating the mean of the squares of the sample
    var <- mean((sample)^2)

    # Initialize a vector to hold the empirical Frechet means for each iteration
    mu_n <- rep(0, R)

    # For each iteration up to R,
    for (i in 1:R) {
        # Generate a sample of size 'n' from the specified distribution
        # Compute the Frechet mean of the sample and store it in 'mu_n'
        mu_n[i] <- fre.Mean(generate.Sample(distr = distr, parameter = parameters, n = n, hole = hole))
    }

    # Calculate the variance of the Frechet means
    var_n <- mean((mu_n)^2)

    # Compute the variance modulation by multiplying 'n' with the ratio of 'var_n' to 'var'
    # Return the variance modulation
    return(n * var_n / var)
}


### Simulation for Table 3

# List of parameter sets for the von Mises distribution
parameter.list <- list(
    # Each set includes the upper parameter, a concentration parameter, and a location parameter
    c(3, 1, 0),
    c(3, 1 / 2, 0),
    c(3, 1 / 2, 1 / 2)
)

# Create a vector of distribution names, all set to "VonMises"
distr <- rep("VonMises", length(parameter.list))

# Initialize an array to hold the variance modulation values
variance_modulation_table3 <- array(0, dim = c(length(N.vector), length(parameter.list), length(hole.vec)))

# For each combination of sample size 'n' and parameter,
for (i in 1:length(N.vector)) {
    for (j in 1:length(parameter.list)) {
        for (k in 1:length(hole.vec)) {
            # Compute the variance modulation for the given sample size, parameter, and hole size
            variance_modulation_table3[i, j, k] <- computeVarianceModulation(n = N.vector[i], distr = distr[j], parameters = parameter.list[[j]], hole = hole.vec[k])
        }
    }
}


### Simulation for Table 4

# List of parameter sets for the von Mises distribution
parameter.list <- list(
    # Each set includes the upper parameter, a concentration parameter, and a location parameter
    c(3, 1),
    c(3, 1.25),
    c(3, 1.5)
)

# Create a vector of distribution names, all set to "VonMises"
distr <- rep("BvMm", length(parameter.list))

# Initialize an array to hold the variance modulation values
variance_modulation_table4 <- array(0, dim = c(length(N.vector), length(parameter.list), length(hole.vec)))

# For each combination of sample size 'n' and parameter,
for (i in 1:length(N.vector)) {
    for (j in 1:length(parameter.list)) {
        for (k in 1:length(hole.vec)) {
            # Compute the variance modulation for the given sample size, parameter, and hole size
            variance_modulation_table4[i, j, k] <- computeVarianceModulation(n = N.vector[i], distr = distr[j], parameters = parameter.list[[j]], hole = hole.vec[k])
        }
    }
}


### Output of simulation

# Return the variance modulation matrix for table 3
round(variance_modulation_table3,1)
# , , 1
#
#      [,1] [,2] [,3]
# [1,]    1  3.7  7.6
# [2,]    1  4.0 11.8
# [3,]    1  4.1 15.6
#
# , , 2
#
#      [,1] [,2] [,3]
# [1,]    1  2.7  5.5
# [2,]    1  2.1  5.4
# [3,]    1  1.4  3.2

# Return the variance modulation matrix  for table 4
round(variance_modulation_table4,1)
# , , 1
#
# [,1] [,2] [,3]
# [1,]  1.1  1.7 26.1
# [2,]  1.1  1.2 47.6
# [3,]  1.1  1.2 32.3
#
# , , 2
#
# [,1] [,2] [,3]
# [1,]  1.1  1.6 22.7
# [2,]  1.0  1.1 32.0
# [3,]  1.0  1.0 11.7

