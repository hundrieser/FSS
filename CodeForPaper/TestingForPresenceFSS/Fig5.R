##### Summary: Generates Figure 5 of manuscript, to illustrate the performance of the test for the presence of finite sample smeariness in case of von Mises mixtures with disks of multiple radii cut out at around the antipode of the mean.


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
N.vector <- unique(ceiling(10^(seq(0.1, 4, 0.2))))

# Number of bootstrap resamples
B <- 100

# Significance level for the test
alpha <- 0.05

# Number of repetitions for the simulation
R <- 100

# Initialize an empty string for INDEX
INDEX <- ""

# Set the parameter for the von Mises distribution to three
upperpar <- 3

# List of parameter sets for the von Mises distribution
parameter.list <- list(
  # Each set includes the upper parameter, a concentration parameter, and a location parameter
  c(upperpar, 1, 0),
  c(upperpar, 1 / 2, 0),
  c(upperpar, 1 / 2, 1 / 2),
  # The last set includes a location parameter that was computed using the find.Lower.Parameter.Smeary.VMM function
  c(upperpar, 1 / 2, 0.8683279)
)

# Create a vector of distribution names, all set to "VonMises"
distr <- rep("VonMises", length(parameter.list))

# Create a vector of hole values
hole.vec <- c(
  0,
  0.1,
  0.2,
  pi / 2
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
  function(distr = c("Uniform", "BimodalNormal", "VonMises"),
           parameter,
           n,
           hole = HOLE) {
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
      # If the third parameter is not -1, generate the sample from a VMM distribution with a FSS.
      else {
        sample <- generate.Sample.FSS.VMM(
          n = n,
          upperParameter = parameter[1],
          alpha = parameter[2],
          lowerParameter = parameter[3]
        )
      }
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

### Calculate Variance

# Rotate data so that Frechet sample mean is at zero
rotateData <- function(x) {
  if (!is.vector(x)) {
    stop("Error: x must be a vector.")
  }

  sample.mean <- fre.Mean(x)
  rotated.x <- modulo(x - sample.mean)

  return(rotated.x)
}

# Calculate Variance of a sample
calculateVariance <- function(x) {
  n <- length(x)
  return(var(rotateData(x)) * (n - 1) / n)
}

### Bootstrap resampling function

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

### Test for the presence of finite sample smeariness

# Test for presence finite sample smeariness, i.e, if the variance modulation is significantly larger than one (which would amount to no presence of finite sample smeariness). The test takes three arguments: the sample x, the number of bootstrap repetitions b, and the significance level alpha.
FSS.Test <- function(x, b, alpha = 0.05) {

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

  # If BVM - 1 is greater than h, return 1; otherwise, return 0
  # This is the result of the Finite Sample Smeariness (FSS) test
  if (BVM - 1 > h) {
    return(1)
  } else {
    return(0)
  }
}


### Simulation of performance of test for presence of FSS

# Initialize a list to hold the results of the simulation
FSS_testResults <- vector(mode = "list", length = length(parameter.list) * length(hole.vec) * length(N.vector))
titleRData <- "FSS_Test_Result"

# Initialize an empty data frame with the following columns: n, rej_prob, hole, par_index, par1, par2, par3
FSS_testResults_dataframe <- data.frame(
  n = numeric(),
  rej_prob = numeric(),
  hole = factor(),
  par_index = numeric(),
  par1 = numeric(),
  par2 = numeric(),
  par3 = numeric()
)



# Initialize an index counter
index <- 1

# For each value in N.vector,
for (i in 1:length(N.vector)) {

  # For each value in hole.vec,
  for (j in 1:length(hole.vec)) {

    # For each value in parameter.list,
    for (k in 1:length(parameter.list)) {

      # Print a message indicating the start of a simulation with the current parameters
      cat("Staring simulation for N = ", N.vector[i], ", Radius = ", hole.vec[j], ", Parameters =", parameter.list[[k]], "\n")

      # Initialize a variable to hold the test result
      testResult <- 0

      # For each iteration up to R,
      for (r in 1:R) {

        # Add the result of the FSS test on a sample generated with the current parameters to testResult
        testResult <- testResult + FSS.Test(generate.Sample(distr[k], parameter.list[[k]], N.vector[i], hole = hole.vec[j]), b = B, alpha = alpha)

        # If r is a multiple of 20,
        if (r %% 20 == 0) {

          # Print a message indicating the progress of the simulation
          cat("   Finished: ", r, " - ", R, " - Current ratio: ", testResult / r, "\n")
        }
      }

      # Store the test result, normalized by R, in FSS_testResults at the current index
      # Also store the current parameters and the current value of hole.vec[j]
      FSS_testResults[[index]] <- list(n = as.integer(N.vector[i]), testResult = testResult / R, index = as.character(index), parameter = c(parameter.list[[k]]), r = as.character(hole.vec[j]))

      FSS_testResults_dataframe <- rbind(FSS_testResults_dataframe, data.frame(n = as.integer(N.vector[i]), rej_prob = testResult / R, hole = as.integer(j), par_index = as.integer(k), par1 = as.numeric(parameter.list[[k]][1]), par2 = as.numeric(parameter.list[[k]][2]), par3 = as.numeric(parameter.list[[k]][3])))

      # Uncomment the following line to save the workspace to a .RData file after each simulation
      # save.image(paste(titleRData,"par",k,"_r",100*hole.vec[j],"_n",N.vector[i],".RData",sep=""))

      # Increment the index counter
      index <- index + 1

      # Print a message indicating the end of a simulation with the current parameters
      cat("Finished simulation for N = ", N.vector[i], ", Radius = ", hole.vec[j], ", Parameters =", parameter.list[[k]], "\n")
    }

    # Uncomment the following line to save the workspace to a .RData file after each set of simulations with a given value of hole.vec[j]
    # save.image(paste(titleRData,"_r",100*hole.vec[j],"_n",N.vector[i],".RData",sep=""))
  }

  # Uncomment the following line to save the workspace to a .RData file after each set of simulations with a given value of N.vector[i]
  # save.image(paste(titleRData,"_n",N.vector[i],".RData",sep=""))
}
print("Finished")
# Save the workspace to a .RData file
save.image(paste(titleRData, ".RData", sep = ""))


### Generate plots

# Define a function to extract the legend from a ggplot object
get_legend <- function(myggplot) {
  # Convert the ggplot object to a gtable object
  tmp <- ggplot_gtable(ggplot_build(myggplot))

  # Find the index of the legend (named "guide-box") in the list of grobs (graphical objects)
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")

  #  Extract the legend from the list of grobs and return it
  return(tmp$grobs[[leg]])
}


# Initialize a vector to hold the titles for each plot
title.vec <- rep(0, length(parameter.list))

# For each set of parameters in parameter.list,
for (i in 1:length(parameter.list)) {

  # Create a title for the plot using the parameters
  # The title is formatted using TeX syntax for mathematical expressions
  # The parameters are inserted into the expression using paste()
  # The expression represents a probability density function with the parameters as subscripts
  # The round() function is used to limit the third parameter to 4 decimal places
  title.vec[i] <- TeX(paste("   $P^{(", parameter.list[[i]][1], ",",
    parameter.list[[i]][2], ",",
    round(parameter.list[[i]][3], 4), ",r)}_{vMm}$",
    sep = ""
  ))
}

# Initialize a list named 'testResult.plot.list' with the same length as 'parameter.list'.
# This list will be used to store the ggplot objects for each set of parameters in 'parameter.list'.
testResult.plot.list <- vector(mode = "list", length = length(parameter.list))

# Loop over the length of the parameter list
for (j in 1:length(parameter.list)) {
  # For each iteration, create a ggplot object and store it in the corresponding index of 'testResult.plot.list'
  testResult.plot.list[[j]] <-
    # Initialize the ggplot with data filtered by 'par_index' equal to the current iteration index
    ggplot(data = FSS_testResults_dataframe[FSS_testResults_dataframe$par_index == j, ], aes(x = n, y = rej_prob, group = hole)) +
    # Add a line geom, with color and group aesthetics mapped to the factor of 'hole'
    geom_line(aes(col = as.factor(hole), group = as.factor(hole))) +
    # Set the title of the plot to the corresponding element in 'title.vec'
    ggtitle(title.vec[j]) +
    # Set the labels of the x and y axes
    labs(y = "Rejection probability", x = TeX("Sample size n")) +
    # Remove minor grid lines
    theme(panel.grid.minor = element_blank()) +
    # Set the limits of the y-axis
    ylim(-0.05, 1) +
    # Manually set the colors for the 'hole' variable
    scale_color_manual(
      name = "r",
      labels = c(
        "0",
        "0.1",
        "0.2",
        expression(pi / 2)
      ),
      values = c(
        "4" = "grey70",
        "3" = "grey60",
        "2" = "grey40",
        "1" = "grey0"
      )
    ) +
    # Manually set the line types
    scale_linetype_manual(values = c("solid", "83")) +
    # Add a horizontal line at y = 0.05
    geom_hline(yintercept = 0.05, linetype = "dotted", col = "black") +
    # Use a black and white theme
    theme_bw() +
    # Set the position of the legend to the bottom and align the title to the left
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0)) +
    # Set the x-axis to a log10 scale and format the breaks and labels accordingly
    scale_x_log10("Sample size n",
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    # Remove minor grid lines (again)
    theme(panel.grid.minor = element_blank()) +
    # Add log ticks on the bottom side of the plot
    annotation_logticks(sides = "b")
}

# Begin outputting plot to a PDF file named "Fig1.pdf" with specified height and width
pdf("Fig5.pdf", height = 3, width = 9 * 3 / 4)
# Use the 'grid.arrange' function from the 'gridExtra' package to arrange multiple plots in a grid
grid.arrange(
  # Create a grid of plots using 'arrangeGrob'
  arrangeGrob(
    # Add the first, third, and fourth plots from 'testResult.plot.list' to the grid
    # The 'theme' function is used to remove the legend from each plot
    testResult.plot.list[[1]] + theme(legend.position = "none"),
    # testResult.plot.list[[2]]+ theme(legend.position = "none"),
    testResult.plot.list[[3]] + theme(legend.position = "none"),
    testResult.plot.list[[4]] + theme(legend.position = "none"),
    # Arrange the plots in a single row
    nrow = 1
  ),
  # Add the legend from the first plot in 'testResult.plot.list' to the grid
  get_legend(testResult.plot.list[[1]]),
  # Arrange the grid and the legend in two rows
  nrow = 2,
  # Specify the heights of the rows
  heights = c(3.7, 0.75)
)

# Close the PDF device, finalizing the plot and writing it to the file
dev.off()


# Save the current R session to a file named "Fig1.RData"
save.image("Fig5.RData")
