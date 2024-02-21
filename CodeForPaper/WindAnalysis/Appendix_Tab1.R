##### Summary: Generates Table 1 of appendix of manuscript, depicting the empirical rejection probability for the test for finite sample smeariness applied on n-out-of-n bootstrap resamples of the collections of daily Frechet means for the cities Basel and Goettingen during the years 2000 until 2019

# Load the 'lubridate' package for date-time manipulation
library(lubridate)

# Load the 'rstudioapi' package for accessing the RStudio API
library(rstudioapi)


# Set the working directory to the directory of the currently active RStudio document
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Define the duration in hours for the wind data analysis
DURATION <- c(24)

# Define the number of data points for the wind data analysis
DATAPOINTS <- c(365)

# Generate a sequence of dates from 2000 to 2020 at the start of each year
current.time <- ymd_hms(paste(as.character(seq(2000, 2019, 1)), "-01-01 00:00:00"))

# Number of Bootstrap repetitions
B <- 1000

# Number of repetitions for the FSS test
R <- 1000

# Significance level for statistical tests
alpha <- 0.05

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

# Define a function to compute the variance of the mean of bootstrap samples
getVarianceOfBootstrapSampleMean <- function(sample, k = B) {
  # Initialize a vector to store the means of the bootstrap samples
  bootstrapped.sample.means <- rep(0, k)

  # Loop over the number of bootstrap samples
  for (i in 1:k) {
    # Generate a bootstrap sample from the input sample and compute its intrinsic sample mean
    bootstrapped.sample.means[i] <- fre.Mean(generateBootstrapSample(x = sample, m = length(sample)))
  }

  # Compute and return the variance of the means of the bootstrap samples
  return(calculateVariance(bootstrapped.sample.means))
}

### Test for the presence of finite sample smeariness

# Test for presence finite sample smeariness, i.e, if the variance modulation is significantly larger than one (which would amount to no presence of finite sample smeariness). The test takes three arguments: the sample x, the number of bootstrap repetitions b, and the significance level alpha. This function also outputs the respective p-value.
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
  pvalue <- 1 - pnorm((BVM - 1) * var_n / sqrt(W - var_n_bootstrap^2) * sqrt(B) / n, 0, 1)

  # If BVM - 1 is greater than h, return 1; otherwise, return 0
  # This is the result of the Finite Sample Smeariness (FSS) test
  if (BVM - 1 > h) {
    return(list(testresult = 1, pvalue = pvalue))
  } else {
    return(list(testresult = 0, pvalue = pvalue))
  }
}



# Obtain daily Frechet means for a specific city and time based on the wind data in the repository
get.wind.means <- function(city = c("Basel", "Goettingen"),
                           current.time,
                           DURATION, DATAPOINTS) {
  # Define the start date
  start.date <- ymd_hms("1985-01-01 00:00:00")

  # Calculate the difference in time between the start date and the current time in hours
  differenceInTime <- as.numeric(start.date %--% current.time, "hours")

  # Read the wind data from a CSV file and convert it to a matrix
  winddata <- as.matrix(read.table(paste(city, "_Winddata.csv", sep = ""), sep = ","))

  # Get the length of the second column of the wind data
  L <- length(winddata[, 2])

  # Convert the third column of the wind data to radians and take the modulo
  angles <- modulo(as.numeric(winddata[11:L, 3]) * 2 * pi / 360)

  # Initialize a vector to store the mean angles
  mean.angles <- rep(0, DATAPOINTS)

  # Initialize a matrix to store the angles
  angles.matrix <- matrix(NA, nrow = DURATION, ncol = DATAPOINTS)

  # Loop over the data points
  for (k in 1:DATAPOINTS) {
    # Loop over the duration
    for (j in (1:DURATION)) {
      # Get the angle for the current duration and data point
      angles.matrix[j, k] <- angles[(k - 1) * DURATION + j + differenceInTime]
    }
    # Compute the mean angle for the current data point
    mean.angles[k] <- fre.Mean(angles.matrix[, k])
  }

  # Return the mean angles
  return(mean.angles)
}


# Define a function to analyse wind data using the Fluctuation Scaling (FSS) test
analyse.wind.data.FSS <- function(r = R, cities, current.time, DURATION, DATAPOINTS) {
  # Set the seed for random number generation to 1 for reproducibility
  set.seed(1)

  # Initialize an array to store the mean angles of the wind data
  mean.angles <- array(dim = c(DATAPOINTS, length(current.time)))

  # Print a message indicating that the wind data is being loaded
  cat("\nLoading wind data for", cities, "...")

  # Loop over the time periods
  for (i in 1:length(current.time)) {
    # Compute the mean angles of the wind data for the current city and time period
    mean.angles[, i] <- get.wind.means(cities, current.time[i], DURATION, DATAPOINTS)
  }

  # Print a message indicating that the wind data has been loaded
  cat(" Finished \n\n")

  # Print a message indicating that the FSS test is being performed
  cat("FSS-Testing for bootstrap samples (with replacement) begins\n")

  # Initialize an array to store the FSS test results
  testresult <- array(rep(0, length(current.time)), dim = c(length(current.time), 1))

  # Loop over the time periods
  for (i in 1:(length(current.time))) {
    # Loop over the bootstrap replications
    for (k in 1:r) {
      # Randomly sample 365 indices with replacement
      indices <- sample(365, 365, TRUE)

      # Compute the mean angles of the bootstrap sample
      bootstrap.mean.angles <- mean.angles[indices, i]

      # Print a message indicating the current replication and time period
      if (k %% 20 == 0) {
        cat(k, " - ", r, " - Performing FSS-test for ", year(current.time[i]), " : ")
      }

      # Perform the FSS test on the bootstrap sample
      output <- FSS.Test(bootstrap.mean.angles, b = B)

      # Update the FSS test result for the current time period
      testresult[i, 1] <- testresult[i, 1] + output$testresult

      # Print the current FSS test result
      if (k %% 20 == 0) {
        cat(testresult[i, 1] / k, "\n")
      }
    }

    # Compute the average FSS test result for the current time period
    testresult[i, 1] <- testresult[i, 1] / r
  }

  # Print the FSS test results
  print(testresult)

  # Return the FSS test results
  return(list(FSS_Tests = testresult[, 1]))
}

cities <- c("Basel", "Goettingen")
testResults.DataFrame <- vector(mode = "list", length = 2)

for (i in 1:2) {
  testResults.DataFrame[[i]] <- analyse.wind.data.FSS(r = R, cities[i], current.time, DURATION, DATAPOINTS)
  cat("Empirical rejection probabilitys for FSS test on bootstrap n-out-of-n resamples for ", cities[i], ":\n")
  cat(testResults.DataFrame[[i]])
}
