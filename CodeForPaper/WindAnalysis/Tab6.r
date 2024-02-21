##### Summary: Generates Table 6 of manuscript, depicting the rejection probability of the quantile and bootstrap based test under sample splitting for the daily Frechet means of wind directions for the cities Basel and Goettingen for the year 2000 until 2019


library(lubridate)
library(rstudioapi)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Define the duration in hours for the wind data analysis
DURATION <- c(24)

# Define the number of data points for the wind data analysis
DATAPOINTS <- c(365)

# Generate a sequence of dates from 2000 to 2020 at the start of each year
current.time <- ymd_hms(paste(as.character(seq(2000, 2020, 1)), "-01-01 00:00:00"))

# Significance level for statistical tests
alpha <- 0.05

# Number of bootstrap samples for the bootstrap test
M <- 1000

# Number of repetitions for sample splitting procedure
R <- 500


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


### Testing methods
# Circular/Toroidal Hotelling-test using asymptotic CLT with p.value as output
clt.hotelling.two.sample.test <- function(x.sample, y.sample, alpha = 0.05) {
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

  if (t2 > qchisq(1 - alpha, m)) {
    return(list(test.result = 1, p.value = p.value))
  } else {
    return(list(test.result = 0, p.value = p.value))
  }
}


# Define a function to generate a bootstrap sample
generateBootstrapSample <- function(x, m = M) {
  # If x is a vector
  if (is.vector(x)) {
    # Generate a bootstrap sample by sampling with replacement from x
    bootstrap.sample <- sample(x, size = m, replace = TRUE)
    # Return the bootstrap sample
    return(bootstrap.sample)
  } else {
    # If x is not a vector (i.e., it is a matrix or data frame)
    # Generate a bootstrap sample by sampling with replacement from the rows of x
    bootstrap.sample <- x[sample(seq(1, length(x[, 1])), m, replace = TRUE), ]
    # Return the bootstrap sample
    return(bootstrap.sample)
  }
}

# Define a function to compute the variance of the mean of bootstrap samples
getVarianceOfBootstrapSampleMean <- function(sample, k = M) {
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

# Circular/Toroidal Hotelling-test using Bootstrap methods with p.value as output
# Define a function to perform a bootstrap Hotelling-test with p.value as output
bootstrap.hotelling.two.sample.test <- function(x.sample, y.sample, alpha = 0.05, bootstrap.rep = M) {
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
      fre.Mean.Multivariate(generateBootstrapSample(x.sample, length(x.sample[, 1]))),
      x.mean
    )
    sigma.x <- sigma.x + d.vec.x[, i] %*% t(d.vec.x[, i])

    # Compute the distance and variance for the y sample
    d.vec.y[, i] <- distance.vec(
      fre.Mean.Multivariate(generateBootstrapSample(y.sample, length(y.sample[, 1]))),
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
  c_boundary <- quantile(t2.bootstrap, probs = 1 - alpha)

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


### Wind analysis

# Obtain daily Frechet means for a specific city and time based on the wind data in the repository
get.wind.means <- function(city = c("Basel", "Goettingen"), # Default cities are Basel and Goettingen
                           current.time, # Current time
                           DURATION = 24, # Duration in hours
                           DATAPOINTS = 365) { # Number of data points
  # Define the start date
  start.date <- ymd_hms("1985-01-01 00:00:00")

  # Calculate the difference in time from the start date to the current time
  differenceInTime <- as.numeric(start.date %--% current.time, "hours")

  # Read the wind data from a CSV file
  winddata <- as.matrix(read.table(paste(city, "_Winddata.csv", sep = ""), sep = ","))

  # Get the length of the second column of the wind data
  L <- length(winddata[, 2])

  # Convert the wind direction to radians
  angles <- modulo(as.numeric(winddata[11:L, 3]) * 2 * pi / 360)

  # Create a vector of wind directions for the wind rose plot
  angles.ForWindrosePlot <- as.numeric(winddata[11 + differenceInTime + 0:(DATAPOINTS * DURATION), 3])

  # Initialize a vector to store the mean angles
  mean.angles <- rep(0, DATAPOINTS)

  # Initialize a matrix to store the angles
  angles.matrix <- matrix(NA, nrow = DURATION, ncol = DATAPOINTS)

  # Loop over the data points
  for (k in 1:DATAPOINTS) {
    # Loop over the duration
    for (j in (1:DURATION)) {
      # Store the angle in the matrix
      angles.matrix[j, k] <- angles[(k - 1) * DURATION + j + differenceInTime]
    }
    # Calculate the mean angle for the current data point
    mean.angles[k] <- fre.Mean(angles.matrix[, k])
  }


  # Return the frechet means of the angles
  return(mean.angles)
}

# Perform the testing wind analysis to compare the wind data for different years in the different cities
sample.splitting.wind.data.testing <- function(r = R, cities, current.time, DURATION, DATAPOINTS) {
  # Initialize an array to store the mean angles for each city and time period
  mean.angles <- array(dim = c(DATAPOINTS, length(cities), 20))

  # Initialize a matrix to store the variance modulation for each city and time period
  varMod <- matrix(nrow = length(cities), ncol = 20)

  # Print a message indicating that the wind data is being loaded
  cat("Loading Wind data ...\n")

  # Loop over each city
  for (i in 1:length(cities)) {
    # Loop over each time period
    for (j in 1:20) {
      # Print a message indicating the city and time period for which the daily Frechet means are being computed
      cat(
        "Computing daily Frechet means for ", cities[i], " during ", month(current.time[j]),
        "/", year(current.time[j]),
        " - ", month(current.time[j + 1] - 1),
        "/", year(current.time[j + 1] - 1), " ... "
      )

      # Compute the daily Frechet means for the current city and time period and store them in the array
      mean.angles[, i, j] <- get.wind.means(cities[i], current.time[j], DURATION, DATAPOINTS)

      # Print a message indicating that the daily Frechet means have been computed
      cat(" Finished.\n")
    }
  }

  # Print a message indicating that the wind data has been loaded
  cat(" Finished \n\n")

  # Print a message indicating that the testing is beginning
  cat("Testing between randomly chosen split samples begins\n")

  # Initialize an array to store the test results
  testResults <- array(rep(0, 20 * 2), dim = c(20, 2))

  # Loop over 20 iterates
  for (i in 1:20) {
    set.seed(1)
    # Print a message indicating the time periods for which the tests are being performed
    cat("Performing both tests for ", year(current.time[i]), " : ")

    # Loop over r iterations
    for (k in 1:r) {
      # Randomly sample 365 indices without replacement
      indices <- sample(365, 365, FALSE)
      # Split the mean angles into two halves based on the sampled indices
      mean.angles1 <- mean.angles[indices[1:182], , i]
      mean.angles2 <- mean.angles[indices[183:365], , i]

      # Perform the naive Hotelling test on the two halves and increment the corresponding test result
      output <- clt.hotelling.two.sample.test(mean.angles1, mean.angles2)
      testResults[i, 1] <- testResults[i, 1] + output$test.result

      # Perform the bootstrap Hotelling test on the two halves and increment the corresponding test result
      output <- bootstrap.hotelling.two.sample.test(mean.angles1, mean.angles2, bootstrap.rep = M)
      testResults[i, 2] <- testResults[i, 2] + output$test.result
    }

    # Compute the average test result for each test
    testResults[i, 1] <- testResults[i, 1] / r
    testResults[i, 2] <- testResults[i, 2] / r

    # Print the average test results
    cat("Quantile ", testResults[i, 1], "  -  Bootstrap ", testResults[i, 2], "\n")
  }

  # Return the test results
  return(testResults)
}


# Initialize a list to store the test results for three different settings: Basel and Goettingen, Basel, and Goettingen
output <- vector(mode = "list", length = 1)

### Test consecutive years for Basel and Goettingen if they are significantly different
cities <- c("Basel", "Goettingen")
output[[1]] <- sample.splitting.wind.data.testing(r = R, cities, current.time, DURATION, DATAPOINTS)

# Print the test results
# Define a vector of labels for the different tests
label <- c(
  "Empirical Rejection probabilities - Quantile Test",
  "Empirical Rejection probabilities - Bootstrap Test"
)

# Loop over the indices 1 and 2
for (i in c(1:2)) {
  # Print the label for the current index
  cat(label[i], "\n")
  # Print the variance modulation for the current index
  cat(round(output[[1]][, i], 3), "\n")
}

# Empirical Rejection probabilities - Quantile Test 
# 0.312 0.394 0.198 0.202 0.424 0.37 0.126 0.584 0.198 0.624 0.544 0.148 0.332 0.694 0.222 0.234 0.406 0.374 0.268 0.254 
# Empirical Rejection probabilities - Bootstrap Test 
# 0.042 0.028 0.016 0.028 0.026 0.024 0.026 0.042 0.012 0.04 0.03 0.026 0.03 0.024 0.018 0.028 0.044 0.038 0.016 0.042 