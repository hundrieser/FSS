####### Preliminary functions

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

####### Compute Frechet sample mean

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
get.Variance.Bootstrap.Sample.Mean <- function(sample, b = 1000) {
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



####### Testing methods

# Circular/Toroidal one-sample Hotelling-test using asymptotic CLT with p.value as output
clt.hotelling.one.sample.test <- function(x.sample, x.mu, significance.level = 0.05) {
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
        return(list(test.statistic = t1, test.result = 1, p.value = p.value))
    } else {
        return(list(test.statistic = t1, test.result = 0, p.value = p.value))
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
        return(list(test.statistic = t2, test.result = 1, p.value = p.value))
    } else {
        return(list(test.statistic = t2, test.result = 0, p.value = p.value))
    }
}


# Circular/Toroidal one-sample Hotelling-test using Bootstrap methods with p.value as output
bootstrap.hotelling.one.sample.test <- function(x.sample, x.mu, significance.level = 0.05, bootstrap.rep = 1000) {
    # Convert the sample to a matrix
    x.sample <- as.matrix(x.sample)

    # Calculate the multivariate mean of the sample
    x.mean <- fre.Mean.Multivariate(x.sample)

    # Initialize the distance vector and the covariance matrix for the bootstrap samples
    d.vec.x <- matrix(0, nrow = length(x.sample[1, ]), ncol = bootstrap.rep)
    sigma.x <- matrix(0, nrow = length(x.sample[1, ]), ncol = length(x.sample[1, ]))

    # For each bootstrap replication
    for (i in 1:bootstrap.rep) {
        # Calculate the distance vector between the multivariate mean of the bootstrap sample and the multivariate mean of the original sample
        d.vec.x[, i] <- distance.vec(
            fre.Mean.Multivariate(generate.Bootstrap.Sample(x.sample, length(x.sample[, 1]))),
            x.mean
        )
        # Update the covariance matrix
        sigma.x <- sigma.x + d.vec.x[, i] %*% t(d.vec.x[, i])
    }

    # Divide the covariance matrix by the number of bootstrap replications to get the average
    sigma.x <- sigma.x / bootstrap.rep

    # Calculate the inverse of the covariance matrix
    sigma.inv <- solve(sigma.x)

    # Initialize the T^2 statistic for the bootstrap samples
    t1.bootstrap <- rep(0, bootstrap.rep)
    for (i in 1:bootstrap.rep) {
        # Calculate the T^2 statistic for each bootstrap sample
        t1.bootstrap[i] <- t(d.vec.x[, i]) %*% sigma.inv %*% (d.vec.x[, i])
    }

    # Calculate the critical value for the T^2 statistic
    c_boundary <- quantile(t1.bootstrap, probs = 1 - significance.level, type = 2)

    # Calculate the T^2 statistic for the original sample
    t1 <- t(distance.vec(x.mean, x.mu)) %*% sigma.inv %*% (distance.vec(x.mean, x.mu))

    # Compute the p-value
    p.value <- (1 - ecdf(t1.bootstrap)(t1))

    # If the T^2 statistic is greater than the critical value, return 1 and the p-value; otherwise, return 0 and the p-value
    if (t1 > c_boundary) {
        return(list(test.statistic = t1, test.result = 1, p.value = p.value))
    } else {
        return(list(test.statistic = t1, test.result = 0, p.value = p.value))
    }
}


# Circular/Toroidal two-sample Hotelling-test using Bootstrap methods with p.value as output
bootstrap.hotelling.two.sample.test <- function(x.sample, y.sample, significance.level = 0.05, bootstrap.rep = 1000) {
    # Convert the samples to matrices
    x.sample <- as.matrix(x.sample)
    y.sample <- as.matrix(y.sample)


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
        return(list(test.statistic = t2, test.result = 1, p.value = p.value))
    } else {
        return(list(test.statistic = t2, test.result = 0, p.value = p.value))
    }
}


# Test for presence finite sample smeariness
FSS.test <- function(x.sample, significance.level = 0.05, bootstrap.rep = 1000) {
    # Initialize a vector to hold bootstrap means
    mu_Bootstrap <- rep(0, bootstrap.rep)

    # Get the length of the input vector x
    n <- length(x.sample)

    # For each iteration up to bootstrap.rep,
    for (i in 1:bootstrap.rep) {
        # Generate a bootstrap sample from x and calculate its Frechet mean
        # Store the mean in mu_Bootstrap
        mu_Bootstrap[i] <- fre.Mean(generate.Bootstrap.Sample(x.sample, n))
    }

    # Calculate the Frechet mean of collection mu_Bootstrap
    mu_mu_Bootstrap <- fre.Mean(mu_Bootstrap)

    # Calculate the variance of sample x
    var_n <- calculateVariance(x.sample)

    # Calculate the variance of mu_Bootstrap
    var_n_bootstrap <- calculateVariance(mu_Bootstrap)

    # Calculate the mean of the fourth power of the distance between mu_mu_Bootstrap and each value in mu_Bootstrap
    W <- mean(distance.vec(rep(mu_mu_Bootstrap, bootstrap.rep), mu_Bootstrap)^4)

    # Calculate the Bootstrap Variance Modulation (BVM)
    BVM <- n * var_n_bootstrap / var_n

    # Calculate the quantile of the standard normal distribution at 1-significance.level
    phi_alpha <- qnorm(1 - significance.level, 0, 1)

    # Calculate the threshold value h
    h <- n * phi_alpha / sqrt(bootstrap.rep) * sqrt(W - var_n_bootstrap^2) / var_n

    # Calculate the p-value for the test
    p.value <- 1 - pnorm((BVM - 1) * var_n / sqrt(W - var_n_bootstrap^2) * sqrt(bootstrap.rep) / n, 0, 1)

    # If BVM - 1 is greater than h, return 1; otherwise, return 0
    # This is the result of the Finite Sample Smeariness (FSS) test
    if (BVM - 1 > h) {
        return(list(test.statistic = BVM, test.result = 1, p.value = p.value))
    } else {
        return(list(test.statistic = BVM, test.result = 0, p.value = p.value))
    }
}

# Define a function to test whether the full is contained in a semi-circle
contained.in.semi.circle <- function(sample) {
    # Calculate the mean of the sample
    sample.mean <- fre.Mean.Multivariate(sample)

    # Shift the sample by subtracting the mean and applying the modulo operation
    sample.shifted <- modulo(sample - sample.mean)

    # If the sample is a vector
    if (is.vector(sample)) {
        # If the range of the shifted sample is greater than pi, return TRUE
        if (max(sample.shifted) - min(sample.shifted) > pi) {
            return(FALSE)
        }
    } else {
        # If the sample is a matrix, iterate over the columns
        m <- length(sample[1, ])
        for (i in 1:m) {
            # If the range of the shifted sample for a column is greater than pi, return TRUE
            if (max(sample.shifted[, m]) - min(sample.shifted[, m]) > pi) {
                return(TRUE)
            }
        }
    }
    # If none of the above conditions are met, return FALSE
    return(FALSE)
}


# Function to decide which test to use in the one-sample case
guideline.one.sample.test <- function(x.sample, x.mu, significance.level = 0.05, bootstrap.rep = 1000, comments = TRUE) {
    # If the sample lies in a semi-circle, use the CLT based Hotelling test
    if (contained.in.semi.circle(x.sample)) {
        if (comments) {
            cat("Sample is not affected by the presence of finite sample smeariness.\n CLT based one-sample Hotelling test is performed.\n\n")
        }
        return(clt.hotelling.one.sample.test(x.sample, x.mu, significance.level))
    } else {
        # If the sample is not affected by the presence of finite sample smeariness, use the CLT based Hotelling test, otherwise use the bootstrap based Hotelling test
        if (FSS.test(x.sample, significance.level, bootstrap.rep)$test.result == 0) {
            if (comments) {
                cat("Sample is not affected by the presence of finite sample smeariness.\nCLT based one-sample Hotelling test is performed.\n\n")
            }
            return(clt.hotelling.one.sample.test(x.sample, x.mu, significance.level))
        } else {
            if (comments) {
                cat("Sample is affected by the presence of finite sample smeariness.\nBootstrap based two-sample Hotelling test is performed.\n\n")
            }
            return(bootstrap.hotelling.one.sample.test(x.sample, x.mu, significance.level, bootstrap.rep))
        }
    }
}


# Guideline for practitioners to perform the right test use in the two-sample case
guideline.two.sample.test <- function(x.sample, y.sample, significance.level = 0.05, bootstrap.rep = 1000, comments = TRUE) {
    # If both sample lies in a semi-circle, use the CLT based Hotelling test
    if (contained.in.semi.circle(x.sample) && contained.in.semi.circle(y.sample)) {
        if (comments) {
            cat("Both samples are not affected by the presence of finite sample smeariness.\n CLT based two-sample Hotelling test is performed.\n\n")
        }
        return(clt.hotelling.two.sample.test(x.sample, y.sample, significance.level))
    } else {
        # If both sample is not affected by the presence of finite sample smeariness, use the CLT based Hotelling test, otherwise use the bootstrap based Hotelling test
        if (FSS.test(x.sample, significance.level, bootstrap.rep)$test.result == 0 && FSS.test(y.sample, significance.level, bootstrap.rep)$test.result == 0) {
            if (comments) {
                cat("Both samples are not affected by the presence of finite sample smeariness.\nCLT based two-sample Hotelling test is performed.\n\n")
            }
            return(clt.hotelling.two.sample.test(x.sample, y.sample, significance.level))
        } else {
            if (comments) {
                cat("At least one sample is affected by the presence of finite sample smeariness.\nBootstrap based two-sample Hotelling test is performed.\n\n")
            }
            return(bootstrap.hotelling.two.sample.test(x.sample, y.sample, significance.level, bootstrap.rep))
        }
    }
}
