##### Summary: Generates Figure 2 of manuscript, depicting the behavior of the variance modulation for varying sample size in relation to the density around the antipode of the mean.


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


### Parameters for initialisation of simulation
# Create a vector of sample sizes
N.vector <- unique(ceiling(10^(seq(0, 4, 0.5))))

# Define the number of repetitions for the simulation
R <- 1000

# Initialize an empty string for INDEX
INDEX <- ""

# Create a vector of x-values for the density plot
x.pos <- seq(-pi, pi, 0.005)

# Create a vector of hole values
hole.vec <- c(0, 0.1, 0.2, round(pi / 2, 2))


# Define a custom modulo function to handle 2*pi as period
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

# Generate samples at the antipode according to a specified cumulative distribution function (cdf)
generate.Sample.Antipode.Rect <- function(N, xcorTotal, cdf) {
  # Generate N uniform random numbers between 0 and 1
  x <- runif(N, 0, 1)

  # Initialize a vector of zeros with length N
  y <- rep(0, N)

  # For each generated random number,
  for (j in 1:N) {
    # For each interval in the cdf (except the last one),
    for (i in 1:(length(cdf) - 1)) {
      # If the random number falls within the current interval,
      if (x[j] >= cdf[i] && x[j] < cdf[i + 1]) {
        # Calculate the corresponding value on the y-axis using linear interpolation
        y[j] <- (x[j] - cdf[i]) / (cdf[i + 1] - cdf[i]) * (xcorTotal[i + 1] - xcorTotal[i]) + xcorTotal[i]
      }
    }
  }

  # Return the generated samples
  return(y)
}


# Generate samples according to a specified weight distribution
generate.Sample.Rect <- function(N, xcor, weight) {

  # Create a sorted vector of x coordinates, including both positive and negative values
  xcorTotal <- sort(c(-xcor, xcor))

  # Create a vector of weights, including both positive and negative values
  weightTotal <- c(weight[length(weight):2], weight)

  # Initialize a vector of zeros to hold the mass of each block
  massBlocks <- rep(0, 2 * length(xcor) - 1)

  # For each interval in the x coordinates,
  for (i in 1:(length(xcorTotal) - 1)) {
    # Calculate the mass of the block as the width of the interval times the weight
    massBlocks[i] <- (xcorTotal[i + 1] - xcorTotal[i]) * weightTotal[i]
  }

  # Calculate the total mass as the sum of the block masses
  totalMass <- sum(massBlocks)

  # Calculate the cumulative distribution function (cdf) as the cumulative sum of the block masses divided by the total mass
  cdf <- c(0, cumsum(massBlocks) / totalMass)

  # Calculate the total mass at the antipode as the total mass divided by 2*pi
  totalMassAtAntipode <- totalMass / (2 * pi)

  # Generate N uniform random numbers between 0 and 1
  topBottom <- runif(N)

  # Calculate the number of samples that should be generated at the bottom (antipode)
  NumberBottom <- sum((topBottom <= totalMassAtAntipode))

  # Initialize a vector of zeros to hold the samples
  sample <- rep(0, N)

  # If no samples should be generated at the bottom,
  if (NumberBottom == 0) {
    # Return the sample vector as is
    return(sample)
  } else {
    # Generate the samples at the antipode
    antipodeSample <- generate.Sample.Antipode.Rect(NumberBottom, xcorTotal, cdf) + pi

    # Adjust the samples to fall within the range [-pi, pi]
    antipodeSample <- antipodeSample + (antipodeSample > pi) * (-2 * pi)

    # Store the antipode samples in the sample vector
    sample[1:NumberBottom] <- antipodeSample

    # Return the sample vector
    return(sample)
  }
}

# Compute the population variance for a rectangular distribution
pop.Variance.Rect <- function(xcor, weight) {
  # Initialize a vector to store the masses of the blocks
  massBlocks <- rep(0, 2 * length(xcor) - 1)
  
  # Create a vector of the x-coordinates sorted in ascending order
  xcorTotal <- sort(c(-xcor, xcor))
  
  # If the weight is a scalar, repeat it for all x-coordinates
  # Otherwise, create a vector of weights for the negative and positive x-coordinates
  if (length(weight) == 1) {
    weightTotal <- weight
  } else {
    weightTotal <- c(weight[length(weight):2], weight)
  }

  # Calculate the mass of each block
  for (i in 1:(length(xcorTotal) - 1)) {
    massBlocks[i] <- (xcorTotal[i + 1] - xcorTotal[i]) * weightTotal[i]
  }
  
  # Calculate the total mass
  totalMass <- sum(massBlocks)

  # Create a vector of the x-coordinates sorted in ascending order, including zero
  xcorTotal <- sort(c(-xcor, 0, xcor))
  
  # If the weight is a scalar, repeat it for all x-coordinates
  # Otherwise, create a vector of weights for the negative and positive x-coordinates
  if (length(weight) == 1) {
    weightTotal <- c(weight, weight)
  } else {
    weightTotal <- c(weight[length(weight):1], weight)
  }

  # Initialize a vector to store the differences
  differences <- rep(0, length(weightTotal) / 2)
  
  # Calculate the differences
  for (i in 1:length(weightTotal) / 2) {
    differences[i] <- ((xcorTotal[i + 1] + pi)^3 - (xcorTotal[i] + pi)^3) / 3 * (weightTotal[i] / (2 * pi))
  }

  # Calculate the true variance
  trueVariance <- sum(differences) * 2
  
  # Return the true variance
  return(trueVariance)
}

# Calculate the density of a rectangular distribution at specified points
density.Rect <- function(x, xcor, weight) {

  # Create a sorted vector of x coordinates, including both positive and negative values
  xcorTotal <- sort(c(-xcor, xcor))

  # If the weight vector has only one element,
  if (length(weight) == 1) {
    # Use the single weight value for all x coordinates
    weightTotal <- weight
  } else {
    # Create a vector of weights, including both positive and negative values
    weightTotal <- c(weight[length(weight):2], weight)
  }

  # Initialize a vector of zeros to hold the density values
  y <- rep(0, length(x))

  # For each x coordinate,
  for (i in 1:length(x)) {
    # For each weight,
    for (w in 1:length(weightTotal)) {
      # If the x coordinate falls within the current interval,
      if (x[i] >= xcorTotal[w] && x[i] < xcorTotal[w + 1]) {
        # Set the density value to the corresponding weight
        y[i] <- weightTotal[w]
      }
    }
  }

  # Return the density values
  return(y)
}


### Simulation of variance modulation

# Computes the variance modulation for varying sample sizes (N.vector), repetitions (r), given parameters for distributional family (xcor, weight)
simulate.Var.Mod <- function(N.vector, r = R, xcor, weight) {

  # Initialize a vector of zeros to hold the variance values
  varianceValues <- rep(0, length(N.vector))

  # Initialize a vector of zeros to hold the sample means
  sampleMeans <- rep(0, r)

  # For each sample size in N.vector,
  for (N in 1:length(N.vector)) {

    # Print a message indicating the current sample size
    cat(paste("   Calculating Variance  for N = ", N.vector[N], " ... "))

    # For each repetition,
    for (n in 1:r) {

      # Calculate the mean of the generated sample and store it in sampleMeans
      sampleMeans[n] <- fre.Mean(generate.Sample.Rect(N.vector[N], xcor, weight))
    }

    # Print a message indicating that the calculations for the current sample size are finished
    cat(paste("Finished \n"))

    # Calculate the variance for the current sample size and store it in varianceValues
    varianceValues[N] <- N.vector[N] * mean(sampleMeans^2)
  }

  # Calculate the true variance
  trueVariance <- pop.Variance.Rect(xcor, weight)

  # Write the variance values to a file
  write.table(varianceValues,
    file = paste("VarianceData_Run_", INDEX, ".txt", sep = ""),
    row.names = FALSE,
    col.names = FALSE
  )

  # Calculate the ratio of the variance values to the true variance
  varianceValuesRatio <- varianceValues / trueVariance

  # Write the variance values ratio to a file
  write.table(varianceValuesRatio,
    file = paste("VarianceDataRatio_Run_", INDEX, ".txt", sep = ""),
    row.names = FALSE,
    col.names = FALSE
  )

  # Print a newline
  cat("\n")

  # Return the variance values
  return(varianceValues)
}



# Define a list of x coordinates for each distribution
XCOR.list <- list(c(1.5), c(0.8, 2), c(0.1, 0.2, 0.5, 2), c(0.8, 2), c(0.1, 0.2, 0.8, 2))

# Define a list of weights for each distribution
WEIGHT.list <- list(c(1 / 2), c(1 / 2, 1), c(1 / 2, 0.8, 0, 1), c(0, 1), c(0, 0.85, 0, 1))

# Define a list of indices for each distribution
INDEX.list <- paste(c(1:5), sep = "")

# Define a list of titles for each distribution
TITLE.list <- paste("Example (", letters[1:length(XCOR.list)], ")", sep = "")

# Initialize a list to hold the variance values for each distribution
varianceValues.list <- vector(mode = "list", length = length(XCOR.list))

# Initialize a list to hold the density values for each distribution
densityValues.list <- vector(mode = "list", length = length(XCOR.list))

# Initialize a list to hold the variance plots for each distribution
variance.plot.list <- vector(mode = "list", length = length(XCOR.list))

# Initialize a list to hold the density plots for each distribution
density.plot.list <- vector(mode = "list", length = length(XCOR.list))


# Loop over each distribution in XCOR.list
for (i in 1:length(XCOR.list)) {

  # Set the seed for the random number generator to ensure reproducibility
  set.seed(5)

  # Set the current index
  INDEX <- INDEX.list[[i]]

  # Compute variance modulation for the current distribution and store the results in varianceValues.list
  varianceValues.list[[i]] <- simulate.Var.Mod(N.vector, r = R, xcor = XCOR.list[[i]], weight = WEIGHT.list[[i]])

  # Compute the density of the current distribution at the specified points and store the results in densityValues.list
  densityValues.list[[i]] <- density.Rect(x.pos, xcor = XCOR.list[[i]], weight = WEIGHT.list[[i]])

  # Uncomment the following line to save the workspace after each iteration
  # save.image(paste("Var_Mod_Data",i,".RData",sep=""))
}

# Print a message indicating that the loop has finished
print("Finished")

# Save the workspace to a file
save.image("Var_Mod_Data.RData")

# Define a vector of limits for each distribution plot
Limit <- c(4, 4, 4, 1, 1)

# Define a vector of scale factors for each distribution plot
scala <- c(6, 6, 6, 2.3, 2.3)

# Loop over each distribution in XCOR.list
for (i in 1:length(XCOR.list)) {

  # Create a data frame with the sample sizes and their corresponding variance values
  variance.data <- data.frame(n = N.vector, variances = varianceValues.list[[i]] / pop.Variance.Rect(xcor = XCOR.list[[i]], weight = WEIGHT.list[[i]]))

  # Create a plot of the variance values and add it to variance.plot.list
  variance.plot.list[[i]] <- ggplot(data = variance.data, aes(x = n, y = variances, group = 1)) +
    ggtitle(TITLE.list[i]) +
    geom_line(col = "black") +
    geom_hline(yintercept = Limit[i], linetype = "dashed") +
    scale_y_continuous(TeX("$n$ E\\[$\\mu_n^2$\\]/$\\sigma$^2"), limits = c(0.6, scala[i])) +
    scale_x_log10("Sample size n",
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    annotation_logticks(sides = "b") +
    theme(panel.grid.minor = element_blank()) +
    theme_bw()

  # Create a data frame with the x positions and their corresponding density values
  density.data <- data.frame(x = x.pos, density = as.numeric(densityValues.list[[i]]))

  # Create a plot of the density values and add it to density.plot.list
  density.plot.list[[i]] <- ggplot(data = density.data, aes(x = x, y = density, group = 1)) +
    geom_line() +
    theme(panel.grid.minor = element_blank()) +
    theme(axis.text.x = element_text(vjust = 1)) +
    scale_y_continuous(TeX("Density $\\cdot 2\\pi$"), limits = c(0, 1)) +
    scale_x_discrete(
      limit = c(0, 2, pi, -pi, -2),
      labels = c(
        expression(paste(pi, "=", -pi, sep = "")),
        expression(paste(-pi, "+", 2, sep = "")),
        expression(paste(0, sep = "")),
        expression(paste(0, sep = "")),
        expression(paste(pi, "-", 2, sep = ""))
      )
    ) +
    theme_bw()
}

# Open a PDF device to save the plots
pdf("Fig2.pdf", height = 4.1, width = 10)

# Arrange the plots in a grid and output them to the PDF device
grid.arrange(variance.plot.list[[1]],
  variance.plot.list[[2]],
  variance.plot.list[[3]],
  variance.plot.list[[4]],
  variance.plot.list[[5]],
  density.plot.list[[1]],
  density.plot.list[[2]],
  density.plot.list[[3]],
  density.plot.list[[4]],
  density.plot.list[[5]],
  nrow = 2, heights = c(2.4, 1.7)
)

# Close the PDF device
dev.off()

# Save the workspace to a file
save.image("Fig2.RData")
