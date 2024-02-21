##### Summary: Generates Figure 1 of manuscript, depicting the variance modulation for multiple von Mises mixture distributions with disks of multiple radii cut out at around the antipode of the mean.


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
# Vector of sample sizes
N.vector <- unique(ceiling(10^(seq(0, 4, 1))))

# Number of repetitions for the simulation
R <- 100

# Initialize an empty string for INDEX
INDEX <- ""

# Set the parameter for the von Mises distribution to three
upperpar <- 3

# Create a list of parameter sets for the von Mises distribution
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
hole.vec <- c(0, 0.1, 0.2, round(pi / 2, 2))


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

### Simulation of variance modulation

# Compute the variance modulation for varying sample sizes (N.vector), repetitions (r), given parameters for distributional family, and for some hole size.
simulate.Var.Mod <-
  function(N.vector, r = R, parameter, distr, hole) {
    # Initialize vectors to hold the variance values and sample means
    varianceValues <- rep(0, length(N.vector))
    sampleMeans <- rep(0, r)

    # For each sample size in N.vector,
    for (N in 1:length(N.vector)) {
      # Print a message to indicate the current sample size
      cat(paste("   Calculating Variance  for N = ", N.vector[N], " ... "))

      # For each repetition,
      for (n in 1:r) {
        # Generate a sample from the specified distribution and compute its mean
        # The mean is computed using the fre.Mean function, which computes the mean of a circular variable
        sampleMeans[n] <- fre.Mean(generate.Sample(distr, parameter, N.vector[N], hole))
      }

      # Print a message to indicate the completion of the current sample size
      cat(paste("Finished \n"))

      # Compute the variance of the sample means and multiply it by the sample size
      # This is done to normalize the variance
      varianceValues[N] <- N.vector[N] * mean(sampleMeans^2)
    }

    # Compute the true variance as the variance of the first sample size
    trueVariance <- varianceValues[1]

    # Compute the variance modulation as the ratio of each variance to the true variance
    varianceModulation <- varianceValues / trueVariance

    # Export the computed variances to a separate file
    write.table(
      varianceValues,
      file = paste("VarianceData_Run_", INDEX, ".txt", sep = ""),
      row.names = FALSE,
      col.names = FALSE
    )

    # Export the computed variance modulations to a separate file
    write.table(
      varianceModulation,
      file = paste("VarianceDataRatio_Run_", INDEX, ".txt", sep = ""),
      row.names = FALSE,
      col.names = FALSE
    )

    # Print a newline character to separate the output for different sample sizes
    cat("\n")

    # Return a list containing the computed variances and variance modulations
    return(list(varianceValues = varianceValues, varianceModulation = varianceModulation))
  }


# Create a list of index names by concatenating the distribution type and hole size
INDEX.list <- paste(rep(distr, each = length(hole.vec)), "hole_", hole.vec, sep = "")

# Initialize a list to hold the variance values for each combination of parameters and hole sizes
varianceValues.list <- vector(mode = "list", length = length(parameter.list) * length(hole.vec))

# Initialize a list to hold the variance values in a dataframe format for each combination of parameters and hole sizes
varianceValues.dataframe.list <- vector(mode = "list", length = length(parameter.list) * length(hole.vec))

# Initialize a list to hold the combined variance values for each set of parameters
variance.combined <- vector(mode = "list", length = length(parameter.list))

# Initialize a list to hold the plot objects for the variance values for each set of parameters
variance.plot.list <- vector(mode = "list", length = length(parameter.list))


# Define a sequence of indices
index.set <- seq(1, 4)

# Run simulation on variance modulation
for (i in index.set) {
  # Initialize a data frame to hold the combined variance values for different hole sizes
  variance.combined.holes <- data.frame(n = integer(), variances = double(), index = as.character(), h = as.character())
  # Set the seed for the random number generator to ensure reproducibility
  set.seed(1)

  # For each hole size,
  for (j in 1:length(hole.vec)) {
    # Determine the index for the current combination of parameters and hole size
    INDEX <- INDEX.list[[i + (j - 1) * length(parameter.list)]]

    # Run the variance modulation simulation for the current combination of parameters and hole size
    varianceValues.list[[i + (j - 1) * length(parameter.list)]] <-
      simulate.Var.Mod(N.vector, r = R, parameter = parameter.list[[i]], distr = distr[i], h = hole.vec[j])

    # Convert the variance values to a data frame format
    varianceValues.dataframe.list[[i + (j - 1) * length(parameter.list)]] <- data.frame(
      n = N.vector,
      variances = varianceValues.list[[i + (j - 1) * length(parameter.list)]]$varianceModulation,
      index = paste(i + (j - 1) * length(parameter.list)),
      r = paste(hole.vec[j])
    )

    # Combine the variance values for different hole sizes
    variance.combined.holes <- rbind(
      variance.combined.holes,
      varianceValues.dataframe.list[[i + (j - 1) * length(parameter.list)]]
    )

    # Uncomment the following line to save the data before the simulation is finished
    # save.image(paste("Var_Mod_Data",i,"_hole",j,".RData",sep=""))
  }

  # Store the combined variance values for the current set of parameters
  variance.combined[[i]] <- variance.combined.holes
}
print("Finished")
# Save the workspace to a .RData file
save.image("Var_Mod_Data.RData")


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

# Define the lower and upper limits for the y-axis of the plots
# These are specified in a logarithmic scale
scalaLower <- c(10^-0.25, 10^-0.25, 10^-0.25, 10^-0.45, 10^-0.45) # Lower limits
scalaUpper <- c(40, 40, 40, 600, 1200) # Upper limits

# Define a list of titles for the plots
# Each title is a LaTeX expression that includes the parameters for the von Mises distribution
# The TeX function is used to convert the expression to a format that can be displayed in a plot
TITLE.list <- c(
  TeX(paste("   $P^{(", parameter.list[[1]][1], ",", parameter.list[[1]][2], ",", parameter.list[[1]][3], ",r)}_{vMm}$", sep = "")),
  TeX(paste("   $P^{(", parameter.list[[2]][1], ",", parameter.list[[2]][2], ",", parameter.list[[2]][3], ",r)}_{vMm}$", sep = "")),
  TeX(paste("   $P^{(", parameter.list[[3]][1], ",", parameter.list[[3]][2], ",", parameter.list[[3]][3], ",r)}_{vMm}$", sep = "")),
  TeX(paste("   $P^{(", parameter.list[[4]][1], ",", parameter.list[[4]][2], ",", round(parameter.list[[4]][3], 4), ",r)}_{vMm}$", sep = ""))
)
# For each index in the set,
for (i in index.set) {
  # Create a ggplot object for the variance values
  variance.plot.list[[i]] <- ggplot(
    data = variance.combined[[i]],
    aes(
      x = n, y = variances,
      group = r
    )
  ) + # Define aesthetics
    geom_line(aes(col = r, group = r)) + # Add lines with color based on 'r'
    ggtitle(TITLE.list[i]) + # Add title from the list of titles
    theme(panel.grid.minor = element_blank()) + # Remove minor grid lines
    scale_y_log10(TeX("$n$ E\\[$\\mu_n^2$\\]/$\\sigma$^2"),
      limits = c(scalaLower[i], scalaUpper[i]), # Set y-axis to log scale with specified limits
      breaks = scales::trans_breaks("log10", function(x) 10^x), # Define breaks for y-axis
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) + # Define labels for y-axis
    scale_x_log10("Sample size n", # Set x-axis to log scale
      breaks = scales::trans_breaks("log10", function(x) 10^x), # Define breaks for x-axis
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) + # Define labels for x-axis
    scale_colour_grey(start = 0, end = 0.75) + # Set color scale to grayscale
    annotation_logticks(sides = "b") + # Add log ticks on the bottom side
    theme_bw() + # Use a theme with a white background
    theme(legend.position = "bottom") # Position the legend at the bottom
}

# Initialize a vector to hold the asymptote values for each set of parameters
asymptotes <- rep(0, 4)

# For each index in the set,
for (i in index.set) {
  # Calculate the asymptote value for the current set of parameters
  # The vmm.Antipodal.Density function is used to calculate the asymptote value
  # The parameters for the function are extracted from the current set of parameters
  # The result is multiplied by 2*pi to convert from a density to a rate
  asymptotes[i] <- vmm.Antipodal.Density(
    parameter.list[[i]][1],
    parameter.list[[i]][3],
    parameter.list[[i]][2], 0
  ) * 2 * pi
}



# Begin outputting plot to a PDF file named "Fig1.pdf" with specified height and width
pdf("Fig1.pdf", height = 3, width = 9 * 3 / 4)

# Arrange multiple ggplot objects in a grid
# The arrangeGrob function is used to create a grob (graphical object) that arranges the plots in a grid
# The get_legend function is used to extract the legend from the first plot
# The grid.arrange function is used to arrange the grob and the legend in a grid with 2 rows
# The heights argument specifies the relative heights of the rows
# The scale_color_discrete function is used to specify the breaks for the color scale
grid.arrange(arrangeGrob(variance.plot.list[[1]] + theme(legend.position = "none") + geom_hline(yintercept = 1 / (1 - asymptotes[1])^2, linetype = "dotted"),
  # variance.plot.list[[2]]+ theme(legend.position = "none") + geom_hline(yintercept = 1/(1-asymptotes[2])^2, linetype = "dotted"),
  variance.plot.list[[3]] + theme(legend.position = "none") + geom_hline(yintercept = 1 / (1 - asymptotes[3])^2, linetype = "dotted"),
  variance.plot.list[[4]] + theme(legend.position = "none") + geom_abline(slope = 2 / 3 + 0.015, intercept = 0, linetype = "dotted"),
  nrow = 1
),
get_legend(variance.plot.list[[1]]),
nrow = 2, heights = c(3.7, 0.75)
) + scale_color_discrete(breaks = paste(hole.vec[length(hole.vec):1]))

# Close the PDF device, finalizing the plot
dev.off()

# Save the current R session to a file named "Fig1.RData"
save.image("Fig1.RData")
