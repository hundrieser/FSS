##### Summary: Generates Figure 2 from appendix of manuscript, depicting the histograms of bootstrapped versionso of daily Frechet means of wind directions for the cities Basel and Goettingen during the year 2000 until 2019


# Load the 'lubridate' package for date-time manipulation
library(lubridate)

# Load the 'rstudioapi' package for accessing the RStudio API
library(rstudioapi)

# Load the 'ggplot2' package for creating graphics
library(ggplot2)

# Load the 'RColorBrewer' package for generating color palettes
library(RColorBrewer)

# Load the 'gridExtra' package for arranging multiple grid-based plots on a page
library(gridExtra)

# Load the 'gtable' package for arranging 'grobs' in tables
library(gtable)

# Load the 'grid' package for low-level graphics functions
library(grid)

# Load the 'gridExtra' package again for arranging multiple grid-based plots on a page
library(gridExtra)

# Load the 'viridis' package for colorblind-friendly color maps
library(viridis)

# Set the working directory to the directory of the currently active RStudio document
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Define the duration in hours for the wind data analysis
DURATION <- c(24)

# Define the number of data points for the wind data analysis
DATAPOINTS <- c(365)

# Generate a sequence of dates from 2000 to 2020 at the start of each year
current.time <- ymd_hms(paste(as.character(seq(2000, 2019, 1)), "-01-01 00:00:00"))

# Number of bootstrap samples for the bootstrap test
M <- 10000

# Initialize a list to store the test results for three different settings:
output <- vector(mode = "list", length = 2)

# Vector with the relevant city names
cities <- c("Goettingen", "Basel")
cities2 <- c("Göttingen", "Basel")


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

modulo360 <- function(x, mod = 360) {
  output <- x
  for (i in 1:length(x)) {
    while (output[i] >= mod) {
      output[i] <- output[i] - mod
    }
    while (output[i] < 0) {
      output[i] <- output[i] + mod
    }
  }
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



# Define a function to generate bootstrap sample means of wind data
generate.wind.bootstrap.sample.means <- function(m = M, label, cities, current.time, DURATION, DATAPOINTS) {
  # Set the seed for random number generation to 0 for reproducibility
  set.seed(0)

  # Initialize an array to store the mean angles of the wind data
  mean.angles <- array(dim = c(DATAPOINTS, length(cities), length(current.time)))

  # Print a message indicating that the wind data is being loaded
  cat("\nLoading wind data for",cities,"...")

  # Loop over the time periods and cities
  for (j in 1:length(current.time)) {
    for (i in 1:length(cities)) {
      # Compute the mean angles of the wind data for the current city and time period
      mean.angles[, i, j] <- get.wind.means(cities[i], current.time[j], DURATION, DATAPOINTS)
    }
  }

  # Print a message indicating that the wind data has been loaded
  cat(" Finished \n")

  # Print a message indicating that the bootstrap sample means are being generated
  cat("Generating bootstrap based means with replacement begins\n")

  # Initialize an array to store the bootstrap sample means
  bootstrap.means <- array(rep(0, length(current.time) * m), dim = c(length(current.time), m))

  # Loop over the time periods
  for (i in 1:(length(current.time))) {
    # Print a message indicating the current time period
    cat("Generating bootstrap sample means for year", year(current.time[i]), "\n")

    # Loop over the bootstrap replications
    for (k in 1:m) {
      # Randomly sample 365 indices with replacement
      indices <- sample(365, 365, TRUE)

      # Compute the bootstrap sample mean for the current time period and replication
      bootstrap.means[i, k] <- fre.Mean(mean.angles[indices, , i])
    }
  }

  # Return the bootstrap sample means
  return(bootstrap.means)
}

# Function to create a wind histogram
# The following code has been adapted from the original windrose function detailed in the following blog post: https://rstudio-pubs-static.s3.amazonaws.com/284981_2e1c4c62b74446008d7ca449b744f7ab.html
wind.histogram <- function(data,
                           spd,
                           dir,
                           spdres = 5,
                           dirres = 11.25,
                           spdmin = 0,
                           spdmax = 65,
                           spdseq = NULL,
                           palette = "YlGnBu",
                           countmax = NA,
                           debug = 0) {
  # Check the type of input data
  if (is.numeric(spd) & is.numeric(dir)) {
    # If speed and direction are numeric, create a data frame
    data <- data.frame(
      spd = spd,
      dir = dir
    )
    spd <- "spd"
    dir <- "dir"
  } else if (exists("data")) {
    # If a data frame is provided, use it as is
  }

  # Clean the input data
  n.in <- NROW(data)
  dnu <- (is.na(data[[spd]]) | is.na(data[[dir]]))
  data[[spd]][dnu] <- NA
  data[[dir]][dnu] <- NA

  # Determine the wind speed bins
  if (missing(spdseq)) {
    spdseq <- seq(spdmin, spdmax, spdres)
  } else {
    if (debug > 0) {
      cat("Using custom speed bins \n")
    }
  }

  # Get information about the number of bins
  n.spd.seq <- length(spdseq)
  n.colors.in.range <- n.spd.seq - 1

  # Create the color map
  spd.colors <- colorRampPalette(brewer.pal(
    min(
      max(
        3,
        n.colors.in.range
      ),
      min(
        9,
        n.colors.in.range
      )
    ),
    palette
  ))(n.colors.in.range)

  # Adjust the speed breaks and labels if the max speed is greater than the specified max
  if (max(data[[spd]], na.rm = TRUE) > spdmax) {
    spd.breaks <- c(
      spdseq,
      max(data[[spd]], na.rm = TRUE)
    )
    spd.labels <- c(
      paste(
        c(spdseq[1:n.spd.seq - 1]),
        "-",
        c(spdseq[2:n.spd.seq])
      ),
      paste(
        spdmax,
        "-",
        max(data[[spd]], na.rm = TRUE)
      )
    )
    spd.colors <- c(spd.colors, "darkblue")
  } else {
    spd.breaks <- spdseq
    spd.labels <- paste(
      c(spdseq[1:n.spd.seq - 1]),
      "-",
      c(spdseq[2:n.spd.seq])
    )
  }

  # Bin the wind speed data
  data$spd.binned <- cut(
    x = data[[spd]],
    breaks = spd.breaks,
    labels = spd.labels,
    ordered_result = TRUE
  )

  # Determine the wind direction bins
  dir.breaks <- c(
    -dirres / 2,
    seq(dirres / 2, 360 - dirres / 2, by = dirres),
    360 + dirres / 2
  )

  # Create labels for the direction bins
  dir.labels <- c(
    paste(360 - dirres / 2, "-", dirres / 2, "  "),
    paste(
      seq(dirres / 2, 360 - 3 * dirres / 2, by = dirres),
      "-",
      seq(3 * dirres / 2, 360 - dirres / 2, by = dirres), "  "
    ),
    paste(360 - dirres / 2, "-", dirres / 2, "  ")
  )

  # Bin the wind direction data
  dir.binned <- cut(data[[dir]],
    breaks = dir.breaks,
    ordered_result = TRUE
  )
  levels(dir.binned) <- dir.labels
  data$dir.binned <- dir.binned

  # Adjust the labels for the direction bins
  dir.labels <- c("N", " ", " ", " ", "NE", " ", " ", " ", "E", " ", " ", " ", "SE", " ", " ", " ", "S", " ", " ", " ", "SW", " ", " ", " ", "W", " ", " ", " ", "NW", " ", " ", " ", "N")

  # Reverse the order of the speed bins and colors
  data$spd.binned <- with(data, factor(spd.binned, levels = rev(levels(spd.binned))))
  spd.colors <- rev(spd.colors)

  # Create the plot
  p.windrose <- ggplot(
    data = data,
    aes(
      x = dir.binned,
      fill = spd.binned
    )
  ) +
    geom_bar() +
    scale_x_discrete(
      drop = FALSE,
      labels = dir.labels
    ) +
    scale_fill_manual(
      name = "Average Wind Speed (km/h)",
      values = spd.colors,
      drop = FALSE
    ) +
    ylim(c(0, NA)) +
    labs(
      x = "",
      y = "",
      title = ""
    ) +
    theme(
      legend.position = "bottom", legend.title = element_text(color = "black", size = 14),
      legend.text = element_text(color = "black", size = 12)
    )

  # Adjust the y-axis if required
  if (!is.na(countmax)) {
    p.windrose <- p.windrose +
      ylim(c(0, NA))
  }

  # Return the plot
  return(p.windrose)
}


# Define a function to generate a bootstrap wind histogram
get.bootstrap.wind.histogram <- function(m = M, yearIndex, cityIndex) {
  # Create a vector of wind speeds, all set to 5
  speed <- rep(65, m)
  
  # Convert the bootstrap means from radians to degrees and ensure they are within the range 0-360
  angles <- modulo360(as.numeric(output[[cityIndex]][yearIndex, ] * 360 / (2 * pi)), 360)

  # Generate a wind histogram using the 'wind.histogram' function from the 'openair' package
  # Remove the legend and set the y-axis limit to 9000
  hist <- wind.histogram(
    spd = speed,
    dir = angles
  ) + theme(legend.position = "none") + ylim(0, 9000)

  # Return the wind histogram
  return(hist)
}


# ### Generate bootstrap sample means from wind data
for (i in 1:2) {
  output[[i]] <- generate.wind.bootstrap.sample.means(m = M, label[1], cities[i], current.time, DURATION, DATAPOINTS)
}

# Open a PDF file to save the plots
pdf(paste("Appendix_Fig2.pdf", sep = ""), height = 9, width = 14)

# Loop over the two cities
for (cityIndex in 1:2) {
  # Arrange multiple plots in a grid using 'grid.arrange' from the 'gridExtra' package
  # Each plot is a bootstrap wind histogram for a specific year from 2000 to 2019
  # The title of each plot is the corresponding year
  # The plots are arranged in 4 rows
  # The left label of the grid is "Frequency"
  # The top label of the grid is the city name and the years
  grid.arrange(arrangeGrob(
    get.bootstrap.wind.histogram(m=M, 1, cityIndex) + ggtitle("2000") + theme_bw() + theme(panel.border = element_blank(), plot.margin = margin(t = 6), legend.position = "none"),
    get.bootstrap.wind.histogram(m=M, 2, cityIndex) + ggtitle("2001") + theme_bw() + theme(panel.border = element_blank(), plot.margin = margin(t = 6), legend.position = "none"),
    get.bootstrap.wind.histogram(m=M, 3, cityIndex) + ggtitle("2002") + theme_bw() + theme(panel.border = element_blank(), plot.margin = margin(t = 6), legend.position = "none"),
    get.bootstrap.wind.histogram(m=M, 4, cityIndex) + ggtitle("2003") + theme_bw() + theme(panel.border = element_blank(), plot.margin = margin(t = 6), legend.position = "none"),
    get.bootstrap.wind.histogram(m=M, 5, cityIndex) + ggtitle("2004") + theme_bw() + theme(panel.border = element_blank(), plot.margin = margin(t = 6), legend.position = "none"),
    get.bootstrap.wind.histogram(m=M, 6, cityIndex) + ggtitle("2005") + theme_bw() + theme(panel.border = element_blank(), plot.margin = margin(t = 6), legend.position = "none"),
    get.bootstrap.wind.histogram(m=M, 7, cityIndex) + ggtitle("2006") + theme_bw() + theme(panel.border = element_blank(), plot.margin = margin(t = 6), legend.position = "none"),
    get.bootstrap.wind.histogram(m=M, 8, cityIndex) + ggtitle("2007") + theme_bw() + theme(panel.border = element_blank(), plot.margin = margin(t = 6), legend.position = "none"),
    get.bootstrap.wind.histogram(m=M, 9, cityIndex) + ggtitle("2008") + theme_bw() + theme(panel.border = element_blank(), plot.margin = margin(t = 6), legend.position = "none"),
    get.bootstrap.wind.histogram(m=M, 10, cityIndex) + ggtitle("2009") + theme_bw() + theme(panel.border = element_blank(), plot.margin = margin(t = 6), legend.position = "none"),
    get.bootstrap.wind.histogram(m=M, 11, cityIndex) + ggtitle("2010") + theme_bw() + theme(panel.border = element_blank(), plot.margin = margin(t = 6), legend.position = "none"),
    get.bootstrap.wind.histogram(m=M, 12, cityIndex) + ggtitle("2011") + theme_bw() + theme(panel.border = element_blank(), plot.margin = margin(t = 6), legend.position = "none"),
    get.bootstrap.wind.histogram(m=M, 13, cityIndex) + ggtitle("2012") + theme_bw() + theme(panel.border = element_blank(), plot.margin = margin(t = 6), legend.position = "none"),
    get.bootstrap.wind.histogram(m=M, 14, cityIndex) + ggtitle("2013") + theme_bw() + theme(panel.border = element_blank(), plot.margin = margin(t = 6), legend.position = "none"),
    get.bootstrap.wind.histogram(m=M, 15, cityIndex) + ggtitle("2014") + theme_bw() + theme(panel.border = element_blank(), plot.margin = margin(t = 6), legend.position = "none"),
    get.bootstrap.wind.histogram(m=M, 16, cityIndex) + ggtitle("2015") + theme_bw() + theme(panel.border = element_blank(), plot.margin = margin(t = 6), legend.position = "none"),
    get.bootstrap.wind.histogram(m=M, 17, cityIndex) + ggtitle("2016") + theme_bw() + theme(panel.border = element_blank(), plot.margin = margin(t = 6), legend.position = "none"),
    get.bootstrap.wind.histogram(m=M, 18, cityIndex) + ggtitle("2017") + theme_bw() + theme(panel.border = element_blank(), plot.margin = margin(t = 6), legend.position = "none"),
    get.bootstrap.wind.histogram(m=M, 19, cityIndex) + ggtitle("2018") + theme_bw() + theme(panel.border = element_blank(), plot.margin = margin(t = 6), legend.position = "none"),
    get.bootstrap.wind.histogram(m=M, 20, cityIndex) + ggtitle("2019") + theme_bw() + theme(panel.border = element_blank(), plot.margin = margin(t = 6), legend.position = "none"),
    nrow = 4,
    left = textGrob("\nFrequency", rot = 90, vjust = 0.5, gp = gpar(fontsize = 14)),
    top = textGrob(paste("Histogram of Bootstrapped Sample Means Wind Directions for ", cities2[cityIndex], " during 2000 to 2019", sep = ""), vjust = 0.5, gp = gpar(fontsize = 17))
  ))
}

# Close the PDF device
dev.off()

# save.image("simulation_UniquenessPopulationMean.RData")
