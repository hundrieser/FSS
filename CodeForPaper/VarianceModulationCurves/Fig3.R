##### Summary: Generates Figure 3 of manuscript, depicting the behavior of the variance modulation for varying sample size in the regime nearby population measures with non-unique means.


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
N.vector <- unique(ceiling(10^(c(seq(0, 2, 0.1), seq(2.5, 4,0.5)))))

# Define the number of repetitions for the simulation
R <- 1000

# Initialize an empty string for INDEX
INDEX <- ""

# Initialize of epsilon values for the distribution
epsilon = c(-0.2,0,0.2,0.5)

# Initialize a vector of weight values for the distribution
weight = c(0.05,0.1,0.25)

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


  # Find the candidates which admits nearly the smallest Frechet value
  # The following is a modification of the standard algorithm from Fig1.R and Fig2.R to sidestep numerical issues related to non-unique Frechet means
   intrinsicSampleMeanAndFrechetValue = as.vector(localMinimisersAndFrechetValues[(localMinimisersAndFrechetValues[, 
    2] <= (min(localMinimisersAndFrechetValues[, 2])) + 
    0.02), 1])

  # Among those candidates, return the one with the largest absolute value
  return(intrinsicSampleMeanAndFrechetValue[which.max(abs(intrinsicSampleMeanAndFrechetValue))])
  # return(intrinsicSampleMeanAndFrechetValue[1])
}


### Sample generating mechanisms

# Generates a sample from a distribution that is concentrated at the North polt and near the equator
generate.Sample.Equator.Dist <- function(N = 10000, epsilon=0.1, weightNearEquator=0.1){
  
  # Generate N uniform random numbers between 0 and 1
  x =  runif(N,0,1)
  
  # Create a vector z that is -1 for x values in the lower tail, 1 for x values in the upper tail, and 0 otherwise
  z =(x >= 1 - 2*weightNearEquator)- 2*(x>= 1 - weightNearEquator)
  
  # Calculate y based on z: for z = -1 or 1, y is (pi/2 + epsilon) times z; for z = 0, y is a scaled and shifted version of x
  y = z*(pi/2 +epsilon) + (z==0)*(x/(1 - 2*weightNearEquator) - 0.5)
  
  # Return the vector y
  return(y)
}


### Simulation of variance modulation

# Computes variance modulation for varying sample sizes (N.vector), repetitions (r), and parameters (epsilon, weightNearEquator)
simulate.Var.Mod <- function(N.vector, r = R,  epsilon=0.1, weightNearEquator=0.1){
  
  # Initialize vectors to hold the variance values and sample means
  varianceValues = rep(0,length(N.vector))
  sampleMeans = rep(0, r)
  
  # Calculate the true variance of the distribution
  trueVariance = 1/12*(1- 2*weightNearEquator) + (pi/2 + epsilon)^2*2*weightNearEquator
  
  # For each sample size in N.vector,
  for(N in 1:length(N.vector)){
    
    # Print a message indicating that the variance calculation is starting for this sample size
    cat(paste("   Calculating Variance  for N = ", N.vector[N], " ... ")) 
    
    # For each trial,
    for(n in 1:r){
      
      # Generate a sample from the distribution and compute its mean
      sampleMeans[n] = fre.Mean(generate.Sample.Equator.Dist(N.vector[N], epsilon, weightNearEquator)) 
    }
    
    # Print a message indicating that the variance calculation is finished for this sample size
    cat(paste("Finished \n")) 
    
    # Compute the variance of the sample means and store it in varianceValues
    varianceValues[N] = N.vector[N]*mean(sampleMeans^2) 
  }
  
  # Write the variance values to a file
  write.table(varianceValues, 
              file=paste("VarianceData_Run_Eps",epsilon,"_Weight",weightNearEquator,".txt",sep=""),
              row.names = FALSE, 
              col.names = FALSE)
  
  # Compute the ratio of the variance values to the true variance
  varianceValuesRatio = varianceValues/trueVariance
  
  # Write the variance ratios to a file
  write.table(varianceValuesRatio, 
              file=paste("VarianceDataRatio_Run_Eps",epsilon,"_Weight",weightNearEquator,".txt",sep=""),
              row.names = FALSE, 
              col.names = FALSE)
  
  # Print a newline
  cat("\n")
  
  # Return a list containing the variance values and variance ratios
  return(list(varianceValues = varianceValues, varianceValuesRatio = varianceValuesRatio))
}

# Define a function to extract the legend from a ggplot object
get_legend <- function(myggplot) {
  # Convert the ggplot object to a gtable object
  tmp <- ggplot_gtable(ggplot_build(myggplot))

  # Find the index of the legend (named "guide-box") in the list of grobs (graphical objects)
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")

  #  Extract the legend from the list of grobs and return it
  return(tmp$grobs[[leg]])
}



# Create a vector epsilon.vec by repeating the epsilon vector as many times as there are elements in weight
epsilon.vec = rep(epsilon,length(weight))

# Create a vector weight.vec by repeating each element of weight as many times as there are elements in epsilon
weight.vec =  rep(weight,each = length(epsilon))

# Create a list of index strings by concatenating the numbers from 1 to the length of epsilon.vec with "_3point"
INDEX.list = paste(c(1:length(epsilon.vec)), "_3point",sep="")

# Initialize a list of titles with zeros, with the same length as weight
TITLE.list= rep(0,length(weight))

# For each index in the length of weight,
for(i in 1:length(weight)){
  
  # Create a title string using TeX formatting and the current weight value, and store it in TITLE.list
  TITLE.list[i] = TeX(paste("   P$^{(","\\epsilon",",",weight[i],")}_{E}$",sep=""))
}

# Initialize a list to hold variance values, with a length equal to the number of elements in epsilon.vec
varianceValues.list = vector(mode = "list", length = length(epsilon.vec))

# Initialize a list to hold dataframes of variance values, with a length equal to the number of elements in epsilon.vec
varianceValues.dataframe.list = vector(mode = "list", length = length(epsilon.vec))

# Initialize a list to hold combined data, with a length equal to the number of elements in weight
combinedList.list = vector(mode = "list", length = length(weight))

# Initialize a list to hold plots of variance, with a length equal to the number of elements in weight
variance.plot.list = vector(mode = "list", length = length(weight))

# For each weight value,
for(j in 1:length(weight)){
  
  # Initialize an empty dataframe to hold combined data
  combinedList = data.frame(n = integer(), variances = double(), index = as.character(), wght = as.character(),  epsilon= as.character())
  
  # For each epsilon value,
  for(i in 1:length(epsilon)){
    
    # Set the seed for random number generation
    set.seed(1+i+(j-1)*length(epsilon)+1)
    
    # Simulate variance modulation and store the ratio of variance values in varianceValues.list
    varianceValues.list[[i+(j-1)*length(epsilon)]] = simulate.Var.Mod(N.vector, r = R, epsilon = epsilon.vec[i+(j-1)*length(epsilon)], weightNearEquator = weight.vec[i+(j-1)*length(epsilon)])$varianceValuesRatio
    
    # Create a dataframe of variance values and store it in varianceValues.dataframe.list
    varianceValues.dataframe.list[[i+(j-1)*length(epsilon)]] <- data.frame(n = N.vector, 
                                                                       variances = varianceValues.list[[i+(j-1)*length(epsilon)]], 
                                                                       index = paste(i), 
                                                                       wght = paste(weight.vec[i+(j-1)*length(epsilon)]), 
                                                                       epsilon= paste(epsilon.vec[i+(j-1)*length(epsilon)]) )
    
    # Add the dataframe to combinedList
    combinedList = rbind(combinedList, 
                         varianceValues.dataframe.list[[i+(j-1)*length(epsilon)]])
  }
  
  # Store the combined data in combinedList.list
  combinedList.list[[j]] = combinedList
  
  # Uncomment the following line to save the workspace after each iteration
  # save.image(paste("Var_Mod_Data",j,".RData",sep=""))
}

# Print a message indicating that the simulation is finished
print("Finished")

# Save the workspace
save.image("Var_Mod_Data.RData")


# Define a vector of scale values
scala = c(2,2,10^3)

# For each weight value,
for(i in 1:length(weight)){
  
  # Create a ggplot of the variance values from the combined data
  variance.plot.list[[i]]<- ggplot(data=combinedList.list[[i]], aes(x=n, y=variances, group = epsilon))+
  
    # Add a line for each epsilon value, with color indicating the epsilon value
    geom_line(aes(col = epsilon,  group=epsilon))+
    
    # Add a title to the plot from TITLE.list
    ggtitle(TITLE.list[i])+
    
    # Remove minor grid lines
    theme(panel.grid.minor = element_blank())+
    
    # Set the y-axis to a log10 scale, with labels formatted as powers of 10
    scale_y_log10(TeX("$n$ E\\[$\\mu_n^2$\\]/$\\sigma$^2"),limit = c(0.6,7),
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    
    # Set the x-axis to a log10 scale, with labels formatted as powers of 10
    scale_x_log10("Sample size n",
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)) ) +
    
    # Set the color scale manually, with specific colors for each epsilon value
    scale_color_manual(name="eps", 
                       labels = c("-0.2", 
                                  "0", 
                                  "0.2", 
                                  "0.5" ), 
                       values = c( "-0.2"="black",
                                   "0"="blue",
                                   "0.2"="purple",
                                   "0.5"="red") )+                               
    
    # Set the color scale to grey, with start and end values
    scale_colour_grey(start =0.75, end = 0 )+                             
    
    # Add log ticks to the bottom of the plot
    annotation_logticks(sides = "b")+
    
    # Set the theme to black and white
    theme_bw()+
    
    # Position the legend at the bottom of the plot
    theme(legend.position = "bottom")
}


# Open a PDF device to save the plots
pdf("Fig3.pdf", height=3, width = 8)

# Arrange the plots from variance.plot.list in a grid, with one row and no legends
# Then add the legend from the third plot below the grid
# The grid of plots is given more height than the legend
grid.arrange(arrangeGrob(variance.plot.list[[1]]+ theme(legend.position = "none"),
             variance.plot.list[[2]]+ theme(legend.position = "none"),
             variance.plot.list[[3]]+ theme(legend.position = "none"),
             nrow=1 ),
       get_legend(variance.plot.list[[3]]), nrow=2,heights=c(3.7, 0.75))

# Close the PDF device, finalizing the plot
dev.off()

# Save the current R session to a file named "Fig3.RData"
save.image("Fig3.RData")


