library(kdensity)
library(GoFKernel)

inFile <- "~/HIVtree/examples/CAP257_ENV_C2C3_W2_P2_QVOA_1_3921.txt"
args = commandArgs(trailingOnly = TRUE)
if(length(args) != 6) {
  stop("Incorrect number of arguments")
} else {
  
  # MCMC file
  inFile <- args[1]
  
  # Integration bounds
  lowerInt <- as.numeric(args[2])
  upperInt <- as.numeric(args[3])

  # Used to convert to calendar time, HIVtree reports in backwards time
  # with a different time scale
  timeUnit <- as.numeric(args[4])
  timeZero <- as.numeric(args[5])
  
  # Number of genes in the input file
  numGenes <- as.numeric(args[6])
}

times<- read.csv(inFile)

if ((dim(times)[2] / 2) != numGenes) {
  error("The number of genes reported does not match the number of columns in the input file. There should be twice as many columns at genes, one column for the posterior and one for the prior of each gene.")
}

maxGenes <- 10 
post <- data.frame(nrow = maxGenes)
prior <- data.frame(nrow = maxGenes)

# Give error message if too many genes
if (numGenes > 10) {
  error("Too many sequences. There must be 10 or fewer genes.")
}

# Calculates a KDE for the priors and posteriors 
for (i in 1:numGenes) {
  if (! all(is.na(times[ , i]))) {
    post[i] <-  kdensity(times[, i], na.rm = TRUE, normalize = FALSE)
    prior[i] <-  kdensity(times[, i + numGenes], na.rm = TRUE, normalize = FALSE)
    
  } else {
    # The data is missing
    post[i] <-  function(x) 1
    prior[i] <-  function(x) 1
  }

}

# There are less than 10 genes being analyzed
for (i in (numGenes+1):maxGenes) {
  post[i] <-  function(x) 1
  prior[i] <-  function(x) 1
}

# Converts between the time scale used in HIVtree to calendar time
adjTime <- function (time) {
  timeZero - time * timeUnit
}

# Calculates normalization constant
constant <- integrate(function(x) post[[1]](x) * post[[2]](x) * post[[3]](x) * post[[4]](x) * post[[5]](x) 
                      * post[[6]](x) * post[[7]](x) * post[[8]](x) * post[[9]](x) * post[[10]](x) / prior[[1]](x) / 
                        prior[[2]](x) / prior[[3]](x) / prior[[4]](x) / prior[[5]](x) / prior[[6]](x) /
                        prior[[7]](x) / prior[[8]](x) / prior[[9]](x) / prior[[10]](x), lower = lowerInt, upper =
                        upperInt)$value
dist <- function (x) {
  integrate(function (x) post[[1]](x) * post[[2]](x) * post[[3]](x) * post[[4]](x) * post[[5]](x) 
            * post[[6]](x) * post[[7]](x) * post[[8]](x) * post[[9]](x) * post[[10]](x) /  prior[[1]](x) / 
              prior[[2]](x) / prior[[3]](x) / prior[[4]](x) / prior[[5]](x) / prior[[6]](x) /
              prior[[7]](x) / prior[[8]](x) / prior[[9]](x) / prior[[10]](x) /constant, lower = lowerInt, upper = x)[[1]]
}

invFunc <- inverse(dist, lower = lowerInt, upper = upperInt)

# Finds the 95% credible set (equal tail probability)
lowerEst <- invFunc(.025)
upperEst <- invFunc(.975)

# Finds mean of combined posterior
meanEst <- integrate(function(x) x * post[[1]](x) * post[[2]](x) * post[[3]](x) * post[[4]](x) * post[[5]](x) 
                     * post[[6]](x) * post[[7]](x) * post[[8]](x) * post[[9]](x) * post[[10]](x) /  prior[[1]](x) / 
                       prior[[2]](x) / prior[[3]](x) / prior[[4]](x) / prior[[5]](x) / prior[[6]](x) /
                       prior[[7]](x) / prior[[8]](x) / prior[[9]](x) / prior[[10]](x) /constant, lower = lowerInt, upper = upperInt)[[1]]

print(paste("Mean:", adjTime(meanEst), ", 95% Credible Interval:", adjTime(upperEst),"-" , adjTime(lowerEst), sep = " "))
#plot(function (x) post[[1]](x) * post[[2]](x) * post[[3]](x) * post[[4]](x) * post[[5]](x) 
#     * post[[6]](x) * post[[7]](x) * post[[8]](x) * post[[9]](x) * post[[10]](x) /  prior[[1]](x) / 
#       prior[[2]](x) / prior[[3]](x) / prior[[4]](x) / prior[[5]](x) / prior[[6]](x) /
#       prior[[7]](x) / prior[[8]](x) / prior[[9]](x) / prior[[10]](x) /constant, xlim = c(lowerInt, upperInt), ylab = "Posterior Density", xlab = "Integration Time")
