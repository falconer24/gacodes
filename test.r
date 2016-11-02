#!/usr/bin/Rscript
source('ga2.r')

### Test code ###
# This code test the genetic algorithm procedures written in file ga2.r


# A test function to optimize
ftarget <- function(x)
{
 x*sin(10*pi*x)+1
}

# A sample GA problem

nPop<-25		

nBit<-22

pCrossover <- 0.25

pMutation <- 0.02

nGeneration <- 20 
# next 2 statements generate the example data and define it in the form of binary arrays
xraw <- as.integer(runif(nPop*nBit)<=0.5)	# Randomly generated population
x <- matrix(xraw, nrow=nPop, ncol=nBit, byrow=T)	

# the following statement defines the implementation specific data structure details, puts it in a
# list structure 
probspec=list(domain=c(-1,2), prec=6, desc="Real nos. encoded as binary strings. Domain and precision information required for encoding and decoding")

population <- list(n=nPop,b=nBit, detail=probspec, chr=x)	# Defines the generic data structure
								# and populates it

# call the genetic algorithm procedure to optimise the function

ga(ftarget,population, nPop, nBit, pCrossover, pMutation, nGeneration,verbose=TRUE)

### End of test code ###
