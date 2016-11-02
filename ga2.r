#!/usr/bin/Rscript


# Decode function
# The interface is problem independent, the internals are problem specific
# This function uses information from sublist `detail' to do problem specific work
# Arguments -
#	1) population object name
#	2) number of chromosome to encode


decode <- function(pop, k)	
{
 # The following code is specific to this problem - i.e. binary string to real number

 d <- pop$detail$domain		# Get domain end points
 r <- d[2]-d[1]			# r is size of real domain
 prec <- pop$detail$prec	# prec contains number od decimal digits in the decoded real no.
 
 l <- pop$b			# l stores bit length
 p <- 2^seq((l-1),0,-1)		# p stores vector [2^(Bit length-1) ... 4, 2, 1] - required for conversion to decimal
 
 s <- pop$chr[k,]		# The bit vector itself is stored in s

 idecimal <- sum( s*p )		# convert bit vector to decimal integer
 
 res <- d[1] + (idecimal*r)/((2^l)-1)	# convert decimal integer to real number using domain and precision info

 return (res)
}

# Encode function
# The interface is problem independent, the internals are problem specific
# This function uses information from sublist `detail' to do problem specific work
# Arguments -
#	1) real number to encode to binary string
#	2) population object name

encode <- function(x,pop)
{
 # The following code is specific to this problem - i.e. real number to binary string

 d <- pop$detail$domain		# Get domain end points
 r <- d[2]-d[1]			# r is size of real domain
 prec <- pop$detail$prec	# prec contains number od decimal digits in the decoded real no.
 l <- pop$b			# l stores bit length

 idecimal <- ((x-d[1])*((2^l)-1))/r
 k<- rev(as.integer(intToBits(idecimal))[1:l])
 return(k)
}

# Mutation Operator function
# The interface is problem independent, the internals are problem specific
# This function uses information from sublist `detail' to do problem specific work
# Arguments -
#	1) population object name
#	2) probability of mutation

mutation <- function(pop, pMutation)
{
 # method 1: more stochastic
 mutants <- which(runif(pop$n*pop$b)<pMutation)		# generate as many random numbers from
 							# U(0,1) as there are bits in all the
							# chromosomes taken together, the indices
							# whose value is less than or equal
							# to prob of mutation are selected
 
 # method 2: more pseudo stochastic but computationally better (?)
 #mutants <- as.integer( runif(pop$n*pop$b*pMutation)*(pop$n*pop$b) )	
 							# generate as many random numbers as
 							# pMutation multiplied by total number of
							# bits and use the numbers themselves as
							# indices after scaling to the range
							# 1:(totalno of bits)
 
 
 if(length(mutants)==0)			# if randomness decides on no mutation, quietly leave
 	{
	 return(pop$chr)
	}
 
 for(i in 1:length(mutants))		# else, flip the bits denoted by the mutants vector
 	{
	 chrno <- (mutants[i] %/% pop$b) + 1	# get chromosome number from index
	 bno <- (mutants[i] %% pop$b) 		# get bit number from index
	 if(bno==0)				# for boundary bits
	 	{
		 bno <- pop$b
		 chrno <- chrno - 1
		}
	 #print(paste(mutants[i],"Bit",bno, "from Chr",chrno))	
	 #print(pop$chr[chrno,bno])  
	 pop$chr[chrno,bno] <- as.integer(xor(1,pop$chr[chrno,bno])) 
	 #print(pop$chr[chrno,bno])  
	}
 return(pop$chr)
}

# Crossover function
# The interface is problem independent, the internals are problem specific
# This function uses information from sublist `detail' to do problem specific work
# Arguments -

crossover <- function(pop,pCrossover)
{
 # the first part involves randomly selecting pairs - common to all data representations

 selected <- which ( runif(pop$n)<pCrossover )	# randomly select which chromosomes will pair
 
 if(length(selected)==0)			# if randomness decides on no crossover, quietly leave
 	{
	 #return(pop$chr)
	 return (-1)
	}

 if(length(selected)%%2==1)			# if odd number selected, convert to even for pairing
 	{
	 selected <- as.integer( c(selected, (1+pop$n*runif(1)) ))	# add random index at the end
	}

 pairs <- matrix(selected,ncol=2)
 
 # this part is specific to the binary string representation
 for(i in 1:nrow(pairs))				# carry out the binary string crossover for each pair
 	{
	 pos <- as.integer(1+pop$b * runif(1))		# random crossover point for each pair

	# print(paste(i,":",pairs[i,1],"X",pairs[i,2],"@",pos))
	 v1 <- pop$chr[pairs[i,1],]
	 v2 <- pop$chr[pairs[i,2],]
	 
	 v3 <- c( v1[1:pos-1], v2[pos:pop$b] )
	 v4 <- c( v2[1:pos-1], v1[pos:pop$b] )
	 
	 pop$chr[pairs[i,1],] <- v3
	 pop$chr[pairs[i,2],] <- v4
	}
 return(pop$chr)
}

# Sampling function
# This function takes the current population X and creates X' - the candidate next generation,
# based on some function of the fitness f(X) where f is the target function to be optimised
# That `some function' is the core of the algorithm.
# This is independent of the internal representation of the chromosomes
# Arguments
#	x : current population
#	y : result of applying target function to current population, y = f(x)

selection_linear <- function(x,y)
{
 l <- length(x)
 rankx <- l-rank(y)+1
 #print(x)
 #print(y)
 #print(rankx)
 
 # Using a linear function on `rankx' to get prob(x is selected)
 # q and r are programmer defined parameters
 # For explanation see [Michealewicz,Genetic Algorithms + Data Structures = Evolution Programs,pg
 # 60, in Ch 4]
 q <- 1.5/l
 r <- 1/(l^2-l)
 
 #print(paste("q:",q,"r:",r))
 
 probRank <- q - (rankx-1)*r		# Compute probability for each rank 
 #print(probRank)
 cumulProbRank <- rep(0,l)		# compute cumulative probability
 for(i in 1:l)
 	cumulProbRank[i] <- sum(probRank[1:i])
 #print(cumulProbRank)
 
 xDash <- NULL		
 				
 s <- runif(l)
 for(i in 1:l)		# Draw using Roulette wheel with slots proportional to probability
 	{		# because we are using cumulative probability
	 #print(s)
	 if(s[i]<cumulProbRank[1])
		 {
		 #print(1)
		 #print(paste("selecting",x[1]))
		 xDash <- c(xDash, x[1])
		 }
	 else
	 for(j in 2:l)
	 	{
		 #print(paste(j-1,j))
		 if(s[i]>= cumulProbRank[j-1] && s[i]<cumulProbRank[j])
		 	{
			 #print(paste("selecting",x[j]))
			 xDash <- c(xDash, x[j])
			 break
			}
		}
	}
 #print(xDash)
 #for(i in 1:l)
 #	 print(paste(x[i],":",  "rank=",rankx[i],sum(as.integer(xDash==x[i])),"times" ))

 # The candidate next generation is kept in xDash, return it
 return(xDash)
}

# Main GA function algorithm
# The main limitation of this code is that it is completely dependent on the next few statements
# generate the inpu data and define it in the form of binary arrays and to put it into the data
# structure for which the various functions have been designed. 
# The reason is that my knowledge of R is limited - I don't know how to make user defined types and
# classes. Thats next on the agenda, but for submission's sake I am keeping it at this.
#
############# MOST IMP SAMPLE CODE STARTS  ##########
# xraw <- as.integer(runif(nPop*nBit)<=0.5)	#// Randomly generated population
# x <- matrix(xraw, nrow=nPop, ncol=nBit, byrow=T)	
# probspec=list(domain=c(-1,2), prec=6, desc="Real nos. encoded as binary strings. Domain and precision information required for encoding and decoding")

# population <- list(n=nPop,b=nBit, detail=probspec, chr=x)	# Defines the generic data structure
############# MOST IMP SAMPLE CODE ENDS  ############
#
# Arguments -
# 1)	ftarget:	function to be optimised
# 2)	population:	population list in the form given in sample code
# 3)	nPop:		number of individuals in the population
# 4)	nBit:		how many bits do we use to represent the chromosome - depends on input data
# 5)	pCrossover:	probability of crossover in each generation
# 6)	pMutation:	probability of mutation in each generation
# 7)	nGeneration:	how many generations should the GA run
# 8)	verbose:	verbosity -T/F default F
# )
# )
# )

ga <- function(ftarget, population, nPop, nBit, pCrossover, pMutation, nGeneration,verbose=FALSE)
{
#this is the start of the GA process
x1 <- rep(0,nPop)	# population in real number form
xBest <- rep(0,nGeneration)
yBest <- rep(0,nGeneration)

if(verbose==TRUE)
	{
	 print("best")
	 print(paste("Gen","ymax","            x of ymax",sep="      "))
	}
for(k in 1:nGeneration)
	{
	for(i in 1:nPop)
		{
		 x1[i] <- decode(population,i) 		# Convert from string to real
		}
	#print("x1")
	#print(x1)

	y1 <- sapply(x1,ftarget)	# Evaluate the target function for this 
					# generation, on this population
	#print("y1")
	#print(y1)
	yBest[k] <- max(y1)
	xBest[k] <- x1[ which(y1==max(y1))[1] ]
	#print("best")
	if(verbose==TRUE)
		print(paste(k,yBest[k],xBest[k],sep="       "))
	#return()
	x1Dash <- selection_linear(x1,y1)	# Run the selection operator to get the next generation

	for(i in 1:nPop)
		{
		 population$chr[i,] <- encode(x1Dash[i], population)	# store the new generation in the
									# main structure
		}

	population$chr <- crossover(population, pCrossover)	# apply crossover operator
	population$chr <- mutation(population, pMutation)	# apply mutation operator


	}
if(verbose==TRUE)
	{print("best ever")
	yBestEver <- max(yBest)
	Ever <- which(yBest==yBestEver)[1]
	xBestEver <- xBest[Ever]

	print(paste("Gen","ymax","            x of ymax",sep="      "))
	print(paste(Ever,yBestEver,xBestEver,sep="      "))
	}
# this is the end of the GA process
}


