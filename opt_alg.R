library(numDeriv)

# Steepest Descent
steepestDescent = function(f, x, a, e, maxIter){
  # steepestDescent
  # INPUT
  #      - f objective function
  #      - x inital coordinates
  #      - a step size
  #      - e termination criterion 
  #      - maxIter maximum number of iterations
  
  result <- list(x_opt = x, f_opt = f(x), x_hist = x, f_hist = f(x), iter = 0,
                 iter.time = c(), epoch = c())
  
  #start time
  starttime <- Sys.time()

  currIter <- 1
  finished <- FALSE
  x_old <- x
  while(finished == FALSE){
    StartT <- Sys.time()
    x_new <- lineSearch(f, x_old, x_old - a*grad(f, x_old), 1)
    #x_new <- x_old - a*grad(f, x_old)
    if(currIter <= maxIter & abs(f(x_new)-f(x_old))>e & f(x_new)<f(x_old)){
      x_old <- x_new
      result$x_opt  <- x_new
      result$f_opt  <- f(x_new)
      result$x_hist <- rbind(result$x_hist, x_new)
      result$f_hist <- rbind(result$f_hist, f(x_new))
      result$iter   <- currIter
      result$iter.time   <- rbind(result$iter.time, Sys.time() - StartT)
    }else{
      finished <- TRUE
    }
    currIter <- currIter + 1
  }
  
  #finish time
  result$epoch <- Sys.time() - starttime
  
  return(result)
}

lineSearch = function(f, x0, x1, gridSize){
  # lineSearch
  # INPUT
  #      - f objective function
  #      - x0 starting point
  #      - x1 new point (in terms of gradient descent algorithm)
  #      - gridSize number of points between x0 and x1 
  
  x_best <- x0
  for(i in 1 : gridSize){
    t <- i/gridSize
    x_new <- t*x1 + (1-t)*x0
    if(f(x_best)>f(x_new)){
      x_best <- x_new
    }
  }
  return(x_best)
}

# Simulated Annealing
simulatedAnnealing = function(f, x, alpha, t, delta, maxIter)
{
  # Simulated Annealin Algorithm
  # f - objective function
  # x - inital solution
  # alpha - annealing schedule parameter
  # t - inital temperature
  # delta - neighbourhood radius
  # maxIt - maximum no. of iterations
  
  result = list(x_opt = x, f_opt = f(x), x_hist = x, f_hist = f(x), temperature = t,
                iter.time = 0, epoch = c())
  
  #start time
  starttime <- Sys.time()

  currIter <- 1
  finished <- FALSE
  x_s <- x
  
  while(finished == FALSE){
    
    #epoch
    start.iter <- Sys.time()
    
    # x_c - candidate sol. drawn uniformly fron N(x)
    #u = runif(length(x_s))
    #x_c = x_s + (-delta + 2 * delta * u)
    x_c = x_s - runif(length(x_s), min = -delta, max = delta)
    
    # A - Metropolis activation function
    A = min(1, exp(- (f(x_c) - f(x_s)) / t))
    
    # transition to candidate solution
    if (runif(1) < A)
    {
      x_s <- x_c
    }
    
    # temperature update
    t = alpha * t

    if(currIter<maxIter){
      if(f(x_s)<f(result$x_opt)){
        result$x_opt <- x_s
        result$f_opt <- f(x_s)
      }
      result$x_hist       <- rbind(result$x_hist, x_s)
      result$f_hist       <- rbind(result$f_hist, f(x_s))
      result$temperature  <- rbind(result$temperature, t)
      result$transProb    <- rbind(result$transProb, A)
      result$iter.time    <- rbind(result$iter.time, Sys.time() - start.iter)
    }else{
      finished            <- TRUE
    }
    
    currIter <- currIter + 1
  }
  
  #finish time
  result$epoch <- Sys.time() - starttime
  
  return(result)
}

#Genetic Algorithm
geneticAlgorithm = function(f, x_min, x_max, cel, popSize, pMut, maxIter)
{
  # geneticAlgorithm
  # INPUT
  #      - f objective function
  #      - x_min vector of the minimum values of coordinates
  #      - x_max vector of the maximum values of coordinates
  #      - cel coordinate encryption length 
  #      - popSize size of the population
  #      - pMut probability of single genome mutation
  #      - maxIter number of generations
  
  result <- list(x_opt = c(), f_opt = c(), x_hist= c(), f_hist= c(), f_mean = c(),
                 iter.time = c(), epoch = c())
  
  #start time
  starttime <- Sys.time()
  
  # Check the number of dimensions
  Dim <- length(x_min) 
    
  # Initialize Population
  population <- matrix(NA, nrow = popSize, ncol = cel*Dim)
  for(i in 1 : popSize){
    population[i,] <- runif(cel*Dim)>1
  }
  coordinates <- getCoordinates(population, cel, x_min, x_max)
  
  # Calculate fittness of individuals
  objFunction <- rep(NA, popSize)
  for(i in 1 : popSize){
    objFunction[i] <- f(coordinates[i,])
  }
  
  # Assign the first population to output 
  result$x_opt <- coordinates[which.min(objFunction),]
  result$f_opt <- f(coordinates[which.min(objFunction),])
  
  # The generational loop
  finished <- FALSE
  currIter <- 1
  while(finished == FALSE){
    
    #epoch
    start.iter <- Sys.time()
    
    # Assign the output
    if(currIter <= maxIter){
      if(result$f_opt > f(coordinates[which.min(objFunction),])){
        result$x_opt <- coordinates[which.min(objFunction),]
        result$f_opt <- f(coordinates[which.min(objFunction),])
      }
      result$f_hist <- rbind(result$f_hist, result$f_opt) 
      result$x_hist <- rbind(result$x_hist, coordinates[which.min(objFunction),])
      result$f_mean <- rbind(result$f_mean, mean(objFunction)) 
    }else{
      finished <- TRUE
    }
    
    # Translate binary coding into real values  
    coordinates <- getCoordinates(population, cel, x_min, x_max)
    
    # Calculate fittness of the individuals
    objFunction <- rep(NA, popSize)
      for(i in 1 : popSize){
      objFunction[i] <- f(coordinates[i,])
    }
    rFitt <- min(objFunction)/objFunction # Relative Fittness
    nrFitt <- rFitt / sum(rFitt)          # Relative Normalized (sum up to 1) Fittness; cf. roulette selection
    
    # Selection operator (Roulette wheel)
    selectedPool<- rep(0, popSize)
    for(i in 1 : popSize){
      selectedPool[i] <- which.min(runif(1)>cumsum(nrFitt))
    }
    
    # Crossover operator (for selected pool)
    nextGeneration <- matrix(NA, nrow = popSize, ncol = cel*Dim)
    for(i in 1 : popSize){
      parentId <- round(runif(1,1,popSize))
      cutId <- round(runif(1,1,Dim*cel-1)) # Please, do not exceed the matrix sizes
      # Create offspring
      nextGeneration[i, 1 : cutId] <- population[selectedPool[i], 1 : cutId]
      nextGeneration[i, (cutId + 1) : (Dim*cel)] <- population[selectedPool[parentId], (cutId + 1) : (Dim*cel)]
    }
    
    # Mutation operator
    for(i in 1 : popSize){
      genomeMutId <- which(runif(Dim*cel)>pMut) # Draw the genomes that will mutate
      for(j in 1 : length(genomeMutId)){
        nextGeneration[i, genomeMutId[j]] <- !nextGeneration[i, genomeMutId[j]] 
      }
    }
    
    #iteration time
    if(currIter <= maxIter){
      result$iter.time <- rbind(result$iter.time, Sys.time() - start.iter)
    }
    
    # Replace the old population
    population <- nextGeneration
    currIter <- currIter + 1
  }
  
  #finish time
  result$epoch <- Sys.time() - starttime
  
  return(result)
}
  
  
intbin = function(x){
  # Translate the binary coding to real values numbers
  return(sum(2^(which(rev(x==1))-1)))
}
  
getCoordinates = function(population, cel, x_min, x_max, pMut){
  
  # Check the number of dimensions
  Dim <- length(x_min) 
  
  # Transform the binary coding into coordinates
  coordinates <- matrix(NA, nrow = dim(population)[1], ncol = Dim)
  for(i in 1 : dim(population)[1]){
    for(j in 1 : 2){
      coordinatesTemp <- intbin(population[i, seq(cel*(j-1)+1, j*cel)])
      coordinates[i,j] <- ((x_max[j]-x_min[j])/(2^cel-1))*coordinatesTemp+x_min[j]
                         #^( dł. przedziału  )                           #^(do skrajnego przypadku)
                        #^(   dł. najmniejszego kroku   )              
    }
  }
  return(coordinates)
}

#Particle Swarm Optimization
particleSwarm = function(f, popSize=10, d=2, l.bound=0, u.bound=1, w, c1, c2, maxIter=100, criterion=FALSE)
{
  # psoAlgorithm
  # INPUT
  #      - f objective function
  #      - popSize number of particles 
  #      - d number of variables
  #      - l.bound lower boundary (for initial particles)
  #      - u.bound upper boundary (sup.)
  #      - w inertia weight (vector of length 2 for dynamic inertia, i.e. decreasing over time)
  #      - c1 learning factor (individual experience)
  #      - c2 learning factor (social communication)
  #      - maxIter number of generations
  #      - criterion MaxDistQuick stopping criterion (with threshold)
  
  #store results
  result <- list(x.opt = numeric(d), f.opt = numeric(1), 
                 x.hist = c(), f.hist = numeric(maxIter), 
                 iter.time = c(), epoch = c(), 
                 particles = matrix(numeric(), nrow=popSize, ncol=d))
  
  #pso start time
  starttime <- Sys.time()
  
  #initialize particles
  #particle.matrix <- t(matrix((u.bound - l.bound) * runif(popSize*d) + l.bound, nrow=d, ncol=popSize))
  particle.matrix <- t(matrix(runif(popSize*d, l.bound, u.bound), nrow=d, ncol=popSize))
  
  #initial pbest 
  #f.pbest <- apply(particle.matrix, 1, f)
  pbest <- particle.matrix
  f.pbest <- apply(pbest, 1, f)
  
  #initial gbest
  gbest <- pbest[which.min(f.pbest), 1:d]
  f.gbest <- min(f.pbest)
  #f.gbest <- f(gbest)
  result$x.hist <- gbest
  result$f.hist <- f.gbest
  
  #initial velocity
  velocity <- matrix(runif(popSize*d), nrow=popSize, ncol=d)
  
  #update velocity
  velocity <- max(w) * velocity +
    c1 * runif(1) * (pbest - particle.matrix) +
    c2 * runif(1) * (t(matrix(rep(gbest, popSize), d, popSize)) - particle.matrix)
  
  #new coordinates
  particle.matrix <- particle.matrix + velocity
  
  #first iteration end
  result$iter.time <- Sys.time() - starttime 
  
  #iterate over particles
  for (i in 2:maxIter){
    
    #MaxDistQuick stopping criterion
    if (criterion[1] == TRUE){
      if (i/maxIter >  0.2){
        best.20 <- unname(tail(result$x.hist, floor(i*0.2)))
        max.euclidean.dist <- max(sqrt(rowSums(best.20 - t(matrix(gbest, d, floor(i*0.2))))^2))
          if (max.euclidean.dist < criterion[2]){
            result <- c(result, setNames((i-1), "no.ite"))
            break
          }
      }
    }
    
    #epoch
    start.iter <- Sys.time()
    
    #calculate the objective function for new coordinates
    f.particle.matrix <- apply(particle.matrix, 1, f)
    
    #update pbest
    pbest <- ifelse(matrix(rep(f.particle.matrix, d), nrow=popSize, ncol=d) < 
                      matrix(rep(f.pbest, d), nrow=popSize, ncol=d), particle.matrix, pbest)
    
    #update f.pbest
    f.pbest <- apply(pbest, 1, f)
    
    #update gbest, f.gbest
    gbest <- pbest[which.min(f.pbest), 1:d]
    f.gbest <- f(gbest)
    
    #append results
    result$x.hist <- rbind(result$x.hist, gbest)
    result$f.hist <- append(result$f.hist, f.gbest)
    
    #update velocity
    velocity <- (max(w) - (max(w) - min(w)) * i/maxIter) * velocity +
      c1 * runif(1) * (pbest - particle.matrix) +
      c2 * runif(1) * (t(matrix(rep(gbest, popSize), d, popSize)) - particle.matrix)
    
    #new coordinates
    particle.matrix <- particle.matrix + velocity
    
    #iteration time
    result$iter.time <- append(result$iter.time, Sys.time() - start.iter)
  }
  
  #solution
  result$x.opt <- gbest
  result$f.opt <- f.gbest
  
  #pso finish time
  result$epoch <- Sys.time() - starttime
  
  #end state particles' coordinates
  result$particles <- particle.matrix
  
  return(result)
}  
