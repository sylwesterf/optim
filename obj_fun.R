# Objective functions:

# RASTRIGIN FUNCTION
rastr <- function(xx)
{
  d <- length(xx)
	
  sum <- sum(xx^2 - 10*cos(2*pi*xx))
	
  y <- 10*d + sum
  return(y)
}

# SCHAFFER FUNCTION 2
schaffer2 <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
	
  fact1 <- (sin(x1^2-x2^2))^2 - 0.5
  fact2 <- (1 + 0.001*(x1^2+x2^2))^2
	
  y <- 0.5 + fact1/fact2
  return(y)
}

# ACKLEY FUNCTION
ackley <- function(xx, a=20, b=0.2, c=2*pi)
{
  d <- length(xx)
  
  sum1 <- sum(xx^2)
  sum2 <- sum(cos(c*xx))

  term1 <- -a * exp(-b*sqrt(sum1/d))
  term2 <- -exp(sum2/d)

  y <- term1 + term2 + a + exp(1)
  return(y)
}

# SPHERE FUNCTION
spheref <- function(xx)
{
  sum <- sum(xx^2)
	
  y <- sum
  return(y)
}

# ROSENBROCK FUNCTION
rosen <- function(xx)
{
  d <- length(xx)
  xi <- xx[1:(d-1)]
  xnext <- xx[2:d]
	
  sum <- sum(100*(xnext-xi^2)^2 + (xi-1)^2)
	
  y <- sum
  return(y)
}

# BEALE FUNCTION
beale <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
	
  term1 <- (1.5 - x1 + x1*x2)^2
  term2 <- (2.25 - x1 + x1*x2^2)^2
  term3 <- (2.625 - x1 + x1*x2^3)^2
	
  y <- term1 + term2 + term3
  return(y)
}

# GOLDSTEIN-PRICE FUNCTION
goldpr <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
	
  fact1a <- (x1 + x2 + 1)^2
  fact1b <- 19 - 14*x1 + 3*x1^2 - 14*x2 + 6*x1*x2 + 3*x2^2
  fact1 <- 1 + fact1a*fact1b
	
  fact2a <- (2*x1 - 3*x2)^2
  fact2b <- 18 - 32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2
  fact2 <- 30 + fact2a*fact2b
	
  y <- fact1*fact2
  return(y)
}

# BOOTH FUNCTION
booth <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
	
  term1 <- (x1 + 2*x2 - 7)^2
  term2 <- (2*x1 + x2 - 5)^2
	
  y <- term1 + term2
  return(y)
}

# BUKIN FUNCTION N. 6
bukin6 <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
	
  term1 <- 100 * sqrt(abs(x2 - 0.01*x1^2))
  term2 <- 0.01 * abs(x1+10)
	
  y <- term1 + term2
  return(y)
}

# MATYAS FUNCTION
matya <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
	
  term1 <- 0.26 * (x1^2 + x2^2)
  term2 <- -0.48*x1*x2
	
  y <- term1 + term2
  return(y)
}

# LEVY FUNCTION N. 13
levy13 <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
	
  term1 <- (sin(3*pi*x1))^2
  term2 <- (x1-1)^2 * (1+(sin(3*pi*x2))^2)
  term3 <- (x2-1)^2 * (1+(sin(2*pi*x2))^2)
	
  y <- term1 + term2 + term3
  return(y)
}

# HUMP FUNCTION
hump <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
	
  term1 <- 4 * x1^2
  term2 <- -2.1 * x1^4
  term3 <- x1^6 / 3
  term4 <- x1*x2
  term5 <- -4 * x2^2
  term6 <- 4 * x2^4
	
  y <- 1.0316285 + term1 + term2 + term3 + term4 + term5 + term6
  return(y)
}

# EASOM FUNCTION
easom <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
	
  fact1 <- -cos(x1)*cos(x2)
  fact2 <- exp(-(x1-pi)^2-(x2-pi)^2)
	
  y <- fact1*fact2
  return(y)
}

# CROSS-IN-TRAY FUNCTION
crossit <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
	
  fact1 <- sin(x1)*sin(x2)
  fact2 <- exp(abs(100 - sqrt(x1^2+x2^2)/pi))
	
  y <- -0.0001 * (abs(fact1*fact2)+1)^0.1
  return(y)
}

# EGGHOLDER FUNCTION
egg <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
	
  term1 <- -(x2+47) * sin(sqrt(abs(x2+x1/2+47)))
  term2 <- -x1 * sin(sqrt(abs(x1-(x2+47))))
	
  y <- term1 + term2
  return(y)
}

# HOLDER TABLE FUNCTION
holder <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
	
  fact1 <- sin(x1)*cos(x2)
  fact2 <- exp(abs(1 - sqrt(x1^2+x2^2)/pi))
	
  y <- -abs(fact1*fact2)
  return(y)
}

# MCCORMICK FUNCTION
mccorm <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
	
  term1 <- sin(x1 + x2)
  term2 <-(x1 - x2)^2
  term3 <- -1.5*x1
  term4 <- 2.5*x2
	
  y <- term1 + term2 + term3 + term4 + 1
  return(y)
}

# SCHAFFER FUNCTION N. 4
schaffer4 <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
	
  fact1 <- cos(sin(abs(x1^2-x2^2))) - 0.5
  fact2 <- (1 + 0.001*(x1^2+x2^2))^2
	
  y <- 0.5 + fact1/fact2
  return(y)
}

# STYBLINSKI-TANG FUNCTION
stybtang <- function(xx)
{
  sum <- sum(xx^4 - 16*xx^2 + 5*xx)
	
  y <- sum/2
  return(y)
}

# Source:

# Authors: Sonja Surjanovic, Simon Fraser University
#          Derek Bingham, Simon Fraser University
# Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
#
# Copyright 2013. Derek Bingham, Simon Fraser University.
#
# THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
# FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
# derivative works, such modified software should be clearly marked.
# Additionally, this program is free software; you can redistribute it 
# and/or modify it under the terms of the GNU General Public License as 
# published by the Free Software Foundation; version 2.0 of the License. 
# Accordingly, this program is distributed in the hope that it will be 
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# For function details and reference information, see:
# http://www.sfu.ca/~ssurjano/
