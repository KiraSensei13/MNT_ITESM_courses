##############################################################################
# Statistics
# Generate a normal distribution
normal <- rnorm(100000, mean= 0, sd =50)
plot(density(normal), main= "Normal distribution", type="h")
plot(ecdf(normal), main = "Normal cdf")

?rnorm
?plot
?density
?ecdf

# Multinomial distribution
multi <- rmultinom(1000, 1, c(0.1, 0.2, 0.7))
multi[ , 1:10]

?rmultinom
?c

# Sum how many times an event ocurred
multi.appeared <- apply(multi,1,sum)

?apply
?rmultinom
?dimnames

###########################################################################

## Compute row and column sums for a matrix:
x <- cbind(x1 = 3, x2 = c(4:1, 2:5))
dimnames(x)[[1]] <- letters[1:8]
apply(x, 2, mean, trim=0.2)
col.sums <- apply(x, 2, sum)
row.sums <- apply(x, 1, sum)
rbind(cbind(x, Rtot = row.sums), Ctot = c(col.sums, sum(col.sums)))

stopifnot( apply(x, 2, is.vector))

## Sort the columns of a matrix
apply(x, 2, sort)

## keeping named dimnames
names(dimnames(x)) <- c("row", "col")
x3 <- array(x, dim = c(dim(x),3),
            dimnames = c(dimnames(x), list(C = paste0("cop.",1:3))))
identical(x,  apply( x,  2,  identity))
identical(x3, apply(x3, 2:3, identity))

##- function with extra args:
cave <- function(x, c1, c2) c(mean(x[c1]), mean(x[c2]))
apply(x, 1, cave,  c1 = "x1", c2 = c("x1","x2"))

ma <- matrix(c(1:4, 1, 6:8), nrow = 2)
ma
apply(ma, 1, table)  #--> a list of length 2
apply(ma, 1, stats::quantile) # 5 x n matrix with rownames

stopifnot(dim(ma) == dim(apply(ma, 1:2, sum)))

## Example with different lengths for each call
z <- array(1:24, dim = 2:4)
zseq <- apply(z, 1:2, function(x) seq_len(max(x)))
zseq         ## a 2 x 3 matrix
typeof(zseq) ## list
dim(zseq) ## 2 3
zseq[1,]
apply(z, 3, function(x) seq_len(max(x)))
# a list without a dim attribute


###########################################################################



# auxiliary variable to plot the distribution
x <- c ( rep(1,multi.appeared[1]), rep(2,multi.appeared[2]), rep(3,multi.appeared[3]))
hist(x, main = "Multinomial pmf")
plot(ecdf(x), main = "Multinomial cmf")


# Central limit theorem

N <- 10000
M <- 100
p = 0.3

binomial <- rbinom(N, M, p)
hist(binomial, freq=F)
lines(density(binomial))

# Binomial mean for one trial
binomial.mean = p
binomial.variance = p*(1-p)

normal <- rnorm(N, binomial.mean*M , sqrt(binomial.variance*M))
lines(density((normal)), col = "red")

## Statistical tests
?t.test
?wilcox.test #paired = F ~ Mann-Whitney U
?cor.test # method=kendall


# chi-squared test
cool.kids <- matrix(c(49,50,69, 24,35,38,19,22,28), byrow=T, ncol = 3)
chisq.test(cool.kids)

cool.kids.2 <- matrix(c(57,87,24,50,42,6,42,22,5), byrow=T, ncol = 3)
chisq.test(cool.kids.2)


##############################################################################
# Gaussian Mixture Model
# install.packages("mixtools")

library(mixtools)
data(faithful)

wait = faithful$waiting
hist(wait)

# normalmixEM K indica la cantidad de estados ocultos que queremos usar
# EM indica el uso del algoritmo de esperanza maximización
mixmdl <- normalmixEM(wait, k= 2)
summary(mixmdl) # lambda significa la proporción que hay en cada grupo

# mixtools sobreescribe la función plot
# Con el valor 2 indicamos que queremos imprimir las curvas de densidad
plot(mixmdl,whichplots = 2) 
lines(density(wait), lty=2, lwd=2)

# Análisis bidimensional
plot(faithful$eruptions, faithful$waiting)
mixmv <- mvnormalmixEM(faithful, k = 2)
plot(mixmv , which = 2)


### HMM
# install.packages(c("Rcpp","RcppArmadillo", "RcppHMM"))
library(RcppHMM)

## Values for a hidden Markov model with continuous observations                          
# Number of hidden states = 3
# Univariate gaussian mixture model

N = c("Low","Normal", "High")
A <- matrix(c(0.5, 0.3,0.2,
              0.2, 0.6, 0.2,
              0.1, 0.3, 0.6),
            ncol= length(N), byrow = TRUE)

Mu <- matrix(c(0, 50, 100), ncol = length(N))
Sigma <- array(c(144, 400, 100), dim = c(1,1,length(N)))
Pi <- rep(1/length(N), length(N))

HMM.cont.univariate <- verifyModel(list( "Model"="GHMM", 
                                         "StateNames" = N,
                                         "A" = A, 
                                         "Mu" = Mu, 
                                         "Sigma" = Sigma, 
                                         "Pi" = Pi))

# Data simulation
set.seed(100)
length <- 100
seqs <- 50

# Multiple sequences to be evaluated
observationSequences<- array(0, dim = c(1, length, seqs) )
for(i in 1:seqs){
  Y <- generateObservations(HMM.cont.univariate , length)$Y
  observationSequences[,,i] <- Y
}

dim(observationSequences)

hist(observationSequences, breaks = 50)

# New model random initialization
# Model to be trained
set.seed(1000)
newModel <- initGHMM(3) 

newModel <- learnEM(newModel,
                    observationSequences,
                    iter= 50, 
                    delta = 1E-5,
                    print = FALSE)


print(newModel)   

