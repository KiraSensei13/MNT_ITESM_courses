# This is a comment

# initial values
x = 10
y = 20

# equation
magnitude=sqrt(x^2+y^2)
v.x = x/magnitude
v.y = y/magnitude

# use class() to return the data type

######################################

# Exercice - Part1
# Compute. log(((3+2)*5)+6)
a = 3 + 2;
b = a * 5;
c = b + 6;
res = log(c);
print(res)

# Create a vector and print the 4th element
v = c(0,1,2,3,4);
print(v[4])


A = "sin(3),pi,4,2.3"
A
B = unlist(strsplit(A,","), "\\d+", as.numeric)
B
f<-function(x) {
  eval(parse(text=B))
  strapply(B, "\\d+", as.numeric)
}
f(B)

# Create a vector of 5 characters
v = c("a","b","c","d","e")

######################################

# Matrices
x = c(1,3,8)
y = c(1,5,3)
xy = rbind(x,y)
xy
xy[2,2]

######################################

n = 1000
randomValues = rnorm(n, mean = 10, sd = 3)
hist(randomValues)

######################################

# print all vaues from 1 to 100 that are divisible by 3
data = 1:100
for (num in data) {
  # check if it is divisible by 3
  if (num%%3 == 0) {
    print(num)
  }
}

######################################

#print all the Fibo numbers up to its 20th element
#calculate the fibonacci sequence recursively
#fibo(n)=fibo(n-1)+fibo(n-2)
fibonacci <- function(n) {
  if(n <= 1) {
    return(n)
  } else {
    return(fibonacci(n-1) + fibonacci(n-2))
  }
}

# print fibonacci sequence up to n-terms
printFibonacci <- function(nterms) {
  # check if the number of terms is valid
  if(nterms <= 0) {
    print("nterms must be greater than 0")
  } else {
    # print Fibonacci sequence
    for(i in 1:(nterms)) {
      print(fibonacci(i))
    }
  }
}

printFibonacci(20)

######################################

#print the entered number is prime
isPrime <- function(num) {
  # Asume all numbers are not prime
  is_prime = FALSE
  # prime numbers are greater than 1
  if(num > 1) {
    is_prime = TRUE
    # check for factors, from 2 to mun-1
    for(i in 2:(num-1)) {
      # shall not be dividers from 2 to mun-1 to be prime
      if ((num %% i) == 0) {
        is_prime = FALSE
        break
      }
    }
  }
  # 2 is a prime number
  if(num == 2)    is_prime = TRUE
  # print the result
  if(is_prime) {
    #print(paste(num,"is a prime number"))
    return(TRUE)
  } else {
    #print(paste(num,"is not a prime number"))
    return(FALSE)
  }
}

isPrime(1)
isPrime(2)
isPrime(3)
isPrime(5)
isPrime(6)
isPrime(7)
isPrime(8)
isPrime(10)
isPrime(11)
isPrime(12)
isPrime(16)
isPrime(18)
isPrime(19)
isPrime(20)

######################################
