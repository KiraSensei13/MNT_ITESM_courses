
#### Exercises

### 1. Use control flow to print all the values from 1 to 100 that are divisible by 3
## Programming Logic:
# I know that all the numbers that are divisible by 3 are also the multiples of 3.
# How do I get a multiple of 3?
# Multiply 3 by any natural number.
# Constraints: Multiples lower or equal than 100 and greater or equal to 1.

maxValue = 100 
x = 3 # First multiple of 3 bigger than 1

while ( x <= maxValue ){  # Condition to satisfy
  print( x )
  x = x + 3  # A product can be decomposed in sums
}

# All done! 

### 2. Use control flow to print all the numbers in the Fibonacci series till its 20th element
## Programming Logic:
# I know that fibonacci uses three values and one operation: last value, present value, and the sum of these two numbers will represent the future value
last = 0
present = 1
future = last + present

# That is the Fibonacci logic. Next, I want the 20th element. At the beginning I only have 2 values: the last and present.
count = 2
maxCount = 20

# Condition to satisfy: count must be less or equal to maxCount
# Lets start

print(last)

while( count <= maxCount){
  print(present)
  # With the past and present values I compute future value
  future = last + present 
  
  # AS we have computed the new value we increase the count by 1
  count = count + 1
  
  # Now we shift the values for the next step in the loop. 
  # The present turns to be the last value. And the future value is the present
  last = present
  present = future
}

# All done!


### 3. Use control flow to print if a numeric variable is a prime number (Intermediate level)
## Programing logic:
# Each time I divide the variable of interest, the max value to query the variable, changes.
# Example value = 11
# if x = 1 then max Value toquery is 11
# Then, if x = 2, the max Value is 11/2
# If x = 3, the max Value is 11/3 
# If x = 4, then x > max value. Therefore, 11 is prime


# The max value is the value we are querying
value = 15

# Initial value
x = 2

# The modulo operator %% it is used to know the residue of a division. Example. x %% y. If it returns 0, x is a multiple of y.
# 4 %% 2 returns 0
# 5 %% 2 returns 1

# Lets use a logical value to see if it is a prime number
prime = TRUE

while ( x < value){
  # If the modulo returns 0. Then the value is not a prime number
  if(value %% x == 0){
    prime = FALSE
  }
  # We update x
  x = x + 1
}

# Finally we print the answer
if ( prime  ){
  print("It is a prime number")
}else{
  print("Not a prime number")
}


#### Extra exercises
### 1. Compare three variables and print the variable with the biggest value. 
x = 5
y = 20
z = 3

if( x > y){
  if( x > z){
    print(x)
  }else{
    print(z)
  }
}else{
  if( y > z){
    print(y)
  }else{
    print(z)
  }
}

### 2. Receive a vector variable and print how many elements are in the vector. (Do not use the length function)
bigvector = c("a", "b", "c")
count = 0

for( i in bigvector){
  count = count + 1
}

print(count)

### 3. Receive a character vector with unitary letter. eg. c("a", "b", "c"). Print TRUE if at least 1 of the characters is a vocal.

vocals = c("a", "e", "i", "o", "u")
bigvector = c("z", "b", "c")
vocal = FALSE

for ( i in bigvector){
  for( j in vocals){
    if(i == j){
      vocal = TRUE
    }
  }
}

print(vocal)

### 4. Receive a vector and a logical variable, if the logical variable is TRUE print the sum of all the elements inside the vector. If the logical value is FALSE return the product of all the elements inside the vector

# R functions to use:
?sum
?prod

bigvector = 1:5
operationToDo = FALSE

if(operationToDo == TRUE){
  print( sum(bigvector) )
}else{
  print( prod(bigvector) )
}

### 5. Receive a vector and print it from the last element to the first one. E.g. c(1,2,3,4) will print:
# 4
# 3
# 2
# 1

bigvector = c(1,2,3,4)
maxValue = length(bigvector)

for( i in 1:maxValue){
  print( bigvector[ maxValue - i + 1 ]  )
}














