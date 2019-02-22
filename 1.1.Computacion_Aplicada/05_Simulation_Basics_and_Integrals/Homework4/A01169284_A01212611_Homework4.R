#*******************************************************************************
#* AUTHOR(S) :
#*     Bruno González Soria          (A01169284)
#*     Antonio Osamu Katagiri Tanaka (A01212611)
#*
#* FILENAME :
#*     Homework4.R
#*
#* DESCRIPTION :
#*     Simulations (Ene 19 Gpo 1)
#*     Homework 4
#*
#* START DATE :
#*     21 Feb 2019
#*******************************************************************************


#PART 1#####################################################################
# Integrate the function f(x)=exp(x+x^2) from -2 to 2, using Rieman sums.
plot(labx=(-2,2), function(x){exp(x+x^2)}, type="h")
step = 0.001
f.x = expression (exp(x + x^2))
x = 1.01+step
eval(f.x)*step

integrate(function(x){exp(x + x^2)}, lower=-2, upper=2)
