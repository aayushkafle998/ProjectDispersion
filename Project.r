#Dhakal_Rabin
#Project1

library(pracma)
library(matlib)
library(stats)
library(Ryacas)
library(base)

#Given Parameters
H <- 5.67                        #Crest to through height, H
d <- 9                             #water depth, d
T <- 7                             #wave period, T
g <- 9.81
Lo <- g*T^2/(2*pi)

#Initial Values of L and a
L <- 60
a <- 2.5

#other given parameters

c <- cosh(2*pi*d/L)

s <- sinh(2*pi*d/L)

C2 <- (3840*c^12-4096*c^10+2592*c^8-1008*c^6+5944*c^4-1830*c^2+147)/(512*(s^10)*(6*c^2-1))

C1 <- (8*c^4-8*c^2+9)/(8*c^4)

B55 <- (192000*c^16-262720*c^14+83680*c^12+20160*c^10-7280*c^8+7160*c^6-1800*c^4-1050*c^2+225)/(12288*(s^10)*(6*c^2-1)*(8*c^4-11*c^2+3))


B35 <- (88128*c^14-208224*c^12+70848*c^10+5400*c^8-21816*c^6+6264*c^4-54*c^2-81)/(12288*(s^12)*(6*c^2-1))

B33 <- 3*(8*c^6+1)/(64*s^6)

beta <- (2*pi)/L

lambda <- beta*a


tolerence <- 0.0001       

x<-as.matrix(c(L,a))     # putting the two variables in a vector

#given two equation

E1 <- x[2]+(-1/2)*H+(3/16)*x[2]^3*x[1]^(-2)*pi^2*(1+8*cosh(2*d*x[1]^(-1)*pi)^6)*sinh(2*d*x[1]^(-1)*pi)^(-6)+(1/782592)*x[2]^5*x[1]^(-4)*pi^4*((-1)+6*cosh(2*d*x[1]^(-1)*pi)^2)^(-1)*sinh(2*d*x[1]^(-1)*pi)^(-12)*(27648*((-3)+2*cosh(2*d*x[1]^(-1)*pi)^2*((-1)+4*cosh(2*d*x[1]^(-1)*pi)^2*(29+(-101)*cosh(2*d*x[1]^(-1)*pi)^2+25*cosh(2*d*x[1]^(-1)*pi)^4+328*cosh(2*d*x[1]^(-1)*pi)^6+(-964)*cosh(2*d*x[1]^(-1)*pi)^8+408*cosh(2*d*x[1]^(-1)*pi)^10)))+5095*(3+(-11)*cosh(2*d*x[1]^(-1)*pi)^2+8*cosh(2*d*x[1]^(-1)*pi)^4)^(-1)*(45+2*cosh(2*d*x[1]^(-1)*pi)^2*((-105)+4*cosh(2*d*x[1]^(-1)*pi)^2*((-45)+179*cosh(2*d*x[1]^(-1)*pi)^2+(-182)*cosh(2*d*x[1]^(-1)*pi)^4+504*cosh(2*d*x[1]^(-1)*pi)^6+2092*cosh(2*d*x[1]^(-1)*pi)^8+(-6568)*cosh(2*d*x[1]^(-1)*pi)^10+4800*cosh(2*d*x[1]^(-1)*pi)^12)))*sinh(2*d*x[1]^(-1)*pi)^2)

E2 <- (-2)*pi+Lo*(2*x[1]^(-1)*pi+(1/16)*x[2]^4*x[1]^(-5)*pi^5*((-1)+6*cosh(2*d*x[1]^(-1)*pi)^2)^(-1)*(147+2*cosh(2*d*x[1]^(-1)*pi)^2*((-915)+4*cosh(2*d*x[1]^(-1)*pi)^2*(743+(-126)*cosh(2*d*x[1]^(-1)*pi)^2+324*cosh(2*d*x[1]^(-1)*pi)^4+(-512)*cosh(2*d*x[1]^(-1)*pi)^6+480*cosh(2*d*x[1]^(-1)*pi)^8)))*sinh(2*d*x[1]^(-1)*pi)^(-10)+x[2]^2*x[1]^(-3)*pi^3*(9+(-8)*cosh(2*d*x[1]^(-1)*pi)^2+8*cosh(2*d*x[1]^(-1)*pi)^4)*sinh(2*d*x[1]^(-1)*pi)^(-4))*tanh(2*d*x[1]^(-1)*pi)



fun <- function(x)
{
  f<-c(x[2]+(-1/2)*H+(3/16)*x[2]^3*x[1]^(-2)*pi^2*(1+8*cosh(2*d*x[1]^(-1)*pi)^6)*sinh(2*d*x[1]^(-1)*pi)^(-6)+(1/782592)*x[2]^5*x[1]^(-4)*pi^4*((-1)+6*cosh(2*d*x[1]^(-1)*pi)^2)^(-1)*sinh(2*d*x[1]^(-1)*pi)^(-12)*(27648*((-3)+2*cosh(2*d*x[1]^(-1)*pi)^2*((-1)+4*cosh(2*d*x[1]^(-1)*pi)^2*(29+(-101)*cosh(2*d*x[1]^(-1)*pi)^2+25*cosh(2*d*x[1]^(-1)*pi)^4+328*cosh(2*d*x[1]^(-1)*pi)^6+(-964)*cosh(2*d*x[1]^(-1)*pi)^8+408*cosh(2*d*x[1]^(-1)*pi)^10)))+5095*(3+(-11)*cosh(2*d*x[1]^(-1)*pi)^2+8*cosh(2*d*x[1]^(-1)*pi)^4)^(-1)*(45+2*cosh(2*d*x[1]^(-1)*pi)^2*((-105)+4*cosh(2*d*x[1]^(-1)*pi)^2*((-45)+179*cosh(2*d*x[1]^(-1)*pi)^2+(-182)*cosh(2*d*x[1]^(-1)*pi)^4+504*cosh(2*d*x[1]^(-1)*pi)^6+2092*cosh(2*d*x[1]^(-1)*pi)^8+(-6568)*cosh(2*d*x[1]^(-1)*pi)^10+4800*cosh(2*d*x[1]^(-1)*pi)^12)))*sinh(2*d*x[1]^(-1)*pi)^2),
    (-2)*pi+Lo*(2*x[1]^(-1)*pi+(1/16)*x[2]^4*x[1]^(-5)*pi^5*((-1)+6*cosh(2*d*x[1]^(-1)*pi)^2)^(-1)*(147+2*cosh(2*d*x[1]^(-1)*pi)^2*((-915)+4*cosh(2*d*x[1]^(-1)*pi)^2*(743+(-126)*cosh(2*d*x[1]^(-1)*pi)^2+324*cosh(2*d*x[1]^(-1)*pi)^4+(-512)*cosh(2*d*x[1]^(-1)*pi)^6+480*cosh(2*d*x[1]^(-1)*pi)^8)))*sinh(2*d*x[1]^(-1)*pi)^(-10)+x[2]^2*x[1]^(-3)*pi^3*(9+(-8)*cosh(2*d*x[1]^(-1)*pi)^2+8*cosh(2*d*x[1]^(-1)*pi)^4)*sinh(2*d*x[1]^(-1)*pi)^(-4))*tanh(2*d*x[1]^(-1)*pi)
   )
  return(f)
}

J <- jacobian(f,c(L,a))

e <- c(0,0)
f <- c(E1,E2)

for (itr in 1:100)
{
  x[1] <- L
  x[2] <- a

  f<-fun(x)

  J <- jacobian(fun,c(L,a))                   #calculation of Jacobian matrix

  e <- -1*inv(J)%*%f

  L <- L+e[1]                  #updatingn value of L and a
  a <- a+e[2]

  if (abs(L-x[1]) < tolerence & abs(a-x[2])< tolerence)  #stopping criterian
  {
    break
  }
  itr <- itr + 1        #calculation of iteration

}

cat("The solution comes in  ",itr," iterations the results are ", L," ,", a)      #printing solution
