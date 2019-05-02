# An example of the averaing metric and the proof of slope estimates
# 
# XJ
# 4/30/2019

library(TeachingDemos)
char2seed("summing")

# simulate three ecosystem functions
F1 <- rnorm(100, 0, 1)
F2 <- rnorm(100, 0, 1)
F3 <- rnorm(100, 0, 1)

# calcualte the averaging multifunctionality metric
Y <- (F1 + F2 + F3)/3

# inspect the variance of the averaging metric 
# and how it comes from the single functions
var(Y)
(var(F1) + var(F2) + var(F3) + 2*cov(F1, F2) + 2*cov(F1, F3) + 2*cov(F2, F3))/9
(sd(F1)^2 + sd(F2)^2 + sd(F3)^2 + 2*cov(F1, F2) + 2*cov(F1, F3) + 2*cov(F2, F3))/9

# simulate biodiversity
div <- rpois(100, 20)

# inspect the covariance and how it comes from the single functions
cov(Y, div) 
(cov(F1, div) + cov(F2, div) + cov(F3, div))/3

# calcualte the slope between bidiversity and the averaging metric
# by the linear model
summary(lm(Y ~ div))
# by the formula
cov(Y, div)/var(div)
# more ways to calculate the slope
cov(Y, div)/cov(div, div)
cor(Y, div)*sd(Y)/sd(div)
# inspect how the slope comes from single functions
(cov(F1, div)/cov(div, div) + cov(F2, div)/cov(div, div) + cov(F3, div)/cov(div, div))/3
(cor(F1, div)*sd(F1)/sd(div) + cor(F2, div)*sd(F2)/sd(div) + cor(F3, div)*sd(F3)/sd(div))/3
