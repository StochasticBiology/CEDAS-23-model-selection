library(mombf)

# this vignette is basically identical to that in the mombf documentation, just with some additional points added

# set random seed for reproducibility
set.seed(1)

# Simulate data
x <- matrix(rnorm(100*3),nrow=100,ncol=3)
theta <- matrix(c(1,1,0),ncol=1)
y <- x %*% theta + rnorm(100)

#Specify prior parameters
priorCoef <- normalidprior(tau=1)
priorDelta <- modelunifprior()

### Alternative parameter prior
#priorCoef <- momprior(tau=0.348)
### Alternative parameter prior: wider
#priorCoef <- momprior(tau=0.0001)
### Alternative model space prior: 0.5 prior prob for including any covariate
#priorDelta <- modelbinomprior(p=0.5)
### Alternative: Beta-Binomial prior for model space
#priorDelta <- modelbbprior(alpha.p=1,beta.p=1)

#Model selection
fit1 <- modelSelection(y=y, x=x, center=FALSE, scale=FALSE,
                       priorCoef=priorCoef, priorDelta=priorDelta)

# comments removed here -- the challenge is to interpret these outputs!
postProb(fit1) 
fit1$margpp 

# this won't work with the normal prior; other packages like rstanarm and BAS could be used 
# but we can use the momprior here (l.18, for example)
coef(fit1) 

####

# (frequentist) multiple regression
my.lm = lm(y ~ x)
summary(my.lm)

my.lm.1 = lm(y ~ x[,1])
my.lm.2 = lm(y ~ x[,2])
my.lm.3 = lm(y ~ x[,3])
my.lm.12 = lm(y ~ x[,1] + x[,2])
my.lm.23 = lm(y ~ x[,2] + x[,3])
my.lm.13 = lm(y ~ x[,1] + x[,3])
my.lm.123 = lm(y ~ x[,1] + x[,2] + x[,3])

AIC(my.lm.1, my.lm.2, my.lm.3, my.lm.12, my.lm.23, my.lm.13, my.lm.123)


