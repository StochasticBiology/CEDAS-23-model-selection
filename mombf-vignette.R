library(mombf)

# this vignette is basically identical to that in the mombf documentation, just with some additional points added

# set random seed for reproducibility
set.seed(1)

# Simulate data
x <- matrix(rnorm(100*3),nrow=100,ncol=3)
theta <- matrix(c(1,1,0),ncol=1)
y <- x %*% theta + rnorm(100)

#Specify prior parameters
priorCoef <- momprior(tau=0.348)
priorDelta <- modelunifprior()

### Alternative parameter prior: wider
#priorCoef <- momprior(tau=10)
### Alternative parameter prior: thinner
#priorCoef <- momprior(tau=0.001)
### Alternative parameter prior: normal
#priorCoef <- normalidprior(tau=1)
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
coef(fit1) 