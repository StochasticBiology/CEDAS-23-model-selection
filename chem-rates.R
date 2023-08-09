library(mombf)

set.seed(1)

# read data
chem.df = read.csv("chem-data.csv")

# cast data in x, y format for convenience
y = chem.df$rate
x = as.matrix(chem.df[,2:5])

# what does a frequentist approach suggest?
summary(lm(y~x))
# only the x^4 term is particularly highlighted. we could explore e.g. AIC further to compare different model structures

# use default prior choice (but can play with this!)
priorCoef <- momprior(tau=0.348)
priorDelta <- modelbbprior(1,1)

# model selection
fit1 <- modelSelection(y=y, x=x, center=FALSE, scale=FALSE,
                       priorCoef=priorCoef, priorDelta=priorDelta)
postProb(fit1) #posterior model probabilities
# here we see that a model with ALL polynomial terms has the highest posterior probability -- but isn't overwhelmingly favoured over the others
# the true generative model involves a reaction in X (x term) which increases rate, X+X (x^2 term) which decreases rate, and X+X+X+X (x^4 term) which increases rate

fit1$margpp #posterior marginal inclusion prob
# here we see fairly strong evidence for the true terms; x^3 is less convincing

coef(fit1) 
# and this is borne out by the parameter estimates
