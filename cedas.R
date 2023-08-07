library(ggplot2)
library(gridExtra)

### simple model selection quiz

set.seed(1)  # seed for reproducibility
sf = 2       # resolution for PNG outputs

# construct and subsample dataset
npoint = 60
x = (1:npoint)/2
y = 0.01*x**2+rnorm(npoint, sd=0.3)
df.all = data.frame(x=x, y=y)
df = df.all[(npoint/4):(3*npoint/4),]

ggplot(df, aes(x=x,y=y)) + geom_point()

# fit various polynomial models; use summary stats for plot titles
my.lm.1 = lm(y ~ x, data=df)
my.lm.2 = lm(y ~ I(x) + I(x**2), data=df)
my.lm.9 = lm(y ~ I(x) + I(x**2) + I(x**3) + I(x**4) + I(x**5) + I(x**6) + I(x**7) + I(x**8) + I(x**9), data=df)
title.1 = paste(c("y ~ x\nR² = ", signif(summary(my.lm.1)$r.squared, digits=3)), collapse="")
title.2 = paste(c("y ~ x + x²\nR² = ", signif(summary(my.lm.2)$r.squared, digits=3)), collapse="")
title.9 = paste(c("y ~ x + x² + ... + x⁹\nR² = ", signif(summary(my.lm.9)$r.squared, digits=3)), collapse="")

# get predictions of fitted models
df$predict.1 = predict(my.lm.1)
df$predict.2 = predict(my.lm.2)
df$predict.9 = predict(my.lm.9)

# construct plots
g.1 = ggplot(df) + geom_point(aes(x=x,y=y)) + geom_line(aes(x=x,y=predict.1), color="red") + ggtitle(title.1)
g.2 = ggplot(df) + geom_point(aes(x=x,y=y)) + geom_line(aes(x=x,y=predict.2), color="green") + ggtitle(title.2)
g.3 = ggplot(df) + geom_point(aes(x=x,y=y)) + geom_line(aes(x=x,y=predict.9), color="blue") + ggtitle(title.9)

# output this set
png("cedas-model-1.png", width=600*sf, height=300*sf, res=72*sf)
grid.arrange(g.1,g.2,g.3,nrow=1)
dev.off()

# add likelihood values
title.1a = paste(c("y ~ x\nR² = ", signif(summary(my.lm.1)$r.squared, digits=3), "\nlog L = ", signif(logLik(my.lm.1), digits=3)), collapse="")
title.2a = paste(c("y ~ x + x²\nR² = ", signif(summary(my.lm.2)$r.squared, digits=3), "\nlog L = ", signif(logLik(my.lm.2), digits=3)), collapse="")
title.9a = paste(c("y ~ x + x² + ... + x⁹\nR² = ", signif(summary(my.lm.9)$r.squared, digits=3), "\nlog L = ", signif(logLik(my.lm.9), digits=3)), collapse="")

# output with likelihood values
g.1a = ggplot(df) + geom_point(aes(x=x,y=y)) + geom_line(aes(x=x,y=predict.1), color="red") + ggtitle(title.1a)
g.2a = ggplot(df) + geom_point(aes(x=x,y=y)) + geom_line(aes(x=x,y=predict.2), color="green") + ggtitle(title.2a)
g.3a = ggplot(df) + geom_point(aes(x=x,y=y)) + geom_line(aes(x=x,y=predict.9), color="blue") + ggtitle(title.9a)

png("cedas-model-3.png", width=600*sf, height=300*sf, res=72*sf)
grid.arrange(g.1a,g.2a,g.3a,nrow=1)
dev.off()

# compare predictions to withheld data
df.all$predict.1 = predict(my.lm.1, newdata=df.all)
df.all$predict.2 = predict(my.lm.2, newdata=df.all)
df.all$predict.9 = predict(my.lm.9, newdata=df.all)

# add AIC to plot titles
title.1b = paste(c("y ~ x\nR² = ", signif(summary(my.lm.1)$r.squared, digits=3), "\nlog L = ", signif(logLik(my.lm.1), digits=3), "\nAIC = ", signif(AIC(my.lm.1), digits=3)), collapse="")
title.2b = paste(c("y ~ x + x²\nR² = ", signif(summary(my.lm.2)$r.squared, digits=3), "\nlog L = ", signif(logLik(my.lm.2), digits=3), "\nAIC = ", signif(AIC(my.lm.2), digits=3)), collapse="")
title.9b = paste(c("y ~ x + x² + ... + x⁹\nR² = ", signif(summary(my.lm.9)$r.squared, digits=3), "\nlog L = ", signif(logLik(my.lm.9), digits=3), "\nAIC = ", signif(AIC(my.lm.9), digits=3)), collapse="")

# plot this set with withheld data
g.1b = ggplot(df.all) + geom_point(aes(x=x,y=y)) + geom_line(aes(x=x,y=predict.1), color="red") + ggtitle(title.1) + ylim(-2,12)
g.2b = ggplot(df.all) + geom_point(aes(x=x,y=y)) + geom_line(aes(x=x,y=predict.2), color="green") + ggtitle(title.2) + ylim(-2,12)
g.3b = ggplot(df.all) + geom_point(aes(x=x,y=y)) + geom_line(aes(x=x,y=predict.9), color="blue") + ggtitle(title.9) + ylim(-2,12)

png("cedas-model-2.png", width=600*sf, height=300*sf, res=72*sf)
grid.arrange(g.1b,g.2b,g.3b,nrow=1)
dev.off()

# plot without the giveaway
g.1c = ggplot(df.all) + geom_point(aes(x=x,y=y)) + geom_line(aes(x=x,y=predict.1), color="red") + ggtitle(title.1b) + ylim(-2,12)
g.2c = ggplot(df.all) + geom_point(aes(x=x,y=y)) + geom_line(aes(x=x,y=predict.2), color="green") + ggtitle(title.2b) + ylim(-2,12)
g.3c = ggplot(df.all) + geom_point(aes(x=x,y=y)) + geom_line(aes(x=x,y=predict.9), color="blue") + ggtitle(title.9b) + ylim(-2,12)

png("cedas-model-4.png", width=600*sf, height=300*sf, res=72*sf)
grid.arrange(g.1c,g.2c,g.3c,nrow=1)
dev.off()

#########

# just a simple illustrative scatterplot for Bayesian inference
set.seed(1)

x = (1:10)/10
y = 0.5*x+rnorm(10, sd=0.2)

png("scatter.png", width=200*sf, height=200*sf, res=72*sf)
ggplot(data.frame(x=x,y=y), aes(x=x,y=y)) + geom_point()
dev.off()

####################

# bayesian model selection
set.seed(1)

# inefficient for loops and data frame binding -- will take a while
npoint = 10

# initialise data frames
res = data.frame()
src = data.frame()
lik.1.df = lik.2.df = data.frame()

# loop through beta2, the coefficient of x^2 in the synthetic data generation
for(beta2 in (0:5)/5) {
  # construct a synthetic dataset 
  x = (1:npoint)/(npoint/2)
  y = beta2*x**2 + 0.5*x +rnorm(npoint, sd=0.001)
  src = rbind(src, data.frame(beta2=beta2, x=x, y=y))
  
  # initialise evidence integrals
  total.lik1 = total.lik2 = 0
  # loop through prior on beta, the x coefficient
  for(tbeta in (0:100)/100) {
    tlik2 = 0
    # loop through prior on beta2, the x^2 coefficient
    for(tbeta2 in (0:100)/100) {
      # get likelihood for x^2 model
      d2 = y-(tbeta2*(x**2) + tbeta*x)
      lik2 = sum(dnorm(d2, mean=0, sd=0.001))
      # add likelihood*prior
      tlik2 = tlik2 + lik2*0.01*0.01
      # record likelihood at this particular prior point for a subset of examples
      if(beta2 == 0 | beta2 == 1) {
        lik.2.df = rbind(lik.2.df, data.frame(beta2=beta2, tbeta=tbeta, tbeta2=tbeta2, lik=lik2))
      }
    }
    # get likelihood for x model
    d1 = y-(tbeta*x)
    lik1 = sum(dnorm(d1, mean=0, sd=0.001))
    # record likelihood at this particular prior point for a subset of examples
    if(beta2 == 0 | beta2 == 1) {
    lik.1.df = rbind(lik.1.df, data.frame(beta2=beta2, tbeta=tbeta, tbeta2=0, lik=lik1))
    }
    # add likelihood*prior for model x
    total.lik1 = total.lik1 + lik1*0.01
    # add this integrated strip for model x^2
    total.lik2 = total.lik2 + tlik2
  }
 
  # add these results to our data frame
  res = rbind(res, data.frame(beta2=beta2, model=1, evidence=total.lik1))
              res = rbind(res, data.frame(beta2=beta2, model=2, evidence=total.lik2))
}

# plot the behaviour of model evidence as synthetic data changes
r.1 = ggplot(src, aes(x=x,y=y)) + geom_point() + facet_wrap(~beta2, nrow=1)
r.2 = ggplot(res, aes(x=factor(model),y=evidence)) + geom_col()  +  facet_wrap(~beta2, nrow=1)
png("cedas-model-sel-x2.png", width=600*sf, height=200*sf, res=72*sf)
grid.arrange(r.1, r.2, nrow=2)
dev.off()

# plot likelihood surfaces over prior support
lik.1.plot1 = ggplot(lik.1.df[lik.1.df$beta2==0,], aes(x=tbeta, y=tbeta2, fill=lik)) + geom_tile()
lik.2.plot1 = ggplot(lik.2.df[lik.2.df$beta2==0,], aes(x=tbeta, y=tbeta2, fill=lik)) + geom_tile()

lik.1.plot2 = ggplot(lik.1.df[lik.1.df$beta2==1,], aes(x=tbeta, y=tbeta2, fill=lik)) + geom_tile()
lik.2.plot2 = ggplot(lik.2.df[lik.2.df$beta2==1,], aes(x=tbeta, y=tbeta2, fill=lik)) + geom_tile()

grid.arrange(lik.1.plot1, lik.2.plot1, lik.1.plot2, lik.2.plot2, nrow=2)
