# this fairly directly taken from https://github.com/StochasticBiology/odna-loss
# with corresponding paper Giannakis, K., Arrowsmith, S.J., Richards, L., Gasparini, S., Chustecki, J.M., Røyrvik, E.C. and Johnston, I.G., 2022. Evolutionary inference across eukaryotes identifies universal features shaping organelle gene retention. Cell Systems, 13(11), pp.874-884.

library(mombf)
library(ggplot2)
library(ggrepel)

# read in data
mt.df = read.csv("mt-data.csv")

# settings, and data structure, for bayesian model selection
nf = ncol(mt.df)
priorCoef <- momprior(tau=0.348)
priorDelta <- modelbbprior(1,1)
mt.x = as.matrix(mt.df[,4:12])
mt.y = mt.df$Index

# frequentist picture -- always pays to take a quick look
mt.freq <- lm(mt.y ~ mt.x)

# perform model selection and store results
mt.fit <- modelSelection(y=mt.y, x=mt.x, priorCoef=priorCoef, priorDelta=priorDelta)

# grab model and parameter inclusion posteriors as before
mt.pp = postProb(mt.fit)
mt.fit$margpp

## prepare a plot! this can certainly be done more smoothly -- a bit hacky here
# this pulls the names of our predictors
feature.labels = colnames(mt.x)
# give each model structure a string label
mt.pp$modelid = as.character(mt.pp$modelid)
# get the direction (positive or negative) of each parameter in each model
mt.signs = ifelse(coef(mt.fit)[,1] < 0, "-", "+")

# convert model labels to human-readable format using names of features from dataframe, and get top 6
mt.results = data.frame(Model=NULL,pp=NULL)
maxresults = 6 
for(i in 1:maxresults) {
  this.label = NULL
  this.result = strsplit(mt.pp$modelid[i], ",")[[1]]
  # just build a character string reflecting the structure of this model
  if(length(this.result) > 0) {
    for(j in 1:length(this.result)) {
      this.trait = as.numeric(this.result[j])
      this.label = c(this.label, paste(c(feature.labels[this.trait], mt.signs[this.trait]), collapse=""))
    }
  }
  # add to a growing data frame for plotting
  mt.results = rbind(mt.results, data.frame(Model=paste(this.label,collapse="\n"), pp=mt.pp$pp[i]))
}


# produce subplot
mt.model.sel.plot = ggplot(mt.results, aes(x=factor(Model, levels=Model),y=mt.pp$pp[1:maxresults])) +
  geom_col(fill="#FF8888", colour="#000000") +
  theme_light() + #theme(axis.text.x = element_text(size = 10, margin = unit(c(t = -4, r = 0, b = 0, l = 0), "cm"))) +
  xlab("Model structure") + ylab("Posterior probability")

mt.model.sel.plot

## now take a look at the plastid data
pt.df = read.csv("pt-data.csv")

# settings, and data structure, for bayesian model selection
nf = ncol(pt.df)
priorCoef <- momprior(tau=0.348)
priorDelta <- modelbbprior(1,1)
pt.x = as.matrix(pt.df[,4:12])
pt.y = pt.df$Index

# frequentist
pt.freq <- lm(pt.y ~ pt.x)

# perform model selection and grab posterior probabilities and signs of coefficients
pt.fit <- modelSelection(y=pt.y, x=pt.x, priorCoef=priorCoef, priorDelta=priorDelta)

# familiar outputs
pt.pp = postProb(pt.fit)
pt.fit$margpp

# construct labels for models as above
feature.labels = colnames(pt.x)
pt.pp$modelid = as.character(pt.pp$modelid)
pt.signs = ifelse(coef(pt.fit)[,1] < 0, "-", "+")

# convert model labels to human-readable format using names of features from dataframe, and get top 6
pt.results = data.frame(Model=NULL,pp=NULL)
maxresults = 6 
for(i in 1:maxresults) {
  this.label = NULL
  this.result = strsplit(pt.pp$modelid[i], ",")[[1]]
  if(length(this.result) > 0) {
    for(j in 1:length(this.result)) {
      this.trait = as.numeric(this.result[j])
      this.label = c(this.label, paste(c(feature.labels[this.trait], pt.signs[this.trait]), collapse=""))
    }
  }
  pt.results = rbind(pt.results, data.frame(Model=paste(this.label,collapse="\n"), pp=pt.pp$pp[i]))
}


# produce subplot for final figure
pt.model.sel.plot = ggplot(pt.results, aes(x=factor(Model, levels=Model),y=pt.pp$pp[1:maxresults])) +
  geom_col(fill="#8888FF", colour="#000000") +
  theme_light() + #theme(axis.text.x = element_text(size = 10, margin = unit(c(t = -4, r = 0, b = 0, l = 0), "cm"))) +
  xlab("Model structure") + ylab("Posterior probability")



# initialise data frame for validation results
mt.test.acc = mt.training.acc = mt.cross.acc = pt.test.acc = pt.training.acc = pt.cross.acc = NULL
# number of training-test splits
nsamp = 100
# test proportion (50% training, 50% test)
testpropn = 0.5
# what correlation to report?
cor.method = "spearman"
# loop over splits
for(i in 1:nsamp) {
  mt.sample.n = nrow(mt.df)
  # construct mtDNA training set
  mt.training.refs = sample(seq(from=1,to=mt.sample.n), size=mt.sample.n*testpropn)
  mt.training.set = mt.df[mt.training.refs,]
  # construct mtDNA test set
  mt.test.set = mt.df[-mt.training.refs,]
  # train linear model involving our selected structure
  mt.trained.lm = lm(Index ~ Hydro+GC, mt.training.set)
  # get predictions on training and test data
  mt.training.predictions = predict(mt.trained.lm, mt.training.set)
  mt.test.predictions = predict(mt.trained.lm, mt.test.set)
  
  pt.sample.n = nrow(pt.df)
  # construct ptDNA training set
  pt.training.refs = sample(seq(from=1,to=pt.sample.n), size=pt.sample.n*testpropn)
  pt.training.set = pt.df[pt.training.refs,]
  # construct ptDNA test set
  pt.test.set = pt.df[-pt.training.refs,]
  # train linear model involving our selected structure
  pt.trained.lm = lm(Index ~ Hydro+GC, pt.training.set)
  # get predictions on training and test data
  pt.training.predictions = predict(pt.trained.lm, pt.training.set)
  pt.test.predictions = predict(pt.trained.lm, pt.test.set)

  # use ptDNA-trained model to predict mtDNA response
  mt.cross.predictions = predict(pt.trained.lm, mt.df)
  # use mtDNA-trained model to predict ptDNA response
  pt.cross.predictions = predict(mt.trained.lm, pt.df)
  
  # report (spearman) correlations between predictions and data for each test
  mt.cross.accuracy = cor(mt.cross.predictions, mt.df$Index, method=cor.method)
  pt.cross.accuracy = cor(pt.cross.predictions, pt.df$Index, method=cor.method)
  mt.training.accuracy = cor(mt.training.predictions, mt.training.set$Index, method=cor.method)
  mt.test.accuracy = cor(mt.test.predictions, mt.test.set$Index, method=cor.method)
  pt.training.accuracy = cor(pt.training.predictions, pt.training.set$Index, method=cor.method)
  pt.test.accuracy = cor(pt.test.predictions, pt.test.set$Index, method=cor.method)

  # store statistics
  mt.test.acc = c(mt.test.acc, mt.test.accuracy)
  mt.training.acc = c(mt.training.acc, mt.training.accuracy)
  pt.test.acc = c(pt.test.acc, pt.test.accuracy)
  pt.training.acc = c(pt.training.acc, pt.training.accuracy)
  mt.cross.acc = c(mt.cross.acc, mt.cross.accuracy)
  pt.cross.acc = c(pt.cross.acc, pt.cross.accuracy)
}
# look at mean accuracies
mean(mt.training.acc)
mean(mt.test.acc)
mean(pt.training.acc)
mean(pt.test.acc)
mean(mt.cross.acc)
mean(pt.cross.acc)
results = data.frame(method="LM-reduced", mt.training=mean(mt.training.acc), mt.test=mean(mt.test.acc), pt.training=mean(pt.training.acc), pt.test=mean(pt.test.acc), mt.cross=mean(mt.cross.acc), pt.cross=mean(pt.cross.acc))

# dataframes for training-test data plots
mt.plot = rbind(data.frame(Predicted = mt.test.predictions, True = mt.test.set$Index, Class = "Test", GeneLabel = mt.test.set$GeneLabel), data.frame(Predicted = mt.training.predictions, True = mt.training.set$Index, Class = "Training", GeneLabel = mt.training.set$GeneLabel))
pt.plot = rbind(data.frame(Predicted = pt.test.predictions, True = pt.test.set$Index, Class = "Test", GeneLabel = pt.test.set$GeneLabel), data.frame(Predicted = pt.training.predictions, True = pt.training.set$Index, Class = "Training", GeneLabel = pt.training.set$GeneLabel))
mt.plot$Class = factor(mt.plot$Class, levels = c("Training", "Test"))
pt.plot$Class = factor(pt.plot$Class, levels = c("Training", "Test"))

mt.cross.plot = data.frame(Predicted = mt.cross.predictions, True = mt.df$Index, GeneLabel = mt.df$GeneLabel) 
pt.cross.plot = data.frame(Predicted = pt.cross.predictions, True = pt.df$Index, GeneLabel = pt.df$GeneLabel) 

# produce plots
mt.model.test.plot = ggplot(mt.plot, aes(x=Predicted, y=True, col=Class, fill=Class)) +
  geom_smooth(fullrange=TRUE, method="lm") +
  geom_point() +
  geom_text_repel(aes(label=GeneLabel), size = 3, segment.color = "#AAAAAA") +
  theme_light() + theme(legend.position = c(0.15,0.8), legend.background = element_rect(fill=alpha("#FFFFFF", 0.8))) +
  scale_fill_manual(values=c("#AAAAAA", "#FF8888")) +
  scale_color_manual(values=c("#AAAAAA", "#FF8888")) +
  xlab("Predicted retention index") + ylab("True retention index")
pt.model.test.plot = ggplot(pt.plot, aes(x=Predicted, y=True, col=Class, fill=Class)) +
  geom_smooth(fullrange=TRUE, method="lm") +
  geom_text_repel(aes(label=GeneLabel), size = 3, segment.color = "#AAAAAA") +
  geom_point() +
  theme_light() + theme(legend.position = c(0.15,0.8), legend.background = element_rect(fill=alpha("#FFFFFF", 0.8))) +
  scale_fill_manual(values=c("#AAAAAA", "#8888FF")) +
  scale_color_manual(values=c("#AAAAAA", "#8888FF")) +
  xlab("Predicted retention index") + ylab("True retention index")

mt.model.test.cross.plot = ggplot(mt.cross.plot, aes(x=Predicted, y=True)) +
  geom_smooth(method="lm", fullrange=TRUE, color="#FF8888", fill="#FFCCCC") +
  geom_text_repel(aes(label=GeneLabel), size = 3, segment.color = "#AAAAAA", color="#888888") +
  geom_point()  +
  theme_light() + xlab("Predicted retention index (from pt fit)") + ylab("True retention index")
pt.model.test.cross.plot = ggplot(pt.cross.plot, aes(x=Predicted, y=True)) +
  geom_smooth(method="lm", fullrange=TRUE, color="#8888FF", fill="#CCCCFF") +
  geom_text_repel(aes(label=GeneLabel), size = 3, segment.color = "#AAAAAA", color="#888888") +
  geom_point()  +
  theme_light() + xlab("Predicted retention index (from mt fit)") + ylab("True retention index")

# lay out overall figure
grid.arrange(mt.model.sel.plot, mt.model.test.plot, mt.model.test.cross.plot, 
             pt.model.sel.plot, pt.model.test.plot, pt.model.test.cross.plot, nrow=2)
