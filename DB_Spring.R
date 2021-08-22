# journal article

# Priyadarshana, T. S., Lee, M-B., Ascher, J. S., Qiu, L., & Goodale, E. (2021). Crop heterogeneity is positively associated with beneficial insect diversity in subtropical farmlands. Journal of Applied Ecology 									

# data 

# Priyadarshana, T. S., Lee, M.-B., Ascher, J. S., Qiu, L., and Goodale, E. (2021). Data from: Crop heterogeneity is positively associated with beneficial insect diversity in subtropical farmlands.  Dryad Digital Repository. https://doi.org/10.5061/dryad.brv15dv9v

# note 1: most of the following R codes were adopted from Zuur et al., (2009), chapter 4

# Zuur, A. F., Ieno, E. N., Walker, N., Saveliev, A. A., & Smith, G. M. (2009). Mixed Effects Models and Extensions in Ecology with R. Springer. https://doi.org/10.1007/978-0-387-87458-6

# note 2: here is one example for dung beetles sampled during the spring surveys 

# for other analysis, the same procedure will be applied  

rm(list=ls()) # clean the R environment 

# load the packages 
library(pgirmess)
library(gstat)
library(nlme)
library (multcomp) 
library (MuMIn) 
library(car)
library (ape)
library (lme4)
library (visreg)
library(sp)
library(ecodist)
library(gstat)
library(emmeans)
library(ggplot2)
library(caret)
library(randomForest)

# check R version and information about the current R session

version
sessionInfo()
RStudio.Version()

# set the working directory 

setwd("~/Desktop/........") 

# following the site ID, fist, compile the site level and transect level data into one data frame  

DBSpring <- read.csv("DBSpring.csv",header=TRUE, sep=",",na.strings=TRUE)

# response variables 

# DBSWInd = dung beetle Shannon-Wiener diversity

# site-level scale explanatory variables

# CropCompositionalHeterogeneity = #crop compositional heterogeneity / 
                                   #crop Shannon-Wiener diversity 

# FieldMarginLength = # crop configurational heterogeneity / 
                      # field margin length in m

# CropArea100 = croplands % at 100 m radius / 

# FieldMarginType = field margin type

# landscape-level explanatory variables

# CropArea500 = croplands % at 500 m radius / 


# check for correlation between explanatory variables

cor.test(DBSpring$FieldMarginLength , DBSpring$CropCompositionalHeterogeneity, method="pearson")

cor.test(DBSpring$FieldMarginLength , DBSpring$CropArea100, method="pearson")

cor.test(DBSpring$FieldMarginLength , DBSpring$CropArea500, method="pearson")

cor.test(DBSpring$CropCompositionalHeterogeneity , DBSpring$CropArea100, method="pearson")

cor.test(DBSpring$CropCompositionalHeterogeneity , DBSpring$CropArea500, method="pearson")

cor.test(DBSpring$CropArea100,DBSpring$CropArea500, method="pearson")

# start with lm model 

mod1 <- lm(DBSWInd ~ FieldMarginLength + CropCompositionalHeterogeneity + CropArea100 + as.factor (FieldMarginType) + CropArea500, data = DBSpring)

plot(mod1) 
ncvTest(mod1) # heteroscedasticity is problematic 
vif(mod1) 
sqrt(vif(mod1)) > 2 
hist(resid(mod1))

# drop each explanatory variable one by one and re-fit the models to find which variable is causing the violation of heteroscedasticity or normal distribution 

# for example drop CropCompositionalHeterogeneity and re-fit the model

# now lm model assumption are met 

mod2 <- lm(DBSWInd ~ FieldMarginLength  + CropArea100 + as.factor (FieldMarginType) + CropArea500, data = DBSpring)
plot (mod2) 
ncvTest(mod2) 

# check for let's check the spatial autocorrelation 

inf.dists <- as.matrix(dist(cbind(lon=DBSpring$ddlongX, lat=DBSpring$ddlatY)))
inf.dists.inv <- 1/inf.dists
diag(inf.dists.inv) <- 0
Moran.I(resid(mod1), inf.dists.inv) # spatial autocorrelation is problematic

# to deal with heteroscedasticity and spatial autocorrelation now run GLS models 

# dealing with heteroscedasticity 

# put all the  variance structure into different formulas, then no need to re-type

# instead of focusing on the increase in variance with CropCompositionalHeterogeneity, maybe, the focus on the different variance per transect type is more important 

vf1 <- varFixed(~CropCompositionalHeterogeneity)
vf2 <- varIdent(form= ~ 1 | FieldMarginType)
vf3 <- varPower(form =~ CropCompositionalHeterogeneity)
vf4 <- varPower(form =~ CropCompositionalHeterogeneity | FieldMarginType)
vf5 <- varExp(form =~ CropCompositionalHeterogeneity)
vf6 <- varConstPower(form =~ CropCompositionalHeterogeneity)
vf7 <- varConstPower(form=~ CropCompositionalHeterogeneity | FieldMarginType)
vf8 <- varComb(varIdent(form =~ 1 | FieldMarginType), varExp(form =~ CropCompositionalHeterogeneity))

# now add all the above error structures to the GLS model 

mod1gls0 <- gls(DBSWInd ~ FieldMarginLength + CropCompositionalHeterogeneity + CropArea100 + as.factor (FieldMarginType) + CropArea500, data = DBSpring) 

# mod1gls0 is exactly like the mod1 since the error structures are not added yet

# also assign the following formula 'f3' when adding the error structures, so no need to keep retyping it

f3 <- formula(DBSWInd ~ FieldMarginLength + CropCompositionalHeterogeneity + CropArea100 + as.factor (FieldMarginType) + CropArea500)

mod1gls1 <- gls(f3, weights = vf1, data = DBSpring)
mod1gls2 <- gls(f3, weights = vf2, data = DBSpring)
mod1gls3 <- gls(f3, weights = vf3, data = DBSpring)
mod1gls4 <- gls(f3, weights = vf4, data = DBSpring)
mod1gls5 <- gls(f3, weights = vf5, data = DBSpring)
mod1gls6 <- gls(f3, weights = vf6, data = DBSpring)
mod1gls7 <- gls(f3, weights = vf7, data = DBSpring)
mod1gls8 <- gls(f3, weights = vf8, data = DBSpring)

# check AIC 

AIC(mod1gls1, mod1gls2, mod1gls3, mod1gls4, mod1gls5, mod1gls6, mod1gls7, mod1gls8)

# mod1gls7 is the best / optimal (i.e. corrected for heteroscedasticity) 

# graphical validation

# plot the standardized residuals vs. fitted for the original linear model and the model that introduces a variance structure and was chosen as optimal (i.e. mod1gls7)

plot(mod1gls0, which = c(1), main ="linear model")
plot(mod1gls7, which = c(1), main ="after adding the varaince structure")
E1 <- resid(mod1gls7)
coplot(E1 ~ CropCompositionalHeterogeneity | FieldMarginType, ylab = "Ordinary residuals", data = DBSpring)
E2 <- resid(mod1gls7, type = "normalized")
coplot(E2 ~ CropCompositionalHeterogeneity | FieldMarginType, data = DBSpring, ylab = "Normalised residuals")

# deal with violation of spatial autocorrelation

# this is a shortcut way to plot a variogram of the residuals from a gls model; resType='pearson' uses standardized residuals

plot(Variogram(mod1gls0, form =~ddlongX + ddlatY, robust = TRUE, maxDist = 2000, resType = "pearson")) 

# extract the residuals of the initial lm model

E<-rstandard(mod1)

# create a small data frame with the residuals and the coordinates

residAndCoord<-data.frame(E,DBSpring$ddlongX,DBSpring$ddlatY)

coordinates(residAndCoord)<-c("DBSpring.ddlongX","DBSpring.ddlatY")

# plot the residuals in space

bubble(residAndCoord,"E",col=c("black","grey"), main="Residuals",xlab="X-coordinates", ylab="Y-coordinates")

# when there is spatial autocorrelation similar residuals tend to be clustered in some areas, there is a tendency for the big grey bubbles to be located together

# add different spatial autocorrelation structures to the mod1gls7 (the model corrected for heteroscedasticity)

B1A<-gls(f3,correlation=corSpher(form=~ddlongX + ddlatY,nugget=T), weights = vf7, data=DBSpring)
B1B<-gls(f3,correlation=corLin(form=~ddlongX + ddlatY,nugget=T), weights = vf7,data=DBSpring)
B1C<-gls(f3,correlation=corRatio(form=~ddlongX + ddlatY, nugget=T),weights = vf7, data=DBSpring)
B1D<-gls(f3,correlation=corGaus(form=~ddlongX + ddlatY,nugget=T), weights = vf7, data=DBSpring)
B1E<-gls(f3,correlation=corExp(form=~ddlongX + ddlatY,nugget=T), weights = vf7, data=DBSpring)

AIC(mod1gls7,B1A,B1C,B1D,B1E)

# B1E is the best model, and this model has been model corrected for both heteroscedasticity and spatial autocorrelation

# to verify that the problems of independence are solved we plot again the variogram for the new model B1E

Vario1E <- Variogram(B1E, form =~ ddlongX + ddlatY, robust = TRUE, maxDist = 2000, resType = "pearson")
plot(Vario1E,smooth=FALSE)
Vario2E <- Variogram(B1E, form =~ ddlongX + ddlatY, maxDist = 2000, resType = "normalized")
plot(Vario2E, smooth = FALSE)

# extract standardized residuals from B1E (the best model)

gls.resids <- residuals(B1E, type = "normalized")

# create spatial points data frame of residuals and x/y coordinates 

gls.resids.spdf <- data.frame(gls.resids, DBSpring$ddlongX, DBSpring$ddlatY, data =DBSpring)

coordinates(gls.resids.spdf) <- c("DBSpring.ddlongX", "DBSpring.ddlatY")

bubble(gls.resids.spdf, "gls.resids", col = c("black", "grey"), main = "Residuals", xlab = "X-coordinates", ylab = "Y-coordinates")

# create Moran Correlogram

gls.moran <- correlog(coordinates(gls.resids.spdf), gls.resids.spdf$gls.resids, method = "Moran", nbclass = NULL)

plot(gls.moran)

# spatial autocorrelation is greatly reduced, look at the values on the y-axis

# successive difference contrasts

options(contrasts=c("contr.sdif", "contr.sdif"))

B1E<-gls(f3,correlation=corExp(form=~ddlongX+ddlatY,nugget=T), weights = vf7, data=DBSpring)
summary(B1E)
coef(B1E)
Anova(B1E)
confint(B1E)

# partial residual plot showing the effect of crop compositional heterogeneity

visreg(B1E, "CropCompositionalHeterogeneity", xlab="Crop Shannon-Weiner Diversity", ylab="Dung Beetle Shannon-Weiner Diversity", points=list(size=5, pch=1), gg=TRUE) + theme_bw()

# partial residual plot showing the effect of crop compositional heterogeneity

visreg(B1E, "FieldMarginLength", xlab="Field Margin Length (m)", ylab="Dung Beetle Shannon-Weiner Diversity", points=list(size=5, pch=1), gg=TRUE) + theme_bw()


# check the effects of FieldMarginType on dung beetle diversity 

em <- emmeans(B1E, "FieldMarginType", mode = "df.error")
pairs(em)
plot(em, comparisons = TRUE)
plot(em, xlab="Estimated Dung Beetle Shannon-Wiener Marginal Mean", ylab="Field Margin Types", horizontal = FALSE, comparisons = TRUE)+theme_bw()

# ranking the importance of the variables 

# load data

head(DBSpring)

# tuning the random forest model using the "caret" package

# select the explanatory variables / predictors

predictors <- DBSpring[,5:9] 
head(predictors)

# define parameter tuning method

fitControl <- trainControl(method = "repeatedcv", number = 10, repeats = 10, returnResamp ="all")

# tuning for dung beetle diversity 

tuning_total <- train(predictors, DBSpring$DBSWInd, method = "rf", trControl = fitControl)
tuning_total

# further confirmation of "mtry"

# mtry = number of variables randomly chosen at each split

dim(DBSpring)

# training Sample with 100 observations

train=sample(1:nrow(predictors),100)

dat.rf=randomForest(DBSWInd
  ~ FieldMarginLength + CropCompositionalHeterogeneity + CropArea100 + FieldMarginType + CropArea500, data = DBSpring , subset = train)
dat.rf

plot(dat.rf)

oob.err=double(5)
test.err=double(5)

set.seed(2345)
for(mtry in 1:5) 
{
  rf=randomForest(DBSWInd
                  ~ FieldMarginLength + CropCompositionalHeterogeneity + CropArea100 + FieldMarginType + CropArea500, data = DBSpring , subset = train,mtry=mtry,ntree=500) 
  oob.err[mtry] = rf$mse[500] # error of all trees fitted
  
  pred<-predict(rf,DBSpring[-train,]) # predictions on test set for each tree
  test.err[mtry]= with(DBSpring[-train,], mean( (DBSWInd - pred)^2)) # mean squared test error
  
  cat(mtry," ") # printing the output to the console
}

test.err
oob.err

matplot(1:mtry , cbind(oob.err,test.err), pch=19 , col=c("red","blue"),type="b",ylab="Mean Squared Error",xlab="Number of Predictors Considered at each Split")
legend("topright",legend=c("Out of Bag Error","Test Error"),pch=19, col=c("red","blue"))

# now build the final random forest model 

modRF1 <- randomForest(DBSWInd ~ FieldMarginLength + CropCompositionalHeterogeneity + CropArea100 + FieldMarginType + CropArea500, data = DBSpring, ntree=500, mtry=2, importance=F, proximity= TRUE, metric =Accuracy)
print(modRF1, metric =Kappa)
print(modRF1)
VerbImp <- as.data.frame(sort(importance(modRF1)[,1],decreasing = TRUE),optional = T)
VerbImp
varImpPlot(modRF1, pch=2, main="Variable Importance") # plot the importance of the variables 
