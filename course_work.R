#Importing Packages

library(ggplot2)
library(dplyr)
library(tidyverse)
library(corrplot)
library(EasyABC)
library(car)
#Task-1 -Preliminary data analysis

X = read.csv("X.csv")
Y = read.csv("y.csv")
T=read.csv("time.csv")

x = data.matrix(X)
y = data.matrix(Y)
time=data.matrix(T)



ones = matrix(0,2,nrow=100)

x1= x[x[,2]==0,]
x2= x[x[,2]==1,]
y1= y[1:100]
y2= y[101:200]


#neutral_audio plot
plot(time[1:100],x1[,1], type="l", lty=1, xlab="Time", ylab="Neutral Audio", main="Time Series Plot Neutral Audio")

#emotional_audio plot
plot(time[101:200],x2[,1], type="l", lty=1, xlab="Time", ylab="Emotional Audio", main="Time Series Plot Emotional Audio")

#Input Audio
plot(time,x[,1], type="l", lty=1, xlab="Time", ylab="Input Audio", main="Time Series Plot Input",col="red")

#Output Audio
plot(time,y, type="l", lty=1, xlab="Time", ylab="Output Signal", main="Time Series Plot Output",col="blue")


#Histogram of output signal
xin=X$x1[1:100]
xin
hist(y,col="blue",main = "Histogram of output y",xlab = "Output")
hist(y1,col="Red",main = "Histogram of output y when x2=0",xlab = "Output")
hist(y2,col = "Purple",main = "Histogram of output y when x2=1",xlab = "Output")

#Histogram of Input signal
xin1=X$x1[101:200]
xin1
hist(xin,col="blue",main = "Histogram of input x2=0",xlab = "Input")
hist(xin1,col="Orange",main = "Histogram of input x2=1",xlab = "Input")

#Correlation X==0
xin=X$x1[1:100]
xin
data=cbind(xin,y1)
data=as.data.frame(data)
data
corr=cor(data)
corrplot::corrplot.mixed(corr,upper = "pie",lower = "number")

#Correlation X==1
xin1=X$x1[101:200]
xin1
data1=cbind(xin1,y2)
data1=as.data.frame(data1)
data1
corr1=cor(data1)
corrplot::corrplot.mixed(corr1,upper = "pie",lower = "number")


#scatter plot between the audio input and output brain signal
x=as.data.frame(x)
x1=x$x1
x1=x1[1:100]
x2=x$x1
x2=x2[101:200]

plot(x1,y1,xlab="Input",ylab="Output",main="scatter plot when x2=0",col="blue")
plot(x2,y2,xlab="Input",ylab="Output",main="scatter plot when x2=1",col="red")

#Box plot of output signal x2=0
boxplot(y1,main="Box Plot of output signal when x2=0 ",col = "green")

#Box plot of output signal x2=1
boxplot(y2,main="Box Plot of output signal ",col = "Blue")



#Task 2:Regression – modelling the relationship between EEG signals

#Task-2.1:Estimate model parameters

comb=cbind(x,y)
regdata=as.data.frame(comb)
first_input=regdata$x1
second_input=regdata$x1^2
third_input=regdata$x1^3
fourth_input=regdata$x1^4
fifth_input=regdata$x1^5
x2=regdata$x2


#model 1

model1=data.frame(input0=c(1),third_input,fifth_input,x2)
mod1=as.matrix(model1)
md1=solve(t(mod1)%*%mod1)%*%t(mod1)%*%y
md1
length(md1)

#model2

model2=data.frame(input0=c(1),first_input,x2)
mod2=as.matrix(model2)
md2=solve(t(mod2)%*%mod2)%*%t(mod2)%*%y
md2
length(md2)

#model 3

model3=data.frame(input0=c(1),first_input,second_input,fourth_input,x2)
mod3=as.matrix(model3)
md3=solve(t(mod3)%*%mod3)%*%t(mod3)%*%y
md3
length(md3)

#model 4

model4=data.frame(input0=c(1),first_input,second_input,third_input,fifth_input,x2)
mod4=as.matrix(model4)
md4=solve(t(mod4)%*%mod4)%*%t(mod4)%*%y
md4
length(md4)

#model 5

model5=data.frame(input0=c(1),first_input,third_input,fourth_input,x2)
mod5=as.matrix(model5)
md5=solve(t(mod5)%*%mod5)%*%t(mod5)%*%y
md5
length(md5)

#Task 2.2:compute the model residual(error)sum of squared errors 

rss1=norm((y-mod1%*%md1)^2)
rss1

rss2=norm((y-mod2%*%md2)^2)
rss2

rss3=norm((y-mod3%*%md3)^2)
rss3


rss4=norm((y-mod4%*%md4)^2)
rss4


rss5=norm((y-mod5%*%md5)^2)
rss5


#Task 2.3:Compute the log-likelihood function

lik_fn1=-(nrow(Y)/2*log(2*pi))-(nrow(Y)/2*(log(rss1/(nrow(Y)-1))))-rss1*(1/(2*(rss1/(nrow(Y)-1))))
lik_fn1

lik_fn2=-(nrow(Y)/2*log(2*pi))-(nrow(Y)/2*(log(rss2/(nrow(Y)-1))))-rss2*(1/(2*(rss2/(nrow(Y)-1))))
lik_fn2

lik_fn3=-(nrow(Y)/2*log(2*pi))-(nrow(Y)/2*(log(rss3/(nrow(Y)-1))))-rss3*(1/(2*(rss3/(nrow(Y)-1))))
lik_fn3


lik_fn4=-(nrow(Y)/2*log(2*pi))-(nrow(Y)/2*(log(rss4/(nrow(Y)-1))))-rss4*(1/(2*(rss4/(nrow(Y)-1))))
lik_fn4

lik_fn5=-(nrow(Y)/2*log(2*pi))-(nrow(Y)/2*(log(rss5/(nrow(Y)-1))))-rss5*(1/(2*(rss5/(nrow(Y)-1))))
lik_fn5


#Task 2.4:Compute the Akaike information criterion (AIC) and Bayesian information criterion (BIC) 

k1=4
aic1=2*k1-(2*lik_fn1)
aic1
bic=k1*log(nrow(Y))-(2*(lik_fn1))
bic

k2=3
aic2=2*k2-(2*lik_fn2)
aic2
bic2=k2*log(nrow(Y))-(2*(lik_fn2))
bic2


k3=5
aic3=2*k3-(2*lik_fn3)
aic3
bic3=k3*log(nrow(Y))-(2*(lik_fn3))
bic3


k4=6
aic4=2*k4-(2*lik_fn4)
aic4
bic4=k4*log(nrow(Y))-(2*(lik_fn4))
bic4

k5=5
aic5=2*k5-(2*lik_fn5)
aic5
bic5=k5*log(nrow(Y))-(2*(lik_fn5))
bic5

# Task 2.5
error1=(y-(mod1%*%md1))
head(error1)

qqnorm(error1,)
qqline(error1,col="red")
qqnorm(md1)
qqline(md1)

error2=(y-(mod2%*%md2))
head(error2)

qqnorm(error2,)
qqline(error2,col="orange")
qqnorm(md2)
qqline(md2)
error3=(y-(mod3%*%md3))
head(error3)

qqnorm(error3,)
qqline(error3,col="blue")
qqnorm(md3)
qqline(md3)

error4=(y-(mod4%*%md4))
head(error4)

qqnorm(error4,)
qqline(error4,col="green")
qqnorm(md4)
qqline(md4)

error5=(y-(mod5%*%md5))
head(error5)

qqnorm(error5,)
qqline(error5,col="red")

qqnorm(md5)
qqline(md5)
#Task 2.6

#In my perspective, Model 3 is the best fit,because of low AIC and BIC value.

#Task-2.7

#Taking sample of 100 rows
set.seed(100)
train_sample=sample(1:200,140)
head(train_sample)

#Traning_data
training=regdata[train_sample,]
head(training)

#Testing_data
testing=regdata[-train_sample,]
head(training)

#Training_model

model03=data.frame(input0=c(1),training$x1,training$x1^2,training$x1^4,training$x2)
mod03=as.matrix(model03)
head(mod03)

#Testing model
y1=training$y
y1=as.matrix(y1)
md03=solve(t(mod03)%*%mod03)%*%t(mod03)%*%y1
md03

#Prediction of Testing data
estimation=sum(md03)/length(y1)
estimation

prediction=estimation*testing
head(prediction)

#Confidence interval 95% of output signal

sd=sqrt(var(prediction$y))
means=mean(prediction$y)
sd
means
paste("confidence interval",means-1.96*sd/sqrt(nrow(testing)),means+1.96*sd/sqrt(nrow(testing)))

#confidence interval 95% of input signal

sd0=sqrt(var(prediction$x1))
means0=mean(prediction$x1)
sd0
means0
paste("confidence interval",means0-1.96*sd0/sqrt(nrow(testing)),means0+1.96*sd0/sqrt(nrow(testing)))

#error

error=testing$y-prediction$y
error1=testing$y+prediction$y

#plotting confidence interval 

prediction_lower=error
prediction_upper=error1

prediction_test_y=testing$y

ggplot(prediction,aes(x1,y))+geom_errorbar(ymin=prediction_lower,ymax=prediction_upper)+geom_point(aes(x1),col="blue")+geom_point(aes(x1,testing$y))+ggtitle("Error Bar Plot")+theme(plot.title = element_text(hjust = 0.5))

#Task=3-Approximate Bayesian Computation (ABC)

# 2 parameters with largest absolute values in  least squares estimation 

mean1=colSums(abs(mod1))/nrow(mod1)
head(mean1)

mean2=colSums(abs(mod2))/nrow(mod2)
head(mean2)

mean3=colSums(abs(mod3))/nrow(mod3)
head(mean3)

mean4=colSums(abs(mod4))/nrow(mod4)
head(mean4)

mean5=colSums(abs(mod5))/nrow(mod5)
head(mean5)

####The biggest least square estimate is for model 1.

model001=as.data.frame(mod1)
head(model001)


# determine the range of the prior distribution.

p1=runif(min = min(model001$third_input),max = max(model001$third_input),nrow(model001))
p1

p2=runif(min=min(model001$fifth_input),max=max(model001$fifth_input),nrow(model001))
p2

#sample 
n=50


#joint distribution parameter

sample=function(x)
{
  c(x[3]^3+p1,x[5]^5+p2)
}
prior=list(c('unif',0,1))
simulation=ABC_rejection(model = sample,prior = prior,nb_simul = n,tol = 0.2)
densityPlot(simulation$param)

# Marginal  Distribution of input_three parameter

sample1=function(x)
{
  x[3]^3+p1
  }
prior1=list(c('unif',0,1))
simulation1=ABC_rejection(model = sample1,prior = prior1,nb_simul = n,tol = 0.2)
densityPlot(simulation1$param,xlab="Simulated parameter",main=" Marginal  distribution of input_three parameter",col="red")

#plot marginal distribution of input_five

sample2=function(x)
{
  x[5]^5+p2
}
prior2=list(c('unif',0,1))
simulation2=ABC_rejection(model = sample2,prior = prior2,nb_simul = n,tol = 0.2)
densityPlot(simulation2$param,xlab = "simulated parameter",main=" Marginal  distribution of input_five parameter",col="blue")
# Result:- For both the marginal and joint distribution of the two posterior parameters, the density curves appear to be normally distributed.

