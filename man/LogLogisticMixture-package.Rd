\name{LogLogisticMixture-package}
\alias{LogLogisticMixture-package}
\alias{LogLogisticMixture}
\docType{package}
\title{
R library for fitting a Log Logistic Mixture model
~~ LogLogisticMixture ~~
}
\description{
The distribution of human response times, including eye fixation duration and saccadic response time, comprises single mode and a long tail. In the case of oculomotor response times, here are theoretical and empirical basis to argue that the observed distribution may be the result of 3 distinct processes (Feng, 2012, see link below), each following a power-law-like distribution. A 3-component mixture of log-logistic model seem to capture the distributional changes in reading and saccadic tasks.  

This library accompanies the following paper: http://www.academia.edu/1188136/Feng_2012_-_Is_there_a_Common_Control_Mechanism_for_Anti-saccades_and_Reading_Eye_Movements_Evidence_from_Distributional_Analyses. It requires the "muhaz" library.
}
\details{
\tabular{ll}{
Package: \tab LogLogisticMixture\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-11-10\cr
License: \tab What license is it under?\cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Gary Feng <gary.feng@gmail.com>
}
\references{
Please refer to the following paper: http://www.academia.edu/1188136/Feng_2012_-_Is_there_a_Common_Control_Mechanism_for_Anti-saccades_and_Reading_Eye_Movements_Evidence_from_Distributional_Analyses
}
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
# script for Feng (Vis Res revision 2010), Study 1
load('proanti1.rda')
source("fitll2.r")	# MLE estimate and a whole load of junk

# load data
datapro		<-list(pro=proanti$pro, anti_incorr=proanti$anti_incorrect);
dataanti	<-list(anti=proanti$anti, anti_corr=proanti$anti_correct);
s1<-c(datapro, dataanti)

# set the fixed parameters for the eXpress saccade component, which we only need to estimate its weight
lx.default<-80; sx.default<-3; # for the Px component

pdf('study1_mlefit_fixP.pdf', width=11, height=8.5);
# full models for pro-saccades, with the scale(location) parameters of P and C tied
fitprofull<-fitll2(datapro,  "?150/180/130, ?8/12/6, =1, ?6/7/2, ?0.7|0.01|0.98, ?100|150|50, ?0.05|0.15|0.01", maxtime=600,ylim=0.014, maintitle='Feng (Vis Res revision 2010, Study 1) Pro-saccades, Full model')
# full models for anti-saccades, with the scale(location) parameters of P and C tied
fitantifull	<-fitll2(dataanti,"?150/180/130, ?8/12/6, =1, ?6/7/2, ?0.3|0.01|0.9, ?100|150|60, ?0.05|0.15|0.01", maxtime=600,maintitle='Feng (Vis Res revision 2010, Study 1) Anti-saccades, Full model')
# restricted model with P component fixed
fit1<-fitll2(s1,"150, 10, 150, ?6/7/2, ?0.5|0|1, ?100|150|60, ?0.05|0.15|0.01", maxtime=600,ylim=0.014, , maintitle='Feng (Vis Res revision 2010, Study 1) Pro/Anti-saccades, Restricted model')
dev.off()
# and you can compare the AIC/BIC of the two models
}
