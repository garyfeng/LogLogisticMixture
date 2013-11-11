# script for Feng (Vis Res revision 2010), Study 1
load('proanti1.rda')
source("fitll2.r")	# MLE estimate and a whole load of junk

# load data
datapro		<-list(pro=proanti$pro, anti_incorr=proanti$anti_incorrect);
dataanti	<-list(anti=proanti$anti, anti_corr=proanti$anti_correct);
s1<-c(datapro, dataanti)
dataall <-list(pro=proanti$pro,anti=proanti$anti);
dataantionly <-list(anti_incorr=proanti$anti_incorrect, anti_corr=proanti$anti_correct);
# set the fixed parameters for the eXpress saccade component, which we only need to estimate its weight
lx.default<-80; sx.default<-3; # for the Px component


pdf('study1_mlefit_fixP.pdf', width=11, height=8.5);
# full models for pro-saccades, with the scale(location) parameters of P and C tied
fitprofull<-fitll2(datapro,  "?150/180/130, ?8/12/6, =1, ?6/7/2, ?0.7|0|1, ?100|150|40, ?0.05|0.15|0.01", maxtime=600,ylim=0.014, maintitle='Feng (Vis Res revision 2010, Study 1) Pro-saccades, Full model')
# full models for anti-saccades, with the scale(location) parameters of P and C tied
fitantifull	<-fitll2(dataanti,"?150/180/130, ?8/12/6, =1, ?6/7/2, ?0.3|0|1, ?100|150|40, ?0.05|0.15|0.01", maxtime=600,maintitle='Feng (Vis Res revision 2010, Study 1) Anti-saccades, Full model')
# restricted model with P component fixed
fit1<-fitll2(s1,"150, 9, 150, ?6/7/2, ?0.5|0|1, ?100|150|40, ?0.05|0.15|0.01", maxtime=600,ylim=0.014, , maintitle='Feng (Vis Res revision 2010, Study 1) Pro/Anti-saccades, Restricted model')
dev.off()
# and you can compare the AIC/BIC of the two models

pdf('study1_mlefit3.pdf', width=11, height=8.5);
fitall<-fitll2(dataall,  "?150/180/130, ?8/12/6, =1, ?6/7/2, ?0.5|0|1, ?100|150|40, ?0.05|0.15|0.01", maxtime=600,ylim=0.014, maintitle='Study 1, Full model')
fit1<-fitll2(dataall,"150, 9, 150, ?6/7/2, ?0.5|0|1, ?100|150|40, ?0.05|0.15|0.01", maxtime=600,ylim=0.014, , maintitle='Study 1, Restricted model')

fitantionly<-fitll2(dataantionly,  "?150/180/130, ?8/12/6, =1, ?6/7/2, ?0.5|0|1, ?100|150|40, ?0.05|0.15|0.01", maxtime=600,ylim=0.014, maintitle='Study 1, Anti-saccade Full model')
# we now fix the P, C, and d parameters. The new model should be better ANTI-only model above.
fitantionly<-fitll2(dataantionly,  "150, 9, 150, 4.04, ?0.5|0|1, 109, ?0.05|0.15|0.01", maxtime=600,ylim=0.014, maintitle='Study 1, Anti-saccades, Restricted model')
dev.off()

# model for the correct rate:
#	eXpress saccades -- 50% correct rate
#	P component: 100% correct for Pro-saccade and 0% correct for anti-saccades
#	C component: q% correct
times<-c(-10000,seq(0,1000, by=20), 100000);
pro_corr<-hist(as.numeric(proanti$pro_correct), breaks=times, plot=F)
pro_err<-hist(as.numeric(proanti$pro_incorrect), breaks=times, plot=F)
anti_corr<-hist(as.numeric(proanti$anti_correct), breaks=times, plot=F)
anti_err<-hist(as.numeric(proanti$anti_incorrect), breaks=times, plot=F)
pro_rate <-pro_corr$counts/(pro_err$counts+pro_corr$counts)
anti_rate <-anti_corr$counts/(anti_err$counts+anti_corr$counts)
xrate=0.5; prate=1; crate=0.9; lx<-lx.default; sx<-sx.default;
l1=fit1$fit$pro$parameters[1];
s1=fit1$fit$pro$parameters[2];
l2=fit1$fit$pro$parameters[3];
s2=fit1$fit$pro$parameters[4];
mixp=fit1$fit$pro$parameters[5];
d=fit1$fit$pro$parameters[6];
px=fit1$fit$pro$parameters[7];
x<-times[2:52]+10;
pdf_pro<- px * dllogis(x, lx, sx, 0) + (1-px)*(mixp * dllogis(x, l1, s1, 0) + (1-mixp) * dllogis(x, l2, s2, d));
pdf_pro_correct<- xrate*px * dllogis(x, lx, sx, 0) + (1-px)*(prate*mixp * dllogis(x, l1, s1, 0) + crate*(1-mixp) * dllogis(x, l2, s2, d));
model_pro_rate<-pdf_pro_correct/pdf_pro;
l1=fit1$fit$anti$parameters[1];
s1=fit1$fit$anti$parameters[2];
l2=fit1$fit$anti$parameters[3];
s2=fit1$fit$anti$parameters[4];
mixp=fit1$fit$anti$parameters[5];
d=fit1$fit$anti$parameters[6];
px=fit1$fit$anti$parameters[7];
pdf_anti<- px * dllogis(x, lx, sx, 0) + (1-px)*(mixp * dllogis(x, l1, s1, 0) + (1-mixp) * dllogis(x, l2, s2, d));
pdf_anti_correct<- xrate*px * dllogis(x, lx, sx, 0) + (1-px)*((1-prate)*mixp * dllogis(x, l1, s1, 0) + crate*(1-mixp) * dllogis(x, l2, s2, d));
model_anti_rate<-pdf_anti_correct/pdf_anti;

plot(x, pro_rate[1:51],type='n', ylim=c(0,1),frame.plot=FALSE, xlab="Fixation Duration", ylab="Correct Rate")
lines(x-20, pro_rate[1:51], col='Blue', lty='dashed')
lines(x-20, anti_rate[1:51], col='Red', lty='dashed')
lines(x, model_pro_rate, col='Blue', lty='solid')
lines(x, model_anti_rate, col='Red', lty='solid')


