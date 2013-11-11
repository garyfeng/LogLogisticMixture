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
fitprofull<-fitll2(datapro,  "?150/180/130, ?8/12/6, =1, ?6/7/2, ?0.7|0.01|0.98, ?100|150|50, ?0.05|0.15|0.01", ylim=0.014, maintitle='Feng (Vis Res revision 2010, Study 1) Pro-saccades, Full model')
# full models for anti-saccades, with the scale(location) parameters of P and C tied
fitantifull	<-fitll2(dataanti,"?150/180/130, ?8/12/6, =1, ?6/7/2, ?0.3|0.01|0.9, ?100|150|60, ?0.05|0.15|0.01", maintitle='Feng (Vis Res revision 2010, Study 1) Anti-saccades, Full model')
# restricted model with P component fixed
fit1<-fitll2(s1,"150, 10, 150, ?6/7/2, ?0.5|0|1, ?100|150|60, ?0.05|0.15|0.01", ylim=0.014, , maintitle='Feng (Vis Res revision 2010, Study 1) Pro/Anti-saccades, Restricted model')
dev.off()
# and you can compare the AIC/BIC of the two models
