# script for Feng (Vis Res revision 2010), Study 1
load('study2.rda')
source("fitll2.r")	# MLE estimate and a whole load of junk

# data:
study2$overlap<-list(pro=study2$fixdur[(study2$gap=="Gap= -200") & (study2$cond=="Pro-saccade")], anti=study2$fixdur[(study2$gap=="Gap= -200") & (study2$cond=="Anti-saccade")])
study2$gap0<-list(pro=study2$fixdur[(study2$gap=="Gap=0") & (study2$cond=="Pro-saccade")], anti=study2$fixdur[(study2$gap=="Gap=0") & (study2$cond=="Anti-saccade")])
study2$gap200<-list(pro=study2$fixdur[(study2$gap=="Gap=200") & (study2$cond=="Pro-saccade")], anti=study2$fixdur[(study2$gap=="Gap=200") & (study2$cond=="Anti-saccade")])
# correct rate data:
crates<-list();
crates$overlap$pro$correct<-study2$fixdur[(study2$gap=="Gap= -200") & (study2$cond=="Pro-saccade")& (study2$correct==1)]
crates$overlap$pro$incorrect<-study2$fixdur[(study2$gap=="Gap= -200") & (study2$cond=="Pro-saccade")& (study2$correct==0)]
crates$overlap$anti$correct<-study2$fixdur[(study2$gap=="Gap= -200") & (study2$cond=="Anti-saccade")& (study2$correct==1)]
crates$overlap$anti$incorrect<-study2$fixdur[(study2$gap=="Gap= -200") & (study2$cond=="Anti-saccade")& (study2$correct==0)]
crates$gap0$pro$correct<-study2$fixdur[(study2$gap=="Gap=0") & (study2$cond=="Pro-saccade")& (study2$correct==1)]
crates$gap0$pro$incorrect<-study2$fixdur[(study2$gap=="Gap=0") & (study2$cond=="Pro-saccade")& (study2$correct==0)]
crates$gap0$anti$correct<-study2$fixdur[(study2$gap=="Gap=0") & (study2$cond=="Anti-saccade")& (study2$correct==1)]
crates$gap0$anti$incorrect<-study2$fixdur[(study2$gap=="Gap=0") & (study2$cond=="Anti-saccade")& (study2$correct==0)]
crates$gap200$pro$correct<-study2$fixdur[(study2$gap=="Gap=200") & (study2$cond=="Pro-saccade")& (study2$correct==1)]
crates$gap200$pro$incorrect<-study2$fixdur[(study2$gap=="Gap=200") & (study2$cond=="Pro-saccade")& (study2$correct==0)]
crates$gap200$anti$correct<-study2$fixdur[(study2$gap=="Gap=200") & (study2$cond=="Anti-saccade")& (study2$correct==1)]
crates$gap200$anti$incorrect<-study2$fixdur[(study2$gap=="Gap=200") & (study2$cond=="Anti-saccade")& (study2$correct==0)]

times<-c(-10000,seq(0,1000, by=20), 100000);
for(cond in names(crates)) {
	for (pa in names(crates[[cond]])) {
		cat(cond, ' : ', pa, '\n' )
		crates[[cond]][[pa]]$corr<-hist(as.numeric(crates[[cond]][[pa]]$correct),   breaks=times, plot=F)
		crates[[cond]][[pa]]$err <-hist(as.numeric(crates[[cond]][[pa]]$incorrect), breaks=times, plot=F)
		crates[[cond]][[pa]]$rate<-crates[[cond]][[pa]]$corr$counts/(crates[[cond]][[pa]]$corr$counts+crates[[cond]][[pa]]$err$counts)
		cat(crates[[cond]][[pa]]$rate,'\n')
	}
}

# MLE
lx.default<-80; sx.default<-3; # for the Px component
xrate=0.5; 
# Note that the correct rate varies across conditions
prate=list(overlap=1, gap0=0.875, gap200=0.9);
crate=list(overlap=1, gap0=0.875, gap200=0.9); 

lx<-lx.default; sx<-sx.default;
x<-times[2:52]+10;

pdf('study2.pdf', width=11, height=8.5);
fit<-list();
for(cond in c( "overlap", "gap0","gap200" )) {
	fit[[cond]]$full<-fitll2(study2[[cond]],"?150/180/130, ?8/12/6, =1, ?6/7/2, ?0.5|0|1, ?100|150|40, ?0.03|0.15|0.01", maxtime=600, ylim=0.015, maintitle=paste('Study 2 ',cond, ' Full model'))
	# restricted model with P component fixed
	fit[[cond]]$fixp<-fitll2(study2[[cond]],"150, 9, 150, ?6/7/2, ?0.5|0|1, ?100|150|40, ?0.03|0.15|0.01", maxtime=600,ylim=0.015, maintitle=paste('Study 2 ',cond, ' Restricted model'))

	l1	=fit[[cond]]$fixp$fit$pro$parameters[1];
	s1	=fit[[cond]]$fixp$fit$pro$parameters[2];
	l2	=fit[[cond]]$fixp$fit$pro$parameters[3];
	s2	=fit[[cond]]$fixp$fit$pro$parameters[4];
	mixp=fit[[cond]]$fixp$fit$pro$parameters[5];
	d	=fit[[cond]]$fixp$fit$pro$parameters[6];
	px	=fit[[cond]]$fixp$fit$pro$parameters[7];
	pdf_pro<- px * dllogis(x, lx, sx, 0) + (1-px)*(mixp * dllogis(x, l1, s1, 0) + (1-mixp) * dllogis(x, l2, s2, d));
	pdf_pro_correct<- xrate*px * dllogis(x, lx, sx, 0) + (1-px)*(prate[[cond]]*mixp * dllogis(x, l1, s1, 0) + crate[[cond]]*(1-mixp) * dllogis(x, l2, s2, d));
	model_pro_rate<-pdf_pro_correct/pdf_pro;

	l1	=fit[[cond]]$fixp$fit$anti$parameters[1];
	s1	=fit[[cond]]$fixp$fit$anti$parameters[2];
	l2	=fit[[cond]]$fixp$fit$anti$parameters[3];
	s2	=fit[[cond]]$fixp$fit$anti$parameters[4];
	mixp=fit[[cond]]$fixp$fit$anti$parameters[5];
	d	=fit[[cond]]$fixp$fit$anti$parameters[6];
	px	=fit[[cond]]$fixp$fit$anti$parameters[7];
	pdf_anti<- px * dllogis(x, lx, sx, 0) + (1-px)*(mixp * dllogis(x, l1, s1, 0) + (1-mixp) * dllogis(x, l2, s2, d));
	pdf_anti_correct<- xrate*px * dllogis(x, lx, sx, 0) + (1-px)*((1-prate[[cond]])*mixp * dllogis(x, l1, s1, 0) + crate[[cond]]*(1-mixp) * dllogis(x, l2, s2, d));
	model_anti_rate<-pdf_anti_correct/pdf_anti;

	plot(x, x,type='n', xlim=c(0,600), ylim=c(0,1),frame.plot=FALSE, xlab="Fixation Duration", ylab="Correct Rate");
	title(main = paste('Study 2 ',cond, ' Correct rate'))

	lines(x-20, crates[[cond]]$pro$rate[1:51], col='Red', lty='dashed')
	lines(x-20, crates[[cond]]$anti$rate[1:51], col='Blue', lty='dashed')
	lines(x, model_pro_rate, col='Red', lty='solid')
	lines(x, model_anti_rate, col='Blue', lty='solid')
}
#dev.off()

# note that the GAP200 condition doesn't really fit the 150/9 model
# eye-balling suggests 135 is a better number for the P
# we need to redo this condition.
cond="gap200";

# parameters set to 135
# pdf('study2_gap200_newfixp.pdf', width=11, height=8.5);

fit[[cond]]$fixp<-fitll2(study2[[cond]],"135, 9, 135, ?6/7/2, ?0.5|0|1, ?100|150|40, ?0.03|0.15|0.01", maxtime=600,ylim=0.014, maintitle=paste('Study 2 ',cond, ' Restricted model'))

l1	=fit[[cond]]$fixp$fit$pro$parameters[1];
s1	=fit[[cond]]$fixp$fit$pro$parameters[2];
l2	=fit[[cond]]$fixp$fit$pro$parameters[3];
s2	=fit[[cond]]$fixp$fit$pro$parameters[4];
mixp=fit[[cond]]$fixp$fit$pro$parameters[5];
d	=fit[[cond]]$fixp$fit$pro$parameters[6];
px	=fit[[cond]]$fixp$fit$pro$parameters[7];
pdf_pro<- px * dllogis(x, lx, sx, 0) + (1-px)*(mixp * dllogis(x, l1, s1, 0) + (1-mixp) * dllogis(x, l2, s2, d));
pdf_pro_correct<- xrate*px * dllogis(x, lx, sx, 0) + (1-px)*(prate[[cond]]*mixp * dllogis(x, l1, s1, 0) + crate[[cond]]*(1-mixp) * dllogis(x, l2, s2, d));
model_pro_rate<-pdf_pro_correct/pdf_pro;

l1	=fit[[cond]]$fixp$fit$anti$parameters[1];
s1	=fit[[cond]]$fixp$fit$anti$parameters[2];
l2	=fit[[cond]]$fixp$fit$anti$parameters[3];
s2	=fit[[cond]]$fixp$fit$anti$parameters[4];
mixp=fit[[cond]]$fixp$fit$anti$parameters[5];
d	=fit[[cond]]$fixp$fit$anti$parameters[6];
px	=fit[[cond]]$fixp$fit$anti$parameters[7];
pdf_anti<- px * dllogis(x, lx, sx, 0) + (1-px)*(mixp * dllogis(x, l1, s1, 0) + (1-mixp) * dllogis(x, l2, s2, d));
pdf_anti_correct<- xrate*px * dllogis(x, lx, sx, 0) + (1-px)*((1-prate[[cond]])*mixp * dllogis(x, l1, s1, 0) + crate[[cond]]*(1-mixp) * dllogis(x, l2, s2, d));
model_anti_rate<-pdf_anti_correct/pdf_anti;

plot(x, x,type='n', xlim=c(0,600), ylim=c(0,1),frame.plot=FALSE, xlab="Fixation Duration", ylab="Correct Rate");
title(main = paste('Study 2 ',cond, ' Correct rate'))

lines(x-20, crates[[cond]]$pro$rate[1:51], col='Red', lty='dashed')
lines(x-20, crates[[cond]]$anti$rate[1:51], col='Blue', lty='dashed')
lines(x, model_pro_rate, col='Red', lty='solid')
lines(x, model_anti_rate, col='Blue', lty='solid')

dev.off()
