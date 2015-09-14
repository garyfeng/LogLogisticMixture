# need to first change working dir to the following directory before sourcing
source("D:\\SourceForge\\loglogisticmixture\\fitll2.r")

# need to swith workign directory back to the data directory

library(foreign)
study3<-read.spss("D:/My Projects/2007 Dundee Corpus/Dundee WFREQ2 FFD for r.sav")
save(study3, file="study3.rda")
names(study3)


lx.default<-100; sx.default<-4; # for the Px component
dundee<-list(freq1=study3$FDUR[study3$wfreqclass==1], freq2=study3$FDUR[study3$wfreqclass==2], freq3=study3$FDUR[study3$wfreqclass==3], freq4=study3$FDUR[study3$wfreqclass==4], freq5=study3$FDUR[study3$wfreqclass==5])
study3$wl <-factor(study3$WLEN)
study3$wf <-factor(study3$wfreqclass)
levels(study3$wl)

lx.default <- 110; sx.default <- 3; # for the Px component

########## FREQ effect at each WLEN #########
pdf('study3_mlefit_freqXwlen_full_2.pdf', width=11, height=8.5);
for (l in levels(study3$wl)) {
	dundee<-list()
	for (f in levels(study3$wf)) {
		tmp <-study3$FDUR[study3$wf==f & study3$wl==l];
		name<-paste("freq",f, sep=""); #cat(name, " length= ", length(tmp), '\n');
		if (length(tmp)>200) {dundee[[name]]<-tmp;}
		# make sure there are data in it.
	}
	fit3full 	<- fitll2(dundee, "?180/220/160, ?6/9/4, =1, ?4/6/3, ?0.8|0.3|0.98, ?100|150|60, ?0.05|0.2|0.01" )	
	print(paste("Wordlengh = ", l))
	for( f in fit3full$fit) {
		print(paste(f$n, f$value, f$BIC))
	}
}
dev.off();
########## FREQ effect at each WLEN, with fixed P parameter #########
pdf('study3_mlefit_freqXwlen_fixp_2.pdf', width=11, height=8.5);
for (l in levels(study3$wl)) {
	dundee<-list()
	for (f in levels(study3$wf)) {
		tmp <-study3$FDUR[study3$wf==f & study3$wl==l];
		name<-paste("freq",f, sep=""); #cat(name, " length= ", length(tmp), '\n');
		if (length(tmp)>200) {dundee[[name]]<-tmp;}
		# make sure there are data in it.
	}
	fit3res <- fitll2(dundee, paste(as.numeric(l)*3+181,", 7.5, ",as.numeric(l)*3+181,", ?4/6/3, ?0.8|0.3|0.98, ?100|150|60, ?0.05|0.2|0.01", sep='') )	
	print(paste("Wordlengh = ", l))
	for( f in fit3res$fit) {
		print(paste(f$n, f$value, f$BIC, f$par[1], f$par[2], f$par[3], f$par[4]))
	}
}
dev.off();


########## FREQ effect at each WLEN, with fixed P parameter #########
pdf('study3_mlefit_freqXwlen_fixpd95_2.pdf', width=11, height=8.5);
for (l in levels(study3$wl)) {
	dundee<-list()
	for (f in levels(study3$wf)) {
		tmp <-study3$FDUR[study3$wf==f & study3$wl==l];
		name<-paste("freq",f, sep=""); #cat(name, " length= ", length(tmp), '\n');
		if (length(tmp)>200) {dundee[[name]]<-tmp;}
		# make sure there are data in it.
	}
	fit3res <- fitll2(dundee, paste(as.numeric(l)*3+181,", 7.5, ",as.numeric(l)*3+181,", ?4/6/3, ?0.8|0.3|0.98, 95, ?0.05|0.2|0.01", sep='') )	
	print(paste("Wordlengh = ", l))
	for( f in fit3res$fit) {
		print(paste(f$n, f$value, f$BIC, f$par[1], f$par[2], f$par[3]))
	}
}
dev.off();



for (l in levels(study3$wl)) {
	dundee<-list()
	for (f in levels(study3$wf)) {
		tmp <-study3$FDUR[study3$wf==f & study3$wl==l];
		if (length(tmp)>0) print (paste(mean(tmp), sd(tmp), sep="\t"))
	}
}
