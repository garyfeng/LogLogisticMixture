
load("zstring.rda")


##############################
# study 4

## remember to set fitdata() function: ylim in both pdf and haz plots
pdf('study4_mlefit_full_fixp_x60_3.pdf', width=11, height=8.5);

lx.default<-60; sx.default<-3; # for the Px component
fit4full  <- fitll2(zstring, "?180/220/160, ?6/9/4, =1, ?4/6/3, ?0.8|0.3|0.98, ?100|150|60, ?0.1|0.2|0.03" )
fit4fixp <- fitll2(zstring, "185, 5.3, 185, ?4/6/3, ?0.8|0.3|0.98, ?100|150|60, ?0.1|0.2|0.03" )
dev.off();

for( f in fit4full$fit) {
	print(paste(f$n, f$value, f$BIC, f$par[1], f$par[2], f$par[3]))
}

for( f in fit4fixp$fit) {
	print(paste(f$n, f$value, f$BIC, f$par[1], f$par[2], f$par[3], f$par[4]))
}
