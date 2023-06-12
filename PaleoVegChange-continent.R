#################################################
# Code to analyze bidoversity and community change in Neotoma db pollen datasets
# Copyright: 2023 M. Allison Stegner
#################################################

# LOAD DATA #######################################
# LOAD NEOTOMA R OBJECT___
load('pol_dlx_2020-12-01.Rdata')

# LOAD AGE MODELS___
# Bayesian age models used here are from Wang et al. (2019) Bayesian ages for pollen records since the last glaciation. Sci Data 6, 176
# code to generate these age models is available at https://github.com/yuewangpaleo/BaconAgeNeotoma
# Age models here are provided with permission from Y. Wang

path<-"Wang-et-al-Cores/Cores_all" # file path for age model data

# LOAD FUNCTIONS #######################################
source('PaleoVegChange_functions.R', chdir = TRUE)

# ANALYSIS ########################################
# DEFINE BINS
step<- -250
ageolder<-13950
bins<-seq(ageolder,step,step)
bin.means<-(bins[-length(bins)]+bins[-1])/2 

# step<- -100
# ageolder<-13950
# bins<-seq(ageolder,step,step) 

# step<- -500
# ageolder<-13950
# bins<-seq(ageolder,step,step) 

################ CONTINENT SCALE ANALYSIS
# settings
iter.siteresamp=50000 
iter.amresamp=1000 
N.sites<-100 

# # test
# iter.siteresamp=50
# iter.amresamp=5 
# N.sites<-10

# CALCULATE ABRUPT CHANGES________________________________________
# per site______________________
# Using a Neotoma download object, run principle curve on the pollen from each site, and calculate BCP
eco.group<-c("TRSH","UPHE") # identify which ecological groups to include in pollen dataset
threshold<-0.5 # bcp posterior probability threshold

multi.bcp<-list()
for (i in 1:length(pol_dlx)){
	print(i)
	pol_ds<-pol_dlx[[i]]
	am_posteriors<-get_am_posteriors(pol_ds)
	am_posteriors<-am_posteriors[complete.cases(am_posteriors),]
	prep.pol<-prep_pollen(pol_ds,eco.group,type="pct")
	
	# remove rows with NA pollen data
	na.inds<-complete.cases(prep.pol$prep.pol)
	prep.pol$prep.pol<-prep.pol$prep.pol[na.inds,]
	prep.pol$depth<-prep.pol$depth[na.inds]
	am_posteriors<-am_posteriors[na.inds,]
	
	multi.bcp[[i]]<-ac_agejigger(prep.pol,am_posteriors,iters=iter.amresamp,type="bcp",threshold)
}

names(multi.bcp)<-names(pol_dlx)

# Abrupt change iteration and slotting into bins__
# integrate a resampling routine so that the same number of sites (N) is used for each time bin
bcp.resamp<-resample.bin.tally(bins,multi.bcp,iter=iter.siteresamp,N.sites=N.sites)

# # Quick plot
# bin.quants<-bcp.resamp$bin.quants
# bin.ends<-bins[-1] # plot at the end of the bin	
# ccs<-complete.cases(t(bin.quants))
# bin.quants<-bin.quants[,ccs]
# bin.ends<-bin.ends[ccs]

# plot(bin.ends,bin.quants["50%",],type="l",ylim=c(0,20))
# polygon(c(bin.ends,rev(bin.ends)),c(bin.quants["2.5%",],rev(bin.quants["97.5%",])))

# LOSS/GAIN ANALYSIS_________________________
# run pollen loss gain iteration
N.years.prior<-500
eco.group<-c("TRSH","UPHE")

# # test run
# gl.out<-gain_loss(POL_DL=pol_dly,iter.siteresamp=50,iter.amresamp=10,N.years.prior=1000,N.sites=100)

gl.out<-gain_loss(POL_DL=pol_dlx,iter.siteresamp=iter.siteresamp,iter.amresamp=iter.amresamp,N.years.prior=N.years.prior,N.sites=N.sites)
# first counter is age model resampling for every site
# second counter is time bins completed

# # quick plot
# lossk<-gl.out[[6]] # Loss
# gaink<-gl.out[[5]] # Gain
# bin.ends<-bins[-1] # plot at the end of the bin
# plot.ts(lossk,bin.ends,YLAB="",new=F,add=F,line.col="black",poly.col="gray",YLIM=c(0,10),XLIM=c(0,14000),hatched=F)
# plot.ts(gaink,bin.ends,YLAB="",new=F,add=T,line.col="blue",poly.col="lightblue",XLIM=c(0,14000),hatched=F)
# axis(1)
# axis(2)
# legend("topleft",c("Gain","Loss"),col=c("lightblue","gray"),pch=15,bty="n",cex=1.5)

# # DIVERSITY and FIRST/LAST OCCURENCES______________________
# # test run
# div.flad.run<-div.flad.in.bin(pol_dlx,iter.siteresamp=50,am.iter=10,bins=bins,N.sites=100,burnin=3)

div.flad.run<-div.flad.in.bin(pol_dlx,iter.siteresamp=iter.siteresamp,am.iter=iter.amresamp,bins=bins,N.sites=N.sites,burnin=3)
# first counter tracks diversty estimation per site
# second counter tracks FAD/LAD estimation per site

# # quick plot
# bin.ends<-bins[-1] # plot at the end of the bin
# sqsk<-div.flad.run[[3]]
# plot.ts(sqsk,bin.ends,YLAB="",new=F,add=F,line.col="black",poly.col="gray",XLIM=c(14000,0),yaxt="n",bty="n",YLIM=c(0,10),hatched=F)
# axis(1)
# axis(2)

# fadk<-div.flad.run[[1]]
# ladk<-div.flad.run[[2]]

# plot.ts(fadk,bin.ends,YLAB="",new=F,add=F,line.col="red",poly.col="pink",YLIM=c(0,max(fadk,ladk,na.rm=T)),XLIM=c(14000,0),xaxt="n",yaxt="n",bty="n",hatched=T,angle=35)
# plot.ts(ladk,bin.ends,YLAB="",new=F,add=T,line.col="blue",poly.col="lightblue",XLIM=XLIM,hatched=T,angle=-75)
# axis(1)
# axis(2)

# PLOT COMPOSITE FIGURE_________________________

# trim times series to start at the same time bin____
max.age<-11575
cut.inds<-which(bin.means > max.age)
bin.means<-bin.means+step/2 # plot at end of bin
bin.means.sub<-bin.means[-cut.inds]

# SQS
sqsk<-div.flad.run[[3]]
sqsk<-sqsk[-cut.inds,]
		
# BCP
bcpk<-t(bcp.resamp$bin.quants)[,c(1,3,5)]
bcpk<-bcpk[-cut.inds,]

# Loss
lossk<-gl.out[[6]]
lossk<-lossk[-cut.inds,]
	
# Gain
gaink<-gl.out[[5]]
gaink<-gaink[-cut.inds,]

# FAD
fadk<-div.flad.run[[1]]
fadk<-fadk[-cut.inds,]
	
# LAD
ladk<-div.flad.run[[2]]
ladk<-ladk[-cut.inds,]

# Plot___
dev.new(height=6,width=4)
par(mfrow=c(4,1),mar=c(0,4,0,1),oma=c(4,1,4,1))
XLIM<-c(11750,300)
axis.ticks<-seq(13000,0,-1000)

sqs.rgb<-c(170,51,119)
ac.rgb<-c(34,136,51)
gain.rgb<-c(204,187,68)
loss.rgb<-c(102,204,238)
fad.rgb<-c(238,102,119)
lad.rgb<-c(68,119,170)

plot.ts(bcpk,bin.means.sub,YLAB="",new=F,add=F,line.col=rgb(ac.rgb[1],ac.rgb[2],ac.rgb[3],max=255),poly.col=rgb(ac.rgb[1],ac.rgb[2],ac.rgb[3],alpha=150,max=255),XLIM=XLIM,yaxt="n",bty="n",YLIM=c(0,25),hatched=F)
axis(2,at=c(0,5,10,15,20,25),label=F,las=1)
axis(4,at=c(0,5,10,15,20,25),label=F,las=1)
axis(2,at=c(5,15,25),label=T,las=1,cex.axis=1.4)
axis(3,at=axis.ticks,labels=F)
mtext("ACs",2,2.75,cex=1.2)

plot.ts(sqsk,bin.means.sub,YLAB="",new=F,add=F,line.col=rgb(sqs.rgb[1],sqs.rgb[2],sqs.rgb[3],max=255),poly.col=rgb(sqs.rgb[1],sqs.rgb[2],sqs.rgb[3],alpha=150,max=255),XLIM=XLIM,yaxt="n",bty="n",YLIM=c(4,7),hatched=F)
axis(2,at=c(4:7),labels=F,las=1,cex.axis=1.4)
axis(2,at=c(5,7),las=1,cex.axis=1.4)
axis(4,at=c(4:7),labels=F,las=1)
mtext("SQS",2,2.75,cex=1.2)

plot.ts(lossk,bin.means.sub,YLAB="",new=F,add=F,line.col=rgb(loss.rgb[1],loss.rgb[2],loss.rgb[3],max=255),poly.col=rgb(loss.rgb[1],loss.rgb[2],loss.rgb[3],alpha=150,max=255),YLIM=c(0,3),XLIM=XLIM,xaxt="n",yaxt="n",bty="n",hatched=F)
axis(2,at=seq(0,3,1),label=F,las=1)
axis(4,at=seq(0,3,1),label=F,las=1)
axis(2,at=seq(0,3,1),las=1,cex.axis=1.4)
plot.ts(gaink,bin.means.sub,YLAB="",new=F,add=T,line.col=rgb(gain.rgb[1],gain.rgb[2],gain.rgb[3],max=255),poly.col=rgb(gain.rgb[1],gain.rgb[2],gain.rgb[3],alpha=175,max=255),XLIM=XLIM,hatched=F)
mtext("Gain/Loss",2,2.75,cex=1.2)
legend("bottomleft",c("Gain","Loss"),col=c(rgb(gain.rgb[1],gain.rgb[2],gain.rgb[3],alpha=125,max=255),rgb(loss.rgb[1],loss.rgb[2],loss.rgb[3],alpha=175,max=255)),pch=15,bty="n",cex=1.5)

plot.ts(fadk,bin.means.sub,YLAB="",new=F,add=F,line.col=rgb(fad.rgb[1],fad.rgb[2],fad.rgb[3],max=255),poly.col=rgb(fad.rgb[1],fad.rgb[2],fad.rgb[3],alpha=150,max=255),YLIM=c(0,max(fadk,ladk,na.rm=T)),XLIM=XLIM,xaxt="n",yaxt="n",bty="n",hatched=T,angle=35)

axis(2,at=seq(0,4,2),las=1,cex.axis=1.4)
axis(2,at=seq(0,4,1),label=F,las=1)
axis(4,at=seq(0,4,1),label=F,las=1)
axis(1,at=seq(12000,0,-2000),cex.axis=1.4)
axis(1,at=axis.ticks, label=F)
plot.ts(ladk,bin.means.sub,YLAB="",new=F,add=T,line.col=rgb(lad.rgb[1],lad.rgb[2],lad.rgb[3],max=255),poly.col=rgb(lad.rgb[1],lad.rgb[2],lad.rgb[3],alpha=125,max=255),XLIM=XLIM,hatched=T,angle=-75)
mtext("FAD/LAD",2,2.75,cex=1.2)
legend("topleft",c("FAD","LAD"),col=c(rgb(fad.rgb[1],fad.rgb[2],fad.rgb[3],alpha=150,max=255),rgb(lad.rgb[1],lad.rgb[2],lad.rgb[3],alpha=125,max=255)),pch=15,bty="n",cex=1.5)
	
mtext("Cal YBP",1,2.5)

