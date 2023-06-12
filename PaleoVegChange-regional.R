#################################################
# Code to analyze bidoversity and community change in pollen datasets from the Neotoma database
# in regional subsets
# Copyright: 2023 M. Allison Stegner
#################################################

# LIBRARIES #######################################
library(rgdal)
library(geosphere)

# LOAD DATA #######################################
# NEOTOMA R OBJECT___
load('PaleoVegChange/pol_dlx_2020-12-01.Rdata')

# LOAD AGE MODELS___
# Bayesian age models used here are from Wang et al. (2019) Bayesian ages for pollen records since the last glaciation. Sci Data 6, 176
# code to generate these age models is available at https://github.com/yuewangpaleo/BaconAgeNeotoma
# Age models here are provided with permission from Y. Wang

path<-"PaleoVegChange/Wang-et-al-Cores/Cores_all" # file path for age model data

# LOAD ECOREGIONS SHAPEFILE__
# Download the Level I Ecoregions of North America Shapefile available at:
# https://www.epa.gov/eco-research/ecoregions-north-america
ecoregions<-readOGR("na_cec_eco_l1/NA_CEC_Eco_Level1.shp") # replace with appropriate file path

# FUNCTIONS #####################
# LOAD FUNCTIONS__
source('PaleoVegChange/PaleoVegChange_functions.R', chdir = TRUE)

# ADDITIONAL FUNCTIONS FOR DEFIINING REGIONS___
extract.coords<-function(pol_dl_obj){
	coords<-matrix(NA,nc=3,nr=length(pol_dl_obj))
	for (i in 1:length(pol_dl_obj)){
		coords[i,]<-c(pol_dl_obj[[i]]$dataset$dataset.meta$dataset.id,pol_dl_obj[[i]]$dataset$site.data$long,pol_dl_obj[[i]]$dataset$site.data$lat)	
	}
	return(coords)
}

# point.in.region________________________________
point.in.region<-function(point.long,point.lat,regionk){
	in.region<-c()
	for (i in 1:length(regionk@polygons)){
		in.polygon<-c()
		in.hole<-c()
		for(j in 1:length(regionk@polygons[[i]]@Polygons)){
			coords.ij<-regionk@polygons[[i]]@Polygons[[j]]@coords
			in.polygon[j]<-point.in.polygon(point.long,point.lat,coords.ij[,1],coords.ij[,2])
			if (regionk@polygons[[i]]@Polygons[[j]]@hole==TRUE){
				in.hole[j]<-1
			} else {
				in.hole[j]<-0
			}	
		}
		temp<-cbind(in.hole,in.polygon)
		if (sum(c(temp[,2]-temp[,1]) %in% 1)>0 && sum(c(temp[,2]-temp[,1]) %in% 0)==0){
			in.region[i]<-1
		} else {
			in.region[i]<-0
		}	
	}

	if (sum(in.region)>0){ in.ecoregion<-1
	} else { in.ecoregion<-0 }	
	return(in.ecoregion)
}

# assign.point________________________________________________________
assign.point<-function(point.long,point.lat,ecoregions){
	
	region.list<-unique(ecoregions$NA_L1NAME)
	region.list<-region.list[-1]
	classify.pt<-matrix(NA,nr=length(region.list),nc=2)
	for (i in 1:length(region.list)){
		regionk<-ecoregions[ecoregions$NA_L1NAME==region.list[i],]
		inYN<-point.in.region(point.long,point.lat,regionk)
		classify.pt[i,]<-c(as.character(region.list[i]),inYN)
	}
	classify.pt<-as.data.frame(classify.pt)
	
	if (sum(as.numeric(as.vector(classify.pt[,2])))==0){
		
	# some sites fall beyond the EPA ecoregion shapefiles
	# those sites are assigned to the nearest polygon in this if statement
	
		min.dist3<-c()
		for (i in 1:length(region.list)){
			regionk<-ecoregions[ecoregions$NA_L1NAME==region.list[i],]
			min.dist2<-c()
			for (j in 1:length(regionk@polygons)){
				min.dist<-c()
				for (k in 1: length(regionk@polygons[[j]]@Polygons)){
					coords.ij<-regionk@polygons[[j]]@Polygons[[k]]@coords
					min.dist[k]<-min(distGeo(c(point.long,point.lat),coords.ij))
				}
				min.dist2[j]<-min(min.dist)	
			}
			min.dist3[i]<-min(min.dist2)
		}	
		out<-as.character(region.list[which.min(min.dist3)])
	} else {
		out<-as.character(classify.pt[which(classify.pt[,2]==1),1])
	}
	return(out)
}

# ANALYSES #######################################
# SORT SITES INTO ECOREGIONS___
ecoregions<-spTransform(ecoregions,"+proj=longlat +datum=WGS84")

site.coords<-extract.coords(pol_dlx)

site.names<-c()
for (i in 1:length(pol_dlx)){
	site.names[i]<-pol_dlx[[i]]$dataset$site.data$site.name	
}

ecoregs<-c()
for (k in 1:nrow(site.coords)){
	# print(k)
	ecoregs[k]<-assign.point(site.coords[k,2],site.coords[k,3],ecoregions)
}
# table(ecoregs)

site.meta<-cbind(site.coords[,1],site.names,site.coords[,2:3],ecoregs,rep(NA,nrow(site.coords)))
colnames(site.meta)<-c("dataset.id","site.name","long","lat","ecoregion","site.group")
site.meta<-as.data.frame(site.meta)

# classify EPA ecoregions into combined regions
groupA<-c("ARCTIC CORDILLERA","TUNDRA","TAIGA","HUDSON PLAIN")
groupB<-"NORTHERN FORESTS"
groupC<-c("NORTHWESTERN FORESTED MOUNTAINS","MARINE WEST COAST FOREST")
groupD<-"EASTERN TEMPERATE FORESTS"
groupE<-"GREAT PLAINS"
groupF<-c("MEDITERRANEAN CALIFORNIA","NORTH AMERICAN DESERTS","SOUTHERN SEMIARID HIGHLANDS","TEMPERATE SIERRAS")
groupG<-c("TROPICAL DRY FORESTS","TROPICAL WET FORESTS")

site.meta$site.group[which(site.meta$ecoregion %in% groupA)]<-"A" # populate site.group column in site.meta table
site.meta$site.group[which(site.meta$ecoregion %in% groupB)]<-"B"
site.meta$site.group[which(site.meta$ecoregion %in% groupC)]<-"C"
site.meta$site.group[which(site.meta$ecoregion %in% groupD)]<-"D"
site.meta$site.group[which(site.meta$ecoregion %in% groupE)]<-"E"
site.meta$site.group[which(site.meta$ecoregion %in% groupF)]<-"F"
site.meta$site.group[which(site.meta$ecoregion %in% groupG)]<-"G"

Aobj<-site.meta[which(site.meta$site.group %in% "A"),1] # list of neotoma db id numbers for sits in each group
Bobj<-site.meta[which(site.meta$site.group %in% "B"),1]
Cobj<-site.meta[which(site.meta$site.group %in% "C"),1]
Dobj<-site.meta[which(site.meta$site.group %in% "D"),1]
Eobj<-site.meta[which(site.meta$site.group %in% "E"),1]
Fobj<-site.meta[which(site.meta$site.group %in% "F"),1]
Gobj<-site.meta[which(site.meta$site.group %in% "G"),1]

sites.inbin<-list(Aobj,Bobj,Cobj,Dobj,Eobj,Fobj,Gobj)
regions.lab<-c("A","B","C","D","E","F","G")
names(sites.inbin)<-regions.lab

# COMPUTE DIVERSITY AND COMMUNITY CHANGE METRICS---
# Settings___
step<- -250
ageolder<-13950
bins<-seq(ageolder,step,step)
bin.means<-(bins[-length(bins)]+bins[-1])/2

# # test
# iter.siteresamp<-50
# am.iter<-10
# N.sites<-20

iter.siteresamp<-5000
am.iter<-1000
# N.sites<-10
# N.sites<-15
N.sites<-20
# N.sites<-25

# Abrupt changes___
threshold<-0.5 # bcp posterior probability threshold
eco.group<-c("TRSH","UPHE") # identify which ecological groups to include in pollen dataset

bcp.resamp.list<-c()
for (m in 1:length(sites.inbin)){
	print(m)
	geog.sub<-sites.inbin[[m]]	
	if (length(geog.sub)<N.sites) {
		bcp.resamp.list[[m]]<-NA
		next
	} else {
		# pull sites for region m
		pol_geog<-pol_dlx[geog.sub]
		
		# run multi AC
		multi.bcp<-list()
		for (i in 1:length(pol_geog)){
			# print(i)
			pol_ds<-pol_geog[[i]]
			am_posteriors<-get_am_posteriors(pol_ds)
			am_posteriors<-am_posteriors[complete.cases(am_posteriors),]
			prep.pol<-prep_pollen(pol_ds,eco.group,type="pct")
			
			na.inds<-complete.cases(prep.pol) # remove rows where pollen abundances are NA
			am_posteriors<-am_posteriors[na.inds,]
			prep.pol$prep.pol<-prep.pol$prep.pol[na.inds,]
			prep.pol$depth<-prep.pol$depth[na.inds]
	
			multi.bcp[[i]]<-ac_agejigger(prep.pol,am_posteriors,iters=am.iter,type="bcp",threshold)
		}

		names(multi.bcp)<-names(pol_geog)

		bcp.resamp<-resample.bin.tally(bins,multi.bcp,iter=iter.siteresamp,N.sites=N.sites)
		
		bcp.resamp.list[[m]]<-bcp.resamp	
	} 
}

# Loss/gain___
gl.out.list<-c()
for (m in 1:length(sites.inbin)){
	# print(m)
	geog.sub<-sites.inbin[[m]]
	if (is.na(geog.sub[1])){
		gl.out.list[[m]]<-NA
		next
	} else {
		# pull sites
		pol_geog<-pol_dlx[geog.sub]
		
		# run gain_loss
		gl.out<-gain_loss(pol_geog,iter.siteresamp=iter.siteresamp,iter.amresamp=am.iter,N.years.prior=500,N.sites=N.sites)
		gl.out.list[[m]]<-gl.out
		
		print(paste("gain loss calculation for region",m,"complete"))
	
	} 
}

# richness and FAD/LADs___
div.flad.list<-c()
for (m in 1:length(sites.inbin)){
	geog.sub<-sites.inbin[[m]]
	if (length(geog.sub)==0){
		gl.out.list[[m]]<-NA
		next
	} else {
		# pull sites
		pol_geog<-pol_dlx[geog.sub]
		
		# run diversity metrics		
		div.flad.run<-div.flad.in.bin(pol_geog,iter.siteresamp=iter.siteresamp,am.iter=am.iter,bins=bins,N.sites=N.sites,burnin=3)
			
		div.flad.list[[m]]<-div.flad.run
		
		print(paste("diversity calculation for bin",m,"complete"))
		
	}
}

# COMPOSITE PLOTS #######################################
# trim times series to start at the same time bin____
reg.id<-c("A","B","C","D","E","F","G")
bin.means<-bins[-1]

sqs.rgb<-c(170,51,119)
ac.rgb<-c(34,136,51)
gain.rgb<-c(204,187,68)
loss.rgb<-c(102,204,238) # 0,119,187
fad.rgb<-c(238,102,119) #238,51,119
lad.rgb<-c(68,119,170) #51,187,238

# for (i in 1:length(regions)){
for (i in 1:4){
	print(i)
	
	max.age<-14000
	XLIM<-c(14000,0)
	axis.ticks<-seq(14000,0,-2000)
	
	# Raw Richness
	# richnessk<-div.run[[6]]
	# richnessk<-richnessk[-cut.inds,]
	
	# SQS
	sqsk<-div.flad.list[[i]][[3]]
	sqs.ylim<-c(0,10)
	
	# BCP
	bcp.i<-bcp.resamp.list[[i]]
	bcpk<-t(bcp.i$bin.quants[c(1,3,5),])
	ac.ylim<-c(0,9)
	
	# Loss
	lossk<-gl.out.list[[i]][[6]]
	
	# Gain
	gaink<-gl.out.list[[i]][[5]]
	lossgain.ylim<-c(0,5)

	# FAD
	fadk<-div.flad.list[[i]][[1]]
		
	# LAD
	ladk<-div.flad.list[[i]][[2]]
	flad.ylim<-c(0,6)
	
	# trim to length
	cut.inds<-which(bin.means > max.age)
	
	if (length(cut.inds)==0){
		bin.means.sub<-bin.means
		sqsk<-sqsk
		bcpk<-bcpk
		lossk<-lossk
		gaink<-gaink
		fadk<-fadk
		ladk<-ladk
	} else {
		bin.means.sub<-bin.means[-cut.inds]
		sqsk<-sqsk[-cut.inds,]
		bcpk<-bcpk[-cut.inds,]
		lossk<-lossk[-cut.inds,]
		gaink<-gaink[-cut.inds,]
		fadk<-fadk[-cut.inds,]
		ladk<-ladk[-cut.inds,]
	}

	# PLOT_______
	dev.new(width=4,height=6)
	par(mfrow=c(4,1),mar=c(0.2,4,0.1,1),oma=c(4,1,4,1))
		plot.ts(bcpk,bin.means.sub,YLAB="",new=F,add=F,line.col=rgb(ac.rgb[1],ac.rgb[2],ac.rgb[3],max=255),poly.col=rgb(ac.rgb[1],ac.rgb[2],ac.rgb[3],alpha=150,max=255),YLIM=ac.ylim,XLIM=XLIM,yaxt="n",bty="n",hatched=F)
	axis(2,at=seq(0,8,2),label=T,las=1,cex.axis=1.4)
	axis(4,at=seq(0,8,2),label=F,las=1,cex.axis=1.4)
	axis(3,at=axis.ticks,labels=F)
	mtext("ACs",2,2.75,cex=1.2)
	plot.ts(sqsk,bin.means.sub,YLAB="",new=F,add=F,line.col=rgb(sqs.rgb[1],sqs.rgb[2],sqs.rgb[3],max=255),poly.col=rgb(sqs.rgb[1],sqs.rgb[2],sqs.rgb[3],alpha=150,max=255),YLIM=sqs.ylim,XLIM=XLIM,yaxt="n",bty="n",hatched=F)
	axis(2,at=seq(0,10,2),labels=T,las=1,cex.axis=1.4)
	axis(4,at=seq(0,10,2),labels=F,las=1,cex.axis=1.4)
	mtext("SQS",2,2.75,cex=1.2)

plot.ts(lossk,bin.means.sub,YLAB="",new=F,add=F,line.col=rgb(loss.rgb[1],loss.rgb[2],loss.rgb[3],max=255),poly.col=rgb(loss.rgb[1],loss.rgb[2],loss.rgb[3],alpha=150,max=255),XLIM=XLIM,YLIM=lossgain.ylim,xaxt="n",yaxt="n",bty="n",hatched=F)
	# axis(2,label=F,at=seq(0,0.2,0.1),las=1,las=1,cex.axis=1.4)
	axis(2,label=T,at=seq(0,5,1),las=1,cex.axis=1.4)
	axis(4,label=F,at=seq(0,5,1),las=1)
	plot.ts(gaink,bin.means.sub,YLAB="",new=F,add=T,line.col=rgb(gain.rgb[1],gain.rgb[2],gain.rgb[3],max=255),poly.col=rgb(gain.rgb[1],gain.rgb[2],gain.rgb[3],alpha=175,max=255),XLIM=XLIM,hatched=F)
	mtext("Gain/Loss",2,2.75,cex=1.2)
	
	# legend("topleft",c("Gain","Loss"),col=c(rgb(gain.rgb[1],gain.rgb[2],gain.rgb[3],alpha=125,max=255),rgb(loss.rgb[1],loss.rgb[2],loss.rgb[3],alpha=175,max=255)),pch=15,bty="n",cex=1.5)
	plot.ts(fadk,bin.means.sub,YLAB="",new=F,add=F,line.col=rgb(fad.rgb[1],fad.rgb[2],fad.rgb[3],max=255),poly.col=rgb(fad.rgb[1],fad.rgb[2],fad.rgb[3],alpha=150,max=255),YLIM=flad.ylim,XLIM=XLIM,xaxt="n",yaxt="n",bty="n",hatched=T,angle=35)

	axis(2,at=seq(0,6,2),las=1,cex.axis=1.4)
	axis(4,at=seq(0,6,2),labels=F,las=1,cex.axis=1.4)
	axis(1,at=seq(12000,0,-4000),cex.axis=1.4)
	axis(1,at=seq(14000,0,-1000), label=F)
	plot.ts(ladk,bin.means.sub,YLAB="",new=F,add=T,line.col=rgb(lad.rgb[1],lad.rgb[2],lad.rgb[3],max=255),poly.col=rgb(lad.rgb[1],lad.rgb[2],lad.rgb[3],alpha=125,max=255),XLIM=XLIM,hatched=T,angle=-75)
	mtext("FAD/LAD",2,2.75,cex=1.2)
	# legend("topleft",c("FAD","LAD"),col=c(rgb(fad.rgb[1],fad.rgb[2],fad.rgb[3],alpha=150,max=255),rgb(lad.rgb[1],lad.rgb[2],lad.rgb[3],alpha=125,max=255)),pch=15,bty="n",cex=1.5)
	
	mtext("Cal YBP",1,2.5)
	mtext(reg.id[i],3,2,cex=1.4,outer=T)

}

