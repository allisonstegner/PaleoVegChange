#################################################
# Code to choose Neotoma pollen datasets with Bayesian age models
# Copyright: 2023 M. Allison Stegner
#################################################

# LIBRARIES #######################################
# library(neotoma,quietly=T) # depreciated

# LOAD DATA #######################################
# Load age models___
# Bayesian age models used here are from Wang et al. (2019) Bayesian ages for pollen records since the last glaciation. Sci Data 6, 176
# code to generate these age models is available at https://github.com/yuewangpaleo/BaconAgeNeotoma
# Age models here are provided with permission from Y. Wang

wang.cores.meta<-read.csv("Wang et al Cores/SiteInfo_allsites.csv")
path<-"Wang-et-al-Cores/Cores_all"

# GET DATA #######################################
# set criteria for including sites based on sample number, resolution, etc.

# limit to sites with Wang et al. 2019 bacon age models
min.samp<-20 # minimum number of samples in time series, regardless of resolution
max.resolution<-250 # mean resolution (# samples/time must be this value or lower)
prop.length.cutoff<-0.1 # sites with hiatuses longer that this proportion will be flagged and info about time segments before and after each hiatus will be collected
min.amt.time<-2000 # total minimum duration for the entire time series. 

full.cores<-wang.cores.meta[wang.cores.meta$AgeModel_completeness %in% "full",] # age models for sites where the full dataset was modeled
part.cores<-wang.cores.meta[wang.cores.meta$AgeModel_completeness %in% "part",] # age models for sites where only part of the dataset was modeled

# get datasets for cores with full age models_____
handles<-as.vector(full.cores$handle)

keepers<-c()
for (i in 1:length(handles)){
	# print(i)
	handlei<-handles[i] # get site handle
	handle.path<-paste(path,"/",handlei,sep="") 
	am.file<-dir(handle.path,pattern="ages.txt") 
	am.path<-paste(handle.path,"/",am.file,sep="")
	am<-read.table(am.path,header=T)
	
	minage<-min(am$median,na.rm=T)
	maxage<-max(am$median,na.rm=T)	
	duration<-maxage-minage
	avg.res<-duration/nrow(am)

	sampsize.test<-(nrow(am)>=min.samp) # sufficient sample size?
	dur.test<-(duration>=min.amt.time) # sufficient duration?
	res.test<-(avg.res<=min.amt.time) # sufficient resolution?
	
	if (sum(sampsize.test,dur.test,res.test)==3) {
		keepers[i]<-handlei
	} else {
		keepers[i]<-NA
	}	
}

keepers<-keepers[!is.na(keepers)]

full.ds<-full.cores[full.cores$handle %in% keepers,c("datasetid","handle")]

# get datasets for cores with partial age models_____
handles<-as.vector(part.cores$handle)

keepers.p<-c()
for (i in 1:length(handles)){
	print(i)
	handlei<-handles[i] # get site handle
	handle.path<-paste(path,"/",handlei,sep="") 
	am.file<-dir(handle.path,pattern="ages.txt") 
	am.path<-paste(handle.path,"/",am.file,sep="")
	am<-read.table(am.path,header=T)
	
	minage<-min(am$median,na.rm=T)
	maxage<-max(am$median,na.rm=T)	
	duration<-maxage-minage
	avg.res<-duration/nrow(am)

	sampsize.test<-(nrow(am)>=min.samp) # sufficient sample size?
	dur.test<-(duration>=min.amt.time) # sufficient duration?
	res.test<-(avg.res<=min.amt.time) # sufficient resolution?
	
	if (sum(sampsize.test,dur.test,res.test)==3) {
		keepers.p[i]<-handlei
	} else {
		keepers.p[i]<-NA
	}	
}

keepers.p<-keepers.p[!is.na(keepers.p)]
part.ds<-part.cores[part.cores$handle %in% keepers.p,c("datasetid","handle")]

# combine full and partial cores list of dataset ids, and download_____
all.cores.ds<-c(full.ds[,1],part.ds[,1])

# get_download (from the noetoma package) is depreciated. An Rdata object of pol_dlx has been saved and posted with this github repo
# pol_dly<-get_download(all.cores.ds)   
# pol_dlx<-cull_spp_poor_sites(pol_dly,min.spp=5) # eliminate sites with too few species
