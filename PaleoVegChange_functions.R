#################################################
# Functions for analyzing bidoversity and community change in pollen datasets from the Neotoma database
# Copyright: 2023 M. Allison Stegner
#################################################

# LIBRARIES #######################################
library(vegan,quietly=T)
# library(neotoma,quietly=T)
library(bcp,quietly=T)
library(princurve,quietly=T)

# FUNCTIONS #######################################
# get_am_posteriors___________________________________________________
# imports age model posterior estimates and depths for a Neotoma dataset using the collection handle
# pol_ds is a single Neotoma dataset

get_am_posteriors<-function(pol_ds){
	handlej<-pol_ds$dataset$dataset.meta$collection.handle # get site handle
	handle.path<-paste(path,"/",handlej,sep="") 
	file.name<-dir(handle.path,pattern="posteriorout.csv") # use handle to import bacon iterations
	depths_file<-dir(handle.path,pattern="depths.txt") # use handle to import bacon model depths

	amdepths<-read.table(paste(handle.path,"/",depths_file,sep=""))
	file.path<-paste(handle.path,"/",file.name,sep="")
	am_iters<-read.csv(file.path)
	am_posteriors<-cbind(amdepths,am_iters)
	return(am_posteriors=am_posteriors)
}

# clean_pollen___________________________________________________
# this function compiles pollen using Neotoma-defined pollen lists or selects user-specified ecological groups
# and optionally limits to abundant taxa, then calculates pollen percent or returns counts  
# pol_ds: a single site from a Neotoma download object
# type: If type="pct", return is proportion. Else, counts are returned.
# topN: integer; allows user to limit result to the topN most abundant taxa. Default (NULL) is to return all taxa.
# eco.group: Neotoma ecological group codes. e.g., UPHE=upland herbs
# list.name: Neotoma pre-defined pollen list name 
# returns a matrix of pollen counts or proportions

clean_pollen<-function(pol_ds,type,topN=NULL,eco.group,list.name=NULL){
			
			if (!is.null(list.name)){
				counts.subset <- compile_taxa(pol_ds,list.name)
				sitei.counts <- counts.subset$counts
				keep.taxa<-counts.subset$taxon.list[counts.subset$taxon.list$ecological.group %in% eco.group,"compressed"]
			} else {
				sitei.counts <- pol_ds$counts
				taxa<-pol_ds$taxon.list
				keep.taxa<-taxa[taxa$ecological.group %in% eco.group,"taxon.name"]
			}
			
			sitei.counts<-sitei.counts[,which(colnames(sitei.counts) %in% keep.taxa)]
			
			if (sum(colnames(sitei.counts) %in% "Other")>0){
				sitei.counts<-sitei.counts[,-which(colnames(sitei.counts)=="Other")]
			}
		
			if (is.numeric(topN)==T){
				abundance.order<-order(colSums(sitei.counts),decreasing=T)
				pollen<-sitei.counts[,abundance.order[1:topN]]
			} else {
				pollen<-sitei.counts
			}
				
			if (type=="pct") pollen<-pollen/rowSums(pollen,na.rm=TRUE)
			
			rownames(pollen)<-pol_ds$sample.meta$age
			# pollen<-pollen[which(rowSums(pollen)>0),]
				
			return(pollen)
	}

# prep_pollen___________________________________________________
# prepare pollen datset for analysis by calling the clean_pollen function
# ensure that age model and pollen dataset depths match
# label rownames as sample depths
# returns the final depths and pollen dataset

prep_pollen<-function(pol_ds,eco.group,type){
	# limit pollen to terrestrial: trees, shrubs, upland herbs
	pollen.i<-clean_pollen(pol_ds,type="pct",eco.grou=eco.group)
	
	# bring in the age model posteriors
	am_posteriors<-get_am_posteriors(pol_ds)
	am_posteriors<-am_posteriors[complete.cases(am_posteriors),]
	
	am.depth<-am_posteriors[,1]
	pol.depth<-pol_ds$sample.meta$depth
	# match the pollen count depths to the age model depths
	
	if (nrow(pol_ds$count)!=nrow(am_posteriors)){
		pol.inds<-match(am.depth,pol.depth)
		pol.matx<-pollen.i[pol.inds,] # limit to pollen samples with modeled ages
	} else {
		pol.matx<-pollen.i
	}
	
	depths<-am.depth
	rownames(pol.matx)<-depths
	
	rpollen.i<-pol.matx[c(nrow(pol.matx):1),] # reverse direction so oldest sample is on top
	# rpollen.i<-rpollen.i[complete.cases(rpollen.i),]
	
	depth<-as.numeric(rownames(rpollen.i))

	out<-list(prep.pol=rpollen.i,depth=depth)
	return(out)
}

# # singlesite_coniss ____________________________________________

# singlesite_coniss<-function(prep.pol,age=NULL){
	
	# rpollen.i<-prep.pol$prep.pol
	# depths<-prep.pol$depth
		
	# # if (is.null(age)){
		# # age<-as.numeric(rownames(rpollen.i)) # define age vector as depth
	# # } else {
		# # age=age
	# # }
		
	# # CONISS
	# dist<-dist(sqrt(rpollen.i),method="euclidean") # run CONISS on site i
	# coniss.pol<-chclust(dist,method="coniss") 
	# # bstick.test<-bstick(coniss.pol,plot=T,ng=nrow(rpollen.i)) # test for significant clusters
	# bstick.test<-bstick(coniss.pol,plot=F) # test for significant clusters
	# bstick.test<-bstick.test[complete.cases(bstick.test),] # black line should be above red line
	
	# # if (sum(bstick.test[,2]>bstick.test[,3])==0){
		# # coniss.depths<-NA
		# # coniss.inds<-NA
	# # } else 
	
	# if (bstick.test[1,2]<bstick.test[1,3]) {
		# coniss.depths<-NA
		# coniss.inds<-NA
	# } else {
			
		# sig.matx<-bstick.test[which(bstick.test[,2]>bstick.test[,3]),] # select only the significant clusters
		# k<-max(sig.matx[,1]) # find the max number of significant clusters
	
		# con.groups<-cutree(coniss.pol,k=k) # locate the significant cluster breaks in time
		# breaks<-c()
		# for (j in 1:(length(con.groups)-1)){
			# if (con.groups[j]!=con.groups[j+1]) breaks[j]<-1
			# else breaks[j]<-0
		# }

		# breaks2<-cbind(breaks,depths[-length(depths)])
		# coniss.inds<-which(breaks2[,1]==1)
		# coniss.depths<-breaks2[coniss.inds,2] # time of cluster breaks

	# }
	# out<-list(coniss.depths=coniss.depths,coniss.inds=coniss.inds)
	# return(out)
# }

# singlesite_bcp ____________________________________________
# run Bayesian changepoint analysis on a single pollen dataset and return depths of changpoints
# prep.pol: output of the prep_pollen function (contains two objects: prep.pol, a pollen data table, and dpeth, a depth vector)
# threshold: the posterior probability of a changepoint (0 to 1)

singlesite_bcp<-function(prep.pol,threshold){
	
	rpollen.i<-prep.pol$prep.pol
	depths<-prep.pol$depth
	
	# principal curve
	pc<-principal_curve(sqrt(rpollen.i))
		
	# BCP _________________________
	# bcp.out<-bcp(sqrt(rpollen.i),mcmc=10000,burnin=1000,w0=0.2,p0=0.05) # run bcp without principle curve
	bcp.out<-bcp(pc$lambda,mcmc=10000,burnin=1000,w0=0.2,p0=0.05) # run bcp
	
	bcp.locs<-which(bcp.out$posterior.prob>threshold) # find time points with BCP post prob > threshold
	
	if (length(bcp.locs)==0){
		bcp.inds<-NA
		bc.depths<-NA
	} else {
		bcp.inds<-bcp.locs
		bc.depths<-as.numeric(depths[bcp.inds]) # locate significant BCPs in time
	}
	
	out<-list(bc.depths=bc.depths,bcp.inds=bcp.inds)
	return(out)
}

# ac.agejigger_______________________________________________
# this function runs an abrupt change method (singlesite_bcp) and resamples the age model posterior estimates

ac_agejigger<-function(prep.pol,am_posteriors,iters,type,threshold=NULL){
	# pol_ds is a pollen download object
	# am.iters is a table of bacon posterior estimates, with many iterations
	# iters is number of times to re-sample am_iters and generate new gamls models

	# this step identifies which pollen samples have age estimates
	depths<-prep.pol$depth
		
	if (type=="coniss"){
		locs<-singlesite_coniss(prep.pol)$coniss.inds
	} else if (type=="bcp"){
		if (is.null(threshold)){
			print("threshold must be provided")
		}
		locs<-singlesite_bcp(prep.pol,threshold)$bcp.inds
	} else {
		print("type must be 'coniss' or 'bcp'")
		break
	}
	
	jigger<-matrix(NA,nrow=length(locs),ncol=iters)
	ami<-matrix(NA,nrow=length(depths),ncol=iters)
	am.ind<-c()
	for (j in 1:iters){
		samp.ind<-sample(2:ncol(am_posteriors),1)
		am.ind[j]<-samp.ind
		am<-rev(round(am_posteriors[,samp.ind]))
		ami[,j]<-am
		if (sum(is.na(locs))>0){
			jigger[,j]<-NA
		} else {
			jigger[,j]<-am[locs] # age estimate of the bayesian changepoints identified by singlesite_bcp
		}
	} 	

	out<-list(ami=ami,jigger=jigger,am.ind=am.ind)
	return(out)
}

# bin_tally _____________________________________________________________
# function to 
# acs: an output from ac_agejigger
# bins: time bins in which to tally abrupt changes.

bin_tally<-function(acs,bins){
	all_bin_inds_ac<-c()
	all_bin_inds_am<-c()
	for (k in 1:length(acs)){
		sitek<-acs[[k]] # choose a site
		
		# bin the abrupt changes
		ac_k<-sitek$jigger # choose the jiggered/re-sampled ACs
		iterk_acs<-ac_k[,sample(1:ncol(ac_k),1)] # randomly sample the iteration

		bin_older_ac<-c()
		for (j in 1:length(iterk_acs)){
			binTFac<-which(iterk_acs[j]>bins)
			if (length(binTFac)==0){
				bin_older_ac[j]<-NA
			} else {
				bin_older_ind_ac<-min(binTFac)-1
				bin_older_ac[j]<-unique(bin_older_ind_ac)
			}
		}
		all_bin_inds_ac<-c(all_bin_inds_ac,bin_older_ac)
		
		# bin the full age model
		am_k<-sitek$ami
		iterk_am<-am_k[,sample(1:ncol(am_k),1)]

		bin_older_am<-c()
		for (j in 1:length(iterk_am)){
			binTFam<-which(iterk_am[j]>bins)
			if (length(binTFam)==0){
				bin_older_am[j]<-NA
			} else {
				bin_older_ind_am<-min(binTFam)-1
				bin_older_am[j]<-unique(bin_older_ind_am)
			}
		}
		all_bin_inds_am<-c(all_bin_inds_am,bin_older_am)	
	}
		
	
	bin.tally_ac<-as.matrix(table(all_bin_inds_ac)) #how many ACs are in the bin?
	bin.tally_am<-as.matrix(table(all_bin_inds_am)) #how many sites available in each bin?
	
	out<-list(bins=bins,bin.tally_ac=bin.tally_ac,bin.tally_am=bin.tally_am)
	return(out)
}

# quants _________________________________________________________
# calculate a range of quantiles

quants<-function(x){ 
	return(quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm=T))
}

# resample.bin.tally _________________________________________________________
resample.bin.tally<-function(bins,acs,iter,N.sites){
	bintally.iter<-matrix(NA,nr=(length(bins)-1),nc=iter)
	
	resampled.ids.list<-list()
	for (k in 1:iter){	
		n.in.bin<-c()
		nACs.bin.i<-c()
		
		resampled.ids.tab.perbin<-array(NA,dim=c(N.sites,2,(length(bins)-1)),dimnames=list(c(1:N.sites),c("site.ids","am.inds"),c(1:(length(bins)-1))))
		for (i in 1:(length(bins)-1)) {
			in.bin<-c()
			choose.iterj<-c()
			acsj<-list()
			for (j in 1:length(acs)) {
				# choose the age model to use for each of the siies
				choose.iter<-sample(1:ncol(acs[[j]]$ami),size=1)
				choose.iterj[j]<-choose.iter 				
				sitej.iteri<-acs[[j]]$ami[,choose.iter]
				age.old<-max(sitej.iteri,na.rm=T)
				age.young<-min(sitej.iteri,na.rm=T)
				if ((age.old>bins[i] && age.young>bins[i]) || (age.old<bins[i+1] && age.young<bins[i+1])) {
					in.bin[j]<-NA
				} else {
					in.bin[j]<-names(acs)[j]
				}
			acsj[[j]]<-acs[[j]]$jigger[,choose.iter]
			}
		
			choose.iterj.inbin<-choose.iterj[!is.na(in.bin)]
			names(acsj)<-names(acs)
			in.bin<-in.bin[!is.na(in.bin)]
			n.in.bin[i]<-length(in.bin)
			acsj<-acsj[in.bin]
		
			if (length(in.bin) < N.sites){
				nACs.bin.i[i]<-NA
				
				resampled.ids.tab.perbin[,,i]<-matrix(NA, nrow=N.sites,ncol=2)
			} else {
				sample.ind<-sample(1:length(in.bin),N.sites) #set up which sites willbe sampled
				resampled.ids<-in.bin[sample.ind]
				resampled.ams<-choose.iterj.inbin[sample.ind]
				
				resampled.ids.tab.perbin[,,i]<-cbind(resampled.ids,resampled.ams)
				
				acsj.resamp<-acsj[resampled.ids]
				ac.TF<-c()
				for (l in 1:length(acsj.resamp)){
					acj.l<-acsj.resamp[[l]]
					ac.in.bin<-intersect(which(acj.l<bins[i]),which(acj.l>bins[i+1]))
					if (length(ac.in.bin)==0) ac.TF[l]<-length(ac.in.bin)
					else ac.TF[l]<-1
				}
				nACs.bin.i[i]<-sum(ac.TF)
			}					
		}
		resampled.ids.list[[k]]<-resampled.ids.tab.perbin
		bintally.iter[,k]<-nACs.bin.i
	}
	
	bin.quants<-apply(bintally.iter,1,quants)
	names(bin.quants)
	
	# resamples.sites is a list of arrays with 1 array per iteration
	# x= site id and age model for the nth site
	# y= column for site id and for age model index, specific to the site
	# z= the ith time bin
	out<-list(bin.quants=bin.quants,resamples.sites=resampled.ids.list,bintally.iter=bintally.iter)
	return(out)	
}

# gain_loss_________
# calculate taxa gained/lost in each sample proceeding forward in time
# pol_dl is a neotoma download object
# iter.siteresamp: number of times to resample sites
# iter.amresamp: number of times to resample the age model per site
# N.years.prior: the period of time (in years) over which to calculate gain loss. e.g. if N.years.prior is 500, calculate how many taxa appeared for the first time in the last 500 years (gain) or were lost after having been present for the last 500 years (loss)
# N.sites: number of sites to resample

gain_loss<-function(POL_DL,iter.siteresamp,iter.amresamp,N.years.prior,N.sites){

	eco.group<-c("TRSH","UPHE")

	sum.gain.per.site<-list()
	sum.loss.per.site<-list()
	pct.novel.per.site<-list()
	pct.lost.per.site<-list()
	for (k in 1:length(POL_DL)) { #for each site
		print(k)
		
		pol_ds<-POL_DL[[k]]
	
		# get age model posteriors
		am_posteriors<-get_am_posteriors(pol_ds)
		am_posteriors<-am_posteriors[complete.cases(am_posteriors),]

		# get pollen data
		prep.pol<-prep_pollen(pol_ds,eco.group,type="pct")
		pol.data<-prep.pol$prep.pol
		pol.pres.ab<-pol.data
		pol.pres.ab[which(pol.pres.ab>0)]<-1
	
		sum.gain.iter<-array(NA,dim=c((nrow(pol.pres.ab)-1),2,iter.amresamp),dimnames=list(c(1:(nrow(pol.pres.ab)-1)), c("am","pct.novel"),c(1:iter.amresamp)))
		sum.loss.iter<-array(NA,dim=c((nrow(pol.pres.ab)-1),2,iter.amresamp),dimnames=list(c(1:(nrow(pol.pres.ab)-1)), c("am","pct.novel"),c(1:iter.amresamp)))	

		pct.novel.iter<-array(NA,dim=c((nrow(pol.pres.ab)-1),2,iter.amresamp),dimnames=list(c(1:(nrow(pol.pres.ab)-1)), c("am","pct.novel"),c(1:iter.amresamp))) 
		pct.lost.iter<-array(NA,dim=c((nrow(pol.pres.ab)-1),2,iter.amresamp),dimnames=list(c(1:(nrow(pol.pres.ab)-1)), c("am","pct.lost"),c(1:iter.amresamp))) 


		for (i in 1:iter.amresamp) {
			# print(i)
			# choose an age model
			am.ind<-sample(1:ncol(am_posteriors),1)
			ami<-am_posteriors[,am.ind]
				
			# reverse age model so oldest samples come first; pollen is already in this format
			ami<-rev(ami)

			gain.sum<-c()
			loss.sum<-c()
			novel.pct<-c()
			lost.pct<-c()
			
			for (n in 2:(length(ami)-1)){
				# print(n)
				target.age<-ami[n]
				older.samps<-ami[1:(n-1)]
				samps.in.range<-older.samps[which(older.samps<(target.age+N.years.prior))]
				if (length(samps.in.range)==0){
					gain.sum[n]<-NA
					loss.sum[n]<-NA
					novel.pct[n]<-NA
					lost.pct[n]<-NA
				} else {
					# find the samples that fall within the last X years (N.years.prior)
					# test whether the pollen type has been seen in those samples
			
					ami.inrange.ind<-which(ami %in% samps.in.range)
					# ami[ami.inrange.ind]
					target.pollen<-pol.pres.ab[n,]
					
					#find which types are not present in the target sample
					zeros<-which(target.pollen==0) 					
					
					# combine the target sample with the sample within range
					pol.sub<-pol.pres.ab[c(n,ami.inrange.ind),]
					
					# eliminate cols where pollen type was not present in the target sample
					present.tab<-pol.sub[,-zeros]
					absent.tab<-as.matrix(pol.sub[,zeros])
				
					if ((ncol(pol.sub)-length(zeros))==1){
						if (present.tab[1]==present.tab[2]){
							novel.pct[n]<-0
						} else {
							novel.pct[n]<-1
						}
					} else {
						# add the column sums
						pol.sum<-colSums(present.tab)
						# if the value is greater than 1, then the pollen type is not new
						no.pol.sum<-colSums(absent.tab)
					
						# if 1, the pollen type did not occur in the previous X years
						novel.sum<-length(which(pol.sum==1)) # total number of novel types
						gain.sum[n]<-novel.sum
						novel.pct[n]<-novel.sum/ncol(present.tab) # percent novel types (out of all types found in the target sample)
						
						absent.sum<-length(which(no.pol.sum>1))
						loss.sum[n]<-absent.sum
						lost.pct[n]<-absent.sum/ncol(present.tab) # lost species as a proportion of species that are present
					}		
				}	
			}
			sum.gain.iter[,,i]<-cbind(ami[-1],gain.sum)	
			sum.loss.iter[,,i]<-cbind(ami[-1],loss.sum)	
			pct.novel.iter[,,i]<-cbind(ami[-1],novel.pct)	
			pct.lost.iter[,,i]<-cbind(ami[-1],lost.pct)	
		}
		sum.gain.per.site[[k]]<-sum.gain.iter
		sum.loss.per.site[[k]]<-sum.loss.iter
		pct.novel.per.site[[k]]<-pct.novel.iter
		pct.lost.per.site[[k]]<-pct.lost.iter
	}
	names(sum.gain.per.site)<-names(sum.loss.per.site)<-names(pct.novel.per.site)<-names(pct.lost.per.site)<-names(POL_DL)
	
	print("age model resampling complete")
	
	total.sum.gain.in.bin<-matrix(NA,nrow=(length(bins)-1),ncol=3)
	total.sum.loss.in.bin<-matrix(NA,nrow=(length(bins)-1),ncol=3)
	full.sum.gain.per.bin<-list()
	full.sum.loss.per.bin<-list()
	
	total.pct.nov.in.bin<-matrix(NA,nrow=(length(bins)-1),ncol=3)
	total.pct.lost.in.bin<-matrix(NA,nrow=(length(bins)-1),ncol=3)
	full.pct.nov.per.bin<-list()
	full.pct.lost.per.bin<-list()
	for (i in 1:(length(bins)-1)) { # for each bin
		print(paste("bin",i,"out of",length(bins)-1))
	
		sum.gain.out.list<-list()
		sum.gain.50quant<-c()
		sum.loss.out.list<-list()
		sum.loss.50quant<-c()
		
		pct.nov.out.list<-list()
		pct.nov.50quant<-c()
		pct.lost.out.list<-list()
		pct.lost.50quant<-c()
		for (j in 1:iter.siteresamp) { # for each iteration
			# print(j)
			choose.iterk<-c() # save which age model iteration you used
			in.bin<-c() # track which sites are in the bin
			
			sum.gain.k<-list()
			sum.loss.k<-list()
			
			pct.novel.k<-list()
			pct.lost.k<-list()
			for (l in 1:length(pct.novel.per.site)) { # for each site
				# print(l)
				# choose the age model to use for each of the sites
				choose.iter<-sample(1:dim(pct.novel.per.site[[l]])[3],size=1)
				choose.iterk[l]<-choose.iter 				
				amk<-pct.novel.per.site[[l]][,1,choose.iter]
			
				in.bin.ind<-intersect(which(amk<bins[i]), which(amk>bins[i+1]))
						
				if (length(in.bin.ind)==0){
					in.bin[l]<-NA
				} else {
					in.bin[l]<-names(pct.novel.per.site)[l] # save the id if its in the bin
				}
				sum.gain.k[[l]]<-cbind(sum.gain.per.site[[l]][,,choose.iter]) # save the pct novel data for the site for the given age model
				sum.loss.k[[l]]<-cbind(sum.loss.per.site[[l]][,,choose.iter]) # save the pct novel data for the site for the given age model
				
				pct.novel.k[[l]]<-cbind(pct.novel.per.site[[l]][,,choose.iter]) # save the pct novel data for the site for the given age model
				pct.lost.k[[l]]<-cbind(pct.lost.per.site[[l]][,,choose.iter]) # save the pct novel data for the site for the given age model
			}
		
			choose.iterk.inbin<-choose.iterk[!is.na(in.bin)]
			
			names(sum.gain.k)<-names(sum.gain.per.site)
			names(sum.loss.k)<-names(sum.loss.per.site)
			
			names(pct.novel.k)<-names(pct.novel.per.site) # name with dataset ids
			names(pct.lost.k)<-names(pct.lost.per.site)
			
			in.bin<-in.bin[!is.na(in.bin)]
			# n.in.bin[i]<-length(in.bin)
			
			sum.gaink<-sum.gain.k[in.bin] #trim the sites that weren't in the bin
			sum.lossk<-sum.loss.k[in.bin]
			
			pct.novk<-pct.novel.k[in.bin] #trim the sites that weren't in the bin
			pct.lostk<-pct.lost.k[in.bin]
		
			if (length(pct.novk)<N.sites){
				# total.pct.nov.in.bin[i]<-NA
				
				sum.gain.out.list[[j]]<-NA
				sum.gain.50quant[j]<-NA
				
				sum.loss.out.list[[j]]<-NA
				sum.loss.50quant[j]<-NA
				
				pct.nov.out.list[[j]]<-NA
				pct.nov.50quant[j]<-NA
				
				pct.lost.out.list[[j]]<-NA
				pct.lost.50quant[j]<-NA
				
			} else { 
		
				samp.ind<-sample(1:length(pct.novk),N.sites)
				
				sum.gain.samp<-sum.gaink[samp.ind]
				sum.loss.samp<-sum.lossk[samp.ind]
				
				pct.nov.samp<-pct.novk[samp.ind]
				pct.lost.samp<-pct.lostk[samp.ind]
				# names(pct.nov.samp)<-names(pct.novk)[samp.ind]
		
				sum.gain.out<-c()
				sum.loss.out<-c()
				pct.nov.out<-c()
				pct.lost.out<-c()
				for (m in 1:length(pct.nov.samp)){
					# print(m)
					sum.gain.samp.m<-sum.gain.samp[[m]]
					sum.loss.samp.m<-sum.loss.samp[[m]]
					
					pct.nov.samp.m<-pct.nov.samp[[m]]
					pct.lost.samp.m<-pct.lost.samp[[m]]
					
					am<-pct.nov.samp.m[,1]
			
					in.bin.ind<-intersect(which(am<bins[i]), which(am>bins[i+1]))
			
					sum.gain.temp<-sum.gain.samp.m[in.bin.ind,2]
					sum.loss.temp<-sum.loss.samp.m[in.bin.ind,2]
					
					pct.nov.temp<-pct.nov.samp.m[in.bin.ind,2]
					pct.lost.temp<-pct.lost.samp.m[in.bin.ind,2]
			
					if (length(pct.nov.temp)>1){
						# sum.gain.avg<-median(sum.gain.temp)
						# sum.gain.avg<-sum(sum.gain.temp)/length(sum.gain.temp)
						sum.gain.avg<-quantile(sum.gain.temp,0.5,na.rm=T)
						sum.gain.out[m]<-sum.gain.avg
						
						pct.nov.avg<-sum(pct.nov.temp)/length(pct.nov.temp) #if there is more than one sample in the bin, average the two values
						pct.nov.out[m]<-pct.nov.avg
					} else {
						sum.gain.out[m]<-sum.gain.temp
						pct.nov.out[m]<-pct.nov.temp
					}
				
					if (length(pct.lost.temp)>1){
						# sum.loss.avg<-median(sum.loss.temp)
						# sum.loss.avg<-sum(sum.loss.temp)/length(sum.loss.temp)
						sum.loss.avg<-quantile(sum.loss.temp,0.5,na.rm=T)
						sum.loss.out[m]<-sum.loss.avg
						
						pct.lost.avg<-sum(pct.lost.temp)/length(pct.lost.temp) #if there is more than one sample in the bin, average the values
						pct.lost.out[m]<-pct.lost.avg
					} else {
						sum.loss.out[m]<-sum.loss.temp
						pct.lost.out[m]<-pct.lost.temp
					}
				}
				# print(m,length(sum.gain.out))
				names(sum.gain.out)<-names(sum.gain.samp)
				names(sum.loss.out)<-names(sum.loss.samp)
				
				names(pct.nov.out)<-names(pct.nov.samp)
				names(pct.lost.out)<-names(pct.lost.samp)
			
				sum.gain.out.list[[j]]<-sum.gain.out
				sum.loss.out.list[[j]]<-sum.loss.out
				
				pct.nov.out.list[[j]]<-pct.nov.out
				pct.lost.out.list[[j]]<-pct.lost.out
		
				sum.gain.50quant[j]<-quantile(sum.gain.out,0.5,na.rm=T)
				sum.loss.50quant[j]<-quantile(sum.loss.out,0.5,na.rm=T)
				
				pct.nov.50quant[j]<-quantile(pct.nov.out,0.5,na.rm=T)
				pct.lost.50quant[j]<-quantile(pct.lost.out,0.5,na.rm=T)
			}			
		}	
		
		total.sum.gain.in.bin[i,]<-quantile(sum.gain.50quant,c(0.25,0.5,0.75),na.rm=T)
		full.sum.gain.per.bin[[i]]<-sum.gain.out.list

		total.sum.loss.in.bin[i,]<-quantile(sum.loss.50quant,c(0.25,0.5,0.75),na.rm=T)
		full.sum.loss.per.bin[[i]]<-sum.loss.out.list

		total.pct.nov.in.bin[i,]<-quantile(pct.nov.50quant,c(0.25,0.5,0.75),na.rm=T)
		full.pct.nov.per.bin[[i]]<-pct.nov.out.list
	
		total.pct.lost.in.bin[i,]<-quantile(pct.lost.50quant,c(0.25,0.5,0.75),na.rm=T)
		full.pct.lost.per.bin[[i]]<-pct.lost.out.list
	}
	
	out<-list(total.pct.nov.in.bin=total.pct.nov.in.bin,full.pct.nov.per.bin=full.pct.nov.per.bin,total.pct.lost.in.bin=total.pct.lost.in.bin,full.pct.lost.per.bin=full.pct.lost.per.bin,total.sum.gain.in.bin=total.sum.gain.in.bin,total.sum.loss.in.bin=total.sum.loss.in.bin)
	return(out)
}

# sqs___________________________
# from Alroy, J (2010) Fair sampling of taxonomic richness and unbiased estimation of orgination and extinction rates. The Paleontological Society Papers, 16: 55-80.
# available at https://bio.mq.edu.au/~jalroy/SQS.html

sqs<-function(ab,q,trials,method,ignore.singletons,dominant)	{

	params <- array(data=NA,dim=0) # MAS - edited from original to remove dimnames
	if (missing(trials))	{
		trials <- 100
	}
	if (missing(method))	{
		method <- ""
	} else if (method != "" && method != "rarefaction" && method != "CR")	{
		return(print('If the method is rarefaction enter method="rarefaction" or "CR"',quote=F))
	}
	if ((q <= 0 || q >= 1) && method != "rarefaction" && method != "CR")	{
		
		return(print("If the method is SQS the quota must be greater than zero and less than one",quote=F))
	} else if (q < 1 && (method == "rarefaction" || method == "CR"))	{
		
		return(print("If the method is rarefaction the quota must be an integer",quote=F))
	}
	ignore <- 0
	if (! missing(ignore.singletons) && (ignore.singletons == T || ignore.singletons == "yes" || ignore.singletons == "y"))	{
		ignore <- 1
	}
	if (missing(dominant))	{
		dominant <- 0
	} else if (dominant != "" && dominant != "exclude" && dominant != "no")	{
		return(print('To exclude the dominant taxon, enter dominant="exclude" or "no"',quote=F))
	}

	# compute basic statistics
	specimens <- sum(ab)
	singletons <- 0
	doubletons <- 0
	highest <- 0
	for (i in 1:length(ab))	{
		if (ab[i] == 1)	{
			singletons <- singletons + 1
		} else if (ab[i] == 2)	{
			doubletons <- doubletons + 1
		}
		if (ab[i] > highest)	{
			highest <- ab[i]
			mostfrequent <- i
		}
	}
	u <- 1 - singletons / specimens
	# optionally, exclude the dominant taxon
	if (dominant == "exclude" || dominant == "no")	{
		u <- 1 - singletons / (specimens - highest)
	}

	if (u == 0)	{
		return(print("Coverage is zero because all taxa are singletons",quote=F))
	}

	# compute raw taxon frequencies (temporarily)
	freq <- ab / specimens

	# standard recursive equation for Fisher's alpha
	alpha <- 10
	oldalpha <- 0
	while (abs(alpha - oldalpha) > 0.0000001)	{
		oldalpha <- alpha
		alpha <- length(ab) / log(1 + specimens/alpha)
	}

	params["raw richness"] <- length(ab)
	params["Good's u"] <- u
	params["subsampled richness"] <- NA
	params["subsampled u"] <- NA
	params["Chao 1"] <- length(ab) + singletons**2/(2* doubletons)
	params["subsampled Chao 1"] <- NA
	# governing parameter of the geometric series distribution
	params["k"] <- abs(lm(log(sort(freq)) ~ c(1:length(freq)))$coefficients[2])
	params["Fisher's alpha"] <- alpha
	params["Shannon's H"] <- -1 * sum(freq * log(freq))
	params["Hurlbert's PIE"] <- (1 - sum(freq**2)) * length(ab) / (length(ab) - 1)
	params["dominance"] <- highest / specimens
	params["specimens"] <- specimens
	params["singletons"] <- singletons
	params["doubletons"] <- doubletons
	params["specimens drawn"] <- 0
	

	if (dominant != "exclude" && dominant != "no")	{
		highest <- 0
		mostfrequent <- 0
	}

	# return if the rarefaction quota is equal to or higher than the
	#  specimen count
	if (method == "rarefaction" && q >= specimens - highest)	{
		return(params)
	}

	# compute taxon frequencies (tweak added in version 3.1)
	freq <- ab - (singletons + doubletons / 2) / length(ab)
	freq <- freq / (specimens - highest)

	# return if the quorum target is higher than estimated coverage
	if ((q > sum(freq) && method != "rarefaction" && method != "CR") || (q >= sum(ab)))	{
		return(params)
	}

	# create an array in which each cell corresponds to one specimen
	ids <- array()
	n <- 0
	for (i in 1:length(ab))	{
		for (j in 1:ab[i])	{
			n <- n + 1
			ids[n] <- i
		}
	}

	# subsampling trial loop
	# s will be the subsampled taxon count
	s <- array(rep(0,trials))
	drawn <- array(rep(0,trials))
	mostfrequentdrawn <- array(rep(0,trials))
	subsingle <- array(rep(0,trials))
	subdouble <- array(rep(0,trials))
	subchao <- array(rep(0,trials))

	for (trial in 1:trials)	{
		pool <- ids
		left <- length(pool)
		seen <- array(data=rep(0,length(ab)))
		subfreq <- array(rep(0,length(ab)))
		if (method != "rarefaction" && method != "CR")	{
			udrawn <- 0
			# equation new to version 3.0
			# the exponent corrects for downwards bias
			while (udrawn < q)	{
			# draw a specimen
				x <- floor(runif(1,min=1,max=left+1))
			# add to frequency and taxon sums if species has
			#  not been drawn previously
				subfreq[pool[x]] <- subfreq[pool[x]] + 1
				if (seen[pool[x]] == 0)	{
					if (pool[x] != mostfrequent && (ignore == 0 || ab[pool[x]] > 1))	{
						udrawn <- udrawn + freq[pool[x]]
					}
					seen[pool[x]] <- 1
			# randomly throw back some draws that put the sum over q
			#  (an even better algorithm added in version 3.0)
					if (runif(1) <= freq[pool[x]] || udrawn < q)	{
						s[trial] <- s[trial] + 1
					} else	{
						subfreq[pool[x]] <- subfreq[pool[x]] - 1
					}
				}
			# decrease pool of specimens not yet drawn
				pool[x] <- pool[left]
				left <- left - 1
			}
		} else	{
			i <- 0
			draws <- 0
			while (i < q)	{
				draws <- draws + 1
				x <- floor(runif(1,min=1,max=length(ids)-draws+2))
				subfreq[pool[x]] <- subfreq[pool[x]] + 1
				if (pool[x] != mostfrequent)	{
					i <- i + 1
				}
				if (seen[pool[x]] == 0)	{
					seen[pool[x]] <- 1
					s[trial] <- s[trial] + 1
				}
				pool[x] <- pool[length(ids)-draws+1]
			}
		}
		for (i in 1:length(ab))	{
			if (subfreq[i] == 1 && i != mostfrequent)	{
				subsingle[trial] <- subsingle[trial] + 1
			} else if (subfreq[i] == 2 && i != mostfrequent)	{
				subdouble[trial] <- subdouble[trial] + 1
			}
		}
		if (subsingle[trial] > 0 && subdouble[trial] > 0)	{
			subchao[trial] <- s[trial] + subsingle[trial]**2/(2*subdouble[trial])
		} else	{
			subchao[trial] <- s[trial]
		}
		drawn[trial] <- sum(subfreq)
		if (mostfrequent != 0)	{
			mostfrequentdrawn[trial] <- subfreq[mostfrequent]
		}
	}

	# compute vectors of non-zero counts
	options(warn=-1)
	s2 <- sort(sqrt(s-1))^2+1
	d2 <- sort(sqrt(drawn-1))^2+1
	m2 <- sort(sqrt(mostfrequentdrawn-1))^2+1
	ss2 <- sort(sqrt(subsingle-1))^2+1
	options(warn=0)

	# compute geometric means
	params["subsampled richness"] <- exp(mean(log(s2))) * length(s2)/length(s)
	params["specimens drawn"] <- exp(mean(log(d2))) * length(d2)/length(drawn)
	meanmost <- 0
	if (sum(mostfrequentdrawn) > 0)	{
		meanmost <- exp(mean(log(m2))) * length(m2)/length(mostfrequentdrawn)
	}
	meansubsingle <- exp(mean(log(ss2))) * length(ss2)/length(subsingle)

	params["subsampled u"] <- 1 - meansubsingle / (params["specimens drawn"] - meanmost)
	params["subsampled Chao 1"] <- exp(mean(log(subchao)))
	
	#return(params)
	s.temp<-exp((log(s2))) * length(s2)/length(s)
	std.err<-sd(s.temp)/sqrt(length(s.temp))
	
	#out<-list(params,std.err)
	#return(out)
	return(params)
}

# site.diversity_____________________________________
# calculate raw richness and standardized richness (SQS) for each sample for each site
site.diversity<-function(POL_DL,am.iter){
	eco.group<-c("TRSH","UPHE")

	richness.matx<-list()
	sqs.matx<-list()
	am.list<-list()
	for (k in 1:length(POL_DL)) { #for each site
		print(k)
			
		pol_ds<-POL_DL[[k]]
	
		# get pollen data
		prep.pol<-prep_pollen(pol_ds,eco.group,type="pct")
		
		na.inds<-complete.cases(prep.pol$prep.pol)
		pol.data<-prep.pol$prep.pol[na.inds,]
		pol.pres.ab<-pol.data
		pol.pres.ab[which(pol.pres.ab>0)]<-1
		
		# get age model posteriors
		am_posteriors<-get_am_posteriors(pol_ds)
		am_posteriors<-am_posteriors[complete.cases(am_posteriors),]
		am_posteriors<-am_posteriors[na.inds,]
		
		# calculate raw richness for each time point
		raw.richness<-apply(pol.pres.ab,1,sum)
		richness.matx[[k]]<-raw.richness
		
		# calculate sqs richness for each time point
		sqs.richness<-c()
		for (i in 1:nrow(pol.data)){
			pol.temp<-pol.data[i,]*100
			pol.temp<-pol.temp[which(pol.temp>0)]
			sqs.calc<-sqs(pol.temp,q=0.8)
			sqs.richness[i]<-sqs.calc[["subsampled richness"]]
		}
		sqs.matx[[k]]<-sqs.richness
		
		# resample the age model
		am.matx<-matrix(NA,nrow=nrow(pol.data),ncol=am.iter)
		for (j in 1:am.iter){
			am.ind<-sample(1:ncol(am_posteriors),1)
			ami<-am_posteriors[,am.ind]
			ami<-rev(ami) # reverse age model so oldest samples come first; pollen is already in this format
			am.matx[,j]<-ami
		}
		am.list[[k]]<-am.matx
	}
	
	out<-list(richness=richness.matx,sqs=sqs.matx,am.list=am.list)
}

# site.flad___________________________________________
# calculate FADs and LADs for each sample for each site

site.flad<-function(POL_DL,burnin,am.iter){
	# Note: there should be no first occurences in the first (burnin) sample
	# in the sample after the initial burnin, we check whether the pollen type was already present
	# Note: there should be no last occurences in the last (burnin) sample
	# in the sample before the final burnin, we check whether the pollen type continues to be present

	eco.group<-c("TRSH","UPHE")
	
	FAD.list<-list()
	LAD.list<-list()
	am.list<-list()
	for (k in 1:length(POL_DL)) { #for each site
		print(k)

		pol_ds<-POL_DL[[k]] # choose site k
	
		# get age model posteriors
		am_posteriors<-get_am_posteriors(pol_ds)
		am_posteriors<-am_posteriors[complete.cases(am_posteriors),]
						
		# get pollen data
		prep.pol<-prep_pollen(pol_ds,eco.group,type="pct")
		pol.data<-prep.pol$prep.pol
		pol.pres.ab<-pol.data
		pol.pres.ab[which(pol.pres.ab>0)]<-1
		
		# get first occurences
		FAD<-rep(NA,length=nrow(pol.pres.ab))
		for (i in burnin:nrow(pol.pres.ab)){
			fad<-intersect(which(pol.pres.ab[i,]==1),which(colSums(pol.pres.ab[c(1:i-1),])==0))
			if (length(fad)==0){ FAD[i]<-0
			} else { FAD[i]<-length(fad) }
		}
		
		# get last occurences
		LAD<-rep(NA,length=nrow(pol.pres.ab))
		for (i in 1:(nrow(pol.pres.ab)-burnin)){
			if (i %in% nrow(pol.pres.ab):(nrow(pol.pres.ab)-burnin)) {
				LAD[i]<-NA
			} else {
				lad<-intersect(which(pol.pres.ab[i,]==1),which(colSums(as.matrix(pol.pres.ab[c((i+1):nrow(pol.pres.ab)),]))==0))
			}
				
			if (length(lad)==0){ LAD[i]<-0
			} else { LAD[i]<-length(lad) }		
		}
		
		# generate a list object for FADs and LADs per site
		FAD.list[[k]]<-FAD
		LAD.list[[k]]<-LAD
		
		# resample the age model
		am.matx<-matrix(NA,nrow=nrow(pol.data),ncol=am.iter)
		for (j in 1:am.iter){
			am.ind<-sample(1:ncol(am_posteriors),1)
			ami<-am_posteriors[,am.ind]
			ami<-rev(ami) # reverse age model so oldest samples come first; pollen is already in this format
			am.matx[,j]<-ami
		}
		am.list[[k]]<-am.matx
		
	}	

	out<-list(FAD.list=FAD.list,LAD.list=LAD.list,am.list=am.list)
	return(out)
}

# div.flad.in.bin_____________________
# calculate FADs, LADs, and diversity (raw and SQS) in bins

div.flad.in.bin<-function(POL_DL,iter.siteresamp,am.iter,bins,N.sites,burnin){
	# run the diversity estimates for each site	
	site.div<-site.diversity(POL_DL,am.iter)
	
	s.richness<-site.div$richness
	s.sqs<-site.div$sqs
	s.ams<-site.div$am.list
	names(s.sqs)<-names(POL_DL)
	names(s.ams)<-names(POL_DL)
	
	print("site-level diversity metrics complete")
	
	# run flad
	site.flad.data<-site.flad(POL_DL,burnin=burnin,am.iter=am.iter)
	
	FAD.list<-site.flad.data$FAD.list
	LAD.list<-site.flad.data$LAD.list
	am.list<-site.flad.data$am.list
	names(FAD.list)<-names(POL_DL)
	names(LAD.list)<-names(POL_DL)
	names(am.list)<-names(POL_DL)
	
	print("Site FADs and LADs complete, starting age model iterations and site resampling")
	
	# set up empty objects
	quants.fadk<-matrix(NA,nrow=(length(bins)-1),ncol=3)
	quants.ladk<-matrix(NA,nrow=(length(bins)-1),ncol=3)
	full.fadk.per.bin<-list()
	full.ladk.per.bin<-list()
	
	quants.sqs.in.bin<-matrix(NA,nrow=(length(bins)-1),ncol=3)
	full.sqs.per.bin<-list()
	for (i in 1:(length(bins)-1)) { # for each bin
		# print(i)
		
		fadk.list<-c()
		fadk.50quant<-c()
		ladk.list<-c()
		ladk.50quant<-c()
		
		sqsk.list<-c()
		sqsk.50quant<-c()
		for (j in 1:iter.siteresamp) { # for each iteration
			# choose age model resample
			choose.iterk<-c() # save which age model iteration you used
			in.bin<-c() # track which sites are in the bin
			sqsk<-c()
			amk<-c()
			fadk<-c()
			ladk<-c()

			for (k in 1:length(POL_DL)){
				choose.iter<-sample(1:am.iter,size=1)
				choose.iterk[k]<-choose.iter
				am<-s.ams[[k]][,choose.iter]
				amk[[k]]<-am
				in.bin.ind<-intersect(which(am<bins[i]), which(am>bins[i+1]))
				if (length(in.bin.ind)==0){
					in.bin[k]<-NA
				} else {
					in.bin[k]<-names(POL_DL)[k] # save the id if its in the bin
				}
				sqsk[[k]]<-s.sqs[[k]]
				fadk[[k]]<-FAD.list[[k]]
				ladk[[k]]<-LAD.list[[k]]		
			}
			names(amk)<-names(POL_DL)
					
			in.bin<-in.bin[!is.na(in.bin)]
			
			fadk.sites<-FAD.list[in.bin]
			ladk.sites<-LAD.list[in.bin]
			sqsk.sites<-s.sqs[in.bin]
			amk.sites<-amk[in.bin]
			
			if (length(in.bin)<N.sites){
				fadk.list[[j]]<-NA
				fadk.50quant[j]<-NA
				
				ladk.list[[j]]<-NA
				ladk.50quant[j]<-NA
				
				sqsk.list[[j]]<-NA
				sqsk.50quant[j]<-NA
				
			} else { 
				samp.ind<-sample(1:length(sqsk.sites),N.sites)
				
				fadk.samp<-fadk.sites[samp.ind]
				ladk.samp<-ladk.sites[samp.ind]
				sqsk.samp<-sqsk.sites[samp.ind]
				amk.samp<-amk.sites[samp.ind]

				fadk.out<-c()
				ladk.out<-c()
				sqsk.out<-c()
				for (m in 1:N.sites){
					# print(m)
					fadk.samp.m<-fadk.samp[[m]]
					ladk.samp.m<-ladk.samp[[m]]
					sqsk.samp.m<-sqsk.samp[[m]]
					amk.samp.m<-amk.samp[[m]]
			
					in.bin.ind<-intersect(which(amk.samp.m<bins[i]), which(amk.samp.m>bins[i+1]))
			
					fadk.temp<-fadk.samp.m[in.bin.ind]
					ladk.temp<-ladk.samp.m[in.bin.ind]
					sqsk.temp<-sqsk.samp.m[in.bin.ind]
			
					if (length(in.bin.ind)>1){
						fadk.out1<-sum(fadk.temp) # total per bin		
						ladk.out1<-sum(ladk.temp)
						sqsk.out1<-sum(sqsk.temp)/length(sqsk.temp) # average per bin
						# sqsk.out1<-median(sqsk.temp,na.rm=T) # average per bin
					} else {
						fadk.out1<-fadk.temp # total per bin			
						ladk.out1<-ladk.temp
						sqsk.out1<-sqsk.temp
					}
					
					fadk.out[m]<-fadk.out1		
					ladk.out[m]<-ladk.out1 #
					
					sqsk.out[m]<-sqsk.out1
					# names(sqsk.out)<-names(sqsk.samp)

					fadk.list[[j]]<-fadk.out1
					ladk.list[[j]]<-ladk.out1
		
					fadk.50quant[j]<-quantile(fadk.out1,0.5,na.rm=T)
					ladk.50quant[j]<-quantile(ladk.out1,0.5,na.rm=T)
			
					sqsk.list[[j]]<-sqsk.out
					
					sqsk.50quant[j]<-quantile(sqsk.out,0.5,na.rm=T)
				}
			}			
		
		quants.fadk[i,]<-quantile(fadk.50quant,c(0.25,0.5,0.75),na.rm=T)
		# full.fadk.per.bin[[i]]<-fadk.list	
	
		quants.ladk[i,]<-quantile(ladk.50quant,c(0.25,0.5,0.75),na.rm=T)
		# full.ladk.per.bin[[i]]<-ladk.list
		
		quants.sqs.in.bin[i,]<-quantile(sqsk.50quant,c(0.25,0.5,0.75),na.rm=T)
		# full.sqs.per.bin[[i]]<-sqsk.list	
		
		}
	}
	
	out<-list(quants.fadk=quants.fadk,quants.ladk=quants.ladk,quants.sqs.in.bin=quants.sqs.in.bin)
	return(out)
}


# plot.ts___________________________________________________________________
plot.ts<-function(quantile.tab,bin.means,YLAB,new,add,line.col,poly.col,YLIM=NULL,XLIM=NULL,xaxt=NULL,yaxt=NULL,bty=NULL,hatched,angle=NULL,xaxs=NULL,yaxs=NULL){
	# new: should a new quartz device be opened?
	# add: should a new plot be drawn, or should a line be added to an existing plot?
	
	if (new==T){dev.new(width=10,height=4)}
	
	temp.matx<-cbind(quantile.tab,bin.means)
	temp.cc<-temp.matx[complete.cases(temp.matx),]
	
	if (add==F) {
		if (is.null(YLIM)){ ylim<-c(min(temp.cc[,1:3],na.rm=T),max(temp.cc[,1:3],na.rm=T))
		} else { ylim<-YLIM }
		
		if (is.null(XLIM)){ xlim<-c(min(temp.cc[,4],na.rm=T),max(temp.cc[,4],na.rm=T))
		} else { xlim<-XLIM }
		
		if (is.null(xaxt)){ xaxt="n"
		} else { xaxt=xaxt }
		
		if (is.null(yaxt)){ yaxt="n"
		} else { yaxt=yaxt }
		
		if (is.null(bty)){ bty="n"
		} else { bty=bty }
		
		if (is.null(xaxs)){ xaxt="n"
		} else { xaxs="i" }
	
		if (is.null(yaxs)){ yaxt="n"
		} else { yaxs="i" }

		
plot(temp.cc[,4],temp.cc[,2],type="l",ylim=ylim,xlim=xlim,xlab="",ylab=YLAB,las=1,xaxt=xaxt,yaxt=yaxt,bty=bty,lwd=2,xaxs=xaxs,yaxs=yaxs) # richness
	}
	if (hatched==T){
		polygon(c(temp.cc[,4],rev(temp.cc[,4])),c(temp.cc[,1],rev(temp.cc[,3])),col=poly.col,border=F,density=50,angle=angle)
	} else {
		polygon(c(temp.cc[,4],rev(temp.cc[,4])),c(temp.cc[,1],rev(temp.cc[,3])),col=poly.col,border=F)
	}
	
	lines(temp.cc[,4],temp.cc[,2],col=line.col,lwd=2)
	
	abline(v=seq(0,30000,1000),lwd=0.1,col="gray")
}

# END FUNCTIONS #######################################
