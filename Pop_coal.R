
ped.gen<-function(N=1000,gens=10){
	these.dup<-numeric()
	all.ancs<-list()
	prev.anc<-sample(1:N,size=1,replace=FALSE)
	all.ancs[[1]]<-prev.anc
	for(k in 1:gens){
		all.k.ancs<-numeric()
		for(i in 1:2^(k-1)){
			
			if(i %in% these.dup){ 
				these.overlap<-which(all.ancs[[k]][i]==all.ancs[[k]]) ##find the inds that overlap the ith ancs
				first.anc<-min(these.overlap)   ##find first of these
				ancs<-c(all.k.ancs[2*first.anc-1],all.k.ancs[2*first.anc])  #copy their info for next gen
			}else{
				ancs<-sample(1:N,size=2,replace=FALSE)							
			}

			all.k.ancs<-c(all.k.ancs,ancs)
		}
		
		these.dup<-which(duplicated(all.k.ancs))
		print(length(these.dup))
		#if(any(duplicated(all.k.ancs))){recover(); all.k.ancs[duplicated(all.k.ancs)]<-NA}
		prev.anc<-all.k.ancs
		all.ancs[[k+1]]<-all.k.ancs

	}
	all.ancs
}


get.num.blocks<-function(family.chunks.gen){
	sapply(family.chunks.gen,function(ind){ sum(unlist(lapply(ind,function(x){if(is.null(x)){return(0)}; nrow(x)})))})
}

ped.plot<-function(all.ancs,family.chunks=NULL,my.col="black",plot.new=TRUE,adj.col=0.1,gens=10,N=1000){
	
	if(plot.new) plot(y=c(0,gens),x=c(0,1.25),type="n",axes=FALSE,ylab="generation",xlab="")
	axis(side=2,at=0:gens,tick=FALSE,las=1)

	if(!is.null(family.chunks)){
		tmp.ancs<-lapply(1:gens,function(k){
			tmp<-all.ancs[[k]]
			if(k>1){
				num.blocks<-get.num.blocks(family.chunks[[k-1]])
				tmp[num.blocks==0]<-NA
			}
			tmp
		})	
		all.ancs<-tmp.ancs
	}
	anc.points<-list()
	for(k in 1:(gens-1)){

		for(i in 1:2^(k-1)){
			
			 lines(c(all.ancs[[k]][i],all.ancs[[k+1]][2*i-1])/N,c(k-1,k),col=adjustcolor(my.col,adj.col))
			 lines(c(all.ancs[[k]][i],all.ancs[[k+1]][2*i])/N,c(k-1,k),col=adjustcolor(my.col,adj.col))
		}
		points(all.ancs[[k]]/N,rep(k-1,2^(k-1)),pch=19,col=adjustcolor(my.col,adj.col),cex=4*(1/2^k))		
		these.overlap<-which(!is.na(all.ancs[[k+1]]) & duplicated(all.ancs[[k+1]]))
		anc.points[[k]]<-all.ancs[[k+1]][these.overlap]/N
		
		all.ancs.table<-table(all.ancs[[k+1]][!is.na(all.ancs[[k+1]])])
		text(x=1.15,y=k,paste(2^(k),length(all.ancs.table),length(unique(all.ancs[[k+1]][these.overlap])),max(all.ancs.table),sep=", "),cex=.7)
	}
	
	lapply(1:(gens-1),function(k){anc.points.k<-anc.points[[k]];points(anc.points.k,rep(k,length(anc.points.k)))})
}


ped.plot.two<-function(all.ancs,all.ancs.2,my.col=rep("black",2),plot.new=TRUE,adj.col=0.1,gens=10,N=1000){
#	recover()
	if(plot.new) plot(y=c(0,gens),x=c(0,1.25),type="n",axes=FALSE,ylab="generation",xlab="")
	axis(side=2,at=0:gens,tick=FALSE,las=1)

	for(k in 1:(gens-1)){
		for(i in 1:2^(k-1)){
			lines(c(all.ancs[[k]][i],all.ancs[[k+1]][2*i-1])/N,c(k-1,k),col=adjustcolor(my.col[1],adj.col))
			lines(c(all.ancs[[k]][i],all.ancs[[k+1]][2*i])/N,c(k-1,k),col=adjustcolor(my.col[1],adj.col))
		}
		points(all.ancs[[k]]/N,rep(k-1,2^(k-1)),pch=19,col=adjustcolor(my.col[1],adj.col),cex=4*(1/2^k))
	}

	for(k in 1:(gens-1)){
		overlaps<-which((all.ancs.2[[k]] %in% all.ancs[[k]]))
		new.ancs<- all.ancs.2[[k]][overlaps]
		#if(length(overlaps)){print("here 1"); recover(); }  
		
		cat(length(new.ancs)," ")
		for(i in 1:2^(k-1)){
			if(i %in% overlaps){ 
				
				these.overlap<-which(all.ancs[[k]]==all.ancs.2[[k]][i]) ## find out which ancs in ped 1 match the ith 
				first.overlap<-min(these.overlap)  #take first of these overlapping ancs
				all.ancs.2[[k+1]][2*i]<-all.ancs[[k+1]][2*first.overlap]; all.ancs.2[[k+1]][2*i-1]<-all.ancs[[k+1]][2*first.overlap-1]
				}
			lines(c(all.ancs.2[[k]][i],all.ancs.2[[k+1]][2*i-1])/N,c(k-1,k),col=adjustcolor(my.col[2],adj.col))
			lines(c(all.ancs.2[[k]][i],all.ancs.2[[k+1]][2*i])/N,c(k-1,k),col=adjustcolor(my.col[2],adj.col))
		}
		
		cat(sum(all.ancs.2[[k+1]] %in% all.ancs[[k]]),"\n")	
		#if(length(new.ancs)) recover()
		#recover()
		text(x=1.15,y=k-1,paste(
		2^(k-1),mean(c(length(unique(all.ancs[[k]])),length(unique(all.ancs.2[[k]])))),
		length(unique(new.ancs)),sep=", "),cex=0.7)	
		points(all.ancs.2[[k]]/N,rep(k-1,2^(k-1)),pch=19,col=adjustcolor(my.col[2],adj.col),cex=4*(1/2^k))
		points(new.ancs/N,rep(k-1,length(new.ancs)))
	}
		
}


plot.lineage<-function(all.ancs,my.col="black",plot.new=FALSE,adj.col=1,gens=10,N=1000,sex){
#	recover()
	if(plot.new){ 
		plot(y=c(0,gens),x=c(0,1.25),type="n",axes=FALSE,ylab="generation")
		axis(side=2,at=0:gens,tick=FALSE,las=1)
	}
	for(k in 1:(gens-1)){
		if(sex=="female"){ 
			lines(c(all.ancs[[k]][1],all.ancs[[k+1]][1])/N,c(k-1,k),col=adjustcolor(my.col,adj.col))
			points(all.ancs[[k]][1]/N,k-1,pch=19,col=adjustcolor(my.col,adj.col),cex=.3)
		}
		if(sex=="male"){ 
			lines(c(all.ancs[[k]][2^(k-1)],all.ancs[[k+1]][2^k])/N,c(k-1,k),col=adjustcolor(my.col,adj.col))
			points(all.ancs[[k]][2^(k-1)]/N,k-1,pch=19,col=adjustcolor(my.col,adj.col),cex=.3)
		}
		
#		these.overlap<-which(duplicated(all.ancs[[k+1]]))
#		points(all.ancs[[k+1]][these.overlap]/N,rep(k,length(these.overlap)))
#		all.ancs.table<-table(all.ancs[[k+1]])
#		text(x=1.15,y=k,paste(2^(k),length(all.ancs.table),length(unique(all.ancs[[k+1]][these.overlap])),max(all.ancs.table),sep=", "))
	}
}


plot.mtDNA.Y<-function(all.ancs,my.col="black",plot.new=FALSE,adj.col=1,gens=10,N=1000,sex){
#	recover()
	if(plot.new){ 
		plot(y=c(0,gens),x=c(0,1.25),type="n",axes=FALSE,ylab="generation")
		axis(side=2,at=0:gens,tick=FALSE,las=1)
	}
	for(k in 1:(gens-1)){
		if(sex=="female"){ 
			lines(c(all.ancs[[k]][1],all.ancs[[k+1]][1])/N,c(k-1,k),col=adjustcolor(my.col,adj.col))
			points(all.ancs[[k]][1]/N,k-1,pch=19,col=adjustcolor(my.col,adj.col),cex=.3)
		}
		if(sex=="male"){ 
			lines(c(all.ancs[[k]][2^(k-1)],all.ancs[[k+1]][2^k])/N,c(k-1,k),col=adjustcolor(my.col,adj.col))
			points(all.ancs[[k]][2^(k-1)]/N,k-1,pch=19,col=adjustcolor(my.col,adj.col),cex=.3)
		}
		
	}
}

layout(t(1:2)); par(mar=c(.5,4,.5,.5))
ped.1<-ped.gen(N=20,gens=10)
ped.2<-ped.gen(N=20,gens=10)

 ped.plot(ped.1,my.col="blue",adj.col=.2,N=20,gens=10)
ped.plot.two(ped.1,ped.2,my.col=c("blue","red"),adj.col=.2,N=20,gens=10)

ped.1<-ped.gen(N=1e5,gens=16)
ped.2<-ped.gen(N=1e5,gens=16)


layout(t(1:2)); par(mar=c(.5,4,.5,.5))

ped.plot(ped.1,my.col="red",adj.col=.2,N=1e5,gens=16)
#ped.plot(ped.2,my.col="red",adj.col=.2,N=1e5,gens=12,plot.new=FALSE)

ped.plot.two(ped.2,ped.1,my.col=c("blue","red"),adj.col=.2,N=1e5,gens=16)

#show(load("~/Dropbox/Ideas/Genealogy_review/figs/family_tree_w_13_gens.Robj"))
#layout(t(1:2));
#ped.plot(ped.1,my.col="red",adj.col=.2,N=200e3,gens=13)
#ped.plot(ped.1,family.chunks=family.chunks,my.col="red",adj.col=.2,N=200e3,gens=13)
