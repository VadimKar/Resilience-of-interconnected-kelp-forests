#Note: we ran most analyses on a computer cluster (typically 100-200 cores per analysis, simulation times ~24h).


#Some global functions to calculate statistics

#Indentifying states
stateID=function(SpatiotempSeries){
	#SpatiotempSeries is a npatch x Len x 3*nruns array of densities; output is a npatch x Len x nruns T/F array of whether community is forested
	if(length(dim(SpatiotempSeries))==3) statePairs=cbind(as.vector(SpatiotempSeries[,,c(F,T,F)]),as.vector(SpatiotempSeries[,,c(T,F,F)]))
	if(length(dim(SpatiotempSeries))==2) statePairs=cbind(as.vector(SpatiotempSeries[,c(F,T,F)]),as.vector(SpatiotempSeries[,c(T,F,F)]))

	#This is a shortcut to reduce computation time: categorize communities we kno must be forested/barren based on kelp densites
	states=rep(NA,nrow(statePairs)); states[statePairs[,2]>1e3]=TRUE; states[statePairs[,2]<4e2]=FALSE; statesUnclear=unique(c(1,which(is.na(states))));

	#for remaining states, classify state based on probability from fitted distributions
	#Active code gives distribution fits used in analysis
	pBar=dmst(statePairs[statesUnclear,], mu=c(6500,100), sigma=matrix(c(1752523,-100540,-100540,26262),2), delta=rbind(-283,37), dof=5.5, known=NULL, tmethod=1)
	pFor=dmst(statePairs[statesUnclear,], mu=c(3800,1400), sigma=matrix(c(1825358,-730412,-730412,309244),2), delta=rbind(64,44), dof=9.9, known=NULL, tmethod=1)
	#pBar=dmst(statePairs[statesUnclear,], mu=fit$mu[[1]], sigma=fit$sigma[[1]], delta=fit$delta[[1]], dof=fit$dof[1], known=NULL, tmethod=1)
	#pFor=dmst(statePairs[statesUnclear,], mu=fit$mu[[2]], sigma=fit$sigma[[2]], delta=fit$delta[[2]], dof=fit$dof[2], known=NULL, tmethod=1)
	states[statesUnclear]=pBar<pFor
	if(length(dim(SpatiotempSeries))==3) return(array(states,dim=dim(SpatiotempSeries)*c(1,1,1/3)))
	if(length(dim(SpatiotempSeries))==2) return(matrix(states,dim(SpatiotempSeries)*c(1,1/3)))
}


#Summary statistics functions
Pcts=c(0.25,0.75)
perCentiles=function(x,Pcts=c(0.25,0.05,0.75),avg=FALSE){
	if(avg){  out=c(mean(x,na.rm=T),quantile(x,Pcts));  names(out)=c("mean",as.character(Pcts));  return(out);  }
	if(!avg){  out=quantile(x,Pcts,type=5);  names(out)=as.character(Pcts);  return(out);  }
}
stays3=function(x,agg=TRUE,avg=FALSE){ temp=rle(x); last=length(temp$values);
	ForestStays=(temp$lengths[-last][temp$values[-last]==TRUE]); if(is.na(mean(ForestStays))) ForestStays=length(x)*mean(x);
	BarrenStays=(temp$lengths[-last][temp$values[-last]==FALSE]); if(is.na(mean(BarrenStays))) BarrenStays=length(x)*(1-mean(x));
	SwitchFreq=1/mean(temp$lengths[-last]); if(is.na(SwitchFreq)) SwitchFreq=0;
	if(agg) return(c(perCentiles(ForestStays,Pcts,avg), perCentiles(BarrenStays,Pcts,avg), SwitchFreq))
	if(!agg){ m=length(ForestStays); n=length(BarrenStays); out=matrix(NA,max(m,n),2); out[1:m,1]=ForestStays; out[1:n,2]=BarrenStays; return(out); }
}
#Analyze state distributions, durations, and extents
# Trun sets simulation length 
# Len sets period analyzed at the end
SpatiotempSumm=function(temp,Len=1e3,Trun=4e3){
	#Identify all forested states in simulations
	ForestedFin=array(dim=c(dim(temp)[1],Len,0)); ForestedInit=ForestedFin[,50:125,];
	Trials=matrix(1:dim(temp)[3],ncol=3,byrow=TRUE); for(i in 1:nrow(Trials)){
		ForestedInit=abind(ForestedInit, stateID(temp[,50:125,Trials[i,]]), along=3)
		ForestedFin=abind(ForestedFin, stateID(temp[,(Trun-Len+1):Trun,Trials[i,]]), along=3)
	}
	Atrans=apply(temp[,50:125,c(TRUE,FALSE,FALSE)],3,mean); Ftrans=apply(ForestedInit,3,mean);
	Afin=apply(temp[,(Trun-Len):Trun,c(TRUE,FALSE,FALSE)],3,mean); Ffin=apply(ForestedFin,3,mean); Pfin=apply(temp[,(Trun-Len):Trun,c(FALSE,FALSE,TRUE)],3,mean);
	stayLs=t(apply(ForestedFin, 3, function(x) rowMeans(apply(x,1,stays3,avg=TRUE))))
	extents=t(apply(ForestedFin, 3, function(x) rowMeans(apply(x,2,stays3,avg=TRUE))))
	colnames(stayLs)=c("Fd0.5","Fd0.25","Fd0.75","Bd0.5","Bd0.25","Bd0.75","SwitchFreq")
	colnames(extents)=c("Fed0.5","Fe0.25","Fe0.75","Be0.5","Be0.25","Be0.75","SwitchFreq")
	return(round(cbind(Atrans,Ftrans,Afin,Ffin,Pfin,stayLs,extents[,1:6]),2))
}











###PRELIMINARY: Determine community state classification - reproduces stateID in OpenForests_CoreFunctions.R:
StartBarren=StochSim4(3e3,start=c(3e2,4e3,6),parms,stoch=c(3,3,0),giveout=2,npatch=Npatch,startvar=0,SimType=3)
StartForest=StochSim4(3e3,start=c(2e3,2e3,50),parms,stoch=c(3,3,0),giveout=2,npatch=Npatch,startvar=0,SimType=3)
parms["fc"]=0; StartForestNofeedback=StochSim4(3e3,start=c(3e2,4e3,6),parms,stoch=c(3,3,0),giveout=2,npatch=Npatch,startvar=0,SimType=3); parms["fc"]=0.85;

#Fitting bivariate skewed-normal distribution to categorize states in stateID. To accelerate fit, using vissually approximated, fixed means
dat=cbind(as.vector(StartForest[,PlotTimes,2]),as.vector(StartForest[,PlotTimes,1]))
fitTwo=fmmst(2, dat, print=TRUE, eps=1e-2, itmax=10, known=list(mu=list(c(6500,1e2),c(3800,1400))))
#Visualize fitted distributions and attraction basins:
fmmst.contour.2d(dat, model=fitTwo, map="cluster", drawpoints=FALSE, clusters=fitTwo$clusters, cex=0.)
points(c(6500,3800),c(1e2,1400),cex=3,col=6,pch=16)

###FIGURE 4
PlotTimes=2700:2800
SpatiotempPlot(StartForest[,PlotTimes,1])
SpatiotempPlot(StartForestNofeedback[,PlotTimes,1])

hexbinplot(as.vector(StartForest[,PlotTimes,1])/Space~as.vector(StartForest[,PlotTimes,2])/Space,xlim=c(0,14),ylim=c(-0.1,4.1),xbins=30,aspect=1,ylab="Kelp biomass, kg m-2",xlab="Urchin density, m-2")
hexbinplot(as.vector(StartForestNofeedback[,PlotTimes,1])/Space~as.vector(StartForestNofeedback[,PlotTimes,2])/Space,xlim=c(0,14),ylim=c(-0.1,4.1),xbins=10,aspect=1,ylab="Kelp biomass, kg m-2",xlab="Urchin density, m-2")


###FIGURE B.1
#matplot(cbind(0.+colMeans(StartForest[,,1])[c(T,rep(F,10))]/Space,colMeans(StartBarren[,,1])[c(T,rep(F,10))]/Space),type="l",lty=1,pch=1,col=c(1,4))
matplot(seq(1,4e3,length=3e3),cbind(colMeans(StartForest[,,1])/Space,colMeans(StartBarren[,,1])/Space),type="l",lty=1,pch=1,col=c(8,1,2),lwd=c(9,1,1),ylab="Kelp density, kg m-2",xlab="Years")
legend("topright",legend=c("Initially forested","Initially near-barren"),lwd=c(8,1),col=c(8,1),box.col=0,cex=1.25)




###FIGURE 3a
#Simulate community dynamics with feedback and no stochasticity across fishing intensities:
ntests=40; Npatch=2e1; Trun=5e2; parms["psiA"]=0; parms["psiU"]=parms["psiP"]=80; RegionSize=1100;
FPseq=round(seq(0.15,0.3,length=ntests),3); StartsMultipatchDet=as.matrix(expand.grid(FP=FPseq,psi=0,StartState=c(0,1),Pstoch=0));
storeMultipatchDet=array(as.numeric(),dim=c(Npatch,1,0)); for(i in 1:nrow(StartsMultipatchDet)){
	parms["FP"]=StartsMultipatchDet[i,1]; Start=as.numeric(rbind(c(3e2,4e3,4),c(2e3,2e3,40))[1+StartsMultipatchDet[i,3],]);
	storeMultipatchDet=abind(storeMultipatchDet, StochSim4(Trun,start=Start,parms,stoch=c(0,0,0, 0,0,0),giveout=2,npatch=Npatch,startvar=0,SimType=3)[,Trun,], along=3)
}

#Apply state ID function from above and plot results across fishing intensities:
MultipatchDet=apply(stateID(storeMultipatchDet), 3, mean) #Spatial averages of forest frequencies
matplot(FPseq,cbind(MultipatchDet[StartsMultipatchDet[,"StartState"]==0,],MultipatchDet[StartsMultipatchDet[,"StartState"]==1,]),lwd=2:1,col=c(8,1),pch=c(16,1),lty=1:2)



###FIGURE 2:
##a: 1-patch kelp forest resilience under openness:
#Estimating recruit production: in each state
parms["FP"]=0.23; parms[1:3]=0; npatch=8; Trun=2e2;	
RhosForest=round(StochSim4(Trun,start=c(2e3,4e3,40),parms,stoch=c(0,0,0,0,0,0),giveout="Rhos",npatch=npatch,startvar=0,SimType=1)[1,2e2,4:6],2)
RhosBarren=round(StochSim4(Trun,start=c(3e2,4e3,4),parms,stoch=c(0,0,0,0,0,0),giveout="Rhos",npatch=npatch,startvar=0,SimType=1)[1,2e2,4:6],2)
parms[1:3]=0.75*RhosForest+0.25*RhosBarren #Results are RhosForest=c(0,2316.8,24.4) and RhosBarren=c(0,710,0.04)

#Simulate dynamics across initial conditions and openness levels
Trun=50; GamSeq=seq(0,1,length=1e2); GamsExpand=expand.grid(GamSeq,GamSeq); gamUs=GamsExpand[,1]; gamPs=GamsExpand[,2];
StrtSeq=seq(0.01,1,length=10); StartsDetOpen=as.matrix(expand.grid(StrtSeq,StrtSeq,StrtSeq))%*%diag(c(2998,12e3,80));
StoreDetOpen=array(as.numeric(),dim=c(nrow(StartsDetOpen),8,0)); 
for(i in 1:nrow(StartsDetOpen)){ Start=as.numeric(StartsDetOpen[i,]);
	StoreDetOpen=abind(StoreDetOpen, cbind(GamsExpand,t(Start),StochSim4(Trun,start=Start,parms,stoch=c(0,0,0,0,0,0),giveout=2,npatch=length(gamUs),startvar=0,SimType=1)[,Trun,]), along=3); }

#Plot results:
reShape=function(x){ x2=aperm(x, c(1,3,2)); dim(x2)=c(dim(x)[1]*dim(x)[3],dim(x)[2]); return(x2); };
DetOpen=reShape(StoreDetOpen)
DetOpenSumm=aggregate(stateID(DetOpen[,6:8])~DetOpen[,1]+DetOpen[,2], FUN=mean)
SpatiotempPlot(matrix(DetOpenSumm[,3],ncol=length(GamSeq),byrow=FALSE),Zlim=c(min(DetOpenSumm[,3]),1))





##b-c and Figure A.3: 1-patch resilience under openness and environmental stochasticity:
#Estimate external recruit production (rho_i) in the stochastic, 1-patch case:
#Generate stochastic simulation 
parms["FP"]=0.192; npatch=1e2; Trans=50:75; Trun=2e2; gamAs=gamUs=gamPs=rep(0,npatch);
OnepatchRecrDyns=StochSim4(Trun,start=c(2e3,2e3,50),parms,stoch=c(3,3,0,0,0,0),giveout="Rhos",npatch=npatch,startvar=0,SimType=1)

#Identify community states based on dynamics during transient period (Trans; Figure A.2):
Forests=stateID(OnepatchRecrDyns[,,1:3])[,,1];
plot(colMeans(OnepatchRecrDyns[,1:Trun,3])/Space,type="l",lwd=2,ylab="Predator density, m^-2",xlab="Time"); abline(v=c(min(Trans),max(Trans)),col=4,lty=2,lwd=3);
par(new=TRUE); plot(colMeans(Forests),type="l",lwd=2,col=6,ylim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab="")
axis(4); mtext("Forest frequency",side=4,line=3)

#Barren recruit productions:
Forested=Forests; Forested[Forested==TRUE]=NA; Forested=1+Forested
RecrProdBarren=abind(OnepatchRecrDyns[,,4]*Forested,OnepatchRecrDyns[,,5]*Forested,OnepatchRecrDyns[,,6]*Forested, along=3)[,Trans,]
BarrenRhos=apply(RecrProdBarren, 3, mean, na.rm=T)
#Forest recruit productions:
Forested=Forests; Forested[Forested==FALSE]=NA;
RecrProdForested=abind(OnepatchRecrDyns[,,4]*Forested,OnepatchRecrDyns[,,5]*Forested,OnepatchRecrDyns[,,6]*Forested, along=3)[,Trans,]
ForestRhos=apply(RecrProdForested, 3, mean, na.rm=T)


#Simulate dynamics across initial conditions and openness levels
ntests=20; nReps=36; Npatch=ntests^2; Trun=4e3; #NOTE THAT Npatch MUST BE EVEN!! 
GamSeq=seq(0.02,1,length=ntests); GamsExpand=expand.grid(GamSeq,GamSeq); gamUs=double(ntests); gamUs=GamsExpand[,1]; gamPs=GamsExpand[,2];
Starts=expand.grid(Stoch=c(1),ExtForest=c(0.775,1),StartForest=c(0,1),Rep=1:nReps); dim(Starts);

ForestRhos=c(0,1743,11.2); BarrenRhos=c(0,201,12.3); store1pStoch=array(dim=c(Npatch,Trun,0));
for(i in 1:nrow(Starts)){
	parms[1:3]=as.numeric(Starts[i,2]*ForestRhos + (1-Starts[i,2])*BarrenRhos)
	if(Starts[i,3]==1) store1pStoch=abind(store1pStoch, StochSim4(Trun,start=c(2e3,2e3,50),parms,stoch=c(3,3,0),giveout=2,npatch=Npatch,startvar=0,SimType=1), along=3)
	if(Starts[i,3]==0) store1pStoch=abind(store1pStoch, StochSim4(Trun,start=c(3e2,4e3,4),parms,stoch=c(3,3,0),giveout=2,npatch=Npatch,startvar=0,SimType=1), along=3)
}

analys=function(X,stoch=TRUE){
	if(!stoch) return(c(mean(X),mean(X>1e3),rep(NA,4)))
	ForestedFin=as.vector(stateID(X))
	return(c(mean(X[,3]), mean(ForestedFin), stays3(ForestedFin,avg=TRUE)[1:6]))
}
SummStats=t(apply(store1pStoch[,2e3:4e3,], 1, analys))
TrialsAlls=data.frame(cbind(round(GamsExpand,2), ExtForest=rep(Starts[,2],each=ntests^2), SummStats))
names(TrialsAlls)=c("GamU","GamP","ExtForest","meanA","meanF","dF","d25F","d75F","dB","d25B","d75B")
SummStatsAgg=aggregate(cbind(meanA,meanF,dF,d25F,d75F,dB,d25B,d75B)~GamU+GamP+ExtForest, FUN=mean, data=TrialsAlls); head(SummStatsAgg);

par(mfrow=c(1,3)); cont=0; cexAx=1.25; contSpan=2; cexCont=1.5;
#Outside 75% forested - Figure 2 b, c, d
SummStatsTemp=SummStatsAgg[which(SummStatsAgg$ExtForest==unique(Starts[,2])[1]),-(1:3)];
x=matrix(SummStatsTemp[,2],ntests,ntests); SpatiotempPlot(x,cont=cont,Conts=c(0.4,0.6,0.8),contSpan=contSpan,cexCont=cexCont,Zlim=range(x),figtitle="Forest resilience, outside forested",cexAx=cexAx, Xax=GamSeq, Yax=GamSeq, XaxN="Predator Openness",YaxN="Urchin Openness")
x=matrix(SummStatsTemp[,3],ntests,ntests); SpatiotempPlot(x,cont=cont,Conts=c(6,8,10,12),contSpan=contSpan,cexCont=cexCont,Zlim=range(x),figtitle="Forest duration, outside forested",cexAx=cexAx, Xax=GamSeq, Yax=GamSeq, XaxN="Predator Openness",YaxN="Urchin Openness")
x=matrix(SummStatsTemp[,6],ntests,ntests); SpatiotempPlot(log(1+x),cont=cont,Conts=c(2,6,4,8,14),contPlot=x,contSpan=contSpan,cexCont=cexCont,Zlim=range(log(1+x)),figtitle="Barren duration, outside forested",cexAx=cexAx, Xax=GamSeq, Yax=GamSeq, XaxN="Predator Openness",YaxN="Urchin Openness")

#Outside fully forested - Figure A.3
SummStatsTemp=SummStatsAgg[which(SummStatsAgg$ExtForest==unique(Starts[,2])[1]),-(1:3)];
x=matrix(SummStatsTemp[,2],ntests,ntests); SpatiotempPlot(x,cont=cont,Conts=c(0.4,0.6,0.8),contSpan=contSpan,cexCont=cexCont,Zlim=range(x),figtitle="Forest resilience, outside forested",cexAx=cexAx, Xax=GamSeq, Yax=GamSeq, XaxN="Predator Openness",YaxN="Urchin Openness")
x=matrix(SummStatsTemp[,3],ntests,ntests); SpatiotempPlot(x,cont=cont,Conts=c(4,6,8,10),contSpan=contSpan,cexCont=cexCont,Zlim=range(x),figtitle="Forest duration, outside forested",cexAx=cexAx, Xax=GamSeq, Yax=GamSeq, XaxN="Predator Openness",YaxN="Urchin Openness")
x=matrix(SummStatsTemp[,6],ntests,ntests); SpatiotempPlot(x,cont=cont,Conts=c(2,4,10,16,22),contSpan=contSpan,cexCont=cexCont,Zlim=range(x),figtitle="Barren duration, outside forested",cexAx=cexAx, Xax=GamSeq, Yax=GamSeq, XaxN="Predator Openness",YaxN="Urchin Openness")











###Stochastic spatial model with feedbacks - FIGURES 3b, 5:
Trun=4e3; Npatch=2e2;
#Varying fishing intensity:
nRepsFP=10; ntests=50; FPseq=round(seq(0.1,0.25,length=ntests),3);
StartsFP=as.matrix(expand.grid(FP=FPseq,psi=0,StartState=c(0,1),Pstoch=0));
#Varying scale of stochasticity in urchin recruitment, psiSigmaU:
nRepsPsi=25; psiSeq=seq(1,100,length=25)
#To save computation time, zero in on FP levels near shift to barren dominance when varying psiSigmaU:
FPseqVarpsi=FPseq[c(26,28:31,33)] 
StartsVarpsi=as.matrix(expand.grid(FP=FPseqVarpsi,psi=psiSeq,StartState=c(1),Pstoch=0))
#Replicate and assemble trials together:
StartsMp=rbind(do.call("rbind",rep(list(StartsFP),nRepsFP)), do.call("rbind",rep(list(StartsVarpsi),nRepsPsi)))

storeMpStoch=array(dim=c(Npatch,Trun,0))
for(i in 1:nrow(StartsMp)){
	parms["FP"]=StartsMp[i,1]; if(StartsMp[i,2]==0) parmsStoch["sigmaUpsi"]=12; if(StartsMp[i,2]>0) parmsStoch["sigmaUpsi"]=StartsMp[i,2];
	if(StartsMp[i,3]==1) store=abind(store, StochSim4(Trun,start=c(2e3,2e3,40),parms,stoch=c(3,3,0),giveout=2,npatch=Npatch,startvar=0,SimType=3), along=3)
	if(StartsMp[i,3]==0) store=abind(store, StochSim4(Trun,start=c(3e2,4e3,6),parms,stoch=c(3,3,0),giveout=2,npatch=Npatch,startvar=0,SimType=3), along=3)
}
parms["sigmaUpsi"]=12; parms["FP"]=0.192; #reset default parameters for future simulations
StochRunsSumm=apply(storeMpStoch,3,SpatiotempSumm,Trun=Trun,Len=1e3) #Summarizing results
Subselect2=function(Var,IC,psi,Smooth=NA,logT=FALSE,FUN=1){
	StochRunsSumm2=cbind(StochRunsSumm,StartsMp)
	temp=cbind(StartsMp,StochRunsSumm2)[which(StartsMp[,"StartState"]==IC & StartsMp[,"psi"]==psi), c(1,Var+4)]
	if(FUN==1)	temp1=as.matrix(aggregate(temp, by=list(temp[,1]), FUN=mean, na.rm=T)[,-1])
	if(FUN==2)	temp1=as.matrix(aggregate(temp, by=list(temp[,1]), FUN=quantile, 0.75)[,-1])
	if(FUN==3)	temp1=as.matrix(aggregate(temp, by=list(temp[,1]), FUN=quantile, 0.25)[,-1])
	if(logT) temp1[,2]=log(temp1[,2])
	if(is.na(Smooth)) return(temp1)
	if(!is.na(Smooth)){ sm=Smoo(temp1,bass=Smooth); return(cbind(sm$x,sm$y)); }
}


#Plotting parameters
gam=0; FPinclude=10:50; FPxax=seq(.15,.25,by=0.05); Pstoch=0; cexlab=1.5; DurYmax=c(7,9); ExtYmax=c(5.5,80);

#Figure 3b - Forest frequency across fishing levels
Var=4; Smooth=NA; plot(Subselect2(Var,1,gam,Pstoch,Smooth)[FPinclude,],type="b",pch=16,col=1,lty=2,lwd=2,ylim=c(0.3,1),xlab="Fishing intensity",ylab="Proportion forested",cex.lab=cexlab,xaxt="n"); axis(1,at=FPxax);
lines(Subselect2(Var,0,gam,Pstoch,Smooth)[FPinclude,], type="b",pch=16,col=1,lty=1,lwd=2)
polygon(rbind(Subselect2(Var,1,gam,Pstoch,Smooth,FUN=3)[FPinclude,], apply(Subselect2(Var,1,gam,Pstoch,Smooth,FUN=2)[FPinclude,],2,rev)), border=0, col=rgb(red=0,green=0,blue=0,alpha=0.1))
polygon(rbind(Subselect2(Var,0,gam,Pstoch,Smooth,FUN=3)[FPinclude,], apply(Subselect2(Var,0,gam,Pstoch,Smooth,FUN=2)[FPinclude,],2,rev)), border=0, col=rgb(red=0,green=0,blue=0,alpha=0.1))

#Figure 5a
#Forest state duration:
Var=6; Smooth=NA; logT=TRUE; plot((Subselect2(Var,1,gam,Pstoch,Smooth,logT)[FPinclude,]),yaxt="n", type="l",col=1,lwd=2,ylim=c(0,DurYmax[1]),xlab="Fishing intensity",ylab="Forest duration",cex.lab=cexlab,xaxt="n"); axis(1,at=FPxax);
polygon(rbind(Subselect2(Var+2,1,gam,Pstoch,Smooth,logT)[FPinclude,], apply(Subselect2(Var+1,1,gam,Pstoch,Smooth,logT)[FPinclude,],2,rev)), border=0, col=rgb(red=0,green=0,blue=0,alpha=0.1))
axis(2,at=(0:round(DurYmax[1]))[c(F,T)],labels=signif(exp(0:DurYmax[1]),1)[c(F,T)])
#Barren state duration:
Var=9; par(new=TRUE); logT=FALSE; plot((Subselect2(Var,1,gam,Pstoch,Smooth,logT)[FPinclude,]),yaxt="n", type="l",col=1,lty=2,lwd=2,ylim=c(-1,DurYmax[2]),xaxt="n",yaxt="n",xlab="",ylab="");
axis(4); mtext("Barren duration",side=4,line=3,cex=cexlab)
polygon(rbind(Subselect2(Var+2,1,gam,Pstoch,Smooth,logT)[FPinclude,], apply(Subselect2(Var+1,1,gam,Pstoch,Smooth,logT)[FPinclude,],2,rev)), border=0, col=rgb(red=0,green=0,blue=0,alpha=0.1))

#Figure 5b
#Forest state extent:
Var=13; Smooth=NA; logT=TRUE; plot((Subselect2(Var,1,gam,Pstoch,Smooth,logT)[FPinclude,]),yaxt="n", type="l",col=1,lwd=2,ylim=c(0,ExtYmax[1]),xlab="Fishing intensity",ylab="Forest extent",cex.lab=cexlab,xaxt="n"); axis(1,at=FPxax);
axis(2,at=0:round(ExtYmax[1]),labels=signif(5.5*exp(0:round(ExtYmax[1])),1))
polygon(rbind(Subselect2(Var+2,1,gam,Pstoch,Smooth,logT)[FPinclude,], apply(Subselect2(Var+1,1,gam,Pstoch,Smooth,logT)[FPinclude,],2,rev)), border=0, col=rgb(red=0,green=0,blue=0,alpha=0.1))
#Barren state extent:
Var=16; par(new=TRUE); logT=FALSE; plot((Subselect2(Var,1,gam,Pstoch,Smooth,logT)[FPinclude,]%*%diag(c(1,5.5))),type="l",col=1,lty=2,lwd=2,ylim=c(0,ExtYmax[2]),xaxt="n",yaxt="n",xlab="",ylab="");
axis(4); mtext("Barren extent",side=4,line=3,cex=cexlab)
polygon(rbind(Subselect2(Var+2,1,gam,Pstoch,Smooth,logT)[FPinclude,]%*%diag(c(1,5.5)), apply(Subselect2(Var+1,1,gam,Pstoch,Smooth,logT)[FPinclude,]%*%diag(c(1,5.5)),2,rev)), border=0, col=rgb(red=0,green=0,blue=0,alpha=0.1))


#Disturbance scale response plots:
trim=function(x,bnds){ x[x<bnds[1] | x>bnds[2]]=NA; return(x); } 
Toplot=which(StartsMp[,"FP"]==FPseq[30] & StartsMp[,"psi"]>0 & StartsMp[,"StartState"]==1);
Bass=8; Xlim=c(0,80); ptcol=0; Lwd=3; Lty=c(1,2);

#Figure 5c
#Forest state duration:
Var=13; x=StartsMp[Toplot,"psi"]; y=trim(5.5*StochRunsSumm[Toplot,Var],c(25,180)); plot(x,y,type="p",xlim=Xlim,col=ptcol); 
lines(supsmu(x,y,bass=Bass),lwd=Lwd,lty=Lty[1]); par(new=TRUE);
#Barren state duration:
Var=16; x=StartsMp[Toplot,"psi"]; y=trim(5.5*StochRunsSumm[Toplot,Var],c(5,25)); plot(x,y,type="p",xlim=Xlim,col=ptcol,xaxt="n",yaxt="n");
lines(supsmu(x,y,bass=Bass),lwd=Lwd,lty=Lty[2]); axis(4);

#Figure 5d
#Forest state extent:
Var=6; x=StartsMp[Toplot,"psi"]; y=trim(StochRunsSumm[Toplot,Var],c(5,20)); plot(x,y,type="p",ylim=c(10,18),xlim=Xlim,col=ptcol); 
lines(supsmu(x,y,bass=Bass),lwd=Lwd,lty=Lty[1]); par(new=TRUE);
#Barren state extent:
Var=9; x=StartsMp[Toplot,"psi"]; y=trim(StochRunsSumm[Toplot,Var],c(0,5)); plot(x,y,type="p",xlim=Xlim,col=ptcol,xaxt="n",yaxt="n"); 
lines(supsmu(x,y,bass=Bass),lwd=Lwd,lty=Lty[2]); axis(4);





###Spatial model without feedback or without stochasticity - FIGURES 3a, 4b,d:
nReps=10; ntests=80; Npatch=2e2; Trun=2e2; parms["psiA"]=0; parms["psiU"]=parms["psiP"]=80; RegionSize=1100;
FPseq=round(seq(0.15,0.5,length=ntests),3); FPseqVarpsi=FPseq[c(26,28:31,33)]; psiSeq=seq(1,100,length=25);
Starts=as.matrix(expand.grid(FP=FPseq,psi=0,StartState=c(0,1),Pstoch=0)); dim(Starts); head(Starts);
for(i in Torun){
	parms["FP"]=Starts[i,1]
	if(Starts[i,3]==1) store=abind(store, StochSim4(Trun,start=c(2e3,2e3,40),parms,stoch=c(0,0,0, 0,0,0),giveout=2,npatch=Npatch,startvar=0,SimType=3), along=3)
	if(Starts[i,3]==0) store=abind(store, StochSim4(Trun,start=c(3e2,4e3,4),parms,stoch=c(0,0,0, 0,0,0),giveout=2,npatch=Npatch,startvar=0,SimType=3), along=3)
}
#Figure 3a:
DetForestIC=Starts[,3]==1
plot(Starts[!DetForestIC,1],(colMeans(Outs[,20,c(T,F,F)])>620)[!DetForestIC],type="o",pch=16,col=8,lty=1,lwd=6,cex=2,xlim=c(0.16,0.275),cex.lab=cexlab,cex.axis=cexlab)
lines(Starts[DetForestIC,1],(colMeans(Outs[,20,c(T,F,F)])>620)[DetForestIC], type="o",pch=1,col=1,lty=2,lwd=2)








###FIGURE E.1 - Sensitivity to recruitment facilitation:
nReps=6; ntests=50; Npatch=2e2; Trun=4e3; FPseq=round(seq(0.1,0.28,length=ntests),3);
Starts0=as.matrix(expand.grid(FP=FPseq,psi=0,StartState=c(0,1),Pstoch=0,fc=c(0.75,0.8,0.9))); StartsNpSens=do.call("rbind", rep(list(Starts0), nReps));

storeNpSens=array(dim=c(Npatch,Trun,0));
for(i in 1:nrow(StartsNpSens)){
	parms["FP"]=StartsNpSens[i,1]; parms["fc"]=StartsNpSens[i,5];
	if(StartsNpSens[i,3]==1) storeNpSens=abind(storeNpSens, StochSim4(Trun,start=c(2e3,2e3,40),parms,stoch=c(3,3,0),giveout=2,npatch=Npatch,startvar=0,SimType=3), along=3)
	if(StartsNpSens[i,3]==0) storeNpSens=abind(storeNpSens, StochSim4(Trun,start=c(3e2,4e3,6),parms,stoch=c(3,3,0),giveout=2,npatch=Npatch,startvar=0,SimType=3), along=3)
}
StochRunsSumm=apply(storeNpSens,3,SpatiotempSumm,Trun=Trun,Len=1e3) #Summarizing results

#Fishing bifurcation diagram at different facilitation levels:
Subselect3=function(Var,IC,psi,fc,Smooth=NA,logT=FALSE,FUN=1){
	StochRunsSumm2=cbind(StochRunsSumm,StartsNpSens)
	temp=cbind(StartsNpSens,StochRunsSumm2)[which(StartsNpSens[,"StartState"]==IC & StartsNpSens[,"psi"]==psi & StartsNpSens[,"fc"]==fc), c(1,Var+4)]
	if(FUN==1)	temp1=as.matrix(aggregate(temp, by=list(temp[,1]), FUN=mean, na.rm=T)[,-1])
	if(FUN==2)	temp1=as.matrix(aggregate(temp, by=list(temp[,1]), FUN=quantile, 0.75)[,-1])
	if(FUN==3)	temp1=as.matrix(aggregate(temp, by=list(temp[,1]), FUN=quantile, 0.25)[,-1])
	if(logT) temp1[,2]=log(temp1[,2])
	if(is.na(Smooth)) return(temp1)
	if(!is.na(Smooth)){ sm=Smoo(temp1,bass=Smooth); return(cbind(sm$x,sm$y)); }
}
Var=5; psi=0;
fc=unique(StartsNpSens[,"fc"])[2]; plot(rbind(Subselect3(Var,0,psi,fc),c(0.282,0.29)),type="o",pch=16,col=8,lty=1,lwd=6,cex=2,main=fc,ylim=c(0.29,1),xlab="Fishing intensity",ylab="Proportion forested");
lines(rbind(Subselect3(Var,1,psi,fc),c(0.282,0.29)),type="o",pch=1,col=1,lty=2,lwd=2); abline(v=0.28,col=4,lwd=2,lty=2);
fc=unique(StartsNpSens[,"fc"])[3]; plot(Subselect3(Var,0,psi,fc),type="o",pch=16,col=8,lty=1,lwd=6,cex=2,main=fc,ylim=c(0.29,1),xlab="Fishing intensity",ylab="Proportion forested");
lines(Subselect3(Var,1,psi,fc),type="o",pch=1,col=1,lty=2,lwd=2); abline(v=0.232,col=4,lwd=2,lty=2);
fc=unique(StartsNpSens[,"fc"])[5]; plot(Subselect3(Var,1,psi,fc),type="o",pch=16,col=8,lty=1,lwd=6,cex=2,main=fc,ylim=c(0.29,1),xlab="Fishing intensity",ylab="Proportion forested");
lines(Subselect3(Var,0,psi,fc),type="o",pch=1,col=1,lty=2,lwd=2); abline(v=0.140,col=4,lwd=2,lty=2);

#State distributions (urchin densities truncated by Umax for clarity)
StatePlot=function(x){ 
	series=storeNpSens[,(Trun-1e3):Trun,x]
	Include=as.vector(series[,,2])<1e4
	plot(hexbin(as.vector(series[,,2])[Include]/750,as.vector(series[,,1])[Include]/750,xbnds=c(0,Umax/750),ybnds=c(0,3e3/750)),xlab="",ylab="")
}
StatePlot(which(StartsNpSens[,1]==0.28 & StartsNpSens[,5]==unique(Starts2[,"fc"])[2])[1])
StatePlot(which(StartsNpSens[,1]==0.232 & StartsNpSens[,5]==unique(Starts2[,"fc"])[3])[1])
StatePlot(which(StartsNpSens[,1]==0.144 & StartsNpSens[,5]==unique(Starts2[,"fc"])[5])[1])




###FIGURE D.1 - Drivers of state transitions:
Transits=StochSim4(4e3,start=c(2e3,2e3,40),parms,stoch=c(3,3,0),giveout=3,npatch=2e2,startvar=0,SimType=3)[,3e3:4e3,]
Transits[,,c(F,F,F,T,T,T)]=Transits[,,c(F,F,F,T,T,T)]/1e4; apply(Transits,3,mean);

#1 - map out state shifts. Positives and Negs indicate forest and barren formation, respectively
#Note that Last year's environment/pulse relative to what see before pulse at end of this year
states=stateID(Transits[,,1:3])[,,1]; #Remove 1st year bcs no reference environment
N=5; states2=t(apply(states, 1, function(x) filter(x,rep(1/N,N))>0.5)); 
switches=t(apply(states2,1,diff))
storms=Transits[,-dim(Transits)[2],4]
pulses=Transits[,-dim(Transits)[2],5]

Cor=function(x,y) cor(as.vector(x),as.vector(y),use="pairwise.complete.obs")
Cor(switches,pulses); Cor(switches,storms); 

par(mfrow=c(1,2))
plot(as.factor(as.vector(switches)),as.vector(storms),col=8,horizontal=T,xlab="Kelp survival")
plot(as.factor(as.vector(switches)),as.vector(pulses),col=8,horizontal=T,xlab="Urchin recruitment success")









