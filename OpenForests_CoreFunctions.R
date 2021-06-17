#Note: we ran all analyses using R 3.4.3. Package EMMIXuskew appears unavalible for R 3.5 and later.
require(deSolve); require(MASS); require(abind); require(EMMIXuskew); require(sn); 
require(fields); require(hexbin); require(pals);

#Default parameter set. Note that to avoid numerical errors, we estimate all dynamics for a 750m^2 habitat area
#Note that added predator resources beta is referred to as rN here.
Space=750;
parms=c(RhoA=0,RhoU=0,RhoP=0, rA=12,K=4*Space,deltaA=1.09/Space, a=0.41, b=0.0144, deltaU=11.42/Space, deltaN=0.1, muU=0.1,muP=0.05,FP=0.19,rN=(6.67e-5)*Space,gammaU=0,gammaP=0,fc=0.85,muA=0)
parmsStoch=c(sAalpha=0.17, sAbeta=0.3, minmuA=0.01, sApsi=20.8, sigmaUalpha=1.8, sigmaUbeta=20, sigmaUpsi=12, gamma=0); finess=11;
parms["psiA"]=0; parms["psiU"]=parms["psiP"]=80; RegionSize=300;



##Functions for simulation continuous intraannual dynamics:
#Main function for study:
simpatch3d=function(t,start,parms,Fishing=parms["FP"],eat=FALSE){
  with(as.list(parms),{
	A=start[1]; U=start[2]; P=start[3]; Ueat=Peat=0;
	dA=(A*rA)*(1-(A/K)) - A*(deltaA*U+muA)
	dUeat=a*A*deltaA
	dU=-U*(deltaU*P+muU) 
	dPeat=b*U*deltaU
	dP=-P*(muP+Fishing)
	if(eat) return(list(c(dA,dU,dP,dUeat,dPeat)))
	return(list(c(dA,dU,dP)))
  })
}

simpatch3d2p=function(t,start,parms,Fishing,eat=FALSE){
  with(as.list(parms),{
	if(eat){ A1=start[1]; U1=start[2]; P1=start[3]; A2=start[6]; U2=start[7]; P2=start[8]; U1eat=U2eat=P1eat=P2eat=0; }
	if(!eat){ A1=start[1]; U1=start[2]; P1=start[3]; A2=start[4]; U2=start[5]; P2=start[6]; U1eat=U2eat=P1eat=P2eat=0; }
	dA1=((1-psiA)*A1 + psiA*A2)*rA*(1-(A1/K)) - A1*(deltaA*U1+muA)
	dA2=((1-psiA)*A2 + psiA*A1)*rA*(1-(A2/K)) - A2*(deltaA*U2+muA)
	dU1=-U1*(deltaU*P1+muU)
	dU2=-U2*(deltaU*P2+muU)
	dP1=-P1*(muP + Fishing[1])
	dP2=-P2*(muP + Fishing[2])
	dU1eat=a*A1*deltaA
	dU2eat=a*A2*deltaA
	dP1eat=b*U1*deltaU
	dP2eat=b*U2*deltaU
	if(eat) return(list(c(dA1,dU1,dP1,dU1eat,dP1eat,dA2,dU2,dP2,dU2eat,dP2eat)))
	return(list(c(dA1,dU1,dP1,dA2,dU2,dP2)))
  })
}


##Function to implement functions for intraannual dynamics and pulses at the end
#Throughout, SimType 1, 2, and 3 corresponds to models of 1, 2, and n patches respectively. (2 patch not used)
SpatiotempIter2=function(Temp,parms,  KelpSurv,SigmasU,SigmasP,  Da,Du,Dp, SimType){
	#simulate intraannual dynamics - note that for 1-patch case, 2-patch machinery is used in simulating intraannual dynamics to increase speed:
	if(SimType==3){
		Intraannual = t(apply(Temp[,c(1:3,dim(Temp)[2])], 1, function(x) lsoda(c(x[1:3],0,0), seq(1,2,length=finess), simpatch3d, parms, Fishing=x[4], eat=TRUE)[finess,-1]))
	} else {
		TempTemp = matrix(as.numeric(t(Temp[,c(1:3,dim(Temp)[2])])), nrow=dim(Temp)[1]/2, ncol=8, byrow=T)
		TempTemp2 = t(apply(TempTemp, 1, function(x) lsoda(c(x[1:3],0,0,x[5:7],0,0), seq(1,2,length=finess), simpatch3d2p, parms, Fishing=x[c(4,8)], eat=TRUE)[finess,-1]))
		Intraannual = matrix(as.numeric(t(TempTemp2)), nrow=dim(Temp)[1], ncol=5, byrow=T)
	}
	#zero out any negatives, store final abundances, and implement recruit production by urchins and predators surviving till end of year:
	Temp[,c(1:3,5:6)]=t(apply(pmax(Intraannual,0), 1, function(x) c(x[1:3],x[2]*x[4],x[3]*(parms["rN"]+x[5]))))
	
	#recruit dispersal - note that meaning of Rhos depends on simulation type:
	if(SimType!=1){
		Temp[,7]=pmax(Da%*%Temp[,4],0.05)
		Temp[,8]=(Du%*%Temp[,5]) * SigmasU
		Temp[,9]=(Dp%*%Temp[,6]) * SigmasP * (1-parms["fc"]+(parms["fc"]*(KelpSurv*Temp[,1])/parms["K"]))
	} else {
		if(parms["RhoA"]==0) gamAs=0; if(parms["RhoU"]==0) gamUs=0; if(parms["RhoP"]==0) gamPs=0;
		Temp[,7]=pmax(Temp[,4],0.05)*(1-gamAs) + parms["RhoA"]*gamAs + 0.1 #note the small extra kelp recruits to avoid very low densities
		Temp[,8]=(Temp[,5]*(1-gamUs) + parms["RhoU"]*gamUs) * SigmasU
		Temp[,9]=(Temp[,6]*(1-gamPs) + parms["RhoP"]*gamPs) * SigmasP * (1-parms["fc"]+(parms["fc"]*(KelpSurv*Temp[,1])/parms["K"]))
	}
	#integrate recruitment + kelp survival, and transfer abundances to next year. Also, temporarily save end-of-year abundances for output (trimmed off in StochSim4 time loop).
	return(cbind(Temp[,1:3], cbind(KelpSurv*Temp[,1], Temp[,2:3]) +  Temp[,7:9], Temp[,4:10]))
}


##Build dispersal matrices Da, Du, and Dp (with mean displacement = 0, mean dispersal distance = psi, and mean alongshore patch length = 1.26km)
#Includes option to return a spatial autocorrelation matrix for implementing stochasticity
SpatialMatr=function(psi,npatch,patchsize=RegionSize/npatch,Dplot=FALSE,spatcor=FALSE,expon=FALSE){
	if(psi==0 | (!is.na(patchsize) & patchsize==0)) return(diag(npatch))
	if(is.na(patchsize)){ #For 2-patch case
		D=diag(npatch)*(1-psi);  OffD=cbind(double(npatch),diag(rep(c(psi,0),npatch/2))[,-npatch]); D=D+OffD+t(OffD); if(spatcor) diag(D)=1; 
	} else { #For n-patch case
		D=0*diag(npatch); mod=function(x){ x[x==0]=npatch; return(x); }; #Individuals at 0 after modding by npatch are really in the "last" patch
		if(!expon){ 
			maxDisp=round(max(4*psi,npatch),0); dists=-maxDisp:maxDisp; 
			distsPrs=pnorm(patchsize*(dists+0.5), mean=0, sd=psi*(pi/2)^0.5) - pnorm(patchsize*(dists-0.5), mean=0, sd=psi*(pi/2)^0.5); 
			for(i in 1:npatch) D[i,]=aggregate(distsPrs, by=list(mod((dists+i)%%npatch)), sum)[,2]
		}
		if(expon){ #Alternative case of exponential decline in dispersal or correlation with distance:
			maxDisp=npatch/2; dists=-maxDisp:maxDisp; distsPrs=pmax(1 - (1/(2*log(1+psi)))*log(1+patchsize*abs(dists)), 0)
			for(i in 1:npatch) D[i,]=aggregate(distsPrs, by=list(mod((dists+i)%%npatch)), min)[,2]
		}
	}
	#Adjustment for when constructing correlation matrices to ensure they're not singular:
	if(spatcor&!expon) D=(0.95*(D-min(D))/(max(D)-min(D))+diag(npatch)*0.05)
	if(spatcor&expon) D=((1-max(min(distsPrs),0.05))*(D-min(D))/(max(D)-min(D)) + max(min(distsPrs),0.05))
	if(Dplot){ library(fields); image.plot(x=1:npatch,y=1:npatch,z=D,col=rev(rainbow(1000,start=0,end=0.7))); box(); }; return(D);
}



##Main simulation function:
#3s in stoch indicate stochasticity in: kelp survival, urchin recruitment success, and predator recruitment success. 0s indicate constant conditions (no stochasticity)
#spp and predtype deprecated. To conserve memory saves all desired output (as specified in giveout) in integer form.
StochSim4=function(t,start=parms["K"]*c(1,0.067,0.0167),startvar=FALSE,parms,stoch=c(3,3,0),giveout=1,spp=1,npatch=50,predtype=1,SimType){
	#set up storage objects: columns in Temp correspond to [patch state, recruits produced in patch, recruits arrived in patch]
	Temp=matrix(0,npatch,9); spacetimeA=spacetimeU=spacetimeP=matrix(as.integer(0),npatch,t); 
	if(giveout=="Rhos") Rhos=array(as.integer(0),dim=c(npatch,t,3));

	##startvar makes initial conditions (community states) spatially heterogeneous rather than homogeneous as specified by start (not used in study)
	#rand=function(N,Max,BoundScale=1.5,Min=1){ set.seed(4); return(runif(n=N, min=Min, max=BoundScale*Max)); };
	#if(startvar[1]==FALSE) Temp[,1:3]=rep(start,each=npatch)  else  Temp[,1:3]=c(rand(npatch,1400),rand(npatch,7e3),rand(npatch,32))
	#If needed, uniformly space initial conditions among patches
	rand=function(N,Max,Min=1){ set.seed(4); return(sample(seq(Min,Max,length=N))); }
	if(startvar[1]==FALSE) Temp[,1:3]=rep(start,each=npatch)  else  Temp[,1:3]=c(rand(npatch,2998),rand(npatch,15e3),rand(npatch,80))
	
	#Option to simulate spatially heterogeneous fishing intensity (not used in study):
	FP=rep(as.numeric(parms["FP"]),npatch); Temp=cbind(Temp,FP);
	
	#Specifying spatial dimensions and generating dispersal matrices:
	patchsize=ifelse(SimType==1, 0, ifelse(SimType==2, NA, RegionSize/npatch))
	Da=SpatialMatr(parms["psiA"],npatch,patchsize); Du=SpatialMatr(parms["psiU"],npatch,patchsize); Dp=SpatialMatr(parms["psiP"],npatch,patchsize);
	
	#Incorporating stochasticity - calculating spatial autocorrelation matrices
	#To conserve memory, use a limited pool of disturbance history (size maxTstoch). For stochasticity in years i>maxTstoch, re-sample this pool.
	maxTstoch=8e2; if(t>maxTstoch){ Tsim=t; t=maxTstoch; } else {Tsim=t}
	if(stoch[2]==0) SigmasU=matrix(1,npatch,t); if(stoch[3]==0) SigmasP=matrix(1,npatch,t);
	if(stoch[1]>0){   rawvarsA=mvrnorm(n=t, mu=double(npatch), SpatialMatr(parmsStoch["sApsi"],npatch,patchsize,spatcor=TRUE,expon=FALSE)); }
	if(stoch[2]>0){   rawvarsU=mvrnorm(n=t, mu=double(npatch), SpatialMatr(parmsStoch["sigmaUpsi"],npatch,patchsize,spatcor=TRUE,expon=FALSE)); }
	if(stoch[3]>0){   rawvarsP=mvrnorm(n=t, mu=double(npatch), SpatialMatr(parmsStoch["sigmaUpsi"],npatch,patchsize,spatcor=TRUE,expon=FALSE)); }

	#Incorporating stochasticity - generate actual values to be used throughout simulation
	if(stoch[1]==3){
		sAbeta=function(x) pmax(qbeta(pnorm(x), parmsStoch["sAalpha"], parmsStoch["sAbeta"]), parmsStoch["minmuA"]) #Note truncation to avoid numerical errors
		KelpSurv=t(apply(rawvarsA, 2, sAbeta)); parms["muA"]=0;
	} else { KelpSurv=matrix(parmsStoch["sAalpha"]/(parmsStoch["sAalpha"]+parmsStoch["sAbeta"]),npatch,t); parms["muA"]=0; }
	if(stoch[2]==3){ SigmasU=t(qbeta(pnorm(rawvarsU), parmsStoch["sigmaUalpha"], parmsStoch["sigmaUbeta"])); SigmasU=SigmasU/mean(SigmasU); };
	if(stoch[3]==3){ SigmasP=t(qbeta(pnorm(rawvarsP), parmsStoch["sigmaUalpha"], parmsStoch["sigmaUbeta"])); SigmasP=SigmasP/mean(SigmasP); };
	StochHist=c(1:pmin(Tsim,maxTstoch),sample(1:maxTstoch,pmax(Tsim-maxTstoch,0),replace=T))
	
	for(i in 1:Tsim){
		Temp=SpatiotempIter2(Temp,parms, KelpSurv[,StochHist[i]],SigmasU[,StochHist[i]],SigmasP[,StochHist[i]], Da,Du,Dp, SimType)
		spacetimeA[,i]=as.integer(Temp[,spp])
		#Record densities just before the pulse:
		if(giveout>1 | giveout=="Rhos"){ spacetimeA[,i]=as.integer(Temp[,1]); spacetimeU[,i]=as.integer(Temp[,2]); spacetimeP[,i]=as.integer(Temp[,3]); }
		if(giveout=="Rhos") Rhos[,i,]=Temp[,7:9]
		#Remove densities before pulse for next iteration:
		Temp=Temp[,-c(1:3)]
	}
	if(giveout==0){ return(spacetimeA) }; if(giveout==1){ dev.new(); SpatiotempPlot(spacetimeA, parms,parmsStoch); return(spacetimeA); };
	if(giveout==2){ return(abind(spacetimeA,spacetimeU,spacetimeP,along=3)); }; if(giveout=="Rhos") return(abind(spacetimeA,spacetimeU,spacetimeP,Rhos,along=3));
	if(giveout==3){ nconv=function(x) apply(1e4*x,c(1,2),as.integer); return(abind(spacetimeA,spacetimeU,spacetimeP,nconv(KelpSurv[,StochHist]),nconv(SigmasU[,StochHist]),nconv(SigmasP[,StochHist]),along=3)); }
}


###Analysis functions:
#Turn a 3d array into a matrix
# x is a 3-dimensional array; byrow indicates whether the 3rd dimension be added row-wise? (added column-wise if FALSE)
restructure=function(x,byrow=TRUE){
	if(byrow){  x2=aperm(x,c(1,3,2));  dim(x2)=c(dim(x)[1]*dim(x)[3], dim(x)[2]);  } #stack 3rd dimension by rows
	if(!byrow){  x2=aperm(x,c(1,2,3));  dim(x2)=c(dim(x)[1], dim(x)[2]*dim(x)[3]);  } #stack 3rd dimension by columns
	return(x2)
}


#Make heatmap of values in matrix
# spacetime: matrix of numeric entries
# Xax, Yax: x- and y-values; default to index
# XaxN, YaxN, figtitl: names of axes and overall figure
# Zlim: range of z-values; defaults to (0,max(spacetime))
# cont: color of contour lines; defaults to 0 (not plotted)
# cexAx: font size of axes and axes titles
# contPlot: matrix of numeric entries to plot contours; defaults to spacetime
# cexCont: countour number sizes and line widths
# Conts: specific contour levels to draw; defaults to pretty(zlim, nlevels=10)
# contSpan: for smoothing contours, number of adjacent cells in contPlot over which to smooth values; default is non-smoothed
SpatiotempPlot=function(spacetime,Xax=1:dim(spacetime)[2],Yax=1:dim(spacetime)[1], XaxN="X",YaxN="Y",figtitle="Title",Zlim=c(0,max(spacetime)),
	cont=NULL,cexAx=1,contPlot=spacetime,cexCont=1.5*cexAx,Conts=NULL,contSpan=1){
			require(fields); require(pals);
			if(length(Zlim)==1){
				Zlim=quantile(spacetime,c(Zlim,1-Zlim),na.rm=T) #can provide just a quantile percentage to set limits
				Rng=range(spacetime); if(Zlim[1]<Rng[1]) Zlim[1]=Rng[1]; if(Zlim[2]>Rng[2]) Zlim[2]=Rng[2];
			}
			spacetime[which(is.na(spacetime),arr.ind=TRUE)]=max(Zlim)+1 #NAs painted white
			image.plot(x=Xax,y=Yax,z=t(spacetime), zlim=Zlim, xlab=XaxN, ylab=YaxN,
			  cex.axis=cexAx, cex.lab=cexAx, legend.cex=cexAx, main=figtitle, 
			  col=parula(1e3));  box(); 
			  #Add contour lines to plot
			if(!is.null(cont)){
				if(abs(log(max(contPlot,na.rm=T),10))>2)  options(scipen=-10)
				if(contSpan>1){
					smoo1=t(apply(contPlot, 1, function(x) supsmu(1:ncol(contPlot),x,span=contSpan/ncol(contPlot))$y))
					smoo2=t(apply(smoo1, 1, function(x) supsmu(1:ncol(smoo1),x,span=contSpan/ncol(smoo1))$y))
					contPlot=smoo2
				}
				if(is.null(Conts)) contour(x=Xax,y=Yax,z=t(contPlot),add=T,col=cont,lwd=1.5*cexCont,labcex=cexCont);
				if(!is.null(Conts)) contour(x=Xax,y=Yax,z=t(contPlot),levels=Conts,add=T,col=cont,lwd=1.5*cexCont,labcex=cexCont);
				if(abs(log(max(contPlot,na.rm=T),10))>2)  options(scipen=0)
			}
}




#Sample simulation
x=StochSim4(50,start=c(3e2,4e3,6),parms,stoch=c(3,3,0),giveout=2,npatch=200,startvar=0,SimType=3)
SpatiotempPlot(x[,,1]) #Kelp time series
SpatiotempPlot(x[,,2]) #Urchins
SpatiotempPlot(x[,,3]) #Predators

