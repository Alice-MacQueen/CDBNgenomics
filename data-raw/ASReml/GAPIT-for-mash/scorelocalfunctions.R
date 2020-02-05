
# computation of the autocorrelation
autocor=function(x){abs(cor(x[-1],x[-length(x)]))} 


### computation of the lindley process from scores

lindley=function(scores){
  L=length(scores)
	sl=rep(0,L+1)
	for (i in 1:L){
		sl[i+1]= max(0, (sl[i]+scores[i]))
		}
	return(sl[-1])
	}


# computation of the local score from a lindley process and the region under 
# the maximum peak

scorelocal=function(lind){
	nC=max(lind[,1])
	sloc=matrix(0, nrow=nC, ncol=3)
	colnames(sloc)=c('chr','beg','end')
	for (i in 1:nC){
		tmp=which(lind[,1]==i)		
		list=lind[tmp,3]
		M_loc=which.max(list)
		if(length(which(list[1:M_loc]==0))>0){
			m_loc=max(which(list[1:M_loc]==0))+1			
			}else{m_loc=1}
	sloc[i,]=c(i,lind[tmp[m_loc],2],lind[tmp[M_loc],2])	
	}
	return(data.frame(zone=seq(1,nrow(sloc)),sloc))
}

#computation of the significance threshold -ksi= 1,2,3 or 4) if the distribution of the p-values is uniform

thresUnif=function(L, cor, xi, alpha = 0.05){
  coefs=list(list('a'=c(-5.5,6.76,-5.66,-2.51),'b'=c(-1.22,3.17,-1.99)),
             list('a'=c(2.47,-4.16,-1.82,-4.58),'b'=c(0.37,2.14,-2.35)),
             list('a'=c(2.04,-5.76,1.04,-6.95),'b'=c(2.55,-0.02,-2.31)),
             list('a'=c(0.22,-4.08,1.16,-9.16),'b'=c(3.45,-0.98,-2.33))
             )
  cors=c(cor^3,cor^2,cor,1)
  if (missing(xi) | !(xi %in% 1:4)){
    print('xi should be 1, 2, 3 or 4') 
    }else{
    a=log(L)+ coefs[[xi]]$a%*%cors
    b=coefs[[xi]]$b %*%cors[-1]
    #then compute the threshold:
  thres = ( log(-log(1-alpha)) - a ) / b
  return(thres)
  }
}

# computation of the significative regions from a lindley process given a significance threshold

sig_sl=function(lind,pos, th){
	zones=c(0,0,0)	
	list=lind
	auxpos=pos
	while(max(list)>=th){
	  M_loc=which.max(list)
		if(length(which(list[1:M_loc]==0))==0){ #the peak is at the beginning of the chrom
			m_loc=1
			zones=rbind(zones, c(auxpos[m_loc],auxpos[M_loc],max(list)))
			tmp=which.min[which(list[M_loc+1:length(list)]==0)] #first 0 score after peak
			list=list[tmp:length(list)]
			auxpos=pos[tmp:length(list)]
			}else{	
				m_loc=max(which(list[1:M_loc]==0))			
				max=max(list)
				zones=rbind(zones, c(auxpos[m_loc+1],auxpos[M_loc],max))
				tmp=which(list[M_loc:length(list)]==0) #first 0 score after peak
				if (length(tmp)>0){
				  auxpos=auxpos[c(1:m_loc,(min(tmp)+M_loc):length(list))]
				  list=list[c(1:m_loc, (min(tmp)+M_loc):length(list))]
				  }else{ #if the peak is at the end of the chromosome
				    auxpos=auxpos[1:m_loc]
				    list=list[1:m_loc]
				    }				
				}
	  }
	zones=matrix(zones, ncol=3)
	zones=data.table(beg=zones[,1],end=zones[,2],peak=zones[,3])
	if (nrow(zones)>1){zones=zones[-1,]}
	return(zones)
	}

### Estimation of Gumble coefficients.

# Estimation of the coeficients of the Gumbel distributions

lineGumb=function(x){
  x1tmp=x
  if (length(table(x1tmp)) >  5){
    Fn = ecdf(x1tmp)
    Fnx=Fn(x1tmp)
    filtre = ( Fnx < 1 - 10**(-1) )  & ( x1tmp > 1 ) #sinon pbs numq
    lm0 = lm(log(-log(Fnx[filtre])) ~ x1tmp[filtre])$coefficients
  }else{lm0=c(0,0)}  
  return(lm0)
}

# Estimation of the coeficients of the polynomes to compute the Gumbel coefficients 
# depending on the length of the chromosomes and the chromosome autocorrelation

coefsGumb=function(mydata, Ls=seq(10000,80000,10000), nSeq=5000, degMod=2){
  bins=seq(0.025,0.975,0.05) 
  coefs=array(0, dim=c(length(bins),5, length(Ls)))  
  coefs[,5,]=matrix(rep(Ls,length(bins)), ncol=length(Ls),byrow=TRUE)
  as=matrix(0, ncol=length(Ls), nrow=length(bins))
  bs=matrix(0, ncol=length(Ls), nrow=length(bins))
  xx=seq(0,max(mydata$lindley), 1)
  for (j in 1:length(Ls)){
    coefs[,1,j]=bins  
    len=Ls[j]
    tmp=sample(nrow(mydata)-len,nSeq,replace=F) #samples the sequences' beginnings
    DTL=data.table(LocScore = vector('numeric', length=nSeq),cors= vector('numeric', length=nSeq)) 
    for (l in 1:nSeq){
      DTL[l]=mydata[seq(tmp[l],length.out=len),.(max(lindley(score)),autocor(pval))]
    }
    DTL[,bin:=which.min(abs(bins-cors)),cors]
    binNE=sort((unique(DTL$bin)))
    coefs[binNE,4,j]=table(DTL$bin)
    DTCoef=DTL[,.(lineGumb(LocScore)),bin]
    coefs[unique(DTL$bin),c(2,3),j]=matrix(DTCoef$V1, ncol=2, byrow=T)
    ys=coefs[,c(2,3),j]%*%rbind(rep(1,length(xx)), xx)
    pdf(paste('GumbelinesAllxi',xi,'M',len,'.pdf', sep=''))
    for (i in binNE){
      Fn=ecdf(DTL[bin==i,LocScore]) 
      tmp=DTL[bin==i,.(LocScore,Fn=Fn(LocScore)),][(Fn< 1 - 10**(-5)) & (LocScore > 1)]
      if(i==min(binNE)){
        plot(tmp$LocScore, log(-log(tmp$Fn)), col=i,main=paste("M = ",len,sep=""), xlim=quantile(DTL$LocScore,c(0,0.9)),ylim=range(ys[binNE,xx<quantile(DTL$LocScore,0.8)], na.rm=T), pch=20, xlab='Local Score', ylab='log(-log(Fn(LS)))')
        }else{points(tmp$LocScore, log(-log(tmp$Fn)), col=i, pch=20)
      }
      if (sum(ys[i,]!=0)){lines(xx, ys[i,], col=i, lwd=1.5)}
    }
    legend("topright",legend=bins[binNE], pch=16,col=seq(1:length(binNE)), title='rho')
    dev.off()
    as[,j]=as.numeric(coefs[,2,j])
    bs[,j]=as.numeric(coefs[,3,j])		
  }
  
  auxWhich=which(coefs[,4,]> 150)
  rhos=coefs[,1,][auxWhich]
  auxAs=as[auxWhich]-log(coefs[,5,][auxWhich])
  auxBs=bs[auxWhich]
  
  fitA=lm((auxAs)~I(rhos))  
  fitB=lm((auxBs)~I(rhos))  
  
  pdf('FitAandB.pdf')
  par(mfrow=c(1,2))
  plot(rhos, auxAs, pch=16, col=coefs[,5,][auxWhich]/10000, xlab='rho', ylab='a-log(M)')
  xslm=seq(0,1,0.01)
  lines(xslm, fitA$coefficients%*%rbind(rep(1,length(xslm)), xslm))
  plot(rhos, auxBs, pch=16, col=coefs[,5,][auxWhich]/10000, xlab='rho', ylab='b')
  lines(xslm, fitB$coefficients%*%rbind(rep(1,length(xslm)), xslm))
  legend(min(rhos), max(auxBs), legend=sort(unique(coefs[,5,][auxWhich])), col=unique(coefs[,5,][auxWhich]/10000), pch=16, title='M values')
  
  dev.off()
  
  return(list(aCoef=fitA$coefficients, bCoef=fitB$coefficients))
  }

# Computation of the significance threshold given the computed polynomes for computing
# the Gumbel coefficients depending of length and autocorrelation

threshold=function(L, cor, aCoef, bCoef, alpha = 0.05){
  degA=length(aCoef)
  degB=length(bCoef)
  a=log(L)
  b=0
  for (i in 1:degA){
    a=a+aCoef[i]*cor^(i-1)  
    }
  for (i in 1:degB){
    b=b+bCoef[i]*cor^(i-1)  
  }
  #then compute the threshold:
  thres = ( log(-log(1-alpha)) - a ) / b
  return(thres)
  }


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
