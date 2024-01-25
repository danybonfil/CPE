#### Examines the effect of GET on popularity of presidents in Brazil and Mexicos ###
#### The code reads in the data,  estimates many different specifications for     ###
#### each countries, impelments diagnostics. Results are then "collected" and most###
#### tables and graphs are assmelbed at the end. 								  ### 
#### Code replicates Table 3 in the JOP table as well as tables and figures in    ###
#### the appendix. Requires data files "data_popularityBR_AMELIA.RData" and 	  ###
#### "data_popularityMX_AMELIA.Rdata"  with monlhtly (imputed) data 			  ###
#### The code, itself, starts in line 600. Prior do that, functions are declared  ###

rm(list=ls(all=TRUE))
library(car)
library(lmtest)
library(xtable)
library(Hmisc) #this is for Lag command
library(Amelia)
library(fUnitRoots) #for formal tests of unit Roots
library(Zelig)
library(forecast)
library(reshape)
library(mediation)

if(grep("zucco",setwd("~"))==1){
	the.path <- "~/Dropbox/Data/Paper-VotoEconomico/"
	}else{
	the.path <- "~/Dropbox/ECOVOTING/"
	}
setwd(the.path)

##### Declare some functions #####################################

plot.pred <- function(x,vd="pop",crisis=T,mim=1,data=imp,newmar=F){ 
	#plot predicted values and actual values
	#x is a list of MI regressions or a stand alone regression; 
	#vd is the dep var; crisis whether to plot crisis;
	#mim number of mi set to use for "actual" data, popularity is mean across imps
    if(is.element("ls-mi",class(x))){#list of zelig mi-ls
   	cat("Plotting average predicted and actual values across imputations\n")
    	#With Zelig < 4, omit $result
    	used <-  union(c(vd,"date","president"),
    				attr(x[[mim]]$result$terms,"term.labels"))
        the.set <- data[[mim]][,used]
    	the.set[,vd] <- apply(#average dep variable across imputations
    						sapply(data,function(x){x[,"pop"]}),1,mean)
    	#With Zelig < 4, omit $result
    	fitted.values <-apply(sapply(x,function(x)fitted(x$result)),1,mean)	
    	rsq <- mean(sapply(x,function(x)summary(x)$r.sq))
    	}else{
    if(is.element("Arima",class(x[[1]]))){
     	cat('plot.pred no implemented for Arima\nNeed to copy code from plot.dev\n')
     	break()
    	}
    if(is.element("ls",class(x))&
       is.element("zelig",class(x))){#a single  zelig ols (v>=4.0)
		used <-  union(c(vd,"date","president"),attr(x$result$terms,"term.labels"))
    	the.set <- na.omit(x$data[,used]) #stand alone "might"have missing data
    	fitted.values <- x$result$fitted.values
    	rsq <- summary(x)$r.squared}
     }#endelse
    min.y <- 0
    if(mean(the.set[,vd],na.rm=T)<1){max.y<-1}else{max.y<-100}
     if(newmar==F){par(mfrow=c(1,1),mar=c(3,3.2,1,1))}else{par(mfrow=c(1,1))}
     plot(as.Date(the.set$date), the.set[,vd],type="n",xaxt="n",ylab="",xlab="",
    					ylim=c(min.y,max.y))
    all.presidents <- unique(the.set$president)
    for(i in 1:length(all.presidents)){
    		if(i/2==round(i/2)){#for "even" presidents
    	polygon(x=c(min(as.Date(the.set$date)[the.set$president==
    				all.presidents[i]],na.rm=T),
    						max(as.Date(the.set$date)[the.set$president==
    				all.presidents[i]],na.rm=T),
    						max(as.Date(the.set$date)[the.set$president==
    				all.presidents[i]],na.rm=T),
    						min(as.Date(the.set$date)[the.set$president==
    				all.presidents[i]],na.rm=T)),
						y=c(min.y,min.y,max.y,max.y),border=NA,col=gray(0.9))
    			
    		}
    		actual.pop <- the.set[which(the.set$president==all.presidents[i]),vd]
			pred.pop <- fitted.values[which(the.set$president==all.presidents[i])]
    		lines(as.Date(the.set$date)[which(the.set$president==all.presidents[i])],
                  actual.pop,col=gray(0)) 
            lines(as.Date(the.set$date)[which(the.set$president==all.presidents[i])],
                  pred.pop,col=gray(.65),lwd=2) 
              }#end i loop
   	#Add president's names to graph (find "central" point of term to plot, then plot)
   	write.pres <- do.call("c",lapply(
    				by(the.set$date,the.set$president,mean),
    				function(x){as.Date(x,origin="1970-01-01")}))
    lab.names <- gsub(".*\\s(.*)","\\1",names(write.pres)) #get only last names
    text(write.pres,95,labels=lab.names,cex=0.7)
    lim.obs <- as.Date(paste(gsub("-.*","",range(
    			as.Date(the.set$date),na.rm=T)),"-01-01",sep=""))
    ax.points <- seq(lim.obs[1],lim.obs[2],by="years")[-1]
    axis(side=1,at=ax.points,labels=substr(ax.points,1,4),las=2,cex=1)
    text(max(as.Date(the.set$date),na.rm=T),
    		min.y+2,label=paste("R2=",
    		round(rsq,2)),pos=2,cex=1.1)
    mtext("Popularity",side=2,line=2.2,cex=1.1)
    if(crisis==T){points(as.Date(the.set$date)[which(the.set$polcrisis==1)],
    		fitted.values[which(the.set[,"polcrisis"]==1)],pch=20)}
    legend(x="bottomright",lwd=c(1,2,NA),lty=c(1,1,NA),col=c(gray(0),gray(0.65)),
    		legend=c("Actual","Predicted",paste("R2=",round(rsq,2),sep="")),
    		bg=gray(0.99),box.col=gray(1),inset=0.01)
    print(summary(x))  
}

plot.coefs <- function(x,v,varname="Variable",xlims=NULL,substantive=NULL,shade=F,leftmar=11){
			#x is a list of regression outputs
			#v is the variable name as it appears in the output
			#varname is the name of the variable to appear in the plot
		models <- length(x)	
	coefs <- sapply(x,function(x){x[v,'coef']})
	ses <- sapply(x,function(x){x[v,'se']})
	ub <- coefs + qnorm(0.95) * ses
	lb <- coefs - qnorm(0.95) * ses
	if(models<=3){offsets <- seq(-.5,.5,len=models)
					divider <- setdiff(seq(-.5,.5,len=(2*models)-1),offsets)
			   }else{offsets <- seq(-.65,.65,len=length(x))
			   	divider <- setdiff(seq(-.65,.65,len=(2*models)-1),offsets)}
	if(is.null(xlims)){
		      xlims <- c(-max(c(abs(ub),abs(lb)),na.rm=T),
			  max(c(abs(ub),abs(lb)),na.rm=T))
			  }else{xlims <- c(-xlims,xlims)}
	ylims<-c(1,-1)
	par(mar=c(2,leftmar,1,1))
	plot(coefs,offsets,xlim=xlims,ylim=ylims,
			ylab="",xlab="",yaxt="n",bty="n",type="n")
	if(shade==T){
		polygon(x=c(xlims[1],rep(xlims[2],2),
		    xlims[1]),y= c(rep(ylims[1],2),
		    rep(ylims[2],2)),col=gray(0.85),border=NA)
		
		abline(h=divider,lty=4,col=gray(1),lwd=4)}#end dark
	abline(v=0,lty=3)
	segments(x0=lb, y0=offsets, x1 = ub, y1 = offsets, lwd=2)

	if(length(x)<=3){
		pchs <- 21:(20+(length(x)))
		bgs <- c(gray(0.9),gray(0.1))
		cols <- c(gray(0.1),gray(0.9))
	}else{
		pchs <- sort(rep(21:(20+ceiling(models/2))))
		bgs <- sort(rep(c(gray(0.9),gray(0.1)),models/2))
		cols <- sort(rep(c(gray(0.9),gray(0.1)),models/2),decreasing=T)
	}
	points(coefs,offsets,pch=pchs,cex=3,bg=bgs,col=cols)
	axis(side=2,at=0,labels=varname,las=2,tick=F,cex.axis=2,hadj=1)
	if(is.null(substantive)==F){
		#y0=offsets[1]
		#y1=offsets[1]
		#x0=coefs[1]
		#x1=coefs[1]+(xlims[2]-coefs[1])/2
		#arrows(x0,y0,x1,y1,code=2,length=0.1,col=2,lwd=2)
		#text(x1,y1,labels=round(substantive[v,],2), cex=1.5,pos=4)
		text(coefs[1],offsets[1]-0.1,labels=round(substantive[v,],2), cex=1.5,pos=3,col=2)
	}
	return(list(pchs=pchs,bgs=bgs,cols=cols))
}

plot.dev <- function(x,vd="pop",crisis=T,mim=1,data=imp,plot.range=NULL){ 
	#plot how presidents 
	#deviate from preidcted values
	#x is a list of regressions; vd is the dep var; crisis whether to plot crisis;
	#mim number of mi set to use for "actual" data, popularity is mean across imps
	#plot.range is integet to override default yaxis limits
    if(is.element("ls-mi",class(x))){#list of zelig mi-ls
    	cat("Using average predicted and actual values across imputations\n")
    	used <-  union(c(vd,"date","president"),
    				attr(x[[mim]]$result$terms,"term.labels"))
        the.set <- data[[mim]][,used]
    	the.set[,vd] <- apply(#average dep variable across imputations
    						sapply(data,function(x){x[,"pop"]}),1,mean)
    	#Worked with Zelig 
    	#fitted.values <- apply(sapply(x,function(x)fitted.values(x)),1,mean)
    	fitted.values <-apply(sapply(x,function(x)fitted(x$result)),1,mean)	
    	rsq <- mean(sapply(x,function(x)summary(x)$r.sq))
    	}else{
    if(is.element("ls",class(x))&
       is.element("zelig",class(x))){#a single  zelig ols (v>=4.0)
		used <-  union(c(vd,"date","president"),attr(x$result$terms,"term.labels"))
    	the.set <- na.omit(x$data[,used])
    	fitted.values <- x$result$fitted.values
    	rsq <- summary(x)$r.squared}
     if(is.element("Arima",class(x[[1]]))){
       used <-  setdiff(
     				union(c(vd,"date","president"),names(x[[mim]]$coef)),
     				c("ar1","intercept"))
       the.set <- data[[mim]][,used]
       my.fitted <- function(y){#function used to correct series "source"
       		 y$series <- "imp[[1]]$pop[-1]" #because it was invoked inside a lapply
       		 fitted(y) #and extract fitted values
       }
       fitted.values <- apply(sapply(x,my.fitted),1,mean)
       #compute R2
       allvar <- sapply(x,function(y){diag(y$var.coef)}) 	
	   tot.var <- sapply(data,function(y)var(y$pop))
	   res.var <- sapply(x,function(y)var(residuals(y)) )
	   rsq <- mean((tot.var-res.var)/tot.var)
       }
       }#end else
    if(is.null(plot.range)){	
    if(mean(the.set[,vd],na.rm=T)<1){max.y<-0.5;min.y=-0.5}else{max.y<-50;min.y<--50}
    }else{max.y <- plot.range;min.y <- -plot.range}
     par(mfrow=c(1,1),mar=c(3,3.2,1,1))
     plot(as.Date(the.set$date), the.set[,vd],type="n",xaxt="n",ylab="",xlab="",
    					ylim=c(min.y,max.y))
    all.presidents <- unique(the.set$president)
    abline(h=0,lwd=1,lty=3,col= gray(0.5))
    for(i in 1:length(all.presidents)){
    		if(i/2==round(i/2)){#for "even" presidents
    	polygon(x=c(min(as.Date(the.set$date)[the.set$president==
    						all.presidents[i]],na.rm=T),
    						max(as.Date(the.set$date)[the.set$president==
    						all.presidents[i]],na.rm=T),
    						max(as.Date(the.set$date)[the.set$president==
    						all.presidents[i]],na.rm=T),
    						min(as.Date(the.set$date)[the.set$president==
    						all.presidents[i]],na.rm=T)),
						y=c(min.y,min.y,max.y,max.y),border=NA,col=gray(0.9))  			
    		}
    	segments(x0=min(as.Date(the.set$date)[the.set$president==
    					all.presidents[i]],na.rm=T), 
				 x1=max(as.Date(the.set$date)[the.set$president==
				 		all.presidents[i]],na.rm=T),
				 y0=0, y1=0, col = gray(0.5),lwd=1,lty=3)
    		dev<-the.set[which(the.set$president==all.presidents[i]),vd]-	
    				fitted.values[which(the.set$president==all.presidents[i])]
    		lines(as.Date(the.set$date)[which(the.set$president==all.presidents[i])],
                  dev,
                  col=gray(0),lwd=2)
              }#end i loop
   	#Add president's names to graph (find "central" point of term to plot, then plot)
   	write.pres <- do.call("c",lapply(
    				by(the.set$date,the.set$president,mean),
    				function(x){as.Date(x,origin="1970-01-01")}))
    lab.names <- gsub(".*\\s(.*)","\\1",names(write.pres)) #get only last names
    text(write.pres,max.y*.95,labels=lab.names,cex=0.7)
    lim.obs <- as.Date(paste(gsub("-.*","",range(
    			as.Date(the.set$date),na.rm=T)),"-01-01",sep=""))
    ax.points <- seq(lim.obs[1],lim.obs[2],by="years")[-1]
    axis(side=1,at=ax.points,labels=substr(ax.points,1,4),las=2,cex=1)
    text(max(as.Date(the.set$date),na.rm=T),
    		min.y*1.05,label=paste("R2=",
    		round(rsq,2)),pos=2,cex=1.1)
    mtext("<-Worse than predicted              Better than predicted->",
    		side=2,line=2.2,cex=1.1)
     if(crisis==T){points(as.Date(the.set$date)[which(the.set$polcrisis==1)],
    		the.set[which(the.set[,"polcrisis"]==1),vd]-
    		fitted.values[which(the.set[,"polcrisis"]==1)],pch=20)}
    print(summary(x))  
}

my.amelia <- function(xx){
		#xx is a list of arima results or of zelig results
		##if object is a list of arima regressions
		##function computes amelia corrected summaries from a list of regression
		##if object is a list from zelig amelia or zelig stand aline
		##simply rearranges the output to make it the same as arima
		if(is.element("Arima",class(xx[[1]]))){
		cat("XX is a list of objects of class ARIMA\n")
		#collect var
		allcoefs <- sapply(xx,coef) #collect coefs
			allvar <- sapply(xx,function(x){diag(x$var.coef)}) 	
			tot.var <- sapply(imp,function(x)var(x$pop))
			res.var <- sapply(xx,function(x)var(residuals(x),na.rm=T) )
			r2 <- (tot.var-res.var)/tot.var	
		ameliar2s <- mean(r2)
		ameliacoefs <- apply(allcoefs,1,mean) #aggregate coefs
		#The variance of the point estimate is 
		#the average of the estimated variances from within each completed data set, (var1) 
		#plus the sample variance in the point estimates across the data sets  (var2)
		#(multiplied by a factor that corrects for the bias because m <1)
		#Checked this against Zelig's calgulations, and it is correct
			var1 <- apply(allvar,1,mean)
			var2 <- apply(sweep(x=allcoefs,MARGIN=1,apply(allcoefs,1,mean),"-")^2,1,sum)/
					(length(xx)-1)
		ameliases <- sqrt(var1 + var2 * (1 + 1/length(xx)))	
		ameliat <- ameliacoefs/ameliases
		ameliap <- 2*pt(abs(ameliat),
					 		length(xx[[1]]$residuals)-length(ameliacoefs),
					 		lower.tail=FALSE)	
		out <- data.frame(coef=round(ameliacoefs,4),se=round(ameliases,4),
					 t=round(ameliat,4),
					 p.value=round(ameliap,4),
					 sig=ifelse(ameliap<0.01,"**",
					 	 ifelse(ameliap>=0.01&ameliap<0.05,"*",
					 	 ifelse(ameliap>=0.05&ameliap<0.1,".",""))))
		out <- rbind(out,R2=c(ameliar2s,rep(NA,ncol(out)-1)))
		}##END PROCESSING ARIMA MI OUTPUT
		if(is.element("zelig",class(xx[[1]]))){#if its MI zelig output
		cat("xx is a list of objects of class zelig (v>4.0)\n")
			out <- summary(xx)$coef;  #put zelig output in same format as others
			rownames(out) <- gsub("\\(Intercept\\)","intercept",rownames(out))
			rownames(out) <- gsub("TRUE","",rownames(out))
			colnames(out) <- c("coef","se","t","p.value")
			rsq <- mean(sapply(xx,function(x){summary(x)$r.sq}))
			out <- data.frame(rbind(out,R2=c(rsq,rep(NA,3))))
		 	out$sig <- ifelse(out$p.value<0.01,"**",
						 ifelse(out$p.value<0.05&out$p.value>=0.01,"*",
						 ifelse(out$p.value>=0.1,"",".")))	
		}
		if(is.element("numeric",class(xx[[1]]))){#ZELIG standalone input
		cat("xx is single objects of class zelig (v>4.0)\n")
			out <- summary(xx)$coef;
			rownames(out) <- gsub("\\(Intercept\\)","intercept",rownames(out))
			rownames(out) <- gsub("TRUE","",rownames(out))
			colnames(out) <- c("coef","se","t","p.value")
			rsq <- summary(xx)$r.sq
			out <- data.frame(rbind(out,R2=c(rsq,rep(NA,3))))
		 	out$sig <- ifelse(out$p.value<0.01,"**",
						 ifelse(out$p.value<0.05&out$p.value>=0.01,"*",
						 ifelse(out$p.value>=0.1,"",".")))	
		}
		return(out)
	} 

diag.ts <- function(the.reg,to.plot=F){
	#the.reg is a list of regression models on MI datasets (5), or zelig on 1 dataset
	#for ARIMA object, function computes diagnostics statistics & plots
	#for zelig list object, fucntion creates diagnostics plots
	#for zelig stand alone, function creates disangostics plots
	#for ARIMA, to.plot will aways be T
	if(is.element("Arima",class(the.reg[[1]]))){#Arima input
			to.plot<-T
		arsq <- NA # mean(sapply(the.reg,function(x)summary(x)$r.sq))	
		aa <- lapply(the.reg,function(x)auto.arima(residuals(x),ic="bic")) 	
		aas <- sapply(aa,function(x){c(p.ar=length(x$model$phi),
	  	                      d=length(x$model$Delta),
							  q.MA=length(x$model$theta))}) 
		aaas <- apply(aas,1,mean) 
				#average corrections needed on autoarima 
				#use with caution because plots are better tools
		urt1 <- lapply(the.reg,function(x)unitrootTest(residuals(x), lags = 1) )
		urt1s <- sapply(urt1,function(x){c(DFstat=x@test$statistic,
									   round(x@test$p.value,4))})
		aurt1s <- apply(urt1s,1,mean)	  #significant rejects unit root				
		dw <- lapply(the.reg,function(x)durbinWatsonTest(as.numeric(residuals(x))))
		dws <-  sapply(dw,function(x){c(#DWautocorr=x$r,
											DWstat=x#$dw#,
											#DWpval=x$p
											)})
		adws <- c(DWstat.DW=mean(dws), DWpval=NA)
     	box <- lapply(the.reg,function(x)Box.test(residuals(x)))
		boxs <-  sapply(box,function(x){c(Box.stat=as.numeric(x$statistic),
										  Box.pval=x$p.value)})
		abox <- c(BPstat=mean(boxs[1,]),BPpval=1-pchisq(mean(boxs[1,]),1))
		bg <-lapply(the.reg,function(x)#tests on the residuals, is this it?
				bgtest(res~1,order=2,data=data.frame(res=x$residuals),fill=0)
				)
		bgs <- sapply(bg,function(x){c(BGstat=x$statistic,
								BGpval=x$p.value)})
		cat("Breush-Godfrey for MI\n")
		print(bgs)
		abgs <- apply(bgs,1,mean) #average BG statistics across 5 sets 
		
		out<- round(c(R2=arsq,#AIC=aa$aic,BIC=aa$bic,AICC=aa$aicc,
			aaas,
			aurt1s,
			avDW=adws,
			abox,
			abgs
			),3)	
		cat("Reporting AVERAGE statistics across imputed sets\n")		
	#plot for ARIMA
	par(mfrow=c(3,5),mar=c(2,4,2,0.5))
	lapply(the.reg,function(x){
		plot(residuals(x)/sd(residuals(x)),type="l",xlab="",
		main=paste("Residuals"),ylab="Standarized Residuals")
		abline(h=0)})
	lapply(the.reg,function(x){acf(residuals(x),main="Residuals")})
	lapply(the.reg,function(x){pacf(residuals(x),main="Residuals")})
	}#end Arima Case (compute & plot)		
	
	if(is.element("zelig",class(the.reg[[1]]))){#ZELIG MI (list) input
	#plot for ZELIG (Zelig 4.0 requires $result to work)
	arsq <- mean(sapply(the.reg,function(x)summary(x$result)$r.sq))	
		aa <- lapply(the.reg,function(x)auto.arima(as.numeric(residuals(x$result))
											,ic="bic")) 	
		aas <- sapply(aa,function(x){c(p.ar=length(x$model$phi),
	  	                      d=length(x$model$Delta),
							  q.MA=length(x$model$theta))}) 
		aaas <- apply(aas,1,mean) 
				#average corrections needed on autoarima 
				#use with caution because plots are better tools
		urt1 <- lapply(the.reg,function(x)unitrootTest(residuals(x$result), lags = 1) )
		urt1s <- sapply(urt1,function(x){c(DFstat=x@test$statistic,
									   round(x@test$p.value,4))})
				###Significant DF rejects unitroot
		aurt1s <- apply(urt1s,1,mean)	  #significant rejects unit root				
		dw <- lapply(the.reg,function(x)dwtest(lm(x$formula,data=x$data),
										alternative="two.sided"))
		dws <-  sapply(dw,function(x){c(DWstat=x$statistic,
											DWpval=x$p.value
											)})
		adws <- apply(dws,1,mean) #average DW statistics across 5 sets
		box <- lapply(the.reg,function(x)Box.test(residuals(x$result)))
		boxs <-  sapply(box,function(x){c(BPstat=as.numeric(x$statistic),
										  BPpval=x$p.value)})
		abox <- c(BPstat=mean(boxs[1,]),BPpval=1-pchisq(mean(boxs[1,]),1))
		cat("Box–Pierce  for MI\n")
		print(boxs)	
		bg <-lapply(the.reg,function(x)bgtest(x$formula,order=2,data=x$data,fill=0))
		bgs <- sapply(bg,function(x){c(BGstat=x$statistic,
								BGpval=x$p.value)})
		cat("Breush-Godfrey for MI\n")
		print(bgs)
		abgs <- apply(bgs,1,mean) #average BG statistics across 5 sets 
		out<- round(c(R2=arsq,
			aaas,
			aurt1s,
			adws,
			abox,
			abgs
			),3)	
	cat("Reporting AVERAGE statistics across imputed sets\n")	
	par(mfrow=c(3,5),mar=c(2,4,2,0.5))
	lapply(the.reg,function(x){
		plot(residuals(x$result)/sd(residuals(x$result)),type="l",xlab="",
		main=paste("Residuals"),ylab="Standarized Residuals");
		abline(h=0)})
	lapply(the.reg,function(x){acf(residuals(x$result),main="Residuals")})
	lapply(the.reg,function(x){pacf(residuals(x$result),main="Residuals")})}
	
	if(is.element("numeric",class(the.reg[[1]]))){#ZELIG standalone input
	#plot for ZELIG (Zelig 4.0 requires $result to work)
	par(mfrow=c(1,3),mar=c(2,4,2,0.5))
	plot(residuals(the.reg$result)/sd(residuals(the.reg$result)),type="l",xlab="",
		main=paste("Residuals"),ylab="Standarized Residuals");
	abline(h=0)
	acf(residuals(the.reg$result),main="Residuals")
	pacf(residuals(the.reg$result),main="Residuals")
	arsq <- summary(the.reg$result)$r.sq
		aa <- auto.arima(residuals(the.reg$result),ic="bic") 	
		aaas <- c(p.ar=length(aa$model$phi),
	  	                      d=length(aa$model$Delta),
							  q.MA=length(aa$model$theta)) 

		urt1 <- unitrootTest(residuals(the.reg$result), lags = 1) 
		aurt1s <- c(DFstat=urt1 @test$statistic,
					round(urt1 @test$p.value,4))
				###Significant DF rejects unitroot
		        ###Significant rejects unit root				
		dw <- dwtest(lm(the.reg$formula,data=the.reg$data),
										alternative="two.sided")
		adws <-  c(DWstat=dw$statistic,
				DWpval=dw$p.value)
		box <- Box.test(residuals(the.reg$result))
		abox <-  c(BPstat=as.numeric(box$statistic),
										  BPpval=box$p.value)
		bg <-NA#bgtest(the.reg$formula,order=2,data=the.reg$data,fill=0))
		abgs <- c(BGstat=NA,BGpval=NA)# c(BGstat=bg$statistic,
						#		BGpval=bg$p.value)})
		out<- round(c(R2=arsq,
			aaas,
			aurt1s,
			adws,
			abox,
			abgs
			),3)	
	cat("Reporting statistics from stand alone regression\n")
	}
		cat("\nUnit Root: Augmented Dickey-Fuller test",
			"\nSerial correlation tests:",
			"\n\tBreusch–Godfrey Lagrange multiplier test. h0=No serial correlation up to p=2",
			"\n\tBox–Pierce test. ho=independence in a given time series (lag=1)",
			"\n\tDurbin–Watson H statistic: h0=First order autocorrelation is 0\n")
	
	return(t(out))
	#library(tseries)
	#pp.test(residuals(imod1)) #: Tests the null hypothesis that the time series is non-
	#kpss.test(residuals(imod1)): 
	}

stack.regs <- function(x){##TO create a table of coefficients
	require(reshape)
	x$vars <- rownames(x)	
	out <- melt.data.frame(x, "vars", variable_name = "desc")
	out <- subset(out,is.element(desc,c("coef","se","p.value")))
	out <- out[order(out$vars),]
	names(out)[3]<-"v"
	out$v <- round(as.numeric(out$v),3)
	return(out)
	}

edited.table <- function(x){#this function is used for all tables below...
		  tmp <- x[c(
		  	grep("get",x$vars),
		  	grep("irates",x$vars),
			grep("log.commodities",x$vars),
			grep("rel.delta.pib06",x$vars),
			grep("time.pres",x$vars),
			grep("polcrisi",x$vars),
			grep("polcrise",x$vars),
			grep("intercept",x$vars),
			grep("lag.pop",x$vars),
			grep("^ar",x$vars),
			grep("^ma",x$vars),
			grep("R2",x$vars)),]
		 tmp$desc <- ifelse(tmp$desc=="coef",
		 			as.character(tmp$vars),
		 			as.character(tmp$desc))	
		 tmp$vars <- NULL
		 rownames(tmp) <- 1:nrow(tmp)
		 if(sum(apply(is.na(tmp[,-1]),1,sum)==(ncol(tmp)-1))>0){
		 	tmp <- tmp[-which((apply(is.na(tmp[,-1]),1,sum)==(ncol(tmp)-1))),]
		 }
		 return(tmp)
			}

clean.regs <- function(x){
	x$vars<-gsub("(Datafolha)|(Consulta.Mitofsky)"," 1",x$vars)
	x$vars<-gsub("(Ibope)|(Grupo.Reforma)"," 2",x$vars)
	x$vars<-gsub("(Sensus)|(OPRM)"," 3",x$vars)
	x$vars<-gsub("(Vox)"," 4",x$vars)
	x$vars<-gsub("instituteNone","Imputed",x$vars)
	return(x)}
	
unit.resp <- function(x,the.var="geti",lag.var="lag.geti"){
	if(class(x)=="data.frame"){tmp<-x
	b0 <-tmp[the.var,"coef"]
	tryb1 <-try(tmp[lag.var,"coef"],silent=T)
	b1<-if(is.na(tryb1)){0}else{tmp[lag.var,"coef"]}
	rho <- tmp["lag.pop","coef"]	
	}else{
	tmp <-  summary(x)
	b0 <-tmp$coefficients[the.var,"Value"]
	tryb1 <-try(tmp$coefficients[lag.var,"Value"],silent=T)
	b1<-if(class(tryb1)=="try-error"){0}else{tmp$coefficients[lag.var,"Value"]}
	rho <- tmp$coefficients["lag.pop","Value"]
	}
	t1 <- b0
	t2 <- b0+b1+rho*b0
	t3 <- b0+b1+rho*b0+rho*b1+rho^2*b0
	t4 <- b0+b1+rho*b0+rho*b1+rho^2*b0+rho^2*b1+rho^3*b0
	t5 <- b0+b1+rho*b0+rho*b1+rho^2*b0+rho^2*b1+rho^3*b0+rho^3*b1+rho^4*b0
	t6 <- b0+b1+rho*b0+rho*b1+rho^2*b0+rho^2*b1+rho^3*b0+rho^3*b1+rho^4*b0+
		  rho^5*b1+rho^5*b0
	t7 <- b0+b1+rho*b0+rho*b1+rho^2*b0+rho^2*b1+rho^3*b0+rho^3*b1+rho^4*b0+
		  rho^5*b1+rho^5*b0+rho^6*b1+rho^6*b0
	t8 <- b0+b1+rho*b0+rho*b1+rho^2*b0+rho^2*b1+rho^3*b0+rho^3*b1+rho^4*b0+
		  rho^5*b1+rho^5*b0+rho^6*b1+rho^6*b0+rho^7*b1+rho^7*b0
	t9 <- b0+b1+rho*b0+rho*b1+rho^2*b0+rho^2*b1+rho^3*b0+rho^3*b1+rho^4*b0+
		  rho^5*b1+rho^5*b0+rho^6*b1+rho^6*b0+rho^7*b1+rho^7*b0+
		  rho^8*b1+rho^8*b0
	t10<- b0+b1+rho*b0+rho*b1+rho^2*b0+rho^2*b1+rho^3*b0+rho^3*b1+rho^4*b0+
		  rho^5*b1+rho^5*b0+rho^6*b1+rho^6*b0+rho^7*b1+rho^7*b0+
		  rho^8*b1+rho^8*b0+rho^9*b1+rho^9*b0
	t11<- b0+b1+rho*b0+rho*b1+rho^2*b0+rho^2*b1+rho^3*b0+rho^3*b1+rho^4*b0+
		  rho^5*b1+rho^5*b0+rho^6*b1+rho^6*b0+rho^7*b1+rho^7*b0+
		  rho^8*b1+rho^8*b0+rho^9*b1+rho^9*b0+rho^10*b1+rho^10*b0
	t12<- b0+b1+rho*b0+rho*b1+rho^2*b0+rho^2*b1+rho^3*b0+rho^3*b1+rho^4*b0+
		  rho^5*b1+rho^5*b0+rho^6*b1+rho^6*b0+rho^7*b1+rho^7*b0+
		  rho^8*b1+rho^8*b0+rho^9*b1+rho^9*b0+rho^10*b1+rho^10*b0+
		  rho^11*b1+rho^11*b0
	responses <- c(t1=t1,t2=t2,t3=t3,t4=t4,t5=t5,t6=t6,
	t7=t7,t8=t8,t9=t9,t10=t10,t11=t11,t12=t12)
	cat("Unit Reponse Equilibrium Y=",(b0+b1)/(1-rho),"\n")
	
	return(responses)
}

impulse.resp <- function(x,the.var="geti",lag.var="lag.geti"){
	##Effects of increase in 1 unit in get over the next three periods
	if(class(x)=="data.frame"){tmp<-x
	b0 <-tmp[the.var,"coef"]
	tryb1 <-try(tmp[lag.var,"coef"],silent=T)
	b1<-if(is.na(tryb1)){0}else{tmp[lag.var,"coef"]}
	rho <- tmp["lag.pop","coef"]	
	}else{
	tmp <-  summary(x)
	b0 <-tmp$coefficients[the.var,"Value"]
	tryb1 <-try(tmp$coefficients[lag.var,"Value"],silent=T)
	b1<-if(class(tryb1)=="try-error"){0}else{tmp$coefficients[lag.var,"Value"]}
	rho <- tmp$coefficients["lag.pop","Value"]
	}
	t1 <- b0
	t2 <- b1+rho*b0
	t3 <- rho*b1+rho^2*b0
	t4 <- rho^2*b1+rho^3*b0
	t5 <- rho^3*b1+rho^4*b0
	t6 <- rho^4*b1+rho^5*b0
	t7 <- rho^5*b1+rho^6*b0
	t8 <- rho^6*b1+rho^7*b0
	t9 <- rho^7*b1+rho^8*b0
	t10 <- rho^8*b1+rho^9*b0
	t11 <- rho^9*b1+rho^10*b0
	t12 <- rho^10*b1+rho^11*b0
	responses <- c(t1=t1,t2=t2,t3=t3,t4=t4,t5=t5,t6=t6
					,t7=t7,t8=t8,t9=t9,t10=t10,t11=t11,t12=t12)
		cat("Impulse Reponse Equilibrium Y=",0,"\n")
	return(responses)
}

##### BRAZIL #####################################################
# Start by loading imputed data 
setwd(the.path)
load("DATA/data_popularityBR_AMELIA.RData") 
print(comment(imp))

## Estimates on Amelia imputed sets
## All these models are estimated and averaged across five different imputed sets
## Most basic model possible, good for comparisons with domestic model
im.basic <- zelig(pop~irates+log.commodities+institute, data=imp,model="ls",cite=F) 
	plot.pred(im.basic,crisis=F) 

# Domestic model mentioned in paper
dm.basic <- zelig(pop ~ renda+pib+log.inpc06+unemp,data=imp,model="ls",cite=F)
	plot.pred(dm.basic ,crisis=F) 

# Basic GET model shown in appendix Figure A.8a
im.getbasicregs <- zelig(pop ~ geti, data=imp,model="ls",cite=F) 
	plot.pred(im.getbasicregs,crisis=F) 
	plot.dev(im.getbasicregs,crisis=F) 
	#pdf(file="Figures/fig_getbasic_predicted&actual-2014.pdf")
		plot.pred(im.getbasicregs,crisis=F) 
	#dev.off()

im.getbasicregs2 <- zelig(pop ~ geti+institute, data=imp,model="ls",cite=F) 
	plot.pred(im.getbasicregs,crisis=F) 
		
#lagDV international models (commodit+irates)
im.lagdvregs <- zelig(pop ~ institute +lag.pop 
				+irates+log.commodities, data=imp,model="ls",cite=F) 
	diag.ts(im.lagdvregs)
	plot.dev(im.lagdvregs,plot.range=10,crisis=F)

#GET (reported in the paper)		
im.getlagdvregs <- zelig(pop ~ institute +lag.pop
				+geti, data=imp,model="ls",cite=F) 
	diagnosticsBR <- diag.ts(im.getlagdvregs,to.plot=T)
	plot.pred(im.getlagdvregs,crisis=F)
	plot.dev(im.getlagdvregs,plot.range=10,crisis=F)
	#diagnosis suggests lagdv solves time dependencies, maybe AR-1 not necessary

#using GET lagDV, controlling for GDP, as expected, effects goes away
im.getlagdvgdpregs <- zelig(pop ~ institute +lag.pop
				+geti+pib, data=imp,model="ls",cite=F) 
	diag.ts(im.getlagdvgdpregs,to.plot=T)
			
#AR1 model withOUT LAG-DV
ima.AR1regs <- lapply(imp,function(x){
					arima(x$pop, order=c(1,0,0),
					xreg=subset(x,select=c(#time.pres,#polcrisis,
						instituteDatafolha,
						instituteIbope,instituteNone,instituteSensus,
						instituteVox,irates,
						log.commodities)))})
					
#GET AR-1 MODEL 
ima.getAR1regs <- lapply(imp,function(x){
					arima(x$pop, order=c(1,0,0),
					xreg=subset(x,select=c(#time.pres,#polcrisis,
						instituteDatafolha,
						instituteIbope,instituteNone,instituteSensus,
						instituteVox,
						geti)))})
	#Diagnostics of AR-1 Model withOUT LAG-DV
	diag.ts(ima.getAR1regs)

#GET ARMA11 
ima.getARMA11regs <- lapply(imp,function(x){
					arima(x$pop, order=c(1,0,1),
					xreg=subset(x,select=c(#time.pres,#polcrisis,
						instituteDatafolha,
						instituteIbope,instituteNone,instituteSensus,
						instituteVox,
						geti)))})
	diag.ts(ima.getARMA11regs)

ima.ARMA11regs <- lapply(imp,function(x){
					arima(x$pop, order=c(1,0,1),
					xreg=subset(x,select=c(#time.pres,#polcrisis,
						instituteDatafolha,
						instituteIbope,instituteNone,instituteSensus,
						instituteVox,irates,
						log.commodities)))})
	diag.ts(ima.getARMA11regs)

#GET AR-1 MODEL with  GDP
#No result is expected here, as we're controlling for a consequence of GETI			
ima.getAR1gdpregs <- lapply(imp,function(x){
					arima(x$pop, order=c(1,0,0),
					xreg=subset(x,select=c(#time.pres,#polcrisis,
						instituteDatafolha,
						instituteIbope,instituteNone,instituteSensus,
						instituteVox,
						geti,pib)))})

# ADL, to report in appendix
adlregs <- zelig(pop ~ institute+lag.pop+
				geti+lag.geti, data=imp,model="ls",cite=F) 
		diag.ts(adlregs,to.plot=T)
	adl <- my.amelia(adlregs)
	adl.models.br <- list(adl=adl)
	
## Estimates on "raw" (non-imputed) data (Per Chris Achen's request)
the.data <- imp[[1]] #get one version of the Amelia data sets
the.data$lag.pop <- as.numeric(I(Lag(the.data$pop))) #create lag dv
## back out the non-imputed observations (drop imputed ones)
the.data <- imp[[1]][which(round(imp[[1]]$pop,1)==imp[[1]]$pop),]
im.getlagdvregsRaw <- zelig(pop ~ lag.pop+institute
				+geti, data=the.data,model="ls",cite=F) 
	diag.ts(im.getlagdvregsRaw,to.plot=T)
	plot.pred(im.getlagdvregsRaw,crisis=F)
	plot.dev(im.getlagdvregsRaw,plot.range=10,crisis=F)
im.getolsregsRaw <- zelig(pop ~ institute #+time.pres #+polcrisis 
				+geti,  data=the.data,model="ls",cite=F) 
	diag.ts(im.getolsregsRaw,to.plot=T)
	plot.pred(im.getolsregsRaw,crisis=F)
	plot.dev(im.getolsregsRaw,plot.range=10,crisis=F)

im.getlagdvregsRaw <- zelig(pop ~ lag.pop+institute#+time.pres#+polcrisis
				+geti, data=the.data,model="ls",cite=F) 
	diag.ts(im.getlagdvregsRaw,to.plot=T)
	plot.pred(im.getlagdvregsRaw,crisis=F)
	plot.dev(im.getlagdvregsRaw,plot.range=10,crisis=F)

im.lagdvregsRaw <- zelig(pop ~ lag.pop+institute#+time.pres#+polcrisis
				+irates+log.commodities, data=the.data,model="ls",cite=F) 
					
###### Diagnostics #######
### This is part of table A9
tmp0 <- diag.ts(im.getlagdvregs) 
tmp2 <- diag.ts(ima.getAR1regs)
tmp3 <- diag.ts(ima.getARMA11regs)
xtable(t(data.frame(BrazillagDV=tmp0[,-c(1:7)],
						BrazilAR1=tmp2[,-c(1:7)],
						BrazilARMA11=tmp3[,-c(1:7)]) ))
	 		 
###### Collect Estimates #################
im.getlagdv<-my.amelia(im.getlagdvregs)
		diag.ts(im.getlagdvregs,to.plot=T) #looks good enough
ima.getAR1<-my.amelia(ima.getAR1regs)
		diag.ts(ima.getAR1regs) #diagnostics looks good enough
		plot.dev(ima.getAR1regs,plot.range=45,crisis=F) 
im.getlagdvRaw <- my.amelia(im.getlagdvregsRaw)
ima.getARMA11<-my.amelia(ima.getARMA11regs)
im.lagdv <- my.amelia(im.lagdvregs)
ima.AR1 <- my.amelia(ima.AR1regs)
im.lagdvRaw <- my.amelia(im.lagdvregsRaw )
ima.ARMA11<-my.amelia(ima.ARMA11regs)
ima.getAR1gdp<- my.amelia(ima.getAR1gdpregs)
im.getlagdvgdp<- my.amelia(im.getlagdvgdpregs)

### The GET models, to assemble tables, later, with Mexico
get.models.br <- list(ARMA11=ima.getARMA11,
					AR1=ima.getAR1,
					LAGDV=im.getlagdv,LAGDVraw=im.getlagdvRaw)

gdp.models.br <- list(LAGDV=im.getlagdvgdp)

### Table  of coefficients with GET & International models ######
### THe're combined into a single table, below				  ###
im.table <- merge(merge(
		stack.regs(im.lagdvRaw),stack.regs(im.lagdv),
					by=c('vars','desc'),all=T,
					suffixes=c('LagDV-nonimputed','LagDV')),			
		merge(stack.regs(ima.AR1),stack.regs(ima.ARMA11),
					by=c('vars','desc'),all=T,
					suffixes=c('AR1','ARMA11')),
					by=c('vars','desc'),all=T)

get.table <- merge(merge(
		stack.regs(im.getlagdvRaw),stack.regs(im.getlagdv),
					by=c('vars','desc'),all=T,
					suffixes=c('LagDV-nonimputed','LagDV')),			
		merge(stack.regs(ima.getAR1),stack.regs(ima.getARMA11),
					by=c('vars','desc'),all=T,
					suffixes=c('AR1','ARMA21')),
					by=c('vars','desc'),all=T)
		
# Table of coefficeints for GET models 	
print(get.table,include.rownames=F)

# Table combining both GET and IM in a single table
getim.table <- edited.table(
			merge(get.table,im.table[,-7],
			by=c("vars","desc"),suffixes=c("get",""),all=T))
	the.xtable <- xtable(getim.table)
	caption(the.xtable) <- "International Factors and Presidential Popularity in Brazil"
	label(the.xtable) <- "tab-tstable"
	print(the.xtable,include.rownames=F,caption.placement="top")

im.table.br <- im.table #for comparison with mexico, below

######  MEDIATION analysis using OLS LAG-DV model  #################
if(1==2){#commented out because it takes a while, but works!

med <- lapply(imp,function(x){
	model.m  <- lm(pib ~ geti, data=x[-1,])
	model.y  <- lm(pop ~ geti+pib+institute+time.pres+lag.pop,data=x[-1,])
	model.full<-lm(pop ~ geti+institute+time.pres+lag.pop,data=x[-1,])
	med <- mediate(model.m, model.y, sims=1000, boot=T, 
        treat="geti", mediator="pib", 
        control=NULL, conf.level=.95)
	return(med)})
		
	medshare <- mean(sapply(med,function(x)x$d0.sims/x$tau.sims))
	med.effect <- mean(sapply(med,function(x) x$d0.sims ))
	med.effect.se <-  sd(sapply(med,function(x) x$d0.sims ))
    tot.effect <- mean(sapply(med,function(x) x$tau.sims ))
    tot.effect.se <- sd(sapply(med,function(x) x$tau.sims ))
	medshare <- med.effect/tot.effect

cat("In this setup, the average mediator effect across the five imputed datasets was ",round(med.effect,2)," (SE=",round(med.effect.se,2),") while the average total effects were ",round(tot.effect,2)," (SE=",round(tot.effect.se,2),")\n",sep="")

}

#### The LEIGH competence model (Table A.14 in Appendix) ####################
tmp <- lapply(imp,function(x){lm(pib~lag.pib+geti,data=x)})
tmp.res <- lapply(tmp,residuals)
tmp2 <- lapply(imp,function(x){lm(delta.pib~delta.geti,data=x)})
tmp2.res <- lapply(tmp2,residuals)
na.end <- nrow(imp[[1]])-length(tmp.res[[1]])
na.end2 <- nrow(imp[[1]])-length(tmp2.res[[1]])
for(i in 1:length(tmp)){imp[[i]]$merit<-c(tmp.res[[i]],rep(NA,na.end))
						imp[[i]]$meritleigh<-c(tmp2.res[[i]],rep(NA,na.end2))}
#Applying the merit measures in the best time series specification
meritluck_base <- zelig(pop ~ institute +lag.pop
				+geti, data=imp,model="ls",cite=F) 
meritluck <- zelig(pop ~ institute +lag.pop
				+geti+merit, data=imp,model="ls",cite=F) 
meritluck2 <- zelig(pop ~ institute +lag.pop
				+geti+meritleigh, data=imp,model="ls",cite=F) 
#Table in appendix
 merit.table <- merge(stack.regs(my.amelia(meritluck)),stack.regs(my.amelia(meritluck2)),
					by=c('vars','desc'),all=T,
					suffixes=c('lag','delta'))
	merit.table$vars<-ifelse(merit.table$desc=="coef",
						merit.table$vars,
						as.character(merit.table$desc))
	the.merit.table <- xtable(merit.table[,-2])
	caption(the.merit.table) <- "Merit v. Luck and Presidential Popularity in Brazil"
	label(the.merit.table) <- "tab-tsmerit"
	print(the.merit.table,include.rownames=F,caption.placement="top")
	



##### MEXICO #####################################################
setwd(the.path)
load("DATA/data_popularityMX_amelia.RData")
print(comment(imp))

# Expanded ols international model, included in collection later in the code 
im.olsregs <- zelig(pop~irates+log.commodities+
			institute+time.pres, data=imp,model="ls",cite=F) 

dm.olsregs <- zelig(pop ~ institute+time.pres+
					renda+pib+log.inpc06+unemp,data=imp,model="ls",cite=F)
	plot.pred(dm.olsregs,crisis=F) 
	
dm.olsregs <- zelig(pop ~ 
					renda+pib+log.inpc06+unemp,data=imp,model="ls",cite=F)
	plot.pred(dm.olsregs,crisis=F) 

# Basic GET model shown in main paper
im.getbasicregs1 <- zelig(pop ~ geti, data=imp,model="ls",cite=F) 
	#pdf(file="Figures/fig_MXgetbasic_predicted&actual.pdf")
	plot.pred(im.getbasicregs1,crisis=F) 
	#dev.off()

im.getbasicregs2 <- zelig(pop ~ institute   
				+geti, data=imp,model="ls",cite=F) 
	plot.pred(im.getbasicregs2,crisis=F) 
	plot.dev(im.getbasicregs2,crisis=F) 

	
# Complete GET OLS model (includes time), reported later in the code
im.getolsregs <- zelig(pop ~ institute +time.pres  
				+geti, data=imp,model="ls",cite=F) 
	plot.pred(im.getolsregs,crisis=F) 

im.lagdvregs <- zelig(pop ~ institute+time.pres+lag.pop
				+irates+log.commodities, data=imp,model="ls",cite=F) 
	diag.ts(im.lagdvregs)
	plot.dev(im.lagdvregs,plot.range=10,crisis=F)

#GET LAG-DV			
im.getlagdvregs <- zelig(pop ~ institute +lag.pop
				+geti, data=imp,model="ls",cite=F) 
	diag.ts(im.getlagdvregs,to.plot=T)
	plot.pred(im.getlagdvregs,crisis=F)
	plot.dev(im.getlagdvregs,plot.range=10,crisis=F)
	
		
#AR1 model withOUT LAG-DV
ima.AR1regs <- lapply(imp,function(x){
					arima(x$pop, order=c(1,0,0),
					xreg=subset(x,select=c(time.pres,
						instituteConsulta.Mitofsky,
						instituteGrupo.Reforma,instituteNone,instituteOPRM,
						log.commodities,irates)))})
						diag.ts(ima.AR1regs)
#The ideal Commodity,Irates model				
ima.ARMA11regs <- lapply(imp,function(x){#from diagnostics of underlying model
					arima(x$pop, order=c(1,0,1),
					xreg=subset(x,select=c(time.pres,
						instituteConsulta.Mitofsky,
						instituteGrupo.Reforma,instituteNone,instituteOPRM,
						log.commodities,irates)))})
diag.ts(ima.ARMA11regs)

#Domestic AR1 model withOUT LAG-DV				
dma.AR1lagdvregs <- lapply(imp,function(x){
					arima(x$pop, order=c(1,0,0),
					xreg=subset(x,select=c(time.pres,
						instituteConsulta.Mitofsky,
						instituteGrupo.Reforma,instituteNone,instituteOPRM,
						renda,delta.pib06,inpc06,unemp)))})			

#GET AR-1 MODEL withOUT LAG-DV	(per Chris Achen's request)
ima.getAR1regs <- lapply(imp,function(x){
					arima(x$pop, order=c(1,0,0),
					xreg=subset(x,select=c(
						instituteConsulta.Mitofsky,
						instituteGrupo.Reforma,instituteNone,instituteOPRM,
						geti)))})
	#Diagnostics of AR-1 Model withOUT LAG-DV
	diag.ts(ima.getAR1regs)

# The "IDeal" GET Mexico Model
ima.getARMA11regs <- lapply(imp,function(x){#from diagnostics of underlying model
					arima(x$pop, order=c(1,0,1),
					xreg=subset(x,select=c(
						instituteConsulta.Mitofsky,
						instituteGrupo.Reforma,instituteNone,instituteOPRM,
						geti)))})
diagnosticsMX<- diag.ts(ima.getARMA11regs)


### Raw data (no imputation)
the.data <- imp[[1]] #get on verion of the Amelia data sets
the.data$lag.pop <- as.numeric(I(Lag(the.data$pop))) #create lag dv
## back out the non-imputed observations (drop imputed ones)
the.data <- imp[[1]][which(round(imp[[1]]$pop,1)==imp[[1]]$pop),]

im.getolsregsRaw <- zelig(pop ~ institute   
				+geti,  data=the.data,model="ls",cite=F) 
	diag.ts(im.getolsregsRaw,to.plot=T)
	plot.pred(im.getolsregsRaw,crisis=F)
	plot.dev(im.getolsregsRaw,plot.range=10,crisis=F)

im.getlagdvregsRaw <- zelig(pop ~ lag.pop+institute 
				+geti, data=the.data,model="ls",cite=F) 
	diag.ts(im.getlagdvregsRaw,to.plot=T)
	plot.pred(im.getlagdvregsRaw,crisis=F)
	plot.dev(im.getlagdvregsRaw,plot.range=10,crisis=F)

im.lagdvregsRaw <- zelig(pop ~ lag.pop+institute 
				+irates+log.commodities, data=the.data,model="ls",cite=F) 

###### Collect results 
im.getols<-my.amelia(im.getolsregs)
		diag.ts(im.getolsregs,to.plot=T) 
im.getlagdv<-my.amelia(im.getlagdvregs)
		diag.ts(im.getlagdvregs,to.plot=T) 
ima.getAR1<-my.amelia(ima.getAR1regs)
		plot.dev(ima.getAR1regs,plot.range=45,crisis=F) 
im.getlagdvRaw <- my.amelia(im.getlagdvregsRaw)
im.getolsregsRaw <- my.amelia(im.getolsregsRaw)
ima.getARMA11 <- my.amelia(ima.getARMA11regs)
im.ols <- my.amelia(im.olsregs)
im.lagdv <- my.amelia(im.lagdvregs)
ima.AR1 <- my.amelia(ima.AR1regs)
im.lagdvRaw <- my.amelia(im.lagdvregsRaw)
ima.ARMA11 <- my.amelia(ima.ARMA11regs)

###### Diagnostics  IMPORANT
# Diagnostics of the basic international model
diag.ts(im.olsregs)
#looks like the best model is AR=(01),D=(01),MA=1

# Diagnostics on basic model:
diag.ts(im.getolsregs)
#looks like optimal model is AR=0,D=1,MA=1

# Diagnostics for lag-dv (included in appendix)
diag.ts(im.getlagdvregs)
tmp0 <- diag.ts(im.getlagdvregs) 
tmp1 <- diag.ts(ima.getAR1regs)
tmp2 <- diag.ts(ima.getARMA11regs)
# Diagnostics table included in the appendix
# This is part of table A9 in the appendix
	xtable(t(data.frame(MexicoLagDV=tmp0[,-c(1:7)],
						MexicoAR1=tmp1[,-(1:7)],
						MexicoARMA11=tmp2[,-(1:7)])) )

### GET models, for comparison later with Brazil
get.models.mx <- list(ARMA11mx=ima.getARMA11,AR1mx=ima.getAR1,
					LAGDVmx=im.getlagdv,LAGDVmxraw=im.getlagdvRaw)
					
					
### Table  of coefficients with GET & International models ######
### Used in the appendix of the paper (Table A11)			  ###
im.table <- merge(merge(
		stack.regs(im.lagdvRaw),stack.regs(im.lagdv),
					by=c('vars','desc'),all=T,
					suffixes=c('LagDV-nonimputed','LagDV')),			
		merge(stack.regs(ima.AR1),stack.regs(ima.ARMA11),
					by=c('vars','desc'),all=T,
					suffixes=c('AR1','ARMA11')),
					by=c('vars','desc'),all=T)
get.table <- merge(merge(
		stack.regs(im.getlagdvRaw),stack.regs(im.getlagdv),
					by=c('vars','desc'),all=T,
					suffixes=c('LagDV-nonimputed','LagDV')),			
		merge(stack.regs(ima.getAR1),stack.regs(ima.getARMA11),
					by=c('vars','desc'),all=T,
					suffixes=c('AR1','ARMA11')),
					by=c('vars','desc'),all=T)
		
# Table combining both GET and IM in a single table
getim.table <- edited.table(
			merge(get.table,im.table[,-7],
			by=c("vars","desc"),suffixes=c("get",""),all=T))
	the.xtable <- xtable(getim.table)
	caption(the.xtable) <- "International Factors and Presidential Popularity in Mexico"
	label(the.xtable) <- "tab-tstableMX"
	print(the.xtable,include.rownames=F,caption.placement="top")

im.table.mx <- im.table #for comparison with Brazil, below


###### MEXICO & BRAZIL combined, response functions, etc ########

## Compare GET coefficients on ARMA models GRAPHICALLY (not in paper)
getbr <- get.models.br
getmx <- get.models.mx
get.models <- list(getbr[[1]],getmx[[1]])
	par(mfrow=c(1,1),oma=c(2,0,0,0))
	tmp <- plot.coefs(get.models,"geti","",shade=T,leftmar=8)
	axis(side=2,at=c(-.5,.5),labels=c("Brazil","Mexico"),las=2,cex.axis=2,tick=F)
	mtext(side=1,line=2.3,text="Coefficient on GET (w/95% Conf. Interval)")

## Table to with "best models" for each country
## Table 3 in the main body of the paper is a short version of this
## Table A7 in the appendix is the complete version of this
getbr <- get.models.br
getmx <- get.models.mx
br <- clean.regs(stack.regs(getbr$LAGDV))
mx <- clean.regs(stack.regs(getmx$ARMA11))
the.table<-merge(br,mx,by=c("vars","desc"),suffixes=c("BR","MX"),all=T)
the.table<-the.table[c(22:24,4:6,7:21,25:30,1:3,31),]
tests <- data.frame(vars=c(rep("Augmented Dickey-Fuller",3),
					rep("Box–Pierce",2),
				  rep("Breusch–Godfrey",2)),
		   desc=c("coef","t","n","coef","p-value","coef","p-value"),
		   vBR=diagnosticsBR[c(5:7,10:13)],
		   vMX=diagnosticsMX[c(5:7,10:13)])
the.table<- rbind(the.table,tests)
rownames(the.table)<- 1:nrow(the.table)
the.table$vars <- ifelse(the.table$desc=="coef",
					the.table$vars,as.character(the.table$desc))
the.table$desc <- NULL
the.xtable <- xtable(the.table)
caption(the.xtable) <- "Time Series Analysis for the Effect of GET on Popularity"
label(the.xtable) <- "tab-tsresults"
print(the.xtable,hline=30,include.rownames=F)


## UNIT and IMPULSE GRAPHS
## For LAGDV and ADL models for Brazil
## in Appendix

getbr <- get.models.br
the.model <- getbr$LAGDV

#pdf(file="FIGURES/fig-unitimpulseresponse.pdf")
par(mfrow=c(1,1),mar=c(4,4,2,0.5))
max.x <- 12
plot(x=c(1,max.x),y=c(0,max(c(unit.resp(the.model),impulse.resp(the.model))))
		,type="n"
		,xlab="Months"
		,ylab=expression(paste(Delta," Popularity")))
lines(unit.resp(the.model),lty=1)
lines(impulse.resp(the.model),lty=2)
text(max.x-3,unit.resp(the.model)[max.x-2],pos=3,label="Unit Response")
text(max.x-3,impulse.resp(the.model)[max.x-4],pos=3,label="Impulse Response")
#dev.off()

the.model <- adl.models.br$adl
#pdf(file="FIGURES/fig-unitimpulseresponseADL.pdf")
par(mfrow=c(1,1),mar=c(4,4,2,0.5))
max.x <- 12
plot(x=c(1,max.x),y=c(0,max(c(unit.resp(the.model),impulse.resp(the.model))))
		,type="n"
		,xlab="Months"
		,ylab=expression(paste(Delta," Popularity")))
lines(unit.resp(the.model),lty=1)
lines(impulse.resp(the.model),lty=2)
text(max.x-3,unit.resp(the.model)[max.x-2],pos=3,label="Unit Response")
text(max.x-3,impulse.resp(the.model)[max.x-4],pos=3,label="Impulse Response")
#dev.off()


### Table with additional specifications, including controlling for GDP
getbr <- get.models.br
brldvraw <- clean.regs(stack.regs(getbr$LAGDVraw))
bradl <- clean.regs(stack.regs(adl.models.br$adl))
brarma11 <- clean.regs(stack.regs(getbr$ARMA11))
brgdp <- clean.regs(stack.regs(gdp.models.br$LAGDV))
the.table01<-merge(brldvraw,bradl,by=c("vars","desc"),suffixes=c("Non Imputed","ADL"),all=T)
the.table02<-merge(brarma11,brgdp,by=c("vars","desc"),suffixes=c("ARMA(1,1)","GDP"),all=T)
the.table <- merge(the.table01,the.table02,by=c("vars","desc"),all=T)
the.table<-the.table[c(22:24,4:6,34:36,7:21,25:30,1:3,31:33,37),]
rownames(the.table)<- 1:nrow(the.table)
the.table$vars <- ifelse(the.table$desc=="coef",
					the.table$vars,as.character(the.table$desc))
the.table$desc <- NULL
the.xtable <- xtable(the.table)
caption(the.xtable) <- "Additional Specification of Time Series Analysis for Brazil"
label(the.xtable) <- "tab-tsresultsadditional"
print(the.xtable,hline=30,include.rownames=F)


### Table A10 in apendix results with IR and COMMODITIES in Bazil and Mexico

the.table <- merge(im.table.br[-grep("institute",im.table.br$vars),] #remove pollster dummies
		,im.table.mx[-grep("institute",im.table.br$vars),]
		,by=c("vars","desc"),suffixes=c(".br",".mx"),all=T)
the.table <- the.table[c(7:9,13:15,22:24,10:12,1:3,16:21),]
rownames(the.table)<- 1:nrow(the.table)
the.table$vars <- ifelse(the.table$desc=="coef",
					the.table$vars,as.character(the.table$desc))
the.table$desc <- NULL
the.xtable <- xtable(the.table)
caption(the.xtable) <- "Specifications Including Interest Rates and Commodity Prices"
label(the.xtable) <- "tab-tsresultsinternational"
print(the.xtable,include.rownames=F, caption.placement = "top")
