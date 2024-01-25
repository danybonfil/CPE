### Examine the effect of GET on probability of reelection  	 ####
### This code replicates results for the analysis of elections	 ####
### Creates all tables included in main body of paper and in the ####
### appendix. Requires file "data-elections.RData", with 		 ####
### country/election observations. Starts by declaring some 	 ####
### used later in the analysis functions						 ####

fd.badgood <- function(xx,the.data=NULL,cl=T){#reports three first differences....
	#xx is a regression object
	if(is.null(the.data)){the.data<-xx$data}
	base.country <<- the.data$country[1]
	inc.ran <<- inc.allowed <<- 0	   
	#define quantities
	lscebad<-setx(xx,geti12=simgeti[1],lsce=1,country=base.country )
	lscegood<-setx(xx,geti12=simgeti[2],lsce=1,country=base.country)	
	cmprsnbad<-setx(xx,geti12=simgeti[1],lsce=0,country=base.country)
	cmprsngood<-setx(xx,geti12=simgeti[2],lsce=0,country=base.country)	
	#bad to good (for tables)
	if(cl==T){#using my function, using clustered VCOV matrix
		simlscebg <- my.sims(xx,x=lscebad,x1=lscegood)
		simcmprsnbg <- my.sims(xx,x=cmprsnbad,x1=cmprsngood)	
		badgood <- rbind(lsce=simlscebg,comparison=simcmprsnbg)
	}
	if(cl==F){#using Zelig's sim, which uses regular VCOV matrix
		to.select <- grep("First Differences"
					,names(summary(sim(xx,x=lscebad,x1=lscegood))$stats))
		simlscebg <- summary(sim(xx,x=lscebad,x1=lscegood))$stats[[to.select]]
		simcmprsnbg <- summary(sim(xx,x=cmprsnbad,x1=cmprsngood))$stats[[to.select]]
		badgood <- rbind(lsce=simlscebg,comparison=simcmprsnbg)
		rownames(badgood)<-c("lsce","comparison")
	}
	return(badgood)
	rm(inc.ran,inc.allowed,envir=.GlobalEnv) 
	}

fd.badgood.alt <- function(xx,the.subset=NULL,clv=T){
	#for use with regression on the LSCE only
	#xx is a regression object
	#the subset is the data
	if(is.null(the.subset)){the.subset<-xx$data}
	base.country <<- the.subset$country[1]
	inc.ran <<- inc.allowed <<- 0
	if(length(grep("geti",xx$formula))==0){#formula does not contain GET
	#create objects in gloabl envir due to Zeli 4 quirks
	simirates <<- c(mean(the.subset$irates12,na.rm=T)-
					sd(the.subset$irates12,na.rm=T),
				   mean(the.subset$irates12,na.rm=T)+
				   	sd(the.subset$irates12,na.rm=T))
	simcomm <<- c(mean(the.subset$comm12,na.rm=T)-
					sd(the.subset$comm12,na.rm=T),
				   mean(the.subset$comm12,na.rm=T)+
				   	sd(the.subset$comm12,na.rm=T))
	#Bad is high irates and low commodities, Good is the opposite	   
	bad <- setx(xx,irates12=simirates[2],comm12=simcomm[1],country=base.country)
	good <- setx(xx,irates12=simirates[1],comm12=simcomm[2],country=base.country)
	}else{#formula contains GET
	simgeti <<- c(mean(the.data$geti12,na.rm=T)-
					sd(the.data$geti12,na.rm=T),
				   mean(the.data$geti12,na.rm=T)+
				   	sd(the.data$geti12,na.rm=T)) 		
	bad<-setx(xx,geti12=simgeti[1],country=base.country)
	good<-setx(xx,geti12=simgeti[2],country=base.country)	
	}
	if(clv==F){#using Zelig's sim, which uses regular VCOV matrix
				sims <- sim(xx,x=bad,x1=good)	
				tfd <- grep("First",names(sims$stats))#which position is FD
				bg.fd <- sims$stats[[tfd]]}
	if(clv==T){#using my custom my.sims, which uses cl SE VCOV matrix
				bg.fd <- my.sims(xx,x=bad,x1=good)
			}
	out <- bg.fd	
	return(out)
	rm(base.country,inc.ran,inc.allowed,envir=.GlobalEnv) 
	}			

my.sims <- function(xx,x,x1){### performs a one time by hand using cl SE vcov
	require(mvtnorm)
	the.sims <- list()
	the.vcov <- cl(xx$data,xx,"country")$vcov ## use clustered SE
	the.coefs <- as.numeric(na.omit(coef(xx)))#Na omit is for excluded countris
	sims <- rmvnorm(n=1000 #boostrap
				,mean=the.coefs
				,sigma=the.vcov) 				 
	the.draws <- plogis(sims%*%t(x1[[4]]))-plogis(sims%*%t(x[[4]]))
	ci <- quantile(the.draws,prob=c(0.25,0.05,0.975))
	names(ci) <- c('X50.','X2.5.','X97.5.')
	pe <- mean(the.draws)
	sd <- sd(the.draws)
	out <- c(mean=pe,sd=sd,ci)
	return(out)
	}
	
cl <- function(dats,fm, cls){
				##dat is the data (including the clustering variable)
				##fm is the regression
				##cl is the string name of the clustering variable
				##This function computes simple "clustered" standard errors
			   #attach(dats, warn.conflicts = F)
			   cluster <- factor(dats[,cls])
			   ## Zelig objects don~t take estfun, so estimate glm
			   if(is.element("zelig",class(fm))){
			   		the.link <- class(fm)[[2]]
			    	fm <- glm(fm$formula,data=dats,
			    			   family = binomial(link=the.link))
			    	cluster <- factor(fm$data[,cls])}
	           require(sandwich)
               require(lmtest)
	           M <- length(unique(cluster))
	           N <- length(cluster)
	           K <- fm$rank
	           dfc <- (M/(M-1))*((N-1)/(N-K))
	           uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
	           vcovCL <- dfc*sandwich(fm, meat=crossprod(uj)/N)
	           cls <- coeftest(fm, vcovCL) 
             #detach(dats)
             return(list(estimates=cls,vcov=vcovCL))}

reduction.in.error<-function(x,print.all=F){
			#computes reudction in error in a logit/probit
	    	#x is a regression object
	    	if(class(x)[1]=="mer"){#RE MODEL
	    	error <- min(prop.table(table(x@y)))
	        correct <- sum(diag(prop.table(table(x@y,fitted(x)>0.5))))
	        new.error <- 1-sum(diag(prop.table(table(x@y,fitted(x)>0.5))))
	    	print.all <- F}
	    	if(is.element("zelig",class(x))){#Zelig
	    	 if(is.element("logit.mixed",class(x))){#Zelig 4 RE model
	    	 error <- min(prop.table(table(x$result@resp$y)))
	         correct <- sum(diag(prop.table(table(x$result@resp$y,
	         			fitted(x$result)>0.5))))
	         new.error <- 1-sum(diag(prop.table(table(x$result@resp$y,
	         			fitted(x$result)>0.5))))
	    	 print.all <- F
	    	}else{#all other zeligs
	    	 error <- min(prop.table(table(x$result$y)))
	    	 correct <- sum(diag(prop.table(table(x$result$y,
	    	 				fitted(x$result)>0.5))))
			 new.error <- 1-sum(diag(prop.table(table(x$result$y,
			 				fitted(x$result)>0.5))))
	    	}#end else
	    	}#end if zelig
	    	if(is.element("gml",class(x))){		  
	    	error <- min(prop.table(table(x$y)))
	    	correct <- sum(diag(prop.table(table(x$y,fitted(x)>0.5))))
			new.error <- 1-sum(diag(prop.table(table(x$y,fitted(x)>0.5))))
			}#end if/else
		    if(print.all){
	    	detail <- data.frame(year=x$data$year,
	    					   elec=x$data$elec,
	    						actual=x$y,
	    					  predicted=as.numeric(fitted(x)>0.5),
	    					  prob=round(fitted(x),2))
	    					  }else{detail=NULL}
		reduction = 1-(new.error/error)
		out<-list(pre=c(baseline.error=error,
					model.error=new.error,
					pre=reduction),detail=detail)
		return(out)
}

stack.regs <- function(x,data=NULL,clse=T){
	if(class(x)[1]=='mer'|class(x)[1]=='glmerMod'){##RE traditional model
			if(is.null(data)){data<-elecs}
			xx<-data.frame(as.table(summary(x)@coefs)) 
			yy<-data.frame(Var1=c("zN","zCountries","zBaseline Error",
									"zModel Error","zPRE"),
				Var2=rep("Stat",5),
				Freq=c(nrow(x@frame),length(unique(x@frame$country)),
				      reduction.in.error(x)$pre ))
			}
	if(is.element("logit.mixed",class(x))){#RE Zelig 4
			if(is.null(data)){data<-x$data}
			xx<-data.frame(as.table(summary(x)$coef)) 
			yy<-data.frame(Var1=c("zN","zCountries","zBaseline Error",
									"zModel Error","zPRE"),
				Var2=rep("Stat",5),
				Freq=c(nrow(x$result@frame),length(unique(x$result@frame$country)),
				      reduction.in.error(x)$pre ))
			}else{#all others
			if(is.null(data)){data<-x$data}
			if(clse==T){#cluster SE's
			cls <- cl(data,x,"country")$estimates
			xx<-data.frame(as.table(cls))
			}else{#not cluster SE's
			xx<-data.frame(as.table(summary(x)$coef))
			}
			yy<-data.frame(Var1=c("zN","zCountries","zBaseline Error",
									"zModel Error","zPRE"),
						Var2=rep("Stat",5),
						Freq=c(nrow(x$data),length(unique(x$data$country)),
				     reduction.in.error(x)$pre ))
		}
	xx<-xx[-which(xx$Var2=='z value'),]
    xx$Var2 <- car::recode(xx$Var2,"'Pr(>|z|)'='zp-value'")
	xx<-xx[order(xx$Var1),]
		out<-rbind(xx,yy)
		names(out)[3] <- "m."	
	return(out)
}
   
   
   
####### ANALYSIS STARTS HERE ###############################   
   
#### Report Sample Characteristics 
library(Zelig)
library(ZeligMultilevel)
library(xtable)
setwd("~/Dropbox/Data/Paper-VotoEconomico")
load("DATA/data_elections.RData")
print(comment(elecs))
LSCE <- c("arg","bol","bra","chi","col","ecu","nic","per","uru","ven")

print(comment(elecs))
tmp <- subset(elecs,free==T)
cat("Incumbent presidents could have ran for office in only",sum(tmp$inc.allowed),"of the",nrow(tmp),"elections deemed free and fair. They effectively chose to run in",table(tmp$inc.allowed,tmp$inc.ran)[2,2],"cases, having won reelection",table(ran=tmp$inc.ran,won=tmp$reelecb)[2,2],"times.\n")

cat("As the list shows, the rules governing immediate reelection varied across free and fair elections in the sample in",sum(table(tmp$reelec.rule.change==0,tmp$country)[1,]!=0),"countries. The great majority of changes were in the direction of allowing immediate reelection where it had previously not been allowed, but in Peru, Dominican Republic, and Nicaragua there were also changes in the other direction.")

# Typeset description of the sample for  the paper
cat("In order to evaluate the association between GET and the electoral success of presidents we identified a total of ",nrow(elecs)," presidential elections in ",length(unique(elecs$country[elecs$free==T]))," Latin American countries between 1980 and 2012, of which ",sum(elecs$free)," were deemed free of  electoral process and/or franchise violations cite{Main+2010}.  Of these, ",table(elecs$free,elecs$lsce)[2,2]," elections were held in the ",length(LSCE)," countries in the LSCE sample and ",table(elecs$free,elecs$lsce)[2,1]," elections were held in the ",sum(table(elecs$country,elecs$lsce)[,2]==0)," countries in the comparison group. footnote{There are actually 74 elections in the LSCE sample in the period. However, one election (Guatemala 2012)  is coded as missing on the reelection variable due to particularities of the case, and is not included in the analysis. See supplemental materials for details.}",sep="")


#### Create Figure 4 in the paper ####
the.data <- subset(elecs,free==T&is.na(reelec)==F)
the.sample <- subset(elecs,lsce==T)
the.samplenot <- subset(elecs,lsce==F)
reelec<-prop.table(table(the.sample$reelecb,the.sample$decade),2)[c(2,1),]
reelecnot<-prop.table(table(the.samplenot$reelecb,the.samplenot$decade),2)[c(2,1),]
reelecall <- prop.table(table(elecs$reelecb,elecs$decade),2)[c(2,1),]

pdf("Figures/fig_reelection.pdf", width = 10.00, height = 10.00)
par(mar=c(2,6,2,2))
to.plot <- t(rbind('LSCE Countries'=reelec[1,],'Comparison Group'=reelecnot[1,]))
barplot(to.plot,ylim=c(0,0.8),beside=T,legend.text = F,cex=1.5,cex.axis=1.5)
my.colors <- gray.colors(3, start = 0.3, end = 0.9, gamma = 2.2)
mtext("Reelection Rates",side=2,line=3,cex=1.5)
legend("topright",legend=c('1980s','1990s','2000s'),fill=my.colors,bty="n",cex=1.5)
dev.off()

#Typeset paragraph for latex
cat("Figure ref{figreelections}  shows that reelection rates increased markedly from ",round(reelec[1,1]*100,1),"% in the worse period, to ",round(reelec[1,3]*100,1),"% in the best period in the LSCE countries. In contrast, these rates did not follow any noticeable trend and hovered around ",round(mean(reelecnot[1,])*100),"% in the comparison group.",sep="") 


### Replicate Table 2 in the paper, as well as supplemental materials 	   ###
### Regression tables, first differences, confidence intervals for effects ###
the.data <- subset(elecs,free==T&is.na(reelec)==F)

##### a bug in Zelig::setx requires control variables be numeric 
##### and character variables be factors, otherwise setx does not
##### work inside a function
the.data$lsce <- as.numeric(the.data$lsce)
the.data$inc.ran <- as.numeric(the.data$inc.ran)
the.data$country <- factor(the.data$country)

##### Declare values for "badgood" simulations, better done globally 
##### to ensure that all simulations use the same values
simgeti <- c(mean(the.data[,"geti12"],na.rm=T)-
					sd(the.data[,"geti12"],na.rm=T),
				   mean(the.data[,"geti12"],na.rm=T)+
				   	sd(the.data[,"geti12"],na.rm=T)) 

##### Cross sectional models, as reported in the main body of the paper ####
##### There are all models  in which GET is interacted with subsample   ####
##### The use of clustered standard errors for quantities of interest 	####
##### requires some handmade funcitons, as Zelig does not allow use 	####
##### defined variance-covariance matrices								####
   			   	
#Simplest regression 
regg <- zelig(reelecb~geti12*lsce, model="logit",
	data=the.data, cite=F) 
	reggcl <- cl(regg$data,regg,"country")$estimates #clustered SE
fd.g <- fd.badgood(regg) #first differences (bad get -> good get)

#Fixed effects (it is necessary to omit the group dummy)
reggfe <- zelig(reelecb~geti12+geti12:lsce+country-1, model="logit",
	data=the.data, cite=F) 
fd.gfe <- fd.badgood(reggfe,cl=F) #first differences. no clustering

#random effects
reggre <- zelig(formula= reelecb~geti12+geti12*lsce+ tag(1| country), 
				data=the.data, model="logit.mixed", cite=F)
fd.gre <- fd.badgood(reggre,cl=F) 

#with ideology 
reggi <- zelig(reelecb~geti12*lsce+d.gov, model="logit",
	data=the.data, cite=F)
	reggicl <- cl(reggi$data,reggi,"country")$estimates
fd.gi <- fd.badgood(reggi)

#with "quality" (polrisk)
the.data.pr <- subset(the.data,is.na(polrisk12)==F) #deal with NA's
reggpr <- zelig(reelecb~geti12*lsce+polrisk12, model="logit",
	data=the.data.pr, cite=F)
	reggprcl <- cl(the.data.pr,
	reggpr,"country")$estimates
fd.gpr <- fd.badgood(reggpr)

#incumbent ran
regginc1 <- zelig(reelecb~geti12*lsce+inc.ran, model="logit",
	data=the.data,cite=F)
fd.ginc1 <- fd.badgood(regginc1)
	
### Produce main table of coefficients in paper ####
reg.table.get <- merge(
			 merge(
			 merge(stack.regs(regg),stack.regs(reggfe,clse=F),
				   by=c("Var1","Var2"),
				   suffixes=c("cl","fe"),all.x=T),
			 merge(stack.regs(reggre,clse=F),stack.regs(regginc1),
			      by=c("Var1","Var2"),
			      suffixes=c("re","clreelect"),all=T),			      ,
			      by=c("Var1","Var2"),all=T),
			 stack.regs(reggi),by=c("Var1","Var2"),all=T)			      			      		     			
reg.table.get$Var1 <- ifelse(reg.table.get$Var2=="Estimate"|reg.table.get$Var2=="Stat",
				as.character(reg.table.get$Var1),as.character(reg.table.get$Var2))
print(xtable(reg.table.get[,-2],digits=3),include.rownames=F)

# add rows to the table with first differences and confidence intervals
fd.lsce <- rbind('Mod. 1'=fd.g[1,c(1,4,5)],
				'Mod. 2'=fd.gfe[1,c(1,4,5)],
				'Mod. 3'=fd.gre[1,c(1,4,5)],
				'Mod. 4'=fd.ginc1[1,c(1,4,5)],
				'Mod. 5'=fd.gi[1,c(1,4,5)])

fd.comp <- rbind('Mod. 1'=fd.g[2,c(1,4,5)],
				'Mod. 2'=fd.gfe[2,c(1,4,5)],
				'Mod. 3'=fd.gre[2,c(1,4,5)],
				'Mod. 4'=fd.ginc1[2,c(1,4,5)],
				'Mod. 5'=fd.gi[2,c(1,4,5)])
								
xtable(round(rbind(t(fd.lsce),t(fd.comp)),2))

# Typeset paragraph with results for paper
cat("In the LSCE sample, these  effects range from  ",round(fd.lsce["Mod. 4","mean"],2)," in the model that controls for whether the incumbent ran, to ",round(fd.lsce["Mod. 2","mean"],2)," in the model with country fixed-effects.  In contrast, effects in the not determined sample are never larger than ",round(fd.comp["Mod. 5","mean"],2)," and not statistically significant.",sep="") 


#### What follows is currently in the Appendix! ##############

### Variations on running for reelection and political risk (in Appendix)
#incumbent ran
regginc1 <- zelig(reelecb~geti12*lsce+inc.ran, model="logit",
	data=the.data,cite=F)
fd.ginc1 <- fd.badgood(regginc1)
	
#incumbent allowed to run
regginc2 <- zelig(reelecb~geti12*lsce+inc.allowed, model="logit",
	data=the.data,robust=T, cite=F)
fd.ginc2 <- fd.badgood(regginc2)

#incumbent RUNNING, excluding rule changes
regginc3 <- zelig(reelecb~geti12*lsce+inc.ran, model="logit",
	data=subset(the.data,reelec.rule.change==0),robust=T, cite=F)
fd.ginc3 <- fd.badgood(regginc3)

#political risk 
the.data.pr <- subset(the.data,is.na(polrisk12)==F) #deal with NA's
reggpr <- zelig(reelecb~geti12*lsce+polrisk12, model="logit",
	data=the.data.pr,cite=F)
fd.pr <- fd.badgood(reggpr,the.data.pr)

fd.lsceinc <- rbind('Mod. Ran'=fd.ginc1[1,c(1,4,5)],
				'Mod. No Rule Change'=fd.ginc3[1,c(1,4,5)],
				'Mod. Allowed'=fd.ginc2[1,c(1,4,5)],
				'Mol. Polrisk'=fd.pr[1,c(1,4,5)])

fd.compinc <- rbind('Mod. Ran'=fd.ginc1[2,c(1,4,5)],
				'Mod. No Rule Change'=fd.ginc3[2,c(1,4,5)],
				'Mod. Allowed'=fd.ginc2[2,c(1,4,5)],
				'Mol. Polrisk'=fd.pr[2,c(1,4,5)])
				
reg.table.r <-  merge(
			 merge(stack.regs(regginc1),stack.regs(regginc3),
				   by=c("Var1","Var2"),
				   suffixes=c("incran","norulechange"),all=T),
			 merge(stack.regs(regginc2),stack.regs(reggpr),
			 	  by=c("Var1","Var2"),
				   suffixes=c("allowed","polrisk"),all=T),
			      by=c("Var1","Var2"),all=T)		
reg.table.r$Var1 <- ifelse(reg.table.r$Var2=="Estimate"|reg.table.r$Var2=="Stat",
				as.character(reg.table.r$Var1),as.character(reg.table.r$Var2))
print(xtable(reg.table.r[,-2],digits=3),include.rownames=F)
# Create lines to report FD and confidence intervals in table (JOP version, with not graph)
xtable(round(t(fd.lsceinc),2))
xtable(round(t(fd.compinc),2))

cat("Results with the three different operationalizations of personal reelection yield very similar results. Although the magnitude of the coefficient on GET is smaller than in models that do not control for personal reelection, a change from bad to good levels of GET (as defined earlier in the text) still amounts to an increase in the probability of reelection of at least ",round(min(fd.lsceinc[,1]),2)," in the LSCE group.\n",sep="")	

#### Estimate models using IRates and Commodity Prices (in Appendix)
#### To avoid multiple interaction, use only LSCE
the.subset <-subset(the.data,lsce==T) 
the.subset$country <- factor(the.subset$country)

#interest rates only
regi <- zelig(reelecb~irates12, model="logit",
	data=the.subset, cite=F)
	regicl <- cl(the.subset,
	regi,"country")$estimates
fd.i <- fd.badgood.alt(regi,the.subset)

#interest rates and ideology
regiid <- zelig(reelecb~irates12+d.gov, model="logit",
	data=the.subset, cite=F)
	regiidcl <- cl(the.subset,
	regiid,"country")$estimates
fd.iid <- fd.badgood.alt(regiid,the.subset)

#commodities only
regc <- zelig(reelecb~log(comm12), model="logit",
	data=the.subset, cite=F)
	regccl <- cl(the.subset,
	regc,"country")$estimates
fd.c <- fd.badgood.alt(regc,the.subset)

#commodities and ideology
regcid <- zelig(reelecb~log(comm12)+d.gov, model="logit",
	data=the.subset)
	regcidcl <- cl(the.subset,
	regcid,"country")$estimates
	print(regcidcl)
fd.cid <- fd.badgood.alt(regcid,the.subset)

#Irates and commodities
reg <- zelig(reelecb~irates12+log(comm12), model="logit",
	data=the.subset)
	regcl <- cl(the.subset,
	reg,"country")$estimates
	print(regcl)
fd <- fd.badgood.alt(reg,the.subset)

#Irates and commodities and ideology (1 is right)
regid <- zelig(reelecb~irates12+log(comm12)+d.gov, model="logit",
	data=the.subset)
	regidcl <- cl(the.subset,
	regid,"country")$estimates
	print(regidcl)
fd.id <- fd.badgood.alt(regid,the.subset)

#Assemble table (in appendix)
reg.table <- merge(
			 merge(
			 merge(stack.regs(regi),stack.regs(regiid),
				   by=c("Var1","Var2"),
				   suffixes=c("IR","IRideo"),all=T),
			 merge(stack.regs(regc),stack.regs(regcid),
			      by=c("Var1","Var2"),
			      suffixes=c("C","Cideo"),all=T),
			      by=c("Var1","Var2"),all=T),
			 merge(stack.regs(reg),stack.regs(regid),
			      by=c("Var1","Var2"),
			      suffixes=c("both","bothideo"),all=T),
			      by=c("Var1","Var2"),all=T)  
			     			
reg.table$Var1 <- ifelse(reg.table$Var2=="Estimate"|reg.table$Var2=="Stat",
				as.character(reg.table$Var1),as.character(reg.table$Var2))
print(xtable(reg.table[,-2],digits=3),include.rownames=F)

#Organize and add First Differences to the table ####
fd.lscebgalt <- rbind('Mod. Irates'=fd.i[c(1,4,5)],
				'Mod. IratesIdeo'=fd.iid[c(1,4,5)],
				'Mod. Comm'=fd.c[c(1,4,5)],
				'Mol. CommTideo'=fd.cid[c(1,4,5)],
				'Mol. Both'=fd[c(1,4,5)],
				'Mol. BothIdeo'=fd.id[c(1,4,5)])
xtable(round(t(fd.lscebgalt),2))
				

### Control for Leigh's "Merit" indicator 		   ### 
### Impelemnted for two different versions of the  ###
### indicator. See last section of the appendix for###
### details 									   ###

### First version of "merit"
reggcomp <- zelig(reelecb~geti12+merit, model="logit",
	data=the.subset,cite=F) 
fd.comp <- fd.badgood.alt(reggcomp,the.subset)
reggcompfe <- zelig(reelecb~geti12+merit+country-1, model="logit",
	data=the.subset,cite=F)#fixed effects
fd.compfe <- fd.badgood.alt(reggcompfe,the.subset,clv=F)
reggcompre <- zelig(formula= reelecb~geti12+merit+ tag(1| country), 
				data=the.subset, model="logit.mixed", cite=F)
fd.compre <- fd.badgood.alt(reggcompre,the.subset,clv=F)
reggcompinc <- zelig(reelecb~geti12+merit+inc.ran, model="logit",
	data=the.subset, cite=F)#
fd.compinc <- fd.badgood.alt(reggcompinc ,the.subset)	
reggcompi <- zelig(reelecb~geti12+merit+d.gov, model="logit",
	data=the.subset, cite=F)#with ideology
fd.compi <- fd.badgood.alt(reggcompi,the.subset)

## Assemble table (for appendix)
reg.table.comp <-  merge(merge( 
	merge(stack.regs(reggcomp),stack.regs(reggcompfe,clse=F),
				   by=c("Var1","Var2"),
				   suffixes=c("detcl","detfe"),all.x=T),
	merge(stack.regs(reggcompre),stack.regs(reggcompreelec),
				   by=c("Var1","Var2"),
				   suffixes=c("detre","detreelec"),all=T),all=T),
	stack.regs(reggcompi),all=T)
reg.table.comp$Var1 <- ifelse(reg.table.comp$Var2=="Estimate"|reg.table.comp$Var2=="Stat",as.character(reg.table.comp$Var1),as.character(reg.table.comp$Var2))
print(xtable(reg.table.comp[,-2],digits=3),include.rownames=F)

### Add FD for version JOP
effects.merit <- rbind('Mod. 1'=fd.comp[c(1,4,5)],
				'Mod. 2'=fd.compfe[c(1,4,5)],
				'Mod. 3'=fd.compre[c(1,4,5)],
				'Mod. 4'=fd.compinc[c(1,4,5)],
				'Mod. 5'=fd.compi[c(1,4,5)])
				
# Create line to report FD and confidence intervals in table (JOP version, with not graph)
xtable(round(t(effects.merit),2))

# First differences for merit
meritlow <- mean(the.subset$merit)-sd(the.subset$merit)
merithigh <- mean(the.subset$merit)+sd(the.subset$merit)
xlow <- setx(reggcomp,merit=meritlow)
xhigh <- setx(reggcomp,merit=merithigh)
summary(sim(reggcomp,x=xlow,x1=xhigh))

#Second definition of "merit"	 
reggcomp <- zelig(reelecb~geti12+meritleigh, model="logit",
	data=the.subset,cite=F) 
reggcompfe <- zelig(reelecb~geti12+meritleigh+country, model="logit",
	data=the.subset,cite=F)#fixed effects
reggcompre <- zelig(formula= reelecb~geti12+meritleigh+ tag(1| country), 
				data=the.subset, model="logit.mixed", cite=F)
reggcompreelec <- zelig(reelecb~geti12+meritleigh+inc.ran, model="logit",
	data=the.subset, cite=F)#
reggcompi <- zelig(reelecb~geti12+meritleigh+d.gov, model="logit",
	data=the.subset, cite=F)#with ideology 
reg.table.compleigh <-  merge(merge( 
	merge(stack.regs(reggcomp),stack.regs(reggcompfe,clse=F),
				   by=c("Var1","Var2"),
				   suffixes=c("detcl","detfe"),all.x=T),
	merge(stack.regs(reggcompre),stack.regs(reggcompreelec),
				   by=c("Var1","Var2"),
				   suffixes=c("detre","detreelec"),all=T),all=T),
	stack.regs(reggcompi),all=T)
reg.table.compleigh$Var1 <- ifelse(reg.table.compleigh$Var2=="Estimate"|reg.table.compleigh$Var2=="Stat",as.character(reg.table.compleigh$Var1),as.character(reg.table.compleigh$Var2))
print(xtable(reg.table.compleigh[,-2],digits=3),include.rownames=F)
			 			 
		 			 
			 
