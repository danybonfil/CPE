#### Examine the effect of GET on Domestic Economy Indices and Other Variables ###
#### The idea here is to compare the effect of GET on the domestic economy 	   ###
#### Contrasting countries in the sample and out of the sample  		       ###
#### This is the code that replicates Table 1 and Figure 1 in the JOP paper    ###
#### Requires data file "data-panel.RData", with country/year observations     ###

stack.test <- function(x){
		require(reshape)
		require(car)
		tmp <- summary(x)$coef
		tmp <- tmp[grep("geti",rownames(tmp)),-3]
		tmp <- melt(tmp)
		tmp <- tmp[order(tmp$X1),]
		tmp2 <- data.frame(X1="diff",X2="Pr(>|t|)",
						value=linearHypothesis(x,
								"geti:lsceFALSE = geti:lsceTRUE",
								white.adjust=F)$Pr[2])
		tmp3 <- data.frame(X1="N",X2="N",value=nrow(x$model))
		out<-rbind(tmp,tmp2,tmp3)
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

#### Replicate Figure 1 ###########
load("DATA/data_panel.RData")
print(comment(gdpla))
lat80<-gdpla[which(gdpla$country=="lat" & gdpla$year>="1980"),]

pdf("Figures/fig_latcommirates80.pdf", width = 10.00, height = 10.00)
#jpeg("fig_latcommirates80.jpg", width = 10.00, height = 10.00, res=300, units="in")
par(mar=c(7,5,1,6))
plot(lat80$year,lat80$irates,type="l",frame.plot="TRUE",ylab="% (Nominal Rates)", xlab="", lwd=3, lty=2,pch=16,col=gray(0.6),cex.axis=1.8,cex.lab=1.8,bty="n")
polygon(x=c(1975,1990,1990,1975),y=c(0,0,20,20),col=gray(0.7),density=30,border=NA)
polygon(x=c(2003,2015,2015,2003),y=c(0,0,20,20),col=gray(0.7),density=30,border=NA)
lines(lat80$year,lat80$irates,lty=2,lwd=3,col=gray(0.6))
par(new=TRUE)
plot(lat80$year,lat80$comm
    ,type="l"
    ,axes=FALSE
    ,bty='c'
    ,lty = 3
    ,ylab =""
    ,xlab=""
    ,lwd= 3
    ,bty="n"
    ,pch=16
    ,col=gray(0.6)
    ,cex.axis=1.8
    ,cex.lab=2)
axis(4,cex.axis=1.5)
mtext("Index Points (year 2000=100)", 4, line=4,cex=1.8)
par(new=T)
plot(lat80$year,lat80$geti,type="l"
    ,axes=FALSE
    ,bty='c'
    ,lty = 1
    ,ylab =""
    ,xlab=""
    ,lwd= 5
    ,pch=16
    ,bty="n"
    ,col="black"
    ,cex.axis=1.8
    ,cex.lab=2)

legend("bottom",legend=c("US Int. Rates","Commodities","GET Index"), horiz=T,cex=2,col=c(gray(0.6),gray(0.6),"black"),bty="n",lty=c(2,3,1),lwd=c(3,3,4),xpd=NA,inset=c(0,-0.18))
dev.off()


#### Replicate Table 1 ###########
load("DATA/data_panel.RData")
print(comment(gdpla))
library(car)
library(xtable)
gdpla <- subset(gdpla,country!="lat")

#### Using only year >= startyear 
#### Table 1 in the RR version to JOP
Ygdpla <- subset(gdpla,year>=startyear)

reg01 <- lm(log(gdpnorm)~geti:lsce+log(gdpnormlag)+country-1,data=Ygdpla)#intersting
reg02 <- lm(dgdp~geti:lsce+country+dgdplag-1,data=Ygdpla)
reg03 <- lm(inflognorm~geti:lsce+inflognormlag+country-1,data=Ygdpla)
reg04 <- lm(unempnorm~geti:lsce+unempnormlag+country-1,data=Ygdpla)
reg05 <- lm(okun~geti:lsce+okunlag+country-1,data=Ygdpla)
reg06 <- lm(hanke~geti:lsce+hankelag+country-1,data=Ygdpla)
#### Results are the same using non-normalized versions of the variables ####

the.table <- merge(merge(
merge(stack.test(reg01),stack.test(reg02),by=c("X1","X2"),suffixes=c("gdp","growth")),
merge(stack.test(reg03),stack.test(reg04),by=c("X1","X2"),suffixes=c("infl","unemp")),by=c("X1","X2")),
merge(stack.test(reg05),stack.test(reg06),by=c("X1","X2"),suffixes=c("okun","hanke")),by=c("X1","X2"))[c(2,4,3,5,7,6,1,8),]

xtab1 <- xtable(the.table[,-c(2,3)],digits=3,
			caption = "Effects of GET on Domestic Economy Indicators for LSCE and Comparison Group",
			label = "tab-getdomestic")
print(xtab1,include.rownames=F,hline.after=c(-1,0,3,6,7,8),caption.placement = "top")
