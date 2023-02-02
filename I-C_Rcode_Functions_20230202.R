#Trade-off between competition ability and invulnerability to predation in marine microbes – contrasting of protist grazing and viral lysis effects
#Jinny Wu Yang, Feng-Hsun Chang, Yi-Chun Yeh, An-Yi Tsai, Fuh-Kwo Shiah, Kuo-Ping Chiang, Gwo-Ching Gong, Chih-hao Hsieh
#jinnyang@umich.edu

##### Function_TradeOffPlot ####
TradeOffPLot=function(x,y,Xlim,Ylim,XlimLabels,YlimLabels,XLAB,YLAB,label,StationCruise){
  plot(x, y, col="black", pch = 1, cex = 1, ylim = Ylim, xlim = Xlim, xlab="",xaxt='n',yaxt='n',ylab="")
  axis(side = 1, at = XlimLabels, labels = FALSE, tck = -0.03)
  axis(1, at = XlimLabels, line = -1, labels = XLAB, cex.axis = 0.5, tick = F)
  axis(side= 2,at = YlimLabels, labels = F,tck = -0.03)
  axis(side= 2,at = YlimLabels, line=-0.5,labels = YLAB, cex.axis = 0.6, tick = F)
  
  mtext(label, line=-1.25, cex=0.75, at=Xlim[1])
  
  l=lm(y~x)
  if (summary(l)$coef[2,4]<=0.05){
    LTY=1
    lwd=2
  }else{
    LTY=2
    lwd=1
  }
  abline(l,lty=LTY,col="blue",lwd=lwd)
  
  Cor=cor(x,y)
  Corp=cor.test(x,y)$p.value
  if (Corp<=0.001){
    Corp="<0.001"
  }else{
    Corp=Corp
  }
  
  rp = vector('expression',3)
  rp[1] = substitute(expression(italic(n) == MYVALUE), list(MYVALUE = format(length(x), digits = 2)))[2]
  rp[2] = substitute(expression(italic(r) == MYVALUE), list(MYVALUE = format(Cor, digits = 2)))[2]
  rp[3] = substitute(expression(italic(p) == MYOTHERVALUE), list(MYOTHERVALUE = format(Corp, digits = 2)))[2]
  #legend('bottomleft', legend = rp, bty = 'n',cex=0.7)
  legend('bottomleft', legend = rp, bty = 'n',cex=0.7)
  legend('topright',StationCruise, bty = 'n',cex=0.85)
}

#### Fit rank-normalize RAD to zipfbrot model and calculate RAD decay coefficient ####
zipfbrot.gamma=function(table){
  PV=numeric()
  PVName=numeric()
  for(i in 1:ncol(table)){
    mod = rad.zipfbrot(table[,i])
    PV[i] = mod$coefficients[2]
    PVName[i]=colnames(table)[i]
  }
  names(PV)=colnames(table)
  names(PV)=sapply(strsplit(sapply(strsplit(PVName,"F02Or2"), "[", 2) ,"R"), "[", 1)
  return(PV)
}

#### LMM ####
LMMEst=function(x,y,r){
  library(nlme) 
  ms=lme(y~x,random= ~1|r, na.action=na.omit)
  #Goodness_of_fitness_p_value=Chi_test$p.value
  #Goodness_of_fitness_X_square=Chi_test$statistic
  LMM=c(summary(ms)[]$tTable[2,],fixef(ms))
  names(LMM)=c("Value" ,"Std. Error","DF","t_value",
               "p_value","Fix_intercept","Fix_slope")  
  LMM=as.data.frame(t(LMM))
  return(LMM)
}
#### Plot Predation-diversity relationship ####
DiversityPlot=function(l,x,Ylim,LMMTable,Xaxis){
  COL=c("red","tomato2","blue","orange","springgreen4","chartreuse3","maroon","tan4","black")
  PCH=c(1,20,5,18,2,17,8,13)
  
  plot(as.numeric(l[[1]]$Dilution),
       as.numeric(l[[1]][,x]), 
       ylim = Ylim, xlim=c(1,4), xlab="",xaxt='n', 
       cex.axis=0.7, ylab="", pch=PCH[1], col=COL[1],cex=0.7)
  
  ll=lm(as.numeric(l[[1]][,x])~as.numeric(l[[1]]$Dilution))
  
  if (summary(ll)$coeff[2,4]<=0.05){
    LTY=1
  }else{
    LTY=2
  }
  
  abline(ll, lty=LTY,col=COL[1])
  
  for (i in 2:6){
    points(as.numeric(l[[i]]$Dilution),
           as.numeric(l[[i]][,x]), 
           ylim = Ylim, xlim=c(1,4), xlab="",xaxt='n', 
           cex.axis=0.7, ylab="", pch=PCH[i], col=COL[i],cex=0.7)
    
    ll=lm(as.numeric(l[[i]][,x])~as.numeric(l[[i]]$Dilution))
    
    if (summary(ll)$coeff[2,4]<=0.05){
      LTY=1
    }else{
      LTY=2
    }
    abline(ll, lty=LTY,col=COL[i])
  }
  if (Xaxis==TRUE){
    axis(side = 1,cex.axis=0.7, at = c(1,2,3,4), labels = c("25%","50%","75%", "100%"), tck = -0.05)
  }else{
    axis(side = 1,cex.axis=0.7, at = c(1,2,3,4), labels = c("","","", ""), tck = -0.03)
  }

  if (LMMTable$p_value<=0.05){
    LTY=1
  }else{
    LTY=2
  }
  if (LMMTable$p_value<=0.001){
    p="<0.001"
  }else{
    p=LMMTable$p_value
  }
  rp = vector('expression',1)
  rp[1] = substitute(expression(italic(LMM-p) == MYOTHERVALUE), 
                     list(MYOTHERVALUE = format(p, digits = 2)))[2]
  
  
  abline(c(LMMTable$Fix_intercept,LMMTable$Fix_slope), lty=LTY,col="black", lwd=3)
  legend('topright', legend = rp, bty = 'n',cex=0.7)
}
DiversityPlot_Combine_Virus=function(V,x,Ylim,LMMTable){
  COL=c("red","tomato2","blue","midnightblue","springgreen4","chartreuse3","maroon","tan4","black")
  PCH=c(1,20,5,18,2,17,8,13)
  
  plot(as.numeric(l[[1]]$Dilution),
       as.numeric(l[[1]][,x]), 
       ylim = Ylim, xlim=c(1,4), xlab="",xaxt='n', 
       cex.axis=0.7, ylab="", pch=PCH[1], col=COL[1],cex=0.8)
  
  ll=lm(as.numeric(l[[1]][,x])~
          as.numeric(V[[1]]$Dilution))
  
  if (summary(ll)$coeff[2,4]<=0.05){
    LTY=1
  }else{
    LTY=2
  }
  
  abline(ll, lty=LTY,col=COL[1])
  
  for (i in 2:7){
    points(as.numeric(V[[i]]$Dilution),
           as.numeric(l[[i]][,x]), 
           ylim = Ylim, xlim=c(1,4), xlab="",xaxt='n', 
           cex.axis=0.7, ylab="", pch=PCH[i], col=COL[i],cex=0.8)
    
    ll=lm(as.numeric(t(V[[i]][,3]))~
            as.numeric(l[[i]][,x]))
    
    if (summary(ll)$coeff[2,4]<=0.05){
      LTY=1
    }else{
      LTY=2
    }
    abline(ll, lty=LTY,col=COL[i])
  }
  
  axis(side = 1,cex.axis=0.8, at = c(1,2,3,4), labels = c("25%","50%","75%", "100%"), tck = -0.05)
  
  if (LMMTable$p_value<=0.001){
    p="<0.001"
  }else{
    p=LMMTable$p_value
  }
  if (LMMTable$p_value<=0.05){
    LTY=1
  }else{
    LTY=2
  }
  rp = vector('expression',1)
  rp[1] = substitute(expression(italic(LMM-p) == MYOTHERVALUE), 
                     list(MYOTHERVALUE = format(p, digits = 2)))[2]
  
  
  abline(c(LMMTable$Fix_intercept,LMMTable$Fix_slope), lty=LTY,col="black", lwd=3)
  legend('topright', legend = rp, bty = 'n',cex=0.7)
}


#### Permutation and Shapiro-Wilk test ####
Permutation_shapiro_test=function(Slope,IGR){
  library("lme4")
  Observed=summary(lm(Slope~IGR))$coef[2,1] 
  Null_Slope_IGR=matrix(0,1000,3)
  colnames(Null_Slope_IGR)=c("Intercept","Estimates","p-value")
  for (j in 1:1000){
    nSlope=Slope[sample(1:length(Slope))] 
    l=lm(nSlope ~ IGR)
    Null_Slope_IGR[j,1]=summary(l)$coef[1,1] #
    Null_Slope_IGR[j,2]=summary(l)$coef[2,1] #Estimates
    Null_Slope_IGR[j,3]=summary(l)$coef[2,4] # 
  }
  Final=numeric()
  Final[1]=shapiro.test(Null_Slope_IGR[,2])$statistic
  Final[2]=shapiro.test(Null_Slope_IGR[,2])$p.value
  z = (Observed-mean(Null_Slope_IGR[,2]))/sd(Null_Slope_IGR[,2])
  Final[3]=z 
  p = 2*pnorm(-abs(z)) # two-tail p-value
  Final[4]=p
  names(Final)=c("shapiro","shapiro_pvalue","zscore","p-value")
  return(Final)
}
