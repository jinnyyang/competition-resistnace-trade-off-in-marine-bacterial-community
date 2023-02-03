#Trade-off between competition ability and invulnerability to predation in marine microbes – contrasting of protist grazing and viral lysis effects
#Jinny Wu Yang, Feng-Hsun Chang, Yi-Chun Yeh, An-Yi Tsai, Fuh-Kwo Shiah, Kuo-Ping Chiang, Gwo-Ching Gong, Chih-hao Hsieh
#jinnyang@umich.edu
source("I-C_Rcode_Functions_20230202.R")
library(gdata)
library(phyloseq)
library(ggpubr)
library(stats)
library(gdata)
library(lavaan)
library(lme4)
library(nlme)
library(MASS)
library(RADanalysis)
library(reshape2)

######### Part 1: Calculate net growth rate and for each ASV########
##### A. Per-capita net growth rate (PNGR) ~ top-down control dilution factors (TCDF) relationship #########

RawTable=read.table("Final_table_dada2.tsv",header=T, row.names=1) #Input raw ASV table
SampleID=read.xls("SampleName.xlsx", sheet=1,header=F) #Input sample name list
ps.rarefied = rarefy_even_depth(otu_table(t(RawTable),taxa_are_rows=F), 
                                rngseed=1, sample.size=3188, replace=T)
ps.rarefied=t(ps.rarefied)
Sample_all=colnames(ps.rarefied)
StationCruiseNames=c("2014AprSt1","2014OctSt1","2014OctSt9",
                     "2015JulSt1","2015JulSt9","2016MaySt1")
write.csv(ps.rarefied,"Rarefied_table_At3188Reads.csv")

#Input Flow cytometery enumeration data
T0_Flow=read.xls("FlowCytometry_Enumeration_Data_202207.xlsx",sheet=1,header=T,row.names=1)
T0_30_Flow=read.xls("FlowCytometry_Enumeration_Data_202207.xlsx",sheet=3,header=T,row.names=1)
T12_Flow=read.xls("FlowCytometry_Enumeration_Data_202207.xlsx",sheet=2,header=T,row.names=1)
T12_30_Flow=read.xls("FlowCytometry_Enumeration_Data_202207.xlsx",sheet=4,header=T,row.names=1)

GrowthRate_D02_All=list()
GrowthRate_D30_All=list()
for (j in 1:ncol(SampleID)){
  Com_Col=grep(SampleID[1,j], Sample_all, value = TRUE)
  rare.table=ps.rarefied[,Com_Col]  
Flow_Col_T0_D02=T0_Flow[grep(SampleID[2,j], rownames(T0_Flow), value = F),]
Flow_Col_T12_D02=T12_Flow[grep(SampleID[2,j], rownames(T12_Flow), value = F),]
Flow_Col_T0_D30=T12_Flow[grep(SampleID[2,j], rownames(T0_30_Flow), value = F),]
Flow_Col_T12_D30=T12_30_Flow[grep(SampleID[2,j], rownames(T12_30_Flow), value = F),]

T0_Com=rare.table[,grep("T0",colnames(rare.table), value = F)]
T12_D30_Com=rare.table[,grep("D30",colnames(rare.table), value = F)]
T12_D02_Com=rare.table[,grep("DS",colnames(rare.table), value = F)[1:8]]

T0_ForEachASV_D02=list()
T0_ForEachASV_D30=list()
T12_ForEachASV_D02=list()
T12_ForEachASV_D30=list()
for (i in 1:8){
  T0_ForEachASV_D02[[i]]=T0_Com[,1]*Flow_Col_T0_D02$TB_D[i]/3188  
  T0_ForEachASV_D30[[i]]=T0_Com[,1]*Flow_Col_T0_D30$TB_D[i]/3188  
  T12_ForEachASV_D02[[i]]=T12_D02_Com[,i]*Flow_Col_T12_D02$TB_D[i]/3188  
  T12_ForEachASV_D30[[i]]=T12_D30_Com[,i]*Flow_Col_T12_D30$TB_D[i]/3188  
}
NetGrowthRate_D02=matrix(0,nrow(T12_ForEachASV_D02[[1]]),8)
NetGrowthRate_D30=matrix(0,nrow(T12_ForEachASV_D30[[1]]),8)
for (i in 1:8){
  NetGrowthRate_D02[,i]=log(T12_ForEachASV_D02[[i]]/T0_ForEachASV_D02[[i]])/12
  NetGrowthRate_D30[,i]=log(T12_ForEachASV_D30[[i]]/T0_ForEachASV_D30[[i]])/12
}
colnames(NetGrowthRate_D02)=colnames(T12_D02_Com)
rownames(NetGrowthRate_D02)=rownames(T12_D02_Com)
colnames(NetGrowthRate_D30)=colnames(T12_D30_Com)
rownames(NetGrowthRate_D30)=rownames(T12_D30_Com)


write.csv(NetGrowthRate_D02, quote=F,file = paste("NetGrowthRate_D02", 
                                   SampleID[1,j],"csv",sep="."))
write.csv(NetGrowthRate_D30, quote=F,file = paste("NetGrowthRate_D30", 
                                                  SampleID[1,j],"csv",sep="."))

}

##### B. Extract competition ability (intercetp) and resistance to top-down control (regression slope) from PBGR~TCDF relationsihp #########
#Input sample list 
table_list = list.files(pattern="Net*.*csv") 
s = strsplit(table_list, "\\.")
CruiseStation=sapply(s,"[",2)
ss = strsplit(sapply(s,"[",1), "\\_")
Treatment=sapply(ss,"[",2)

for (i in 1:length(table_list)){
  NetGrowth=read.csv(table_list[i],header=T,row.names = 1) #Input 
  is.na(NetGrowth)=sapply(NetGrowth, is.infinite)
  NetGrowth[is.na(NetGrowth)]=0
  
  
  Dilution=sapply(strsplit(colnames(NetGrowth), "D"),"[",3)
  
  EachASV=matrix(0,nrow(NetGrowth),3)
  for (j in 1:nrow(NetGrowth)){
    l=lm(as.numeric(NetGrowth[j,])~as.numeric(Dilution))
    EachASV[j,1]=coef(l)[2] #Slope
    EachASV[j,2]=summary(l)$coefficients[2,4] #p-value
    EachASV[j,3]=coef(l)[1] #intercept
  }
  rownames(EachASV)=rownames(NetGrowth)
  colnames(EachASV)=c("Slope","p-value","intercept")
  
  write.csv(EachASV, quote=F,file = paste("EachASV_Regression_", 
                                          CruiseStation[i],Treatment[i],"csv",sep="."))
}


  
##### C. Plot Figure 2: Plot Competitiveness-Invulnerability trade-off  
  HI_list_D02 = list.files(pattern="EachASV_Regression_.*.D02.csv") 
  HI_list_D30 = list.files(pattern="EachASV_Regression_.*.D30.csv") 

  CruiseStation=sapply(strsplit(HI_list_D02, "\\_|\\.| "), "[", 4)
  Treatment=sapply(strsplit(HI_list_D02, "\\_|\\.| "), "[", 5)
  
  
  par(mfrow = c(4,3)) # 2-by-2 grid of plots
  par(oma = c(4, 4, 2, 2)) # make room for the overall x and y axis titles
  par(mar = c(0, 0, 0, 0)) # make the plots be closer together
  
  ## Protists-caused mortality (Slope_D02) ~ Predation-free growth rate (IGR) 
  Xlim=c(-0.15,0.15)
  XlimLabels=c(-0.15,-0.1,-0.05,0,0.05,0.1,0.15)
  Ylim=c(-0.6,0.6)
  YlimLabels=c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6)
  label=c("a","b","c","d","e","f")
  YLAB=c(T,F,F,T,F,F)
  XLAB=c(F,F,F,F,F,F)
  for (i in 1:6){
    t_D02=read.csv(HI_list_D02[i],header=T,row.names = 1)
    t_D30=read.csv(HI_list_D30[i],header=T,row.names = 1)
    Data=data.frame(t_D02$Slope,t_D30$intercept)
    Data=Data[Data[,1]!=0,]
    Data=Data[Data[,2]!=0,]
    x=Data$t_D02.Slope
    y=Data$t_D30.intercept
    TradeOffPLot(x,y,Xlim,Ylim,XlimLabels,YlimLabels,XLAB[i],YLAB[i],label[i],StationCruiseNames[i])
  }  

  ## Viral+Protists-caused mortality (SlopeD30) ~ Predation-free growth rate (IGR) ####
  label=c("g","h","i","j","k","l")
  YLAB=c(T,F,F,T,F,F)
  XLAB=c(F,F,F,T,T,T)
  for (i in 1:6){
    t_D02=read.csv(HI_list_D02[i],header=T,row.names = 1)
    t_D30=read.csv(HI_list_D30[i],header=T,row.names = 1)
    #Data=data.frame(t_D30$Slope-t_D02$Slope,t_D30$intercept)
    Data=data.frame(t_D30$Slope,t_D30$intercept)
    Data=Data[Data[,1]!=0,]
    Data=Data[Data[,2]!=0,]
    #x=Data$t_D30.Slope...t_D02.Slope
    x=Data$t_D30.Slope
    y=Data$t_D30.intercept
    TradeOffPLot(x,y,Xlim,Ylim,XlimLabels,YlimLabels,
                 XLAB[i],YLAB[i],label[i],StationCruiseNames[i])
  }  
  
  mtext("Resistance to",side=2, outer=T,cex=0.8, line=2.5, at=0.5)
  mtext("Protists",side=2, outer=T,cex=0.7, line=1.25, at=0.875)
  mtext("Protists",side=2, outer=T,cex=0.7, line=1.25, at=0.625)
  mtext("Viruses+protists",side=2, outer=T,cex=0.7, line=1.25, at=0.375)
  mtext("Viruses+protists",side=2, outer=T,cex=0.7, line=1.25, at=0.125)
  mtext("Top-down control-free growth rate",side=1, outer=T,cex=0.8, line=1)
  
  ##### C. Plot Figure S6: Viral-caused mortality (SlopeD30-Slope_D02) ~ Predation-free growth rate (IGR) ####
  par(mfrow = c(2,3)) # 2-by-2 grid of plots
  par(oma = c(4, 4, 2, 2)) # make room for the overall x and y axis titles
  par(mar = c(0, 0, 0, 0)) # make the plots be closer together
  label=c("a","b","c","d","e","f")
  YLAB=c(T,F,F,T,F,F)
  XLAB=c(F,F,F,T,T,T)
  for (i in 1:6){
    t_D02=read.csv(HI_list_D02[i],header=T,row.names = 1)
    t_D30=read.csv(HI_list_D30[i],header=T,row.names = 1)
    Data=data.frame(t_D30$Slope-t_D02$Slope,t_D30$intercept)
    Data=Data[Data[,1]!=0,]
    Data=Data[Data[,2]!=0,]
    x=Data$t_D30.Slope...t_D02.Slope
    y=Data$t_D30.intercept
    TradeOffPLot(x,y,Xlim,Ylim,XlimLabels,YlimLabels,
                 XLAB[i],YLAB[i],label[i],StationCruiseNames[i])
  }  
  mtext("Resistance to viruses",side=2, outer=T,cex=0.8, line=1.5, at=0.5)
  mtext("Top-down control-free growth rate",side=1, outer=T,cex=0.8, line=1)
  
#################################### End of Part 1 ########################################
  
  
  
################# Part 2: Preparing diversity index  ############################
  #### A. Richness and Rank-based normalization #### 
  # Rank-normalized ASV tables are for obtaining RAD and Evenness

  Table=ps.rarefied
  ### Finding minimum ranking among samples (smallest number of ASV among samples)
  TTable=as.data.frame(Table)
  TTable[TTable>=1]=1 #Making ASV number of reads > 1 become 1
  Richness=colSums(TTable)
  ### Normalizes an abundance table to the desired number of ranks 
  nrads=RADnormalization_matrix(input = Table, max_rank = min(Richness),average_over = 1000, sample_in_row = F)
  
  nrad_mat=t(nrads$norm_matrix)
  colnames(nrad_mat)=colnames(Table)
  write.table(nrad_mat,file="RankFixed_table.txt",quote=F)
  write.table(Richness,file="Richness.txt",quote=F)
  
  #### B. Calculate RAD decay coefficient and Evenness from rank-normalized RAD ####
  library(vegan)
  library(reshape)
  library("lme4")
  library(nlme)
  library(viridis)
  ### C. Calculate bacterial RAD and Evenness under protists and protists-viral effect ####
  ## Calculate RAD decay coefficient with rank with rank-normalized table
  gamma=zipfbrot.gamma(nrad_mat)
  gamma=gamma*(-1) # covert to positive value as decay coefficient
  
  ## Calculate Evenness with rank-normalized table
  H = diversity(t(nrad_mat))
  J = H/log(specnumber(t(nrad_mat)))
  
  ##Prepare and organize data
  finalTable=data.frame(gamma,J,H)
  colnames(finalTable)=c("gamma","Evenness","H")
  
  write.table(finalTable,file="RankFixed_RAD_Evenness.txt",quote=F)
  
#################################### End of Part 2 ########################################
  
########### Part 3: LMM analysis ################################################
  
  #### A. Resistance ~ competition ####
  temp = list.files(pattern="HypothesisI_variables.*.txt") 
  l=list()
  for (i in 1:length(temp)){
    table=read.table(temp[i], header=T,row.names=1)
    Sample=sapply(strsplit(temp,"_"), "[", 3)[i]
    Station=sapply(strsplit(Sample,"St"), "[", 2)
    l[[i]]=data.frame(table,rep(Sample,nrow(table)),rep(Station,nrow(table)))
    colnames(l[[i]])=c(colnames(table),"Sample","Station")
    names(l)[i]=Sample
  }
  rl=rbind(l[[1]],l[[2]],l[[3]],l[[4]],l[[5]],l[[6]],l[[7]])
  
  #Do LMM with function "LM_Estimates"
  LMM_Estimates=rbind(LMMEst(rl$IGR,rl$Slope_D02,rl$Sample),
                      LMMEst(rl$IGR,rl$Slope_D30-rl$Slope_D02,rl$Sample),
                      LMMEst(log10(rl$T0_AbsA),rl$Slope_D02,rl$Sample),
                      LMMEst(log10(rl$T0_AbsA),rl$Slope_D30-rl$Slope_D02,rl$Sample))
  
  colnames(LMM_Estimates)=c("Value" ,"Std. Error","DF","t-value","p-value","Fix_intercept","Fix_slope")     
  rownames(LMM_Estimates)=c("Growth_Protists","Growth_Viruses",
                            "Density_Protists","Density_Viruses")
  write.csv(LMM_Estimates,"LMM_R-C.csv", quote=F) # Output LMM result for Hypothesis I
  
  
  #### B. Diversity ~ Top-down control dilution factors #### 
  Dilution_WithT0=sapply(strsplit(names(Richness),"D"),"[",3)
  T0=which(is.na(Dilution_WithT0))
  NoT0_Richness=Richness[-T0]
  Dilution=sapply(strsplit(names(NoT0_Richness),"D"),"[",3)
  Treatment=sapply(strsplit(names(NoT0_Richness),"S"),"[",1)
  CruiseStation=sapply(strsplit(sapply(strsplit(names(NoT0_Richness),"Or"),"[",2),"R"),"[",1)
  
  Richness_D02=data.frame(NoT0_Richness[which(Treatment=="D")],
                          as.numeric(Dilution[which(Treatment=="D")]),
                          CruiseStation[which(Treatment=="D")])
  Richness_D30=data.frame(NoT0_Richness[which(Treatment=="D30")],
                          as.numeric(Dilution[which(Treatment=="D30")]),
                          CruiseStation[which(Treatment=="D30")])
  colnames(Richness_D02)=c("Richness","Dilution","CruiseStation")
  colnames(Richness_D30)=c("Richness","Dilution","CruiseStation")
  
  x=Richness_D02$Dilution
  y=Richness_D02$Richness
  r=Richness_D02$CruiseStation
  Richness_D02_LMM=LMMEst(x,y,r)
  x=Richness_D30$Dilution
  y=Richness_D30$Richness
  r=Richness_D30$CruiseStation
  Richness_D30_LMM=LMMEst(x,y,r)
  
  #For evenness and RAD
  Dilution_WithT0=sapply(strsplit(rownames(finalTable),"D"),"[",3)
  T0=which(is.na(Dilution_WithT0))
  NoT0_finalTable=finalTable[-T0,]
  Dilution=sapply(strsplit(rownames(NoT0_finalTable),"D"),"[",3)
  Treatment=sapply(strsplit(rownames(NoT0_finalTable),"S"),"[",1)
  CruiseStation=sapply(strsplit(sapply(strsplit(rownames(NoT0_finalTable),"Or"),"[",2),"R"),"[",1)
  
  EvennessRAD_D02=data.frame(NoT0_finalTable[which(Treatment=="D"),],
                             as.numeric(Dilution[which(Treatment=="D")]),
                             CruiseStation[which(Treatment=="D")])
  EvennessRAD_D30=data.frame(NoT0_finalTable[which(Treatment=="D30"),],
                             as.numeric(Dilution[which(Treatment=="D30")]),
                             CruiseStation[which(Treatment=="D30")])
  colnames(EvennessRAD_D02)=c("Gamma","Evenness","H","Dilution","CruiseStation")
  colnames(EvennessRAD_D30)=c("Gamma","Evenness","H","Dilution","CruiseStation")
  x=EvennessRAD_D02$Dilution
  y=EvennessRAD_D02$Evenness
  r=EvennessRAD_D02$CruiseStation
  Evenness_D02_LMM=LMMEst(x,y,r)
  x=EvennessRAD_D30$Dilution
  y=EvennessRAD_D30$Evenness
  r=EvennessRAD_D30$CruiseStation
  Evenness_D30_LMM=LMMEst(x,y,r)
  
  x=EvennessRAD_D02$Dilution
  y=EvennessRAD_D02$Gamma
  r=EvennessRAD_D02$CruiseStation
  RAD_D02_LMM=LMMEst(x,y,r)
  x=EvennessRAD_D30$Dilution
  y=EvennessRAD_D30$Gamma
  r=EvennessRAD_D30$CruiseStation
  RAD_D30_LMM=LMMEst(x,y,r)
  
  
  LMM_Diversity=rbind(Richness_D02_LMM,Richness_D30_LMM,
                      Evenness_D02_LMM,Evenness_D30_LMM,
                      RAD_D02_LMM,RAD_D30_LMM)
  rownames(LMM_Diversity)=c("Richness_D02_LMM","Richness_D30_LMM",
                            "Evenness_D02_LMM","Evenness_D30_LMM",
                            "RAD_D02_LMM","RAD_D30_LMM")
  
  write.csv(LMM_Diversity,"LMM_Diversity.csv")
  
#################################### End of Part 3 ########################################
  

########### Part 4: Plot Figure 3:Top-down control effect on diversity ################################################  
  COL=c("red","tomato2","blue","orange","springgreen4","chartreuse3","maroon","tan4","black")
  PCH=c(1,20,5,18,2,17,8,13)
  par(mfrow = c(3, 2)) # 1-by-3 grid of plots
  par(oma = c(3, 3, 2, 2)) # make room axis titles
  par(mar = c(1, 2, 0, 0)) # adjust width between plots 
  
  StationCruiseNames=c("2014AprSt1","2014OctSt1","2014OctSt9",
                       "2015JulSt1","2015JulSt9","2016MaySt1")
  
  RDA_Ev_D02=split(EvennessRAD_D02,EvennessRAD_D02$CruiseStation)
  RDA_Ev_D30=split(EvennessRAD_D30,EvennessRAD_D30$CruiseStation)
  x=1 #gamma
  l=RDA_Ev_D02
  Ylim=c(0.5,3)
  LMMTable=LMM_Diversity[5,]
  Xaxis=FALSE
  DiversityPlot(l,x,Ylim,LMMTable,Xaxis)
  mtext("Protists-diluted", outer=F,side=3, line=0, cex=0.7)
  mtext("a", outer=F,side=3, line=-1, at=1, cex=0.8)
  mtext("RDA decay coefficient", outer=F,side=2, line=2, cex=0.7)
  
  x=1
  l=RDA_Ev_D30
  Ylim=c(0.5,3)
  LMMTable=LMM_Diversity[6,]
  Xaxis=FALSE
  DiversityPlot(l,x,Ylim,LMMTable,Xaxis)
  mtext("Protists+viruses-diluted", outer=F,side=3, line=0, cex=0.7)
  mtext("b", outer=F,side=3, line=-1, at=1, cex=0.8)
  #######################
  ######Evenness#####
  x=2 #Evenness
  l=RDA_Ev_D02
  Ylim=c(0.35,1)
  LMMTable=LMM_Diversity[3,]
  Xaxis=FALSE
  DiversityPlot(l,x,Ylim,LMMTable,Xaxis)
  #mtext("Protists-manipulated", outer=F,side=3, line=0, cex=0.7)
  mtext("c", outer=F,side=3, line=-1, at=1, cex=0.8)
  mtext("Evenness", outer=F,side=2, line=2, cex=0.7)
  
  x=2
  l=RDA_Ev_D30
  Ylim=c(0.35,1)
  LMMTable=LMM_Diversity[4,]
  Xaxis=FALSE
  DiversityPlot(l,x,Ylim,LMMTable,Xaxis)
  mtext("d", outer=F,side=3, line=-1, at=1, cex=0.8)
  mtext("Top-down control dilution factors", outer=T,side=1, line=1, cex=0.7)
  
  #######################
  #####Richness##########
  Ri_D02=split(Richness_D02,Richness_D02$CruiseStation)
  Ri_D30=split(Richness_D30,Richness_D30$CruiseStation)
  
  x=1
  l=Ri_D02
  Ylim=c(40,410)
  LMMTable=LMM_Diversity[1,]
  Xaxis=TRUE
  DiversityPlot(l,x,Ylim,LMMTable,Xaxis)
  #mtext("Protists-diluted", outer=F,side=3, line=0, cex=0.7)
  mtext("e", outer=F,side=3, line=-1, at=1, cex=0.8)
  mtext("Richness", outer=F,side=2, line=2, cex=0.7)
  
  x=1
  l=Ri_D30
  Ylim=c(40,410)
  LMMTable=LMM_Diversity[2,]
  Xaxis=TRUE
  DiversityPlot(l,x,Ylim,LMMTable,Xaxis)
  #mtext("Protists+viruses-diluted", outer=F,side=3, line=0, cex=0.7)
  mtext("f", outer=F,side=3, line=-1, at=1, cex=0.8)
  legend("bottomright",
         legend=StationCruiseNames,
         pch=PCH[1:6],col=COL[1:6],cex=0.6)
  mtext("Top-down control dilution factors", outer=T,side=1, line=1, cex=0.7)
  
#################################### End of Part 4 ########################################
  
############ Part 5: Identify top-down control resistant and susceptible taxa #########################
  
  HI_list_D02 = list.files(pattern="EachASV_Regression_.*.D02.csv") 
  HI_list_D30 = list.files(pattern="EachASV_Regression_.*.D30.csv") 
  Taxonomy=read.table("taxonomy.tsv",header=T,sep="\t")
  ASV_regressionWithTax_D02=list()
  ASV_regressionWithTax_D30=list()
  for (j in 1:6){
    tax_D02=read.csv(HI_list_D02[j],header=T,row.names = 1)
    tax_D30=read.csv(HI_list_D30[j],header=T,row.names = 1)
    tax_all=data.frame(tax_D02$Slope,tax_D30$Slope,tax_D30$intercept,tax_D02$p.value,tax_D30$p.value)
    rownames(tax_all)=rownames(tax_D02)
    tax_all_D02=tax_all[tax_all$tax_D02.Slope!=0,]
    SigOrNot_D02=tax_all_D02$tax_D02.p.value<=0.05
    tax_all_D30=tax_all[tax_all$tax_D30.Slope!=0,]
    SigOrNot_D30=tax_all_D30$tax_D30.p.value<=0.05
    
    TT_D02=numeric()
    TT_D30=numeric()
    for (i in 1:length(rownames(tax_all_D02))){
      TT_D02[i]=grep(rownames(tax_all_D02)[i], Taxonomy[,1], value = F)
    }
    Order_D02=sapply(strsplit(Taxonomy[TT_D02,]$Taxon,";"), "[", 4)
    ASV_regressionWithTax_D02[[j]]=data.frame(tax_all_D02$tax_D02.Slope,
                                              tax_all_D02$tax_D30.intercept,
                                              Taxonomy[TT_D02,]$Feature.ID,
                                              Order_D02,
                                              SigOrNot_D02)
    # if (length(rownames(Sig_D30)==0)){
    
    for (i in 1:length(rownames(tax_all_D30))){
      TT_D30[i]=grep(rownames(tax_all_D30)[i], Taxonomy[,1], value = F) 
      
    } 
    #}
    Order_D30=sapply(strsplit(Taxonomy[TT_D30,]$Taxon,";"), "[", 4)
    ASV_regressionWithTax_D30[[j]]=data.frame(tax_all_D30$tax_D30.Slope,
                                              tax_all_D30$tax_D30.intercept,
                                              Taxonomy[TT_D30,]$Feature.ID,
                                              Order_D30,
                                              SigOrNot_D30)
  }
  
  names(ASV_regressionWithTax_D02)=sapply(strsplit(HI_list_D02, "\\_|\\.| "), "[", 4)
  names(ASV_regressionWithTax_D30)=sapply(strsplit(HI_list_D30, "\\_|\\.| "), "[", 4)
  
  #########Plotting#####
  M_D02=melt(ASV_regressionWithTax_D02,id=1:5)
  Split_Order=sapply(strsplit(M_D02$Order, "__|;"), "[", 2)
  ForPloting_D02=data.frame(M_D02[,1:2],Split_Order,M_D02[,5],M_D02[,6])
  colnames(ForPloting_D02)=c("Slope","Intercept","Order","SigOrNot","Experiments")
  
  M_D30=melt(ASV_regressionWithTax_D30,id=1:5)
  Split_Order=sapply(strsplit(M_D30$Order, "__|;"), "[", 2)
  ForPloting_D30=data.frame(M_D30[,1:2],Split_Order,M_D30[,5],M_D30[,6])
  colnames(ForPloting_D30)=c("Slope","Intercept","Order","SigOrNot","Experiments")
  
  
  PlotFigure5_SigOnly=function(ForPloting_D02,Xlab){
    library(dplyr)
    
    #Remove non-sig data
    ForPloting_D02=filter(ForPloting_D02,ForPloting_D02$SigOrNot==T)
    
    ForPloting_D02$Order[which(is.na(ForPloting_D02$Order))]="uncultured"
    ForPloting_D02=filter(ForPloting_D02,ForPloting_D02$Order!="marine metagenome",
                          ForPloting_D02$Order!="unidentified marine bacterioplankton",
                          ForPloting_D02$Order!="uncultured",
                          ForPloting_D02$Order!="Unclassified",
                          ForPloting_D02$Order!="unidentified marine eubacterium",
                          ForPloting_D02$Order!="uncultured bacterium")
    
    
    
    NewOrderName=rep("a",length(ForPloting_D02$Order))
    
    for (i in 1:length(levels(as.factor(ForPloting_D02$Order)))){
      O=levels(as.factor(ForPloting_D02$Order))[i]
      OO=filter(ForPloting_D02,ForPloting_D02$Order==O)
      OrderWithSig=sprintf("%s (%s/%s)",O,
                           length(which(OO$SigOrNot==T & OO$Slope >= 0)),
                           length(which(OO$SigOrNot==T & OO$Slope <= 0)))
      NewOrderName[which(ForPloting_D02$Order==O)]=OrderWithSig
    }
    
    ForPloting_D02=data.frame(ForPloting_D02,NewOrderName)
    
    #Rank Order based on resistance to predation
    ForPloting_D02$NewOrderName = factor(ForPloting_D02$NewOrderName, 
                                         levels=names(sort(tapply(ForPloting_D02$Slope,ForPloting_D02$NewOrderName, mean))))
    
    
    ggplot(data=ForPloting_D02,aes(y=NewOrderName,x=Slope,horizontal=TRUE)) +
      geom_boxplot(outlier.colour = NA) + theme_bw()+
      geom_jitter(aes(),cex=1,shape=1,
                  position=position_jitter(width = 0,height = 0)) + 
      theme(axis.text=element_text(size=5.5))+
      xlab(Xlab)+
      ylab("")+ theme(legend.text=element_text(size=6),
                      axis.title.x=element_text(size=7),
                      axis.text.y=element_text(face="bold"))+
      geom_vline(xintercept = 0, col="blue")
    
    
    
    
    
  }
  
  
  pD02_tax=PlotFigure5_SigOnly(ForPloting_D02,"Resistance to protists")
  
  pD30_tax=PlotFigure5_SigOnly(ForPloting_D30,"Resistance to protists+viruses")
  
  ggarrange(pD02_tax,pD30_tax,nrow=1, widths = c(1, 1),
            labels = c("A","B"),
            font.label = list(size = 10),
            legend="right", common.legend = TRUE) +
    theme(plot.margin = margin(0.6,0.1,0.1,0.1, "cm"))
  
#################################### End of Part 4 ########################################
  
  
############ Part 6: supplementary table and figures #########################
  
#####Table S2: Permutation test ################################################
  temp = list.files(pattern="EachASV_Regression_.*.D02.csv") #Input variable table list
  #### 1. Permutation test on Predation-free growth rate versus Protists-caused mortality ####
  Final=matrix(0,6,4)
  for (i in 1:6){
    table=read.csv(temp[i],header=T,row.names = 1)
    table=table[!is.na(table[,2]),]
    Final[i,]=Permutation_shapiro_test(table$Slope,table$intercept)
  }
  colnames(Final)=c("Shapiro-Wilk","Shapiro-Wilk_p-value","z-score","p-value")
  rownames(Final)=sapply(strsplit(temp, "_", fixed = TRUE), "[", 3)
  write.csv(Final,"Permutation_shapiro_Protist-IGR_RegresionSlope.csv")
  
  #### 2. Permutation test on Predation-free growth rate versus Protists+viruses-caused mortality ####
  temp = list.files(pattern="EachASV_Regression_.*.D30.csv") #Input variable table list
  
  Final=matrix(0,6,4)
  for (i in 1:6){
    table=read.csv(temp[i],header=T,row.names = 1)
    table=table[!is.na(table[,2]),]
    Final[i,]=Permutation_shapiro_test(table$Slope,table$intercept)
  }
  colnames(Final)=c("Shapiro-Wilk","Shapiro-Wilk_p-value","z-score","p-value")
  rownames(Final)=sapply(strsplit(temp, "_", fixed = TRUE), "[", 3)
  write.csv(Final,"Permutation_shapiro_Protists+viruses-IGR_RegresionSlope.csv")

##########################################################################################
  
############## Figure S3: Total abundance at T0 and T12 plot #####
    
    #Input Flow cytometery enumeration data
    T0_Flow=read.xls("FlowCytometry_Enumeration_Data_202207.xlsx",sheet=1,header=T,row.names=1)
    T0_30_Flow=read.xls("FlowCytometry_Enumeration_Data_202207.xlsx",sheet=3,header=T,row.names=1)
    T12_Flow=read.xls("FlowCytometry_Enumeration_Data_202207.xlsx",sheet=2,header=T,row.names=1)
    T12_30_Flow=read.xls("FlowCytometry_Enumeration_Data_202207.xlsx",sheet=4,header=T,row.names=1)
    Flow_D02=rbind(T0_Flow,T12_Flow)
    Flow_D30=rbind(T0_30_Flow,T12_30_Flow)
    
    Flow_D02$Dilution=as.factor(Flow_D02$Dilution)
    CruiseStation=sapply(strsplit(rownames(Flow_D02),"Or2|D"), "[", 2)
    Flow_D02=data.frame(Flow_D02,CruiseStation)
    
    Flow_D30$Dilution=as.factor(Flow_D30$Dilution)
    CruiseStation=sapply(strsplit(rownames(Flow_D30),"Or2|D"), "[", 2)
    Flow_D30=data.frame(Flow_D30,CruiseStation)
    
    library(ggplot2)
    pD02=ggplot(Flow_D02,aes(x=CruiseStation,y=log(TB_D),col=Incubation,fill=Dilution))+
      geom_boxplot(lwd=0.2)+ theme_bw() + 
      geom_point(aes(col = Incubation,fill=Dilution),
                 alpha = 0.6,size = 1,
                 position = position_jitterdodge(0)) +
      scale_fill_discrete(name = "Dilution factors", labels = c("25%", "50%", "75%","100%")) +
      scale_color_discrete(name = "Collection time", labels = c("T0", "T12")) +
      scale_fill_manual(values=c("pink","orange","cyan4","darkblue"))+
      scale_color_manual(values=c("dimgrey","black"))+
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank()) + ylab("Total density")

    pD30=ggplot(Flow_D30,aes(x=CruiseStation,y=log(TB_D),col=Incubation,fill=Dilution))+geom_boxplot()+ theme_bw() +
      geom_boxplot(lwd=0.2)+ theme_bw() + 
      geom_point(aes(col = Incubation,fill=Dilution),
                 alpha = 0.6,size = 1,
                 position = position_jitterdodge(0)) +
      scale_fill_discrete(name = "Dilution factors", labels = c("25%", "50%", "75%","100%")) +
      scale_color_discrete(name = "Collection time", labels = c("T0", "T12")) +
      scale_fill_manual(values=c("pink","orange","cyan4","darkblue"))+
      scale_color_manual(values=c("dimgrey","black"))+
      theme(axis.title.x = element_blank())+
      ylab("Total density")+scale_x_discrete(labels= StationCruiseNames)

    
     library(ggpubr)
    ggarrange(pD02,pD30,
              labels = c("(a) Protists-diluted",
                         "(b) Protists+viruses-diluted"),
              font.label = list(size = 8, color = "black"),
              label.x=c(0,-0.04), label.y=c(0.2,0.25),
              ncol = 1, nrow = 2, legend="right", common.legend = T) 
    
    
###### Taxonomy plot ######Not sure if use
    Order=sapply(strsplit(Taxonomy[,2], "__|;"), "[", 8)
    Class=sapply(strsplit(Taxonomy[,2], "__|;"), "[", 6)
    Phylum=sapply(strsplit(Taxonomy[,2], "__|;"), "[", 4)
    
    New_tax=data.frame(Order,Class,Phylum)
    rownames(New_tax)=Taxonomy[,1]
    pp=list()
    for (j in 1:6){
      Com_Col=grep(SampleID[1,j], Sample_all, value = TRUE)
      rare.table=ps.rarefied[,Com_Col]    
      p=phyloseq(otu_table(as.matrix(rare.table),taxa_are_rows=F),
                 tax_table(as.matrix(New_tax)))
      pp[[j]]=plot_bar(p,fill="Phylum")+ 
        scale_x_discrete(limits=colnames(rare.table),labels=Label)+
        theme_bw()+
        labs(title=SampleID[1,j], y="# of reads", x="Sample") +
        theme(title=element_text(size=8),
              axis.title.y=element_text(size=8),
              axis.title.x=element_text(size=8),
              axis.text.y = element_text(size = 5),
              axis.text.x = element_text(size = 5),
              legend.text = element_text(size = 6),
              legend.title = element_blank(),
              legend.position="right")
    }
    
    Label=c("D30_0%","D30_25%","D30_50%","D30_75%","D30_100%",
            "D30_0%","D30_25%","D30_50%","D30_75%","D30_100%",
            "D_0%","D_25%","D_50%","D_75%","D_100%",
            "D_0%","D_25%","D_50%","D_75%","D_100%","T0")
    
    ggarrange(pp[[1]],pp[[2]], legend="right", common.legend = TRUE)
    

    
    
    
  
  
  
    
##########################################################################################   
    
############## Figure S4: Top-down control effect on community composition maintenance #####       
    
##########################################################################################   
############## Figure S5: Taxonomy all #####
    
    HI_list_D02 = list.files(pattern="EachASV_Regression_.*.D02.csv") 
    HI_list_D30 = list.files(pattern="EachASV_Regression_.*.D30.csv") 
    Taxonomy=read.table("taxonomy.tsv",header=T,sep="\t")
    ASV_regressionWithTax_D02=list()
    ASV_regressionWithTax_D30=list()
    for (j in 1:6){
      tax_D02=read.csv(HI_list_D02[j],header=T,row.names = 1)
      tax_D30=read.csv(HI_list_D30[j],header=T,row.names = 1)
      tax_all=data.frame(tax_D02$Slope,tax_D30$Slope,tax_D30$intercept,tax_D02$p.value,tax_D30$p.value)
      rownames(tax_all)=rownames(tax_D02)
      tax_all_D02=tax_all[tax_all$tax_D02.Slope!=0,]
      SigOrNot_D02=tax_all_D02$tax_D02.p.value<=0.05
      tax_all_D30=tax_all[tax_all$tax_D30.Slope!=0,]
      SigOrNot_D30=tax_all_D30$tax_D30.p.value<=0.05
      
      TT_D02=numeric()
      TT_D30=numeric()
      for (i in 1:length(rownames(tax_all_D02))){
        TT_D02[i]=grep(rownames(tax_all_D02)[i], Taxonomy[,1], value = F)
      }
      Order_D02=sapply(strsplit(Taxonomy[TT_D02,]$Taxon,";"), "[", 4)
      ASV_regressionWithTax_D02[[j]]=data.frame(tax_all_D02$tax_D02.Slope,
                                                tax_all_D02$tax_D30.intercept,
                                                Taxonomy[TT_D02,]$Feature.ID,
                                                Order_D02,
                                                SigOrNot_D02)
      
      for (i in 1:length(rownames(tax_all_D30))){
        TT_D30[i]=grep(rownames(tax_all_D30)[i], Taxonomy[,1], value = F)
      }
      Order_D30=sapply(strsplit(Taxonomy[TT_D30,]$Taxon,";"), "[", 4)
      ASV_regressionWithTax_D30[[j]]=data.frame(tax_all_D30$tax_D30.Slope,
                                                tax_all_D30$tax_D30.intercept,
                                                Taxonomy[TT_D30,]$Feature.ID,
                                                Order_D30,
                                                SigOrNot_D30)
    }
    
    names(ASV_regressionWithTax_D02)=sapply(strsplit(HI_list_D02, "\\_|\\.| "), "[", 4)
    names(ASV_regressionWithTax_D30)=sapply(strsplit(HI_list_D30, "\\_|\\.| "), "[", 4)
    
    ##Plotting
    M_D02=melt(ASV_regressionWithTax_D02,id=1:5)
    Split_Order=sapply(strsplit(M_D02$Order, "__|;"), "[", 2)
    ForPloting_D02=data.frame(M_D02[,1:2],Split_Order,M_D02[,5],M_D02[,6])
    colnames(ForPloting_D02)=c("Slope","Intercept","Order","SigOrNot","Experiments")
    
    M_D30=melt(ASV_regressionWithTax_D30,id=1:5)
    Split_Order=sapply(strsplit(M_D30$Order, "__|;"), "[", 2)
    ForPloting_D30=data.frame(M_D30[,1:2],Split_Order,M_D30[,5],M_D30[,6])
    colnames(ForPloting_D30)=c("Slope","Intercept","Order","SigOrNot","Experiments")
    
    
    PlotFigure5=function(ForPloting_D02,Xlab){
      library(dplyr)
      ForPloting_D02$Order[which(is.na(ForPloting_D02$Order))]="uncultured"
      ForPloting_D02=filter(ForPloting_D02,ForPloting_D02$Order!="marine metagenome",
                            ForPloting_D02$Order!="unidentified marine bacterioplankton",
                            ForPloting_D02$Order!="uncultured",
                            ForPloting_D02$Order!="Unclassified",
                            ForPloting_D02$Order!="unidentified marine eubacterium",
                            ForPloting_D02$Order!="uncultured bacterium")
      
      NewOrderName=rep("a",length(ForPloting_D02$Order))
      
      for (i in 1:length(levels(as.factor(ForPloting_D02$Order)))){
        O=levels(as.factor(ForPloting_D02$Order))[i]
        OO=filter(ForPloting_D02,ForPloting_D02$Order==O)
        OrderWithSig=sprintf("%s (%s/%s/%s)",O,
                             length(which(OO$SigOrNot==T & OO$Slope >= 0)),
                             length(which(OO$SigOrNot==T & OO$Slope <= 0)),
                             length(OO$SigOrNot))
        NewOrderName[which(ForPloting_D02$Order==O)]=OrderWithSig
      }
      
      ForPloting_D02=data.frame(ForPloting_D02,NewOrderName)
      
      
      #Rank Order based on resistance to predation
      ForPloting_D02$NewOrderName = factor(ForPloting_D02$NewOrderName, 
                                           levels=names(sort(tapply(ForPloting_D02$Slope,ForPloting_D02$NewOrderName, mean))))
      
      ggplot(data=ForPloting_D02,aes(y=NewOrderName,x=Slope,horizontal=TRUE)) +
        geom_boxplot(outlier.colour = NA)+theme_bw()+
        geom_jitter(aes(shape=SigOrNot,color=SigOrNot),cex=0.9,
                    position=position_jitter(width = 0,
                                             height = 0)) + 
        scale_color_manual("",values = c("gray","black"), labels = c("Non-Sig","Sig"))+
        scale_shape_manual("",values = c(1,16), labels = c("Non-Sig","Sig"))+
        theme(axis.text=element_text(size=5.5))+
        xlab(Xlab)+
        ylab("")+ theme(legend.text=element_text(size=6),
                        axis.title.x=element_text(size=7),
                        axis.text.y=element_text(face="bold"))+
        geom_vline(xintercept = 0, col="blue")
      
      
      
      
      
    }
    
    
    pD02_tax=PlotFigure5(ForPloting_D02,"Resistance to protists")
    
    pD30_tax=PlotFigure5(ForPloting_D30,"Resistance to protists+viruses")
    
    
    ggarrange(pD02_tax,pD30_tax,nrow=1, widths = c(1, 1),
              labels = c("A","B"),
              font.label = list(size = 10),
              legend="right", common.legend = TRUE) +
      theme(plot.margin = margin(0.6,0.1,0.1,0.1, "cm"))
    
    
##########################################################################################
############## Figure S7: Initial composition#####   
    
    RawTable=read.table("Final_table_dada2.tsv",header=T, row.names=1) #Input raw ASV table
    colnames(RawTable)
    Taxonomy=read.table("taxonomy.tsv",header=T,sep="\t",row.names=1)
    TT=Taxonomy[,1]
    TTT1=sapply(strsplit(sapply(strsplit(TT,";"),"[",2),"_"),"[",4)
    TTT2=sapply(strsplit(sapply(strsplit(TT,";"),"[",3),"_"),"[",4)
    TTT3=sapply(strsplit(sapply(strsplit(TT,";"),"[",4),"_"),"[",4)
    TAX=data.frame(TTT1,TTT2,TTT3)
    rownames(TAX)=rownames(Taxonomy)
    colnames(TAX)=c("D2","D3","D4")
    T0_number=which(sapply(strsplit(colnames(RawTable),"T"),"[",2)==0)
    T0=RawTable[,T0_number]
    colnames(T0)=c("2014AprSt1","2014OctSt1","2014OctSt9",
                   "2015JulSt1","2015JulSt9","2016MaySt1")
    iht=phyloseq(otu_table(as.matrix(T0),taxa_are_rows=T),tax_table(as.matrix(TAX)))
    
    rT0=rarefy_even_depth(otu_table(as.matrix(T0),taxa_are_rows=T), 
                          rngseed=1, sample.size=3188, replace=T)
    iht=phyloseq(otu_table(as.matrix(rT0),taxa_are_rows=T),tax_table(as.matrix(TAX)))
    
    
    pd = psmelt(iht)
    
    ggplot(pd, mapping=aes(x = Sample, y = Abundance,fill=D3))+
      geom_bar(stat="identity")+
      theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
            legend.key.height = unit(0.5, 'cm'), #change legend key height
            legend.key.width = unit(0.5, 'cm'), #change legend key width
            legend.title = element_text(size=0), #change legend title font size
            legend.text = element_text(size=6)) +
      guides(fill=guide_legend(ncol =3))+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    theme(legend.position="none")
    
    grep(TAX[2,1],rownames(T0), value = F) 
    
    
##########################################################################################
############## Figure S8: Total bacterial community net growth rate plot#####   
    
    StationCruiseNames=c("2014AprSt1","2014OctSt1","2014OctSt9",
                         "2015JulSt1","2015JulSt9","2016MaySt1")
    
    #Input Flow cytometery enumeration data
    T0_Flow=read.xls("FlowCytometry_Enumeration_Data_202207.xlsx",sheet=1,header=T,row.names=1)
    T0_30_Flow=read.xls("FlowCytometry_Enumeration_Data_202207.xlsx",sheet=3,header=T,row.names=1)
    T12_Flow=read.xls("FlowCytometry_Enumeration_Data_202207.xlsx",sheet=2,header=T,row.names=1)
    T12_30_Flow=read.xls("FlowCytometry_Enumeration_Data_202207.xlsx",sheet=4,header=T,row.names=1)
    SampleID=read.xls("SampleName.xlsx", sheet=1,header=F) #Input sample name list
    
    GrowthRate_All=list()
    ps.rarefied=read.csv("Rarefied_table_At3188Reads.csv")
    Sample_all=colnames(ps.rarefied)
    
    for (j in 1:ncol(SampleID)){
      Com_Col=grep(SampleID[1,j], Sample_all, value = TRUE)
      rare.table=ps.rarefied[,Com_Col]  
      Flow_Col_T0_D02=T0_Flow[grep(SampleID[2,j], rownames(T0_Flow), value = F),]
      Flow_Col_T12_D02=T12_Flow[grep(SampleID[2,j], rownames(T12_Flow), value = F),]
      Flow_Col_T0_D30=T0_30_Flow[grep(SampleID[2,j], rownames(T0_30_Flow), value = F),]
      Flow_Col_T12_D30=T12_30_Flow[grep(SampleID[2,j], rownames(T12_30_Flow), value = F),]
      
      
      
      NetGrowthRate_D02=log(Flow_Col_T12_D02[,7]/Flow_Col_T0_D02[,7])/12
      NetGrowthRate_D30=log(Flow_Col_T12_D30[,7]/Flow_Col_T0_D30[,7])/12
      
      GrowthRate_All[[j]]=data.frame(Flow_Col_T0_D02[,1:5],NetGrowthRate_D02,NetGrowthRate_D30)
    }
    
    par(mfrow = c(2, 3)) # 1-by-3 grid of plots
    par(oma = c(3, 3, 2, 2)) # make room axis titles
    par(mar = c(1, 1, 0, 0)) # adjust width between plots 
    Xaxis=c(F,F,F,T,T,T)
    Yaxis=c(T,F,F,T,F,F)
    XlimLabels=c("25%","50%","75%", "100%")
    YLAB=c(1,2,3,4,5)
    for (i in 1:6){
      plot(GrowthRate_All[[i]]$Dilution,GrowthRate_All[[i]]$NetGrowthRate_D02,pch=16,ylim=c(-0.1,0.3),
           xaxt='n',yaxt='n')
      if (Xaxis[i]==TRUE){
        axis(side = 1,cex.axis=0.7, at = c(25,50,75,100), labels = c("25%","50%","75%", "100%"), tck = -0.05)
      }else{
        axis(side = 1,cex.axis=0.7, at = c(25,50,75,100), labels = c("","","", ""), tck = -0.03)
      }
      if (Yaxis[i]==TRUE){
        axis(side = 2,cex.axis=0.7, at = c(-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4), labels = c(-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4), tck = -0.05)
      }else{
        axis(side = 2,cex.axis=0.7, at = c(-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4), labels = c("","","","","","","","",""), tck = -0.03)
      }
      
      l=lm(GrowthRate_All[[i]]$NetGrowthRate_D02~GrowthRate_All[[i]]$Dilution)
      if (summary(l)$coef[2,4]<=0.05){
        LTY=1
      }else{LTY=2}
      abline(l,lty=LTY)
      points(GrowthRate_All[[i]]$Dilution,GrowthRate_All[[i]]$NetGrowthRate_D30,pch=16,col="red") 
      l=lm(GrowthRate_All[[i]]$NetGrowthRate_D30~GrowthRate_All[[i]]$Dilution)
      if (summary(l)$coef[2,4]<=0.05){
        LTY=1
      }else{LTY=2}
      abline(l,col="red",lty=LTY)
      
      mtext(StationCruiseNames[i], outer=F,side=3, line=-1.5, at=43, cex=0.8)
      
    }
    
    mtext("Top-down control dilution factors", outer=T,side=1, line=1.5, at=0.5, cex=0.8)
    mtext("Net growth rate", outer=T,side=2, line=1.5, at=0.5, cex=0.8)
    
    legend("bottomright",legend=c("protists-diluted","protists+viruses-diluted"),pch=16,col=c("black","red"), cex=0.7)
    
##########################################################################################    
    
    
############## Figure S9: Community composition shift with Dilution factor after 12 hours incubation#####  
    
    Taxonomy=read.table("taxonomy.tsv",header=T,sep="\t")
    RawTable=read.table("Final_table_dada2.tsv",header=T, row.names=1) #Input raw ASV table
    
    RarefiedTable=read.csv("Rarefied_table_At3188Reads.csv", row.names=1)
    
    TaxOrder=numeric()
    for (i in 1:nrow(RarefiedTable)){
      TaxOrder[i]=grep(rownames(RarefiedTable)[i], Taxonomy[,1], value = F)
    }
    Tax=Taxonomy[TaxOrder,][,2] #keep the taxonomy name
    Split_taxa=strsplit(Tax, "__|;")
    Genuses=sapply(Split_taxa, "[", 12)
    Families=sapply(Split_taxa, "[", 10)
    Orders=sapply(Split_taxa, "[", 8)
    Classes=sapply(Split_taxa, "[",  6)
    Phylums=sapply(Split_taxa, "[",  4)
    TAX=data.frame(Phylums,Classes,Orders,Families,Genuses)
    rownames(TAX)=rownames(RarefiedTable)
    
    
    
    SampleNames=colnames(RarefiedTable)
    S1=strsplit(SampleNames, "F02|R|D|T")
    Filtration=sapply(S1, "[",  2)
    Cruise=sapply(S1, "[",  3)
    Replicates=sapply(S1, "[",  4)
    Dilution=sapply(S1, "[",  5)
    
    tt=paste(Filtration,Cruise, Dilution)
    
    
    Table=t(rowsum(t(RarefiedTable), group = tt, na.rm = T))
    
    for (i in 1:ncol(Table)){
      if (sum(Table[,i])>3199){
        Table[,i]=Table[,i]/2
      }else{
        Table[,i]=Table[,i]
      }
    }
    
    SS=strsplit(colnames(Table), " ")
    Cruise=sapply(SS, "[",  2)
    CC=names(split(Cruise,Cruise))
    Filtration=sapply(SS, "[",  1)
    FF=names(split(Filtration,Filtration))
    Dilution=sapply(SS, "[",  3)
    #TableT0=Table[,which(Cruise==CC[1] & Dilution=="NA")]
    p_D02=list()
    for (i in 1:6){
      nTable=Table[,which(Cruise==CC[i] & Filtration==FF[1])]  
      colnames(nTable)=c("25%","50%","75%","100%")
      iht=phyloseq(otu_table(nTable,taxa_are_rows=T),tax_table(as.matrix(TAX)))
      pd <- psmelt(iht)
      
      p_D02[[i]]=ggplot(pd, mapping=aes(x = Sample, y = Abundance,fill=Orders))+
        geom_bar(stat="identity")+
        theme(plot.margin = margin(0.1,0,0.5,0, "cm"))+
        theme_minimal()+
        theme(axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())+ 
        scale_x_discrete(limit = c("25%","50%","75%","100%"))+ 
        theme(axis.text = element_text(size = 5)) +
        theme(legend.key.size = unit(0.2, 'cm'))+
        theme(legend.text = element_text(size=3))+
        theme(legend.title = element_text(size=4))
    }
    
    p_D30=list()
    for (i in 1:6){
      nTable=Table[,which(Cruise==CC[i] & Filtration==FF[2] & Dilution!="NA")]  
      colnames(nTable)=c("25%","50%","75%","100%")
      iht=phyloseq(otu_table(nTable,taxa_are_rows=T),tax_table(as.matrix(TAX)))
      pd <- psmelt(iht)
      
      p_D30[[i]]=ggplot(pd, mapping=aes(x = Sample, y = Abundance,fill=Orders))+
        geom_bar(stat="identity")+
        theme(plot.margin = margin(0.1,0,0.5,0, "cm"))+
        theme_minimal()+
        theme(axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())+ 
        scale_x_discrete(limit = c("25%","50%","75%","100%"))+ 
        theme(axis.text = element_text(size = 5)) +
        theme(legend.key.size = unit(0.2, 'cm'))+
        theme(legend.text = element_text(size=3))+
        theme(legend.title = element_text(size=4))
    }
    #+theme(legend.position = "none") +
    library(ggpubr)
    StationCruiseNames=c("2014AprSt1","2014OctSt1","2014OctSt9",
                         "2015JulSt1","2015JulSt9","2016MaySt1")
    figure=ggarrange(p_D02[[1]],p_D02[[2]],p_D02[[3]],p_D02[[4]],p_D02[[5]],p_D02[[6]],
                     p_D30[[1]],p_D30[[2]],p_D30[[3]],p_D30[[4]],p_D30[[5]],p_D30[[6]],
                     ncol = 6, nrow = 2, labels=StationCruiseNames, 
                     font.label = list(size = 7),
                     common.legend = TRUE, legend = "bottom")
    
    annotate_figure(figure,
                    bottom = text_grob("Top-down control dilution factors",
                                       hjust = 1.75, x = 0, size = 6),
                    left = text_grob("Protists+viruses-diluted                                      Prosits-diluted",
                                     hjust = -0.05, x = 1, size=7, rot = 90))
########################################################################################################
