### Analysis of TSE Phenomaster data

### Author: Cristina Vicente Garcia - cvicgar@gmail.com
### Date: June 2022
### Description: This script takes a file 

#!/usr/bin/env Rscript

## Installing and loading the required libraries

if (!("ggplot2" %in% installed.packages())) { 
  BiocManager::install("ggplot2",dependencies=TRUE)
}
if (!("reshape2" %in% installed.packages())) { 
  BiocManager::install("reshape2",dependencies=TRUE)
}
if (!("Rmisc" %in% installed.packages())) { 
  BiocManager::install("Rmisc",dependencies=TRUE)
}
if (!("rlist" %in% installed.packages())) { 
  BiocManager::install("rlist",dependencies=TRUE)
}
if (!("optparse" %in% installed.packages())) { 
  BiocManager::install("optparse",dependencies=TRUE)
}
if (!("dplyr" %in% installed.packages())) { 
  BiocManager::install("dplyr",dependencies=TRUE)
}
if (!("ggsci" %in% installed.packages())) { 
  BiocManager::install("ggsci",dependencies=TRUE)
}

library(ggplot2)
library(reshape2)
library(Rmisc)
library(rlist)
library(optparse)
library(dplyr)
library(ggsci)


## Defining functions: plot_ind(), plot_bxplt(), plot_mean(), processFile()

# plot_ind() takes as inputs the whole dataset and the parameter to plot and returns two types of plots: 
# 1) individual profiles of that parameter over time, highlighting the dark and light phases and grouped by genotype, 
# resulting in as many panels as different genotypes, 2) individual profiles of that parameter over time, 
# highlighting the dark and light phases, resulting in as many panels as animals in the experiment.
plot_ind<-function(parameter,data,plot_labels_temp){
  plots<-list()

  for (i in 1:length(parameter)){
    p<-ggplot(data=data,aes_string(x="TimePoint",y=parameter[i],colour=as.factor(data$Genotype),group="Animal"))+
      facet_wrap(~Genotype)+
      geom_rect(data=dark_phase,mapping=aes(xmin=rep(unlist(dark_phase$dark_start),length(genotype)),xmax=rep(unlist(dark_phase$dark_end),length(genotype)),ymin=-Inf,ymax=+Inf),inherit.aes = F,fill="azure2")+
      geom_point(size=0.1,show.legend=F)+
      geom_line(show.legend=F)+
      scale_color_locuszoom()+
      theme_bw()+
      scale_x_continuous("Time point",breaks=seq(1,length(plot_labels_temp),10),labels=plot_labels_temp[seq(1,length(plot_labels_temp),10)])+
      theme(strip.background = element_rect(fill="white"),strip.text = element_text(face="bold"),panel.grid = element_blank(),axis.text.x = element_text(angle=90))
    g<-ggplot(data=data,aes_string(x="TimePoint",y=parameter[i],colour="Genotype",group="Animal"))+
      geom_rect(data=dark_phase,mapping=aes(xmin=rep(unlist(dark_phase$dark_start),nrow(weight_data)),xmax=rep(unlist(dark_phase$dark_end),nrow(weight_data)),ymin=-Inf,ymax=+Inf),inherit.aes = F,fill="azure2")+
      geom_point(size=0.1)+
      geom_line()+
      scale_color_locuszoom()+
      facet_wrap(~Animal,nrow=length(unique(data$Genotype)))+
      theme_bw()+
      scale_x_continuous("Time point",breaks=seq(1,length(plot_labels_temp),20),labels=plot_labels_temp[seq(1,length(plot_labels_temp),20)])+
      theme(strip.background = element_rect(fill="white"),strip.text = element_text(face="bold"),panel.grid = element_blank(),axis.text.x = element_text(angle=90))
    plots<-list.append(plots,p,g)
  }
  return(plots)
}

# plot_bxplt() 
plot_bxplt<-function(parameter,data,pval_wg,folder,dark_phase){
  write.table("\n######### In boxplots #########","./Results/Warnings.txt",append=T,quote=F,col.names=F,row.names=F)
  write.table("\n######### In boxplots #########","./Results/Warnings_summary.txt",append=T,quote=F,col.names=F,row.names=F)
  special<-c("VO2.3.","VCO2.3.","H.3.")
  plots<-list()
  
  for (i in 1:length(parameter)){
    meandata<-unique(data[c("Genotype","Animal","WgStart")])
    meandata$LIGHT<-rep(0,nrow(meandata))
    meandata$DARK<-rep(0,nrow(meandata))
    
    for(y in 1:length(meandata$Genotype)){
      meandata[y,"DARK"]<-mean(data[data$Animal==as.character(meandata[y,2]) & data$Phase=="DARK",][[parameter[i]]])
      meandata[y,"LIGHT"]<-mean(data[data$Animal==as.character(meandata[y,2]) & data$Phase=="LIGHT",][[parameter[i]]])
    }
    dat.m <- melt(meandata,id.vars=c("Genotype","Animal","WgStart"), measure.vars=c("DARK","LIGHT"))
    colnames(dat.m)<-c("Genotype","Animal","WgStart","Phase",parameter[i])
    
    #Boxplots separated by genotype and phase
    r<-ggplot(data=dat.m,aes_string("Phase",parameter[i]))+geom_boxplot(fill=rep(c("grey36","white"),length(genotype)))+
      facet_wrap(~Genotype)+stat_summary(fun.y=mean,geom="point",color="darkred")+
      stat_summary(fun.y=mean,geom="text",color="black",vjust=-1.5,aes(label=round(..y..,digits=2)))+
      theme_bw()+theme(strip.background = element_rect(fill="white"),strip.text = element_text(face="bold"),panel.grid = element_blank())
    #Individual boxplots
    s<-ggplot(data=data,aes_string("Phase",parameter[i],fill="Genotype"))+geom_boxplot()+
      stat_summary(fun.y=mean,geom="point",color="darkred")+
      stat_summary(fun.y=mean,geom="text",color="black",vjust=-1.5,aes(label=round(..y..,digits=2)))+
      facet_wrap(~Animal,nrow=length(unique(data$Genotype)))+
      theme_bw()+theme(strip.background = element_rect(fill="white"),strip.text = element_text(face="bold"),panel.grid = element_blank())+
      scale_fill_locuszoom()
    plots<-list.append(plots,r,s)
    
    dat.light<-dat.m[dat.m$Phase=="LIGHT",]
    dat.dark<-dat.m[dat.m$Phase=="DARK",]
    dat.ph<-list(dat.light,dat.dark)
    
    if(parameter[i] %in% special){
      v<-ggplot(data=dat.m,aes_string("WgStart",parameter[i],colour="Genotype"))+
        geom_point()+
        facet_wrap(~Phase)+
        theme_bw()+
        theme(strip.background = element_rect(fill="white"),strip.text = element_text(face="bold"))+
        scale_color_locuszoom()
      plots<-list.append(plots,v)
      
      write.table("######################\n### ANCOVA RESULTS ###\n######################",file=paste0(folder,parameter[i],"_ANCOVA_statistics.txt"),quote=F,row.names=F,col.names=F)
      
      for (k in 1:length(dat.ph)){
        if(k==1){
          phase<-"Light"
        } else {
          phase<-"Dark "
        }
        
        write.table(paste0("\n###################\n### ",phase," phase ###\n###################\n"),file=paste0(folder,parameter[i],"_ANCOVA_statistics.txt"),append=T,quote=F,row.names=F,col.names=F)
        
        ancova<-aov(dat.ph[[k]][[parameter[i]]] ~ dat.ph[[k]][["WgStart"]] + dat.ph[[k]][["Genotype"]])
        write.table("\n### WITHOUT INTERACTIONS ###\n",file=paste0(folder,parameter[i],"_ANCOVA_statistics.txt"),append=T,quote=F,row.names=F,col.names=F)
        capture.output(summary(ancova),file=paste0(folder,parameter[i],"_ANCOVA_statistics.txt"),append=T)
        
        if(summary(ancova)[[1]][["Pr(>F)"]][[2]]<0.05){
          write.table(paste(parameter[i],"without interaction (ANCOVA)",sep=" "),"./Results/Significant_5e-2_boxplots.txt",quote=FALSE,append=T,row.names=F,col.names=F)
          if(length(genotype)>2){
            tukey<-TukeyHSD(ancova)
            write.table("\n# Tukey multiple pairwise comparisons #\n",file=paste0(folder,parameter[i],"_ANCOVA_statistics.txt"),append=T,quote=F,row.names=F,col.names=F)
            capture.output(tukey,file=paste0(folder,parameter[i],"_ANCOVA_statistics.txt"),append=T)
          }
        }
        
        aov_residuals <- residuals(object = ancova)
        shap<-shapiro.test(x = aov_residuals )
        
        if(shap$p.value<0.05){
          write.table(paste("\n\n",parameter[i],"- Phase:",phase,"- Warning in Shapiro-Wilk's test without interaction",sep=" "),file="./Results/Warnings.txt",append = T,quote=F,row.names=F,col.names=F)
          write.table(paste(parameter[i],"- Phase:",phase,"- Warning in Shapiro-Wilk's test without interaction",sep=" "),file="./Results/Warnings_summary.txt",append = T,quote=F,row.names=F,col.names=F)
          capture.output(shap,file="./Results/Warnings.txt",append=T)
        }
        
        ancova_int<-aov(dat.ph[[k]][[parameter[i]]] ~ dat.ph[[k]][["WgStart"]] * dat.ph[[k]][["Genotype"]])
        write.table("\n### WITH INTERACTIONS ###\n",file=paste0(folder,parameter[i],"_ANCOVA_statistics.txt"),append=T,quote=F,row.names=F,col.names=F)
        capture.output(summary(ancova_int),file=paste0(folder,parameter[i],"_ANCOVA_statistics.txt"),append=T)
        
        if(summary(ancova_int)[[1]][["Pr(>F)"]][[2]]<0.05){
          write.table(paste(parameter[i],"with interaction (ANCOVA)",sep=" "),"./Results/Significant_5e-2_boxplots.txt",quote=FALSE,append=T,row.names=F,col.names=F)
          if(length(genotype)>2){
            tukey_int<-TukeyHSD(ancova_int)
            write.table("\n# Tukey multiple pairwise comparisons #\n",file=paste0(folder,parameter[i],"_ANCOVA_statistics.txt"),append=T,quote=F,row.names=F,col.names=F)
            capture.output(tukey_int,file=paste0(folder,parameter[i],"_ANCOVA_statistics.txt"),append=T)
          }
        }
        
        aov_residuals_int <- residuals(object = ancova_int)
        shap_int<-shapiro.test(x = aov_residuals_int )
        if(shap_int$p.value<0.05){
          write.table(paste("\n\n",parameter[i],"Phase:",phase,"- Warning in Shapiro-Wilk's test with interaction",sep=" "),file="./Results/Warnings.txt",append = T,quote=F,row.names=F,col.names=F)
          write.table(paste(parameter[i],"Phase:",phase,"- Warning in Shapiro-Wilk's test with interaction",sep=" "),file="./Results/Warnings_summary.txt",append = T,quote=F,row.names=F,col.names=F)
          capture.output(shap_int,file="./Results/Warnings.txt",append=T)
        }
        
      }
      
    }
    
    
    meansddat.m<-summarySE(dat.m,measurevar = parameter[i],groupvars = c("Genotype","Phase"))
    colnames(meansddat.m)<-c("Genotype", "Phase", "N",parameter[i],"sd" , "se","ci")
    meansddat.m<-meansddat.m %>% mutate_if(is.numeric,round,digits=3)
    write.table(as.matrix(meansddat.m),file=paste0(folder,parameter[i],"_boxplots_statistics.txt"),quote=F,col.names=TRUE,row.names=F,sep = "\t")
    
    for (k in 1:length(dat.ph)){
      if(k=="1"){
        phase<-"Light"
      } else {
        phase<-"Dark "
      }
      
      write.table(paste0("\n###################\n### ",phase," phase ###\n###################\n"),file=paste0(folder,parameter[i],"_boxplots_statistics.txt"),append=T,quote=F,row.names=F,col.names=F)
      res.aov<-aov(dat.ph[[k]][[parameter[i]]] ~ dat.ph[[k]][["Genotype"]])
      capture.output(summary(res.aov),file=paste0(folder,parameter[i],"_boxplots_statistics.txt"),append=T)
      if(summary(res.aov)[[1]][["Pr(>F)"]][[1]]<0.05){
        write.table(paste0(parameter[i]," on the ",phase," phase"),"./Results/Significant_5e-2_boxplots.txt",quote=FALSE,append=T,row.names=F,col.names=F,sep="\t")
        if(length(genotype)>2){
          tukey<-TukeyHSD(res.aov)
          write.table("\n# Tukey multiple pairwise comparisons #\n",file=paste0(folder,parameter[i],"_boxplots_statistics.txt"),append=T,quote=F,row.names=F,col.names=F)
          capture.output(tukey,file=paste0(folder,parameter[i],"_boxplots_statistics.txt"),append=T)
        }
      }
      
      bartlett<-bartlett.test(dat.ph[[k]][[parameter[i]]] ~ dat.ph[[k]][["Genotype"]])
      if (bartlett$p.value<0.05){
        write.table(paste("\n\n",parameter[i],"Phase:",phase,"- Warning in Bartlett's test",sep=" "),file="./Results/Warnings.txt",append = T,quote=F,row.names=F,col.names=F)
        write.table(paste(parameter[i],"Phase:",phase,"- Warning in Bartlett's test",sep=" "),file="./Results/Warnings_summary.txt",append = T,quote=F,row.names=F,col.names=F)
        capture.output(bartlett,file="./Results/Warnings.txt",append = T)
      }
      
      aov_residuals <- residuals(object = res.aov)
      shap<-shapiro.test(x = aov_residuals)
      if(shap$p.value<0.05){
        write.table(paste("\n\n",parameter[i]," - Phase:",phase,"- Warning in Shapiro-Wilk's test",sep=" "),file="./Results/Warnings.txt",append = T,quote=F,row.names=F,col.names=F)
        write.table(paste(parameter[i]," - Phase:",phase,"- Warning in Shapiro-Wilk's test",sep=" "),file="./Results/Warnings_summary.txt",append = T,quote=F,row.names=F,col.names=F)
        capture.output(shap,file="./Results/Warnings.txt",append=T)
      }
    }
    
  }
  return(plots)
}

# plot_mean() 
plot_mean<-function(parameter,data,folder,dark_phase,plot_labels_temp){
  write.table("\n######### In meanplots #########","./Results/Warnings.txt",append=T,quote=F,col.names=F,row.names=F)
  write.table("\n######### In meanplots #########","./Results/Warnings_summary.txt",append=T,quote=F,col.names=F,row.names=F)
  plots<-list()
  TPvsP<-unique(data[c("TimePoint","Phase")])
  tpvsp<-rep(TPvsP$Phase,length(unique(data$Genotype)))
  
  for (i in 1:length(parameter)){
    meansddf<-summarySE(data,measurevar=parameter[i],groupvars=c("Genotype","TimePoint"))
    meansddf$Phase<-tpvsp
    colnames(meansddf)<-c("Genotype", "TimePoint", "N","Parameter" ,"sd" , "se","ci","Phase")
    k<-ggplot(meansddf,aes(x=TimePoint,y=Parameter,colour=Genotype,group=Genotype))+
      geom_rect(data=dark_phase,inherit.aes = F,aes(xmin=as.numeric(dark_start),xmax=as.numeric(dark_end),ymin=-Inf,ymax=Inf),fill="azure2")+
      geom_errorbar(aes(ymin=Parameter-se,ymax=Parameter+se), width=.1)+  
      geom_line() +
      geom_point()+
      ylab(parameter[i])+
      scale_x_continuous("Time point",breaks=seq(1,length(plot_labels_temp),5),labels=plot_labels_temp[seq(1,length(plot_labels_temp),5)])+
      scale_color_locuszoom()+
      theme_bw()+
      theme(panel.grid = element_blank(),axis.text.x = element_text(angle=90))
    
    plots<-list.append(plots,k)
    
    meansddf<-summarySE(data,measurevar=parameter[i],groupvars=c("Genotype","TimePoint"))
    meansddf<-cbind(meansddf,tpvsp)
    colnames(meansddf)<-c("Genotype", "TimePoint", "N",parameter[i] ,"sd" , "se","ci","Phase")
    
    
    pval<-list()
    sign<-list()
    tukey<-list()
    signtp<-list()
    n<-0
    
    
    for (j in 1:length(unique(data$TimePoint))){
      data.tp<-data[data$TimePoint==j,c("Genotype", parameter[i])]
      res.aov.data.tp<-aov(data.tp[[parameter[i]]] ~ data.tp$Genotype)
      pval[j]<-summary(res.aov.data.tp)[[1]][["Pr(>F)"]][[1]]
      if(pval[j]=="NaN"){
        pval[j]<-1
      }
      
      if(pval[j]<1e-3){
        sign[j]<-"****"
      } else if(pval[j]<5e-3){
        sign[j]<-"***"
      } else if(pval[j]<1e-2){
        sign[j]<-"**"
      } else if(pval[j]<5e-2){
        sign[j]<-"*"
      } else{
        sign[j]<-"."
      }
      
      bartlett.data.tp<-bartlett.test(data.tp[[parameter[i]]] ~ data.tp$Genotype)
      if(bartlett.data.tp$p.value=="NaN"){
        bartlett.data.tp$p.value<-1
      } else{
        
        if (bartlett.data.tp$p.value<0.05){
          write.table(paste("\n\n",parameter[i],"- Timepoint:",j,"- Warning in Bartlett's test",sep=" "),file="./Results/Warnings.txt",append = T,quote=F,row.names=F,col.names=F)
          write.table(paste(parameter[i],"- Timepoint:",j,"- Warning in Bartlett's test",sep=" "),file="./Results/Warnings_summary.txt",append = T,quote=F,row.names=F,col.names=F)
          capture.output(bartlett.data.tp,file="./Results/Warnings.txt",append = T)
        }
        
        aov_residuals.data.tp <- residuals(object = res.aov.data.tp)
        shap.data.tp<-shapiro.test(x = aov_residuals.data.tp )
        if (shap.data.tp$p.value<0.05){
          write.table(paste("\n\n",parameter[i],"- Timepoint:",j,"- Warning in Shapiro-Wilk's test",sep=" "),file="./Results/Warnings.txt",append = T,quote=F,row.names=F,col.names=F)
          write.table(paste(parameter[i],"- Timepoint:",j,"- Warning in Shapiro-Wilk's test",sep=" "),file="./Results/Warnings_summary.txt",append = T,quote=F,row.names=F,col.names=F)
          capture.output(shap.data.tp,file="./Results/Warnings.txt",append=T)
        }
      }
      
      if(length(genotype)>2 && pval[j]<0.05){
        n<-n+1
        
        tukey[n]<-TukeyHSD(res.aov.data.tp)
        
        signtp[n]<-j
      }
      
      
      
      
    }
    
    
    meansddf$pval<-as.numeric(rep(pval,length(unique(genotype))))
    meansddf$sign<-as.character(rep(sign,length(unique(genotype))))
    colnames(meansddf)<-c("Genotype", "TimePoint", "N",parameter[i],"sd" , "se","ci","Phase","ANOVA_pValue","Sign")
    meansddf<-meansddf %>% mutate_at(c(parameter[[i]],"sd","se","ci"),round,digits=3)
    meansddf<-meansddf[order(meansddf$TimePoint),]
    meansddf$ANOVA_pValue <- formatC(meansddf$ANOVA_pValue, format="e", digits=3)
    write.table(as.matrix(meansddf),file=paste0(folder,parameter[i],"_statistics.txt"),quote=F,col.names=TRUE,row.names=F,sep="\t")
    
    if(length(genotype)>2 && n>0){
      
      for (x in 1:n){
        write.table(paste0("\n### Tukey multiple pairwise comparisons on timepoint ",signtp[x]," ###"),file=paste0(folder,parameter[i],"_statistics.txt"),append=T,quote=F,row.names=F,col.names=F)
        capture.output(tukey[x],file=paste0(folder,parameter[i],"_statistics.txt"),append=T)
        
      }
    }
    
    sign<-meansddf[as.numeric(meansddf$ANOVA_pValue)<=0.05,]
    sign<-sign[order(sign$TimePoint),]
    
    if(nrow(sign)>0){
      write.table(paste0("\n\n### ",parameter[i]," ###\n"),"./Results/Significant_5e-2_meanovertimepoints.txt",quote=FALSE,append=T,row.names=F,col.names=F)
      write.table(sign,"./Results/Significant_5e-2_meanovertimepoints.txt",quote=FALSE,append=T,row.names=F,col.names=T,sep="\t")
    }
    
  }
  return(plots)
}

# processFile() takes the input raw data file and defines how many lines to skip before start reading, 
# which is one line after finding "Date" in the first column
processFile <- function(filepath) {
  con <- file(filepath, "r")
  keep_on<-TRUE
  i<-0
  while(keep_on){
    line <- strsplit(readLines(con, n = 1),",")
    i<-i+1
    if(!is.na(line[[1]][1]) && line[[1]][1]=="Date")  {
      keep_on<-FALSE
      i<-i+1
    }
  }
  close(con)
  return(i)
}

#Reading the input file or exiting the program

option_list<-list(make_option(c("-f","--file"),type="character",default=NULL,help="Input data file",metavar="character"),
                 make_option(c("-w","--weight"),type="character",default=NULL,help="Weight data file",metavar="character"),
                 make_option(c("-s","--start"),type="integer",default=NULL,help="Start day for analysis",metavar="integer"),
                 make_option(c("-e","--end"),type="integer",default=NULL,help="End day for analysis",metavar="integer"),
                 make_option(c("-l","--wheel"),type="character",action="store_true",default=NULL,help="Analyse wheel data",metavar="character"),
                 make_option(c("-d","--dfw"),type="character",action="store_true",default=NULL,help="Analyse drink/feed/weight data",metavar="character"),
                 make_option(c("-a","--activity"),type="character",action="store_true",default=NULL,help="Analyse activity data",metavar="character"),
                 make_option(c("-c","--calorimetry"),type="character",action="store_true",default=NULL,help="Analyse calorimetry data",metavar="character"),
                 make_option(c("-r","--reference"),type="character",action="store_true",default=NULL,help="Analyse reference values",metavar="character"),
                 make_option(c("-t","--width"),type="integer",default=7,help="Width of plots in pdf format",metavar="integer"),
                 make_option(c("-i","--height"),type="integer",default=7,help="Height of plots in pdf format",metavar="integer"),
                 make_option(c("-v","--levels"),type="character",default=NULL,help="Genotype order for results",metavar="character"),
                 make_option(c("-n","--numbers"),type="character",action="store_true",default=NULL,help="Convert timepoints to numbers or use date_time instead",metavar="character"))
opt_parser<-OptionParser(option_list=option_list)
opt<-parse_args(opt_parser)

if(is.null(opt$file) || is.null(opt$weight)){
  print_help(opt_parser)
  stop("Input files missing.\n",call.=FALSE)
}

skip_lines<-processFile(opt$file)
tempdata<-read.csv(opt$file,header=F,skip=skip_lines)
tempdata<-tempdata[,-ncol(tempdata)]
colnames(tempdata)<-c("Date","Time","Animal","Box","TempC","HumC","LightC","S.Flow","Ref.O2","Ref.CO2","Flow","Temp","O2","CO2",
                     "dO2","dCO2","VO2.1.","VO2.2.","VO2.3.","VCO2.1.","VCO2.2.","VCO2.3.","RER","H.1.","H.2.","H.3.",
                     "XT.YT","XT","XA","XF","YT","YA","YF","Z","CenT","CenA","CenF","PerT","PerA","PerF","DistK","DistD","Speed",
                     "Drink","Feed","Weight","Right","Left","SumR.L","SumTime","SumRuns","MaxSpeed","AvgSpeed","MaxLen")
weight_data<-read.table(opt$weight,header=T)
genotype<-unique(weight_data$Genotype)


# Assigning phases
for(i in 1:nrow(tempdata)){
  if(tempdata$LightC[i]>0){
    tempdata$Phase[i]<-"LIGHT"
  } else{
    tempdata$Phase[i]<-"DARK"
  }
}

# Assigning days
dates<-as.data.frame(cbind(as.character(unique(tempdata$Date)),1:length(unique(tempdata$Date))))
colnames(dates)<-c("Date","Day")
tempdata$Day<-as.numeric(unlist(lapply(tempdata$Date,function(x) dates$Day[match(x,dates$Date)])))

# Assigning genotypes and starting weights
tempdata$Genotype<-unlist(lapply(tempdata$Box,function(x) as.character(weight_data$Genotype[match(x,weight_data$Box)])))
tempdata$WgStart<-unlist(lapply(tempdata$Box,function(x) weight_data$WgStart[match(x,weight_data$Box)]))

#Selection of days from input interval and tare of cummulative parameters if necessary.
if(!is.null(opt$start) || !is.null(opt$end)){
  interval<-"yes"
  start<-opt$start
  end<-opt$end
  data<-tempdata[tempdata$Day>=start & tempdata$Day<=end,]
  
  data$TimePoint<-rep(1:nrow(unique(data[c("Date","Time")])),length(unique(data$Box)))
  
  times<-length(unique(data$TimePoint))
  drinktare<-list()
  feedtare<-list()
  distktare<-list()
  
  Box<-unique(data$Box)
  for(i in 1:length(unique(data$Box))){
    drinktare[i]<-data$Drink[data$Box==Box[i]][1]
    feedtare[i]<-data$Feed[data$Box==Box[i]][1]
    distktare[i]<-data$DistK[data$Box==Box[i]][1]
  }
  
  tare<-as.data.frame(cbind(unlist(Box),unlist(drinktare),unlist(feedtare),unlist(distktare)))
  colnames(tare)<-c("Box","drinktare","feedtare","distktare")
  
  data$DrinkT0<-unlist(lapply(data$Box,function(x) tare$drinktare[match(x,tare$Box)]))
  data$TareDrink<-data$Drink-data$DrinkT0
  
  data$FeedT0<-unlist(lapply(data$Box,function(x) tare$feedtare[match(x,tare$Box)]))
  data$TareFeed<-data$Feed-data$FeedT0
  
  data$DistkT0<-unlist(lapply(data$Box,function(x) tare$distktare[match(x,tare$Box)]))
  data$TareDistK<-data$DistK-data$DistkT0
  
} else {
  interval<-"no"
  data<-tempdata
}

# Extracting dark periods for plots
#dark_start <- list()
#dark_end <- list()
#j <- 1
#k <- 1

#simp_data=data[data$Box==1,]
#dark_start[[j]]<-(simp_data$TimePoint[1])
#for(i in 2:nrow(simp_data)){
#  if(simp_data$Phase[i]=="LIGHT"){
#    next()
#  } else if(i>1 && simp_data$Phase[i-1]=="DARK"){
#    if(i<nrow(simp_data) && simp_data$Phase[i+1]=="LIGHT"){
#      dark_end[[k]]<-simp_data$TimePoint[i]
#      k<-k+1
#    } else{
#      next()
#    }
#  } else{
#    j<-j+1
#    dark_start[[j]]<-simp_data$TimePoint[i]
#  }
#}
#dark_end[[k]]<-simp_data$TimePoint[i]
#dark_phase<-as.data.frame(cbind(dark_start,dark_end))



dark <- list()
dark <- unique(split(data,data$Phase)$DARK[,"TimePoint"])

dark_start <- list()
dark_end <- list()
j <- 1

dark_start[[1]] <- dark[[1]]
for(i in 2:length(dark)){
  if(dark[[i]]==dark[[i-1]]+1){
    next() 
  }else{
    dark_end[[j]]<-dark[[i-1]]
    dark_start[[j+1]]<-dark[[i]]
    j<-j+1
  }
}
dark_end[[j]]<-dark[[i]]

dark_phase<-as.data.frame(cbind(dark_start,dark_end))

# Setting the desired genotype order for results if required
if(!is.null(opt$levels)){
  data$Genotype <- factor(data$Genotype, levels=unlist(strsplit(opt$levels,",")))
}

# Adding cummulative data for specific parameters
cum<-c("XT.YT","XT","XA","XF","YT","YA","YF","Z","CenT","CenA","CenF","PerT","PerA","PerF","Right","Left","SumR.L","SumTime","SumRuns")
for(i in 1:length(cum)){
  name<-paste("Cum",cum[i],sep="")
  data[[name]]<-ave(data[[cum[i]]],data$Box,FUN=cumsum)
}

# Replacing non existent data for AvgSpeed by 0.00
data$AvgSpeed<-as.numeric(as.character(data$AvgSpeed))
data$AvgSpeed[is.na(data$AvgSpeed)] <- 0

# Creating the directory to store the results and writing the formatted dataset
dir.create("./Results")
write.table(data,"./Results/Formatted_dataset.txt",quote=F,col.names=TRUE,row.names=F)

#Creating files for only significant pvalues
write.table("PARAMETERS WITH STATISTICALLY SIGNIFICANT DIFFERENCES IN BOXPLOTS (pval<0.05)\n","./Results/Significant_5e-2_boxplots.txt",quote=F,col.names=F,row.names=F)
write.table("PARAMETERS AND TIMEPOINTS WITH STATISTICALLY SIGNIFICANT DIFFERENCES IN ANOVA (pval<0.05)","./Results/Significant_5e-2_meanovertimepoints.txt",quote=F,col.names=F,row.names=F)
write.table("WARNINGS BEFORE CONSIDERING ANOVA RESULS","./Results/Warnings.txt",quote=F,col.names=F,row.names=F)
write.table("WARNINGS BEFORE CONSIDERING ANOVA RESULS - SUMMARY","./Results/Warnings_summary.txt",quote=F,col.names=F,row.names=F)

#Homogeneity of variance test and t.test for weight data
weight_data$Genotype<-factor(weight_data$Genotype)
res.aovwg<-aov(weight_data$WgStart ~ weight_data$Genotype)
write.table("#######################\n#### ANOVA RESULTS ####\n#######################\n\n",file="./Results/ANOVA_weight.txt",quote=F,row.names=F,col.names=F)
capture.output(summary(res.aovwg),file="./Results/ANOVA_weight.txt",append=T)

if(length(genotype)>2){
  tukeywg<-TukeyHSD(res.aovwg)
  write.table("\n\n#############################################\n#### Tukey multiple pairwise comparisons ####\n#############################################\n\n",file="./Results/ANOVA_weight.txt",append=T,quote=F,row.names=F,col.names=F)
  capture.output(tukeywg,file="./Results/ANOVA_weight.txt",append=T)
}

barltestwg<-bartlett.test(weight_data$WgStart ~ weight_data$Genotype)
write.table("\n\n######################################################\n### Barlett's test ###\n######################################################\n\n",file="./Results/ANOVA_weight.txt",append=T,quote=F,row.names=F,col.names=F)
capture.output(barltestwg,file="./Results/ANOVA_weight.txt",append=T)

aov_residualswg <- residuals(object = res.aovwg)
shapwg<-shapiro.test(x = aov_residualswg )
write.table("\n\n######################################################\n### Shapiro-Wilk's test for normality of residuals ###\n######################################################\n\n",file="./Results/ANOVA_weight.txt",append=T,quote=F,row.names=F,col.names=F)
capture.output(shapwg,file="./Results/ANOVA_weight.txt",append=T)

# Assigning labels
if(!is.null(opt$numbers)){
  plot_labels_temp<-unique(data$TimePoint) 
} else{
  plot_labels_temp <- unique(paste0(data$Date,"_",data$Time))
}

# Plotting reference values
if(!is.null(opt$reference)){
  cat("\nGenerating reference plots...\n")
  reference<-c("TempC","HumC","LightC","S.Flow","Ref.O2","Ref.CO2")
  
  pdf("./Results/Reference_plots.pdf",width = opt$width,height = opt$height)
  for (i in 1:length(reference)){
    p<-ggplot(data=data[data$Box=="1",],aes_string(x="TimePoint",y=reference[i],group="Animal"))+
      geom_rect(data=dark_phase,mapping=aes(xmin=as.numeric(dark_phase$dark_start),xmax=as.numeric(dark_phase$dark_end),ymin=-Inf,ymax=+Inf),inherit.aes = F,fill="azure2")+
      geom_point(size=0.1)+
      geom_line()+
      scale_x_continuous("Time point",breaks=seq(1,length(plot_labels_temp),5),labels=plot_labels_temp[seq(1,length(plot_labels_temp),5)])+
      theme_classic()+
      theme(axis.text.x = element_text(angle=90))
    g<-ggplot(data=data[data$Box=="1",],aes_string("Phase",reference[i]))+geom_boxplot(fill=c("grey36","white"))+
      stat_summary(fun.y=mean,geom="point",color="darkred")+
      stat_summary(fun.y=mean,geom="text",color="black",vjust=-0.7,aes(label=round(..y..,digits=2)))+
      theme_classic()
    invisible(print(p))
    invisible(print(g))
  }
  graphics.off()
}

# Plotting calorimetry values
if(!is.null(opt$calorimetry)){
  dir.create("./Results/Calorimetry")
  folder<-"./Results/Calorimetry/"
  cat("\nGenerating calorimetry plots...\n")
  
  write.table("\n\n###############################\n######### CALORIMETRY #########\n###############################","./Results/Warnings.txt",append=T,quote=F,col.names=F,row.names=F)
  write.table("\n\n###############################\n######### CALORIMETRY #########\n###############################","./Results/Warnings_summary.txt",append=T,quote=F,col.names=F,row.names=F)
  calo<-c("Flow","Temp","O2","CO2","dO2","dCO2","VO2.1.","VO2.2.","VO2.3.","VCO2.1.","VCO2.2.","VCO2.3.","RER","H.1.","H.2.","H.3.")
  
  pdf("./Results/Calorimetry/Calorimetry_plots.pdf",width = opt$width,height = opt$height)
  
  indvplots<-plot_ind(calo,data,plot_labels_temp)
  for(i in 1:length(indvplots)){
    invisible(print(indvplots[[i]]))
  }
  
  calo<-calo[2:length(calo)]
  
  bxplts<-plot_bxplt(calo,data,pval_wg,folder,dark_phase)
  for(i in 1:length(bxplts)){
    invisible(print(bxplts[[i]]))
  }
  
  meanplots<-plot_mean(calo,data,folder,dark_phase,plot_labels_temp)
  for(i in 1:length(meanplots)){
    invisible(print(meanplots[[i]]))
  }    
  
  invisible(print(ggplot(data[data$SumR.L>0,],aes_string(x="SumR.L",y="RER",colour="Genotype",group="Animal"))+
                    geom_point(size=1)+
                    facet_wrap(~Phase)+
                    theme_bw()+
                    theme(strip.background = element_rect(fill="white"),strip.text = element_text(face="bold"))+
                    scale_color_locuszoom()))
  invisible(print(ggplot(data[data$SumR.L>0,],aes_string(x="SumR.L",y="RER",colour="Phase",shape="Genotype"))+
                    geom_point(size=1)+
                    facet_wrap(~Animal)+
                    theme_bw()+
                    theme(strip.background = element_rect(fill="white"),strip.text = element_text(face="bold"))+
                    scale_color_locuszoom()))
  
  graphics.off()
  
}

# Plotting activity values
if(!is.null(opt$activity)){
  dir.create("./Results/Activity")
  folder<-"./Results/Activity/"
  cat("\nGenerating activity plots...\n")
  
  write.table("\n\n###############################\n#########   ACTIVITY  #########\n###############################","./Results/Warnings.txt",append=T,quote=F,col.names=F,row.names=F)
  write.table("\n\n###############################\n#########   ACTIVITY  #########\n###############################","./Results/Warnings_summary.txt",append=T,quote=F,col.names=F,row.names=F)
  
  activity<-c("XT.YT","XT","XA","XF","YT","YA","YF","Z","CenT","CenA","CenF","PerT","PerA","PerF")
  if(interval=="no"){
    cumactivity<-c("CumXT.YT","CumXT","CumXA","CumXF","CumYT","CumYA","CumYF","CumZ","CumCenT","CumCenA","CumCenF","CumPerT","CumPerA","CumPerF","DistK")
  } else{
    cumactivity<-c("CumXT.YT","CumXT","CumXA","CumXF","CumYT","CumYA","CumYF","CumZ","CumCenT","CumCenA","CumCenF","CumPerT","CumPerA","CumPerF","TareDistK")
  }
  othersactivity<-c("DistD","Speed")
  
  pdf("./Results/Activity/Activity_plots.pdf",width = opt$width,height = opt$height)
  
  indvplots<-plot_ind(activity,data,plot_labels_temp)
  for(i in 1:length(indvplots)){
    invisible(print(indvplots[[i]]))
  }
  
  bxplts<-plot_bxplt(activity,data,pval_wg,folder,dark_phase)
  for(i in 1:length(bxplts)){
    invisible(print(bxplts[[i]]))
  }
  
  meanplots<-plot_mean(activity,data,folder,dark_phase,plot_labels_temp)
  for(i in 1:length(meanplots)){
    invisible(print(meanplots[[i]]))
  }    
  
  indvplots<-plot_ind(cumactivity,data,plot_labels_temp)
  for(i in 1:length(indvplots)){
    invisible(print(indvplots[[i]]))
  }
  
  meanplots<-plot_mean(cumactivity,data,folder,dark_phase,plot_labels_temp)
  for(i in 1:length(meanplots)){
    invisible(print(meanplots[[i]]))
  }  
  
  indvplots<-plot_ind(othersactivity,data,plot_labels_temp)
  for(i in 1:length(indvplots)){
    invisible(print(indvplots[[i]]))
  }
  
  bxplts<-plot_bxplt(othersactivity,data,pval_wg,folder,dark_phase)
  for(i in 1:length(bxplts)){
    invisible(print(bxplts[[i]]))
  }
  
  meanplots<-plot_mean(othersactivity,data,folder,dark_phase,plot_labels_temp)
  for(i in 1:length(meanplots)){
    invisible(print(meanplots[[i]]))
  } 
  
  graphics.off()
  
}

# Plotting DFW values
if(!is.null(opt$dfw)){
  dir.create("./Results/DFW")
  folder<-"./Results/DFW/"
  cat("\nGenerating DFW plots...\n")
  
  write.table("\n\n###############################\n#########     DFW     #########\n###############################","./Results/Warnings.txt",append=T,quote=F,col.names=F,row.names=F)
  write.table("\n\n###############################\n#########     DFW     #########\n###############################","./Results/Warnings_summary.txt",append=T,quote=F,col.names=F,row.names=F)
  
  if(interval=="no"){
    dfw<-c("Drink","Feed","Weight")
  } else{
    dfw<-c("TareDrink","TareFeed","Weight")
  }
  
  pdf("./Results/DFW/DFW_plots.pdf",width = opt$width,height = opt$height)
  
  indvplots<-plot_ind(dfw,data,plot_labels_temp)
  for(i in 1:length(indvplots)){
    invisible(print(indvplots[[i]]))
  }
  
  meanplots<-plot_mean(dfw,data,folder,dark_phase,plot_labels_temp)
  for(i in 1:length(meanplots)){
    invisible(print(meanplots[[i]]))
  } 
  
  graphics.off()
  
}

# Plotting wheel values
if(!is.null(opt$wheel)){
  dir.create("./Results/Wheels")
  folder<-"./Results/Wheels/"
  cat("\nGenerating wheel plots...\n")
  
  write.table("\n\n###############################\n#########   WHEELS   ##########\n###############################","./Results/Warnings.txt",append=T,quote=F,col.names=F,row.names=F)
  write.table("\n\n###############################\n#########   WHEELS   ##########\n###############################","./Results/Warnings_summary.txt",append=T,quote=F,col.names=F,row.names=F)
  
  wheel<-c("Right","Left","SumR.L","SumTime","SumRuns")
  cumwheel<-c("CumRight","CumLeft","CumSumR.L","CumSumTime","CumSumRuns")
  otherswheel<-c("MaxSpeed","AvgSpeed","MaxLen")
  
  pdf("./Results/Wheels/wheel_plots.pdf",width = opt$width,height = opt$height)
  
  indvplots<-plot_ind(wheel,data,plot_labels_temp)
  for(i in 1:length(indvplots)){
    invisible(print(indvplots[[i]]))
  }
  
  bxplts<-plot_bxplt(wheel,data,pval_wg,folder,dark_phase)
  for(i in 1:length(bxplts)){
    invisible(print(bxplts[[i]]))
  }
  
  meanplots<-plot_mean(wheel,data,folder,dark_phase,plot_labels_temp)
  for(i in 1:length(meanplots)){
    invisible(print(meanplots[[i]]))
  } 
  
  indvplots<-plot_ind(cumwheel,data,plot_labels_temp)
  for(i in 1:length(indvplots)){
    invisible(print(indvplots[[i]]))
  }
  
  meanplots<-plot_mean(cumwheel,data,folder,dark_phase,plot_labels_temp)
  for(i in 1:length(meanplots)){
    invisible(print(meanplots[[i]]))
  } 
  
  indvplots<-plot_ind(otherswheel,data,plot_labels_temp)
  for(i in 1:length(indvplots)){
    invisible(print(indvplots[[i]]))
  }
  
  bxplts<-plot_bxplt(otherswheel,data,pval_wg,folder,dark_phase)
  for(i in 1:length(bxplts)){
    invisible(print(bxplts[[i]]))
  }
  
  meanplots<-plot_mean(otherswheel,data,folder,dark_phase,plot_labels_temp)
  for(i in 1:length(meanplots)){
    invisible(print(meanplots[[i]]))
  } 
  graphics.off()
  
}

cat("\nDone.\n")




