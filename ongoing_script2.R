### Analysis of TSE Phenomaster data

### Author: Cristina Vicente Garcia - cvicgar@gmail.com
### Date: June 2022
### Description: This script takes a the raw output .csv file of a TSE Phenomaster experiment, a metadata file with information 
### regarding the weight of the animals used, and returns plots and statistical analyses of the desired parameters. It also  
### allows to modify plot sizes.

#!/usr/bin/env Rscript

## LIBRARIES
# Installing and loading the required libraries

if (!("ggplot2" %in% installed.packages())) { 
  BiocManager::install("ggplot2", dependencies=TRUE)
}
if (!("ggpubr" %in% installed.packages())) { 
  BiocManager::install("ggpubr", dependencies=TRUE)
}
if (!("reshape2" %in% installed.packages())) { 
  BiocManager::install("reshape2", dependencies=TRUE)
}
if (!("Rmisc" %in% installed.packages())) { 
  BiocManager::install("Rmisc", dependencies=TRUE)
}
if (!("rlist" %in% installed.packages())) { 
  BiocManager::install("rlist", dependencies=TRUE)
}
if (!("optparse" %in% installed.packages())) { 
  BiocManager::install("optparse", dependencies=TRUE)
}
if (!("dplyr" %in% installed.packages())) { 
  BiocManager::install("dplyr", dependencies=TRUE)
}
if (!("ggsci" %in% installed.packages())) { 
  BiocManager::install("ggsci", dependencies=TRUE)
}
if (!("car" %in% installed.packages())) { 
  BiocManager::install("car", dependencies=TRUE)
}
if (!("multcomp" %in% installed.packages())) { 
  BiocManager::install("multcomp", dependencies=TRUE)
}

library("ggplot2")
library("ggpubr")
library("reshape2")
library("Rmisc")
library("rlist")
library("optparse")
library("dplyr")
library("ggsci")
library("car")
library("multcomp")

  
## DEFINING FUNCTIONS

# processFile() takes the input raw data file and returns how many lines to skip before start reading, which is one line 
# after finding "Date" in the first column.
processFile <- function(filepath) {
  con <- file(filepath, "r")
  keep_on <- TRUE
  i <- 0
  while(keep_on){
    line <- strsplit(readLines(con, n=1), ",")
    i <- i+1
    if(!is.na(line[[1]][1]) && line[[1]][1]=="Date")  {
      keep_on <- FALSE
      i <- i+1
    }
  }
  close(con)
  return(i)
}

# assumptions() takes a statistical object (object) and performs the Shapiro-Wilk's test to check for normality of residuals.
# It also takes raw data (dataframe) and the variables to test (form) to perform the Bartlett's test to check for  
# homogeneity of variances. Significant results are included in the Warning files, indicating which statistical analysis
# (caption) needs to be carefully considered.
assumptions <- function(form, dataframe, object, caption, plot_type){
  # Bartlett's test to check for homogeneity of variances.
  bartlett <- bartlett.test(formula=as.formula(form), data=dataframe)
  if (is.na(bartlett$p.value) || bartlett$p.value<0.05){
    if(plot_type=="meanplot"){
      write.table(paste0("\n\n", caption, "- Bartlett's test"), file="./Results/Warnings.txt", append=T, quote=F, 
                  row.names=F, col.names=F)
      write.table(paste0(caption, "- Bartlett's test"), file="./Results/Warnings_summary.txt", append=T, quote=F, 
                  row.names=F, col.names=F)
    } else if(plot_type=="boxplot"){
      write.table(paste0("\n\n", caption, " - Bartlett's test"), file="./Results/Warnings.txt", append=T, quote=F, 
                  row.names=F, col.names=F)
      write.table(paste0(caption, " - Bartlett's test"), file="./Results/Warnings_summary.txt", append=T, quote=F, 
                  row.names=F, col.names=F)
    }
  }
  
  if(!is.na(bartlett$p.value) && bartlett$p.value<0.05){
    capture.output(bartlett, file="./Results/Warnings.txt", append=T)
  } else if (is.na(bartlett$p.value)){
    test <- "Not possible to perform ANOVA testing"
    write.table(test, file="./Results/Warnings.txt", append=T, quote=F, row.names=F, col.names=F)
  }
  
  # Shapiro-Wilk's test to check for normality of residuals.
  aov_residuals <- residuals(object)
  shap <- tryCatch(shapiro.test(x=aov_residuals), error = function(e){invisible(e); shap <- 0})
  if(shap==0 || shap$p.value<0.05){
    if(plot_type=="meanplot"){
      write.table(paste0("\n\n", caption, "- Shapiro-Wilk's test"), file="./Results/Warnings.txt", append=T, quote=F, 
                row.names=F, col.names=F)
      write.table(paste0(caption, "- Shapiro-Wilk's test"), file="./Results/Warnings_summary.txt", append=T, quote=F, 
                row.names=F, col.names=F)
    } else if(plot_type=="boxplot"){
      write.table(paste0("\n\n", caption, " - Shapiro-Wilk's test"), file="./Results/Warnings.txt", append=T, quote=F, 
                row.names=F, col.names=F)
      write.table(paste0(caption, " - Shapiro-Wilk's test"), file="./Results/Warnings_summary.txt", append=T, quote=F, 
                row.names=F, col.names=F)
    }
  }
  
  if(shap!=0 && shap<0.05){
    capture.output(shap, file="./Results/Warnings.txt", append=T)
  } else if (shap==0){
    test <- "Not possible to perform ANOVA testing"
    write.table(test, file="./Results/Warnings.txt", append=T, quote=F, row.names=F, col.names=F)
  }
}

# statistics() takes data in the form of a dataframe and performs a one-way ANOVA statistical test. If more than two
# genotypes are analysed, it also performs Tukey's HSD post-hoc testing. It returns a list with ANOVA and Tukey results 
# as well as a variable, "sign" that indicates if statistically significant results were obtained. This functions calls
# assumptions() to test ANOVA assumptions.
statistics <- function(dataframe, plot_type){
  # Definition of some variables.
  genotypes <- length(unique(dataframe$Genotype))
  if(plot_type=="boxplot"){
    param <- names(dataframe)[5]
    phase <- unique(dataframe$Phase)
    caption <- paste(param, "- Phase:", phase, "- Warning in ANOVA", sep=" ")
  } else if(plot_type=="meanplot"){
    param <- names(dataframe)[2]
    tp <- as.numeric(unique(dataframe$Timepoint))
    caption <- paste(param, "- Timepoint:", tp, "- Warning in ANOVA", sep=" ")
  }
  
  # One-way ANOVA.
  eq <- paste0(param, "~Genotype")
  anova <- glm(eq, data=dataframe, family=gaussian(link="identity"))
  
  # Checking statistical significance and performing Tukey HSD post-hoc testing if applicable.
  sign <- "no"
  tukey <- "no"
  res <- sum(summary(anova)$coeff[-1,ncol(summary(anova)$coeff)]<0.05)
  if(!is.na(res) && res>0){
    sign <- "yes"
    if(genotypes>2){
      tukey <- glht(anova,mcp(Genotype="Tukey"))
    }
  }

  # Testing ANOVA assumptions with the assumptions() function.
  assumptions(eq, dataframe, anova, caption, plot_type)
  
  # Returning results.
  results <- list(anova, tukey, sign)
  return(results)
}

# special_statistics() performs additional statistical tests in special parameters that may be influenced by the animal 
# weight. It functions similarly to the statistics() function.
special_statistics <- function(dataframe, plot_type){
  # Definition of some variables.
  genotypes <- length(unique(dataframe$Genotype))
  parameter <- names(dataframe)[5]
  phase <- unique(dataframe$Phase)
  
  # ANCOVA (NO INTERACTIONS) analysis: this test aims at determining whether the genotype has any statistical
  # significant impact on the parameter under study while controlling for the covariate "weight", that is, if the 
  # potential differences found in a specific parameter are due to their different genotype or to differences in weight.
  eq <- paste0(parameter, "~WgStart+Genotype")
  ancova <- glm(eq, data=dataframe, family=gaussian(link="identity"))

  # Checking statistical significance and performing Tukey HSD post-hoc testing if applicable.
  sign <- "no"
  tukey <- "no"
  res <- sum(summary(ancova)$coeff[-(1:2),ncol(summary(ancova)$coeff)]<0.05)
  if(!is.na(res) && res>0){
    sign <- "yes"
    if(genotypes>2){
      tukey <- glht(ancova,mcp(Genotype="Tukey"))
    }
  }
    
  # Testing ANCOVA assumptions with the assumptions() function.
  caption <- paste(parameter, "- Phase:", phase, "- Warning in ANCOVA (without interaction)", sep=" ")
  form <- paste0(parameter, "~Genotype")
  assumptions(form, dataframe, ancova, caption, plot_type)
    
  # ANOVA WITH INTERACTIONS OR ANCOVA FOR NON-PARALLEL SLOPES analysis: this test determines whether weight and 
  # genotypes interact influencing a particular parameter.
  eq2 <- paste0(parameter, "~WgStart+Genotype+WgStart:Genotype")
  ancova_int <- glm(eq2, data=dataframe, family=gaussian(link="identity"))

  # Checking statistical significance and performing Tukey HSD post-hoc testing if applicable.
  sign_int <- "no"
  tukey_int <- "no"
  res_int <- sum(summary(ancova_int)$coeff[-(1:2),ncol(summary(ancova_int)$coeff)]<0.05)
  if(!is.na(res_int) && res_int>0){
    sign_int <- "yes"
    if(genotypes>2){
      tukey_int <- glht(ancova_int,mcp(Genotype="Tukey"))
    }
  }
    
  # Testing assumptions with the assumptions() function.
  caption_int <- paste(parameter, "- Phase:", phase, "- Warning in ANOVA with interactions", sep=" ")
  assumptions(form, dataframe, ancova_int, caption_int, plot_type)

  # Returning results.
  results_special <- list(ancova, tukey, sign, ancova_int, tukey_int, sign_int)
  return(results_special)
}

# plot_ind() takes as inputs the whole dataset, the parameter to plot, as well as the dark and light periods and plot labels,  
# and returns two types of plots: 1) individual profiles of that parameter over time, highlighting the dark and light phases  
# and grouped by genotype, resulting in as many panels as different genotypes, 2) individual profiles of that parameter over   
# time, highlighting the dark and light phases, resulting in as many panels as animals in the experiment.
plot_ind <- function(parameter, data, dark_phase, plot_labels){
  plots <- list()
  genotype_no <- length(unique(data$Genotype))
  animal_no <- length(unique(data$Animal))
  
  for (i in 1:length(parameter)){
    p <- ggplot(data=data, aes_string(x="TimePoint", y=parameter[i], colour=as.factor(data$Genotype), group="Animal")) +
      facet_wrap(~Genotype) +
      geom_rect(data=dark_phase, mapping=aes(xmin=rep(unlist(dark_phase$dark_start), genotype_no),
            xmax=rep(unlist(dark_phase$dark_end), genotype_no), ymin=-Inf, ymax=+Inf), inherit.aes=F, fill="azure2") +
      geom_point(size=0.1, show.legend=F) +
      geom_line(show.legend=F) +
      scale_color_locuszoom() +
      scale_x_continuous("Time point", breaks=seq(1, length(plot_labels), 10),
            labels=plot_labels[seq(1, length(plot_labels), 10)]) +
      theme_bw() +
      theme(strip.background=element_rect(fill="white"), strip.text=element_text(face="bold"), panel.grid=element_blank(),
            axis.text.x=element_text(angle=90)) +
      ylab(names(parameter[i]))
    g <- ggplot(data=data, aes_string(x="TimePoint", y=parameter[i], colour="Genotype", group="Animal")) +
      geom_rect(data=dark_phase, mapping=aes(xmin=rep(unlist(dark_phase$dark_start), animal_no), 
            xmax=rep(unlist(dark_phase$dark_end), animal_no), ymin=-Inf, ymax=+Inf), inherit.aes=F, fill="azure2") +
      geom_point(size=0.1) +
      geom_line() +
      scale_color_locuszoom() +
      facet_wrap(~Animal, nrow=genotype_no) +
      scale_x_continuous("Time point", breaks=seq(1, length(plot_labels), 20), 
            labels=plot_labels[seq(1, length(plot_labels), 20)]) +
      theme_bw() +
      theme(strip.background=element_rect(fill="white"), strip.text=element_text(face="bold"), panel.grid=element_blank(), 
            axis.text.x=element_text(angle=90)) +
      ylab(names(parameter[i]))
    plots <- list.append(plots, p, g)
  }
  return(plots)
}

# plot_bxplt() takes as inputs the whole dataset and the parameter to plot, as well as the dark and light periods and the 
# folder to save the results, and returns 1) boxplots depicting results grouped by phase, 2) boxplots depicting results 
# grouped by phase and genotype, 3) statistical analyses of all results, including warnings in case assumptions for the 
# corresponding analyses are not met.
plot_bxplt <- function(parameter, data, pval_wg, folder, dark_phase){
  # Appending "In boxplots" in the warning files, which contain information regarding the statistical tests
  # that do not meet the underlying assumptions of homogeneity of variances and normality.
  write.table("\n######### In boxplots #########", "./Results/Warnings.txt", append=T, quote=F, col.names=F, row.names=F)
  write.table("\n######### In boxplots #########", "./Results/Warnings_summary.txt", append=T, quote=F, col.names=F, 
              row.names=F)
  # Definition of some variables.
  special <- c ("VO2.3.", "VCO2.3.", "H.3.")
  plots <- list()
  genotype_no <- length(unique(data$Genotype))
  
  # For each of the parameters to plot: 1) calculation of mean data during the light and dark cycles, 2) plotting.
  for (i in 1:length(parameter)){
    
    # Calculation of mean data during the light and dark cycles for plots.
    meandata <- unique(data[c("Genotype", "Animal", "WgStart")])
    meandata <- cbind(meandata, "LIGHT"=rep(0, nrow(meandata)), "DARK"=rep(0, nrow(meandata)))
    for(y in 1:nrow(meandata)){
      meandata[y,"DARK"] <- mean(data[data$Animal==as.character(meandata[y,2]) & data$Phase=="DARK",][[parameter[i]]])
      meandata[y,"LIGHT"] <- mean(data[data$Animal==as.character(meandata[y,2]) & data$Phase=="LIGHT",][[parameter[i]]])
    }
    dat.m <- melt(meandata, id.vars=c("Genotype", "Animal", "WgStart"), measure.vars=c("DARK", "LIGHT"))
    colnames(dat.m) <- c("Genotype", "Animal", "WgStart", "Phase", parameter[i])
    
    # Boxplots grouped by genotype and phase. One value is plotted per animal and phase: the mean value of the corresponding
    # paramater during each phase. These data will be used for statistical analyses.
    r <- ggplot(data=dat.m, aes_string("Phase", parameter[i], fill="Phase")) +
      geom_boxplot(outlier.shape=NA) +
      scale_fill_manual(values=rep(c("azure2", "white"), genotype_no)) +               
      geom_point(size=1, position=position_jitterdodge((jitter.width=0.3))) +
      facet_wrap(~Genotype) +
      stat_summary(fun.y=mean, geom="point", color="darkred") +
      stat_summary(fun.y=mean, geom="text", color="black", vjust=-1, aes(label=round(..y.., digits=2))) +
      theme_bw() +
      theme(strip.background=element_rect(fill="white"), strip.text=element_text(face="bold"), panel.grid=element_blank(), 
            legend.position="none") +
      ylab(names(parameter[i]))
    
    # Individual boxplots. All values registered at a given phase are plotted.
    s <- ggplot(data=data, aes_string("Phase", parameter[i], fill="Genotype")) +
      geom_boxplot(outlier.shape=NA) +
      geom_point(size=1, color="azure4", position=position_jitterdodge((jitter.width=0.6))) +
      stat_summary(fun.y=mean, geom="point", color="darkred") +
      stat_summary(fun.y=mean, geom="text", color="black", vjust=-1, aes(label=round(..y.., digits=2))) +
      facet_wrap(~Animal, nrow=genotype_no) +
      theme_bw() +
      theme(strip.background=element_rect(fill="white"), strip.text=element_text(face="bold"), panel.grid=element_blank(), 
            legend.position="none") +
      scale_fill_locuszoom() +
      ylab(names(parameter[i]))
    
    plots <- list.append(plots, r, s)
    
    dat.ph <- list("Light"=dat.m[dat.m$Phase=="LIGHT",], "Dark"=dat.m[dat.m$Phase=="DARK",])
    
    # Special parameters VO2.3., VCO2.3., H.3. are first plotted and then subjected to statistical analyses with the 
    # special_statistics() function to check whether potential differences in these parameters are due to differences 
    # in weight and/or in genotype, (ANCOVA), or whether weight and genotypes interact (ANOVA with interactions or 
    # ANCOVA for non-parallel slopes). 
    if(parameter[i] %in% special){
      v <- ggplot(data=dat.m, aes_string("WgStart", parameter[i], colour="Genotype")) +
        geom_point() +
        geom_smooth(method="lm", formula=y~x, se=FALSE) +
        stat_regline_equation(aes(label=paste(..eq.label.., ..rr.label.., sep="-"))) +
        facet_wrap(~Phase) +
        theme_bw() +
        theme(strip.background=element_rect(fill="white"), strip.text=element_text(face="bold")) +
        scale_color_locuszoom() +
        xlab("Starting weight (g)") +
        ylab(names(parameter[i]))
      plots <- list.append(plots, v)
      
      # Statistical results are stored in the corresponding *_ANCOVA_statistics.txt files.
      write.table("######################################################################\n########################### ANCOVA RESULTS ###########################\n######################################################################", 
                  file=paste0(folder, parameter[i], "_ANCOVA_statistics.txt"), quote=F, row.names=F, col.names=F)
      
      for (k in 1:length(dat.ph)){
        phase <- names(dat.ph)[k]
        
        # Calling the special_statistics() function.
        stats_res_special <- special_statistics(dat.ph[[k]], "boxplot")
        names(stats_res_special) <- c("ancova", "tukey_anc", "sign_anc", "ancova_int", "tukey_int", "sign_int")
        
        # ANCOVA (NO INTERACTIONS) results.
        write.table(paste0("\n\n######################################################################\n############################# ", phase, " phase ############################\n######################################################################"), 
                    file=paste0(folder, parameter[i], "_ANCOVA_statistics.txt"), append=T, quote=F, row.names=F, col.names=F)
        write.table("\n###################### ANCOVA (NO INTERACTIONS) ######################\n", file=paste0(folder, parameter[i], "_ANCOVA_statistics.txt"), 
          append=T, quote=F, row.names=F, col.names=F)
        capture.output(summary(stats_res_special$ancova), file=paste0(folder, parameter[i], "_ANCOVA_statistics.txt"), append=T)

        if(stats_res_special$sign_anc=="yes"){
          write.table(paste(parameter[i], "without interaction (ANCOVA)", sep=" "), 
                      "./Results/Significant_5e-2_boxplots.txt", quote=FALSE, append=T, row.names=F, col.names=F)
          if(stats_res_special$tukey_anc!="no"){
            write.table("\n#### Tukey multiple pairwise comparisons ####\n", file=paste0(folder, parameter[i], 
                       "_ANCOVA_statistics.txt"), append=T, quote=F, row.names=F, col.names=F)
            capture.output(summary(stats_res_special$tukey_anc), file=paste0(folder, parameter[i], "_ANCOVA_statistics.txt"), append=T)
          }
        }
        
        # ANOVA WITH INTERACTIONS OR ANCOVA FOR NON-PARALLEL SLOPES results.
        write.table("\n###### ANOVA WITH INTERACTIONS OR ANCOVA FOR NON-PARALLEL SLOPES #####\n", file=paste0(folder, 
                      parameter[i], "_ANCOVA_statistics.txt"), append=T, quote=F, row.names=F, col.names=F)
        capture.output(summary(stats_res_special$ancova_int), file=paste0(folder, parameter[i], "_ANCOVA_statistics.txt"), append=T)
        
        if(stats_res_special$sign_int=="yes"){
          write.table(paste(parameter[i], "with interaction (ANCOVA for non-parallel slopes)", sep=" "), 
                      "./Results/Significant_5e-2_boxplots.txt", quote=FALSE, append=T, row.names=F, col.names=F)
          if(stats_res_special$tukey_int!="no"){
            write.table("\n#### Tukey multiple pairwise comparisons ####\n", file=paste0(folder, parameter[i], 
                       "_ANCOVA_statistics.txt"), append=T, quote=F, row.names=F, col.names=F)
            capture.output(summary(stats_res_special$tukey_int), file=paste0(folder, parameter[i], "_ANCOVA_statistics.txt"), append=T)
          }
        }
        
      }
    }
    
    # For all parameters, including the special ones, one-way ANOVA is performed to compare boxplot averaged data with the
    # statistics() function. 
    # Preparing data.
    meansddat.m <- summarySE(dat.m, measurevar=parameter[[i]], groupvars=c("Genotype", "Phase"))
    colnames(meansddat.m) <- c("Genotype", "Phase", "N", parameter[i], "sd", "se", "ci")
    meansddat.m <- meansddat.m %>% mutate_if(is.numeric, round, digits=3)
    write.table(as.matrix(meansddat.m), file=paste0(folder, parameter[i], "_boxplots_statistics.txt"), 
                quote=F, col.names=TRUE, row.names=F, sep="\t")
    
    # Calling the statistics() function at each phase.
    for (k in 1:length(dat.ph)){
      phase <- names(dat.ph)[k]
      stats_res <- statistics(dat.ph[[k]], "boxplot")
      names(stats_res) <- c("anova","tukey","sign")
      
      write.table(paste0("\n###################\n### ", phase, " phase ###\n###################\n"), 
                  file=paste0(folder, parameter[i], "_boxplots_statistics.txt"), append=T, quote=F, row.names=F, col.names=F)
      capture.output(summary(stats_res$anova), file=paste0(folder, parameter[i], "_boxplots_statistics.txt"), append=T)
      
      if(stats_res$sign=="yes"){
        write.table(paste0(parameter[i], " in the ", phase, " phase"), "./Results/Significant_5e-2_boxplots.txt", 
                    quote=FALSE, append=T, row.names=F, col.names=F, sep="\t")
        if(stats_res$tukey!="no"){
          write.table("\n#### Tukey multiple pairwise comparisons ####\n", file=paste0(folder, parameter[i], 
                      "_boxplots_statistics.txt"), append=T, quote=F, row.names=F, col.names=F)
          capture.output(summary(stats_res$tukey), file=paste0(folder, parameter[i], "_boxplots_statistics.txt"), append=T)
        }
      }
      
    }
    
  }
  return(plots)
}

# plot_mean() takes as inputs the whole dataset and the parameter to plot, as well as the dark and light periods, the 
# folder to save the results and plot labels, and returns 1) mean profiles of that parameter over time for a given genotype, 
# highlighting the dark and light phases, 2) statistical analyses of all results, including warnings in case assumptions 
# for the corresponding analyses are not met.
plot_mean <- function(parameter, data, folder, dark_phase, plot_labels){
  # Appending "In meanplots" in the warning files, which contain information regarding the one-way ANOVA statistical tests
  # that do not meet the underlying assumptions of homogeneity of variances and normality.
  write.table("\n######### In meanplots #########", "./Results/Warnings.txt", append=T, quote=F, col.names=F, row.names=F)
  write.table("\n######### In meanplots #########", "./Results/Warnings_summary.txt", append=T, quote=F, col.names=F, 
              row.names=F)
  
  plots <- list()
  genotype_no <- length(unique(data$Genotype))
  TPvsP <- unique(data[c("TimePoint", "Phase")])

  # For each parameter, data are first formatted for proper analysis. Next, meanplots are plotted.
  for (i in 1:length(parameter)){
    # Formatting data
    meansddf <- summarySE(data, measurevar=parameter[[i]], groupvars=c("Genotype", "TimePoint"))
    meansddf$Phase <- rep(TPvsP$Phase, genotype_no)
    colnames(meansddf) <- c("Genotype", "TimePoint", "N", "Parameter", "sd", "se", "ci", "Phase")
    
    # Plotting data
    k <- ggplot(meansddf, aes(x=TimePoint, y=Parameter, colour=Genotype, group=Genotype)) +
      geom_rect(data=dark_phase, inherit.aes=F, aes(xmin=as.numeric(dark_start), xmax=as.numeric(dark_end), ymin=-Inf, 
          ymax=Inf), fill="azure2") +
      geom_errorbar(aes(ymin=Parameter-se, ymax=Parameter+se), width=.1) +  
      geom_line() +
      geom_point() +
      ylab(parameter[i]) +
      scale_x_continuous("Time point", breaks=seq(1, length(plot_labels), 5), labels=plot_labels[seq(1, length(plot_labels), 5)]) +
      scale_color_locuszoom() +
      theme_bw() +
      theme(panel.grid=element_blank(), axis.text.x=element_text(angle=90)) +
      ylab(names(parameter[i]))
    
    plots <- list.append(plots, k)

    # One-way ANOVA followed by a Tukey's HSD post-hoc test are performed for each timepoint.
    colnames(meansddf) <- c("Genotype", "TimePoint", "N", parameter[[i]], "sd", "se", "ci", "Phase")
    pval <- list()
    signtp <- c()
    tukey <- list()
    n <- 0
    
    for (j in 1:length(unique(data$TimePoint))){
      data.tp <- data[data$TimePoint==j, c("Genotype", parameter[[i]])]
      data.tp$Timepoint <- j
      mean_stats_res <- statistics(data.tp, "meanplot")
      names(mean_stats_res) <- c("anova","tukey","sign")
      
      # Extracting pvalues.
      coeffs <- summary(mean_stats_res$anova)$coeff[-(1),ncol(summary(mean_stats_res$anova)$coeff)]
      pval[[j]] <- data.frame(pvalue=coeffs,significance=rep(".",length(coeffs)))
  
      # Formatting significance levels according to pvalues.
      for(h in 1:nrow(pval[[j]])){
        if(is.na(pval[[j]]$pvalue[h])){
          pval[[j]]$significance[h] <- "NaN"
        } else if(pval[[j]]$pvalue[h]<1e-3){
          pval[[j]]$significance[h] <- "****"
        } else if(pval[[j]]$pvalue[h]<5e-3){
          pval[[j]]$significance[h] <- "***"
        } else if(pval[[j]]$pvalue[h]<1e-2){
          pval[[j]]$significance[h] <- "**"
        } else if(pval[[j]]$pvalue[h]<5e-2){
          pval[[j]]$significance[h] <- "*"
        } 
      }
     
      # Storing Tukey's results in the corresponding timepoints.
      if(mean_stats_res$sign=="yes"){
        n <- n+1
        signtp[n] <- j
        if(mean_stats_res$tukey!="no"){
          tukey[[n]] <- summary(mean_stats_res$tukey)
        } 
      }
    }
      
    # pvalues and significance levels are added to the meansddf dataframe. First, this object is split by genotype. In 
    # the resulting list, elements appear in the same order as the levels of data$Genotype, which can be modified by
    # the user when calling the function. Whether it is modified or not, the reference level is always the first level.
    # The pvalues and significance levels that appear in this dataframe are the result of the one-way ANOVA tests 
    # performed between the reference level and each of the rest of the genotypes. Thus, these columns contain "1"
    # when the genotype in the first column is the reference, because they would be the result of comparing
    # reference vs. reference.
    meansddf <- split(meansddf,meansddf$Genotype)
    meansddf[[1]]$ANOVA_pValue <- as.numeric(1)
    meansddf[[1]]$Sign <- "."
    for (m in 2:length(meansddf)){
      meansddf[[m]]$ANOVA_pValue <-  formatC(sapply(pval, function(x) x$pvalue[m-1]), format="e", digits=3)
      meansddf[[m]]$Sign <- sapply(pval, function(x) x$significance[m-1])
    }
    
    # Turning the modified list into a dataframe again.
    meansddf.df=data.frame()
    for(n in 1:length(meansddf)){
      meansddf.df=rbind(meansddf.df,meansddf[[n]])
    }
    
    # Formatting data and writing results.
    meansddf.df <- meansddf.df %>% mutate_at(c(parameter[[i]], "sd", "se", "ci"), round, digits=3)
    meansddf.df <- meansddf.df[order(meansddf.df$TimePoint),]
    write.table(as.matrix(meansddf.df), file=paste0(folder, parameter[i], "_statistics.txt"), quote=F, col.names=TRUE, 
                row.names=F, sep="\t")
    
    # Adding the results of Tukey's HSD post-hoc tests when applicable.
    if(length(signtp)>0 && genotype_no>2){
      for (x in 1:length(signtp)){
        write.table(paste0("\n### Tukey multiple pairwise comparisons on timepoint ", signtp[x], " ###"), 
                    file=paste0(folder, parameter[i], "_statistics.txt"), append=T, quote=F, row.names=F, col.names=F)
        capture.output(tukey[[x]], file=paste0(folder, parameter[i], "_statistics.txt"), append=T)
      }
    }
    
    # Storing only the statistical significant results in a summary file.
    sign <- meansddf.df[as.numeric(meansddf.df$ANOVA_pValue)<=0.05,]
    if(nrow(sign)>0){
      write.table(paste0("\n\n### ", parameter[i], " ###\n"), "./Results/Significant_5e-2_meanovertimepoints.txt", 
                  quote=FALSE, append=T, row.names=F, col.names=F)
      write.table(sign, "./Results/Significant_5e-2_meanovertimepoints.txt", quote=FALSE, append=T, row.names=F, 
                  col.names=T, sep="\t")
    }
  }
  return(plots)
}


## ARGUMENTS
# Defining and storing the script arguments
option_list <- list(make_option(c("-f", "--file"), type="character", default=NULL, help="Input data file", metavar="character"),
                 make_option(c("-w", "--weight"), type="character", default=NULL, help="Weight data file", metavar="character"), 
                 make_option(c("-s", "--start"), type="integer", default=NULL, help="Start day for analysis", metavar="integer"), 
                 make_option(c("-e", "--end"), type="integer", default=NULL, help="End day for analysis", metavar="integer"), 
                 make_option(c("-l", "--wheel"), type="character", action="store_true", default=NULL, help="Analyse wheel data", 
                             metavar="character"), 
                 make_option(c("-d", "--dfw"), type="character", action="store_true", default=NULL, 
                             help="Analyse drink/feed/weight data", metavar="character"), 
                 make_option(c("-a", "--activity"), type="character", action="store_true", default=NULL, 
                             help="Analyse activity data", metavar="character"), 
                 make_option(c("-c", "--calorimetry"), type="character", action="store_true", default=NULL, 
                             help="Analyse calorimetry data", metavar="character"), 
                 make_option(c("-r", "--reference"), type="character", action="store_true", default=NULL, 
                             help="Analyse reference values", metavar="character"), 
                 make_option(c("-t", "--width"), type="integer", default=7, help="Width of plots in pdf format", 
                             metavar="integer"), 
                 make_option(c("-i", "--height"), type="integer", default=7, help="Height of plots in pdf format", 
                             metavar="integer"), 
                 make_option(c("-v", "--levels"), type="character", default=NULL, help="Genotype order for results, separated by comma", 
                             metavar="character"), 
                 make_option(c("-o", "--omit"), type="character", default=NULL, help="Animals to remove from the analysis, separated by comma", 
                             metavar="character"), 
                 make_option(c("-n", "--numbers"), type="character", action="store_true", default=NULL, 
                             help="Convert timepoints to numbers or use date_time instead", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

## INPUT FILES
# Reading the data (-f, --file) or weight (-w, --weight) files or exiting the program
if(is.null(opt$file) || is.null(opt$weight)){
  print_help(opt_parser)
  stop("Input files missing.\n", call.=FALSE)
}

# TSE output files (-f, --file) have a long header with information regarding animal weights. These rows are skipped.
skip_lines <- processFile(opt$file)

tempdata <- read.csv(opt$file, header=F, skip=skip_lines)
tempdata <- tempdata[,-ncol(tempdata)]
colnames(tempdata) <- c("Date", "Time", "Animal", "Box", "TempC", "HumC", "LightC", "S.Flow", "Ref.O2", "Ref.CO2", "Flow", 
                     "Temp", "O2", "CO2", "dO2", "dCO2", "VO2.1.", "VO2.2.", "VO2.3.", "VCO2.1.", "VCO2.2.", "VCO2.3.", "RER", 
                     "H.1.", "H.2.", "H.3.", "XT.YT", "XT", "XA", "XF", "YT", "YA", "YF", "Z", "CenT", "CenA", "CenF", "PerT", 
                     "PerA", "PerF", "DistK", "DistD", "Speed", "Drink", "Feed", "Weight", "Right", "Left", "SumR.L", "SumTime", 
                     "SumRuns", "MaxSpeed", "AvgSpeed", "MaxLen")

# Reading and storing weight data (-w, --weight). Defining the genotypes to compare
weight_data <- read.table(opt$weight, header=T, colClasses=c("integer", "factor", "numeric"))

## DATA FORMATTING  
# Removing undesired animals
if(!is.null(opt$omit)){
  omit_animals <- as.numeric(unlist(strsplit(opt$omit,",")))
  tempdata <- tempdata[!tempdata$Animal %in% omit_animals,]
  omit_box <- unique(tempdata$Box[tempdata$Animal==omit_animals])
  weight_data <- weight_data[!weight_data$Box %in% omit_box,]
}


# Assigning cycle phases
tempdata$Phase[tempdata$LightC>0] <- "LIGHT"
tempdata$Phase[tempdata$LightC==0] <- "DARK"

# Assigning days
dates <- as.data.frame(cbind(as.character(unique(tempdata$Date)), 1:length(unique(tempdata$Date))))
colnames(dates) <- c("Date", "Day")
tempdata$Day <- as.numeric(unlist(lapply(tempdata$Date, function(x) dates$Day[match(x, dates$Date)])))

# Assigning genotypes and starting weights
tempdata$Genotype <- unlist(lapply(tempdata$Box, function(x) as.character(weight_data$Genotype[match(x, weight_data$Box)])))
genotype_no <- length(unique(tempdata$Genotype))
tempdata$WgStart <- unlist(lapply(tempdata$Box, function(x) weight_data$WgStart[match(x, weight_data$Box)]))

# Replacing non existent data for AvgSpeed by 0.00
suppressWarnings(tempdata$AvgSpeed <- as.numeric(as.character(tempdata$AvgSpeed)))
tempdata$AvgSpeed[is.na(tempdata$AvgSpeed)] <- 0

# Selecting days from input interval, assigning timepoints and taring cumulative parameters if the -s (--start) and -e (--end)
# arguments were indicated
if(!is.null(opt$start) || !is.null(opt$end)){
  interval <- "yes"
  data <- tempdata[tempdata$Day>=opt$start & tempdata$Day<=opt$end,]
  Box <- unique(data$Box)
  data$TimePoint <- rep(1:nrow(unique(data[c("Date", "Time")])), length(Box))
  data[,3:54]=sapply(data[,3:54],function(x) as.numeric(x))
  
  # Getting the baselines for cumulative parameters to subtract them from the rest of the values of the interval 
  drinktare <- c()
  feedtare <- c()
  distktare <- c()
  for(i in Box){
    drinktare[i] <- data$Drink[data$Box==Box[i]][1]
    feedtare[i] <- data$Feed[data$Box==Box[i]][1]
    distktare[i] <- data$DistK[data$Box==Box[i]][1]
  }
  tare <- as.data.frame(cbind(Box, drinktare, feedtare, distktare))

  data$DrinkT0 <- unlist(lapply(data$Box, function(x) tare$drinktare[match(x, tare$Box)]))
  data$TareDrink <- data$Drink-data$DrinkT0
  data$FeedT0 <- unlist(lapply(data$Box, function(x) tare$feedtare[match(x, tare$Box)]))
  data$TareFeed <- data$Feed-data$FeedT0
  data$DistkT0 <- unlist(lapply(data$Box, function(x) tare$distktare[match(x, tare$Box)]))
  data$TareDistK <- data$DistK-data$DistkT0
  
  data <- subset(data,select=-c(DrinkT0,FeedT0,DistkT0))
} else {
  interval <- "no"
  data <- tempdata
}

# Extracting dark periods for shading plots
dark <- c()
dark <- unique(split(data, data$Phase)$DARK[,"TimePoint"])

dark_start <- c()
dark_end <- c()
j <- 1

dark_start[1] <- dark[1]
for(i in 2:length(dark)){
  if(dark[i]==dark[i-1]+1){
    next() 
  }else{
    dark_end[j] <- dark[i-1]
    dark_start[j+1] <- dark[i]
    j <- j+1
  }
}
dark_end[j] <- dark[i]

dark_phase <- as.data.frame(cbind(dark_start, dark_end))

# Setting the desired genotype order for results if the -v (--levels) argument was used
if(!is.null(opt$levels)){
  data$Genotype <- factor(data$Genotype, levels=unlist(strsplit(opt$levels, ",")))
}

# Adding cumulative data for specific parameters
cum <- c("XT.YT", "XT", "XA", "XF", "YT", "YA", "YF", "Z", "CenT", "CenA", "CenF", "PerT", "PerA", "PerF", "Right", "Left", 
         "SumR.L", "SumTime", "SumRuns")
for(i in 1:length(cum)){
  name <- paste("Cum", cum[i], sep="")
  data[[name]] <- ave(data[[cum[i]]], data$Box, FUN=cumsum)
}

## FOLDERS AND STATISTICS FILES
# Creating the directory to store the results and writing the formatted dataset
dir.create("./Results")
write.table(data, "./Results/Formatted_dataset.txt", quote=F, col.names=TRUE, row.names=F)

# Creating files for only significant pvalues and warnings:
# 1) Significant_5e-2_boxplots.txt: Summary file specifying which parameters are statistically significant (ANOVA
# pvalue < 0.05), and in which phase in the analyses of boxplot data. Details of the analysis can be found in the 
# corresponding *_boxplots_statistics.txt files (see plot_bxplt() above)
# 2) Significant_5e-2_meanovertimepoints.txt: Summary file specifying which parameters are statistically significant, 
# (ANOVA pvalue < 0.05) and in which phase and specific time point, as a result of the analyses with mean plots. 
# Statistical details can be found in the corresponding *_statistics.txt files (see plot_mean() above)
# 3) Warnings.txt: File specifying detailed information regarding those one-way ANOVA statistical analyses where care 
# should be taken before reaching conclusions because one or both of the underlying assumptions (homogeneity of 
# variances and normality) are not met
# 4) Warnings_summary.txt: Summary file containing in which ANOVA statistical analysis (parameter, plot type,
# timepoint/phase), the underlying assumptions (homogeneity of variances and normality) are not met
write.table("PARAMETERS WITH STATISTICALLY SIGNIFICANT DIFFERENCES IN BOXPLOTS (pval<0.05)\n", 
            "./Results/Significant_5e-2_boxplots.txt", quote=F, col.names=F, row.names=F)
write.table("PARAMETERS AND TIMEPOINTS WITH STATISTICALLY SIGNIFICANT DIFFERENCES IN ANOVA (pval<0.05)", 
            "./Results/Significant_5e-2_meanovertimepoints.txt", quote=F, col.names=F, row.names=F)
write.table("WARNINGS BEFORE CONSIDERING ANOVA RESULS", "./Results/Warnings.txt", quote=F, col.names=F, row.names=F)
write.table("WARNINGS BEFORE CONSIDERING ANOVA RESULS - SUMMARY", "./Results/Warnings_summary.txt", quote=F, col.names=F, 
            row.names=F)

# One-way ANOVA test for weight data. Tukey HSD post-hoc test is applied if more than 2 genotypes are being compared.
# Homogeneity of variances and normality assumptions are checked with the Barlett's and Shapiro-Wilk tests, respectively.
res.aovwg <- aov(data=weight_data, WgStart ~ Genotype)

write.table("######################################################\n#################### ANOVA RESULTS ###################\n######################################################\n\n", 
            file="./Results/ANOVA_weight.txt", quote=F, row.names=F, col.names=F)
capture.output(summary(res.aovwg), file="./Results/ANOVA_weight.txt", append=T)

if(genotype_no>2){
  tukeywg <- TukeyHSD(res.aovwg)
  write.table("\n\n\n######################################################\n######### Tukey multiple pairwise comparisons ########\n######################################################\n\n", 
              file="./Results/ANOVA_weight.txt", append=T, quote=F, row.names=F, col.names=F)
  capture.output(tukeywg, file="./Results/ANOVA_weight.txt", append=T)
}

barltestwg <- bartlett.test(data=weight_data, WgStart ~ Genotype)
write.table("\n\n######################################################\n################### Barlett's test ###################\n######################################################\n", 
            file="./Results/ANOVA_weight.txt", append=T, quote=F, row.names=F, col.names=F)
capture.output(barltestwg, file="./Results/ANOVA_weight.txt", append=T)

aov_residualswg <- residuals(object=res.aovwg)
shapwg <- shapiro.test(x=aov_residualswg )
write.table("\n\n######################################################\n### Shapiro-Wilk's test for normality of residuals ###\n######################################################\n", 
            file="./Results/ANOVA_weight.txt", append=T, quote=F, row.names=F, col.names=F)
capture.output(shapwg, file="./Results/ANOVA_weight.txt", append=T)

## PLOTS
# The optional -t (--width) and -i (--height) can be used to modify plot sizes

# Assigning labels for plots. If the script was called with the -n (--number) argument, then time points are indicated as 
# numbers. If this argument was not used, then time points are labeled as Date_Hour.
if(!is.null(opt$numbers)){
  plot_labels <- unique(data$TimePoint) 
} else{
  plot_labels <- unique(paste0(data$Date, "_", data$Time))
}

# Plotting reference values if requested with the -r (--reference) argument. Two types of plots are produced: 1) individual 
# profiles of the corresponding parameter over time, highlighting the dark and light phases, and 2) boxplots depicting results 
# grouped by phase. In the latter case, mean values are indicated and marked with a dark red dot. 
if(!is.null(opt$reference)){
  cat("\nGenerating reference plots...\n")
  reference <- c("Chamber temperature (ºC)"="TempC", "Chamber humidity (%)"="HumC", "Brightness (%)"="LightC", 
                 "Sample flow (L/min)"="S.Flow", "Reference O2 value (%)"="Ref.O2", "Reference CO2 value (%)"="Ref.CO2")
  
  pdf("./Results/Reference_plots.pdf", width=opt$width, height=opt$height)
  for (i in 1:length(reference)){
    p <- ggplot(data=data[data$Box=="1",], aes_string(x="TimePoint", y=reference[i], group="Animal")) +
      geom_rect(data=dark_phase, mapping=aes(xmin=as.numeric(dark_phase$dark_start), xmax=as.numeric(dark_phase$dark_end),
          ymin=-Inf, ymax=+Inf), inherit.aes=F, fill="azure2") +
      geom_point(size=0.1) +
      geom_line() +
      scale_x_continuous("Time point", breaks=seq(1, length(plot_labels), 5), 
          labels=plot_labels[seq(1, length(plot_labels), 5)]) +
      theme_classic() +
      theme(axis.text.x=element_text(angle=90)) + 
      ylab(names(reference[i]))
    g <- ggplot(data=data[data$Box=="1",], aes_string("Phase", reference[i], fill="Phase")) +
      geom_boxplot(outlier.shape=NA) + scale_fill_manual(values=c("azure2", "white")) +
      geom_point(size=1, position=position_jitterdodge((jitter.width=0.3))) +
      stat_summary(fun.y=mean, geom="point", color="darkred") +
      stat_summary(fun.y=mean, geom="text", color="black", vjust=-0.7, aes(label=round(..y.., digits=2))) +
      theme_classic() + 
      theme(legend.position="none") +
      ylab(names(reference[i]))
    invisible(print(p))
    invisible(print(g))
  }
  graphics.off()
}

# Plotting calorimetry values if requested with the -c (--calorimetry) argument using the plot_ind(), plot_bxplt() and 
# plot_mean() functions. Results are stored in the Calorimetry folder.
if(!is.null(opt$calorimetry)){
  dir.create("./Results/Calorimetry")
  folder <- "./Results/Calorimetry/"
  cat("\nGenerating calorimetry plots...\n")
  
  # All warnings regarding statistical analyses that are generated with the plot_ind(), plot_bxplt() and plot_mean() 
  # functions, will be stored in the corresponding warning files under the CALORIMETRY subheading. Even if no
  # warnings are produced, these subheadings will appear to inform the user that the program ran succesfully but 
  # without warnings.
  write.table("\n\n###############################\n######### CALORIMETRY #########\n###############################",
              "./Results/Warnings.txt", append=T, quote=F, col.names=F, row.names=F)
  write.table("\n\n###############################\n######### CALORIMETRY #########\n###############################",
              "./Results/Warnings_summary.txt", append=T, quote=F, col.names=F, row.names=F)
  
  # Specifying the calorimetry parameters and units
  calo <- c("Flow per box (L/min)"="Flow", "Temperature (ºC)"="Temp", "O2 concentration (%)"="O2", 
            "CO2 concentration (%)"="CO2", "Difference in O2 concentration relative to the chamber (%)"="dO2", 
            "Difference in CO2 concentration relative to the chamber (%)"="dCO2", 
            "O2 production per total weight (mL/h/kg)"="VO2.1.", "O2 production per lean body mass (mL/h/kg)"="VO2.2.", 
            "O2 production (mL/h)"="VO2.3.", "CO2 consumption per total weight (mL/h/kg)"="VCO2.1.", 
            "CO2 consumption per lean body mass (mL/h/kg)"="VCO2.2.", "CO2 consumption(mL/h)"="VCO2.3.", 
            "Respiratory Exchange Ratio"="RER", "Heat production per total weight (kcal/h/kg)"="H.1.", 
            "Heat production per lean body mass (kcal/h/kg)"="H.2.", "Heat production (kcal/h)"="H.3.")
  
  # Generating individual plots and RER in function of wheel activity
  pdf("./Results/Calorimetry/Calorimetry_individual_plots.pdf", width=opt$width, height=opt$height)
  indvplots <- plot_ind(calo, data, dark_phase, plot_labels)
  for(i in 1:length(indvplots)){
    invisible(print(indvplots[[i]]))
  }

  invisible(print(ggplot(data[data$SumR.L>0,], aes_string(x="SumR.L", y="RER", colour="Genotype", group="Animal")) +
                    geom_point(size=1) +
                    facet_wrap(~Phase) +
                    theme_bw() +
                    theme(strip.background=element_rect(fill="white"), strip.text=element_text(face="bold")) +
                    scale_color_locuszoom() +
                    xlab("Total running distance covered (cm)"))) 
  invisible(print(ggplot(data[data$SumR.L>0,], aes_string(x="SumR.L", y="RER", colour="Phase", shape="Genotype")) +
                    geom_point(size=1) + 
                    facet_wrap(~Animal) +
                    theme_bw() +
                    theme(strip.background=element_rect(fill="white"), strip.text=element_text(face="bold")) +
                    scale_color_locuszoom() +
                    xlab("Total running distance covered (cm)"))) 
  graphics.off()
  
  # Do not generate boxplots or mean plots for the flow within cages, which is constant
  calo <- calo[2:length(calo)]
  
  # Generating boxplots
  pdf("./Results/Calorimetry/Calorimetry_boxplots.pdf", width=opt$width, height=opt$height)
  bxplts <- plot_bxplt(calo, data, pval_wg, folder, dark_phase)
  for(i in 1:length(bxplts)){
    invisible(print(bxplts[[i]]))
  }
  graphics.off()
  
  # Generating mean plots
  pdf("./Results/Calorimetry/Calorimetry_meanplots.pdf", width=opt$width, height=opt$height)
  meanplots <- plot_mean(calo, data, folder, dark_phase, plot_labels)
  for(i in 1:length(meanplots)){
    invisible(print(meanplots[[i]]))
  }    
  graphics.off()
}

# Plotting activity values if requested with the -a (--activity) argument using the plot_ind(), plot_bxplt() and 
# plot_mean() functions. Results are stored in the Activity folder.
if(!is.null(opt$activity)){
  dir.create("./Results/Activity")
  folder <- "./Results/Activity/"
  cat("\nGenerating activity plots...\n")
  
  # All warnings regarding statistical analyses that are generated with the plot_ind(), plot_bxplt() and plot_mean() 
  # functions, will be stored in the corresponding warning files under the ACTIVITY subheading. 
  write.table("\n\n###############################\n#########   ACTIVITY  #########\n###############################",
              "./Results/Warnings.txt", append=T, quote=F, col.names=F, row.names=F)
  write.table("\n\n###############################\n#########   ACTIVITY  #########\n###############################",
              "./Results/Warnings_summary.txt", append=T, quote=F, col.names=F, row.names=F)
  
  # Specifying the activity parameters and units
  activity <- c("Total XY-beam interruptions (counts)"="XT.YT", "Total X-beam interruptions (counts)"="XT", 
                "Ambulatory X-beam interruptions (counts)"="XA", "Fine movement X-beam interruptions (counts)"="XF", 
                "Total Y-beam interruptions (counts)"="YT", "Ambulatory Y-beam interruptions (counts)"="YA", 
                "Fine movement Y-beam interruptions (counts)"="YF", "Rearing Z-beam interruptions (counts)"="Z", 
                "Central total movement (counts)"="CenT", "Central ambulatory movement (counts)"="CenA", 
                "Central fine movement (counts)"="CenF", "Peripheral total movement (counts)"="PerT", 
                "Peripheral ambulatory movement (counts)"="PerA", "Peripheral fine movement (counts)"="PerF")
  cumactivity <- c("Cumulative total XY-beam interruptions (counts)"="CumXT.YT", 
                   "Cumulative total X-beam interruptions (counts)"="CumXT",
                   "Cumulative ambulatory X-beam interruptions (counts)"="CumXA",
                   "Cumulative fine movement X-beam interruptions (counts)"="CumXF", 
                   "Cumulative total Y-beam interruptions (counts)"="CumYT",
                   "Cumulative ambulatory Y-beam interruptions (counts)"="CumYA",
                   "Cumulative fine movement Y-beam interruptions (counts)"="CumYF", 
                   "Cumulative rearing Z-beam interruptions (counts)"="CumZ",
                   "Cumulative central total movement (counts)"="CumCenT",
                   "Cumulative central ambulatory movement (counts)"="CumCenA",
                   "Cumulative central fine movement (counts)"="CumCenF",
                   "Cumulative peripheral total movement (counts)"="CumPerT",
                   "Cumulative peripheral ambulatory movement (counts)"="CumPerA",
                   "Cumulative peripheral fine movement (counts)"="CumPerF")
  if(interval=="no"){
    cumactivity <- c(cumactivity, "Cumulative distance (cm)"="DistK")
  } else{
    cumactivity <- c(cumactivity, "Cumulative distance (cm)"="TareDistK")
  }
  othersactivity <- c("Differential distance (cm)"="DistD", "Speed (cm/s)"="Speed")
  
  all_activity <- list(activity, cumactivity, othersactivity)
  raw_activity <- list(activity, othersactivity)
  
  # Generating individual plots
  pdf("./Results/Activity/Activity_individual_plots.pdf", width=opt$width, height=opt$height)
  for(j in 1:length(all_activity)){
    indvplots <- plot_ind(all_activity[[j]], data, dark_phase, plot_labels)
    for(i in 1:length(indvplots)){
      invisible(print(indvplots[[i]]))
    }
  }
  graphics.off()
  
  # Generating boxplots
  pdf("./Results/Activity/Activity_boxplots.pdf", width=opt$width, height=opt$height)
  for(j in 1:length(raw_activity)){
    bxplts <- plot_bxplt(raw_activity[[j]], data, pval_wg, folder, dark_phase)
    for(i in 1:length(bxplts)){
      invisible(print(bxplts[[i]]))
    }
  }
  graphics.off()
  
  # Generating mean plots
  pdf("./Results/Activity/Activity_meanplots.pdf", width=opt$width, height=opt$height)
  for(j in 1:length(all_activity)){
    meanplots <- plot_mean(all_activity[[j]], data, folder, dark_phase, plot_labels)
    for(i in 1:length(meanplots)){
      invisible(print(meanplots[[i]]))
    }  
  }
  graphics.off()
}

# Plotting DFW values if requested with the -d (--dfw) argument using the plot_ind(), plot_bxplt() and 
# plot_mean() functions. Results are stored in the DFW folder.
if(!is.null(opt$dfw)){
  dir.create("./Results/DFW")
  folder <- "./Results/DFW/"
  cat("\nGenerating DFW plots...\n")
  
  # All warnings regarding statistical analyses that are generated with the plot_ind(), plot_bxplt() and plot_mean() 
  # functions, will be stored in the corresponding warning files under the DFW subheading. 
  write.table("\n\n###############################\n#########     DFW     #########\n###############################", 
              "./Results/Warnings.txt", append=T, quote=F, col.names=F, row.names=F)
  write.table("\n\n###############################\n#########     DFW     #########\n###############################", 
              "./Results/Warnings_summary.txt", append=T, quote=F, col.names=F, row.names=F)
  
  # Specifying the activity parameters and units
  if(interval=="no"){
    dfw <- c("Water consumption (mL)"="Drink", "Food consumption (g)"="Feed", "Weight (g)"="Weight")
  } else{
    dfw <- c("Water consumption (mL)"="TareDrink", "Food consumption (g)"="TareFeed", "Weight (g)"="Weight")
  }
  
  # Generating individual plots
  pdf("./Results/DFW/DFW_individual_plots.pdf", width=opt$width, height=opt$height)
  indvplots <- plot_ind(dfw, data, dark_phase, plot_labels)
  for(i in 1:length(indvplots)){
    invisible(print(indvplots[[i]]))
  }
  graphics.off()
  
  # Generating mean plots
  pdf("./Results/DFW/DFW_meanplots.pdf", width=opt$width, height=opt$height)
  meanplots <- plot_mean(dfw, data, folder, dark_phase, plot_labels)
  for(i in 1:length(meanplots)){
    invisible(print(meanplots[[i]]))
  } 
  graphics.off()
}

# Plotting wheel values if requested with the -l (--wheel) argument using the plot_ind(), plot_bxplt() and 
# plot_mean() functions. Results are stored in the Wheels folder.
if(!is.null(opt$wheel)){
  dir.create("./Results/Wheels")
  folder <- "./Results/Wheels/"
  cat("\nGenerating wheel plots...\n")
  
  # All warnings regarding statistical analyses that are generated with the plot_ind(), plot_bxplt() and plot_mean() 
  # functions, will be stored in the corresponding warning files under the WHEELS subheading. 
  write.table("\n\n###############################\n#########   WHEELS   ##########\n###############################", 
              "./Results/Warnings.txt", append=T, quote=F, col.names=F, row.names=F)
  write.table("\n\n###############################\n#########   WHEELS   ##########\n###############################", 
              "./Results/Warnings_summary.txt", append=T, quote=F, col.names=F, row.names=F)
  
  # Specifying the activity parameters and units
  wheel <- c("Distance covered to the right (cm)"="Right", "Distance covered to the left (cm)"="Left", 
           "Total distance covered (cm)"="SumR.L", "Total running duration (min)"="SumTime", 
           "Total number of running episodes"="SumRuns")
  cumwheel <- c("Cumulative distance covered to the right (cm)"="CumRight", "Cumulative distance covered to the left (cm)"="CumLeft", 
              "Cumulative total distance covered (cm)"="CumSumR.L", "Cumulative running duration (min)"="CumSumTime", 
              "Cumulative number of running episodes"="CumSumRuns")
  otherswheel <- c("Maximum running speed (m/s)"="MaxSpeed", "Average running speed (m/s)"="AvgSpeed", 
                 "Maximum duration of running episodes (s)"="MaxLen")
  
  all_wheel <- list(wheel, cumwheel, otherswheel)
  raw_wheel <- list(wheel, otherswheel)
  
  # Generating individual plots
  pdf("./Results/Wheels/wheel_individual_plots.pdf", width=opt$width, height=opt$height)
  for(j in 1:length(all_wheel)){
    indvplots <- plot_ind(all_wheel[[j]], data, dark_phase, plot_labels)
    for(i in 1:length(indvplots)){
      invisible(print(indvplots[[i]]))
    }
  }
  graphics.off()
  
  # Generating boxplots
  pdf("./Results/Wheels/wheel_boxplots.pdf", width=opt$width, height=opt$height)
  for(j in 1:length(raw_wheel)){
    bxplts <- plot_bxplt(raw_wheel[[j]], data, pval_wg, folder, dark_phase)
    for(i in 1:length(bxplts)){
      invisible(print(bxplts[[i]]))
    }
  }
  graphics.off()
  
  # Generating mean plots
  pdf("./Results/Wheels/wheel_meanplots.pdf", width=opt$width, height=opt$height)
  for(j in 1:length(all_wheel)){
    meanplots <- plot_mean(all_wheel[[j]], data, folder, dark_phase, plot_labels)
    for(i in 1:length(meanplots)){
      invisible(print(meanplots[[i]]))
    } 
  } 
  graphics.off()
}

cat("\nDone.\n")




