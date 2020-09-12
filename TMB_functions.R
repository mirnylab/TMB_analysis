# Figure 1A
plot_all_cancers <- function(data_)
{
  # data_ = data_tmb_wes # For testing
  data_ = data_[c("TMB", "dataset", "response")] 
  data_ = data_[complete.cases(data_),]
  
  data_$dataset = factor(data_$dataset, 
                         levels=c("SCLC", "Sarcoma", "Anal", 'bladder2', 'ccRCC',
                                  'HNSCC', 'uro1', 'Bladder', 'lung1', 'lung2',
                                  'mel1', 'mel2', 'mel3', 'mel4', 'mel5'))
  
  data_$response = factor(data_$response, 
                         levels=c("NR", "R"))
  
  data_ = data_ %>% 
  group_by(dataset) %>% 
  mutate(dataset = paste0(dataset, "\nn=", dplyr::n()))
  
  p <- ggboxplot(data_, x = "response", y = "TMB",
                 color = "response", palette =  c( "darkgoldenrod2", "cadetblue3"), 
                 add = "jitter",
                 facet.by = "dataset", short.panel.labs = T) 
  p = p + yscale("log10", .format = TRUE)  + coord_cartesian(ylim = c(1,10^4.5))
  p = p + stat_compare_means(comparisons = list(c("R", "NR")),
                             label = "p.format", 
                             method = "wilcox.test", 
                             method.args = list(alternative = "t"))+
    facet_wrap(~dataset,  nrow=1)
  return(p)
  
}

# Figure 1B and 1C
plot_stratif <- function(data_, title)
{
  # data_ = mel1  # For testing
  # title = "mel1"  # For testing
  data_ = data.frame(data_[,c("TMB", "response", "Stratif_type")])
  data_ = data_[complete.cases(data_),]
  
  data_$response = factor(data_$response, levels=c('NR', 'R'))
  data_$Stratif_type = factor(data_$Stratif_type, levels=c(
                                      'acral/ mucosal',
                                      'skin/ occult',
                                      'never',
                                      'former/ current',
                                      'all'))
  data_ = data_ %>% 
    group_by(Stratif_type) %>% 
    mutate(Stratif_type = paste0(Stratif_type, "\nn=", dplyr::n()))
  
  p <- ggboxplot(data_, x = "response", y = "TMB",
                 color = "response", palette =  c( "darkgoldenrod2", "cadetblue3"), 
                 add = "jitter",
                 facet.by = "Stratif_type", short.panel.labs = FALSE)
  p = p + stat_compare_means(comparisons = list(c("R", "NR")),
                             label = "p.format", method = "wilcox.test")
  p = p + theme(legend.position = "none")
  p = p + yscale("log10", .format = TRUE)  + coord_cartesian(ylim = c(1,10^4.5)) + ggtitle(title)
  
  return(p)
}

# Figure 2A
plot_raw_data <- function(data_, title)
{
  # data_ = mel1  # For testing
  # title = "mel1"  # For testing
  p = ggscatter(data_, x = "TMB", y = "PFS",
                color = "response",fill = "response",  size = 4, shape=as.numeric(data_$PFS_censorship)+16, alpha = 0.7) + 
    yscale("log10", .format = TRUE) + 
    xscale("log10", .format = TRUE) + coord_cartesian(ylim = c(0.3,10^2), xlim = c(1, 10^5)) +
    scale_colour_manual(values = c("darkgoldenrod2", "cadetblue3"), na.value = "grey")
  p = p + theme(legend.position = "none") + ggtitle(title)
  return(p)
}

# Figure 2C, 2D, S2A and S2B
permutation_analysis <- function(clinical.data_, 
                                 stratification_1, 
                                 stratification_2,
                                 survival_type,
                                 survival_censor)
{
  # clinical.data_ = mel1 # for testing
  # stratification_1 = ("skin/ occult") # for testing
  # stratification_2 = ("acral/ mucosal") # for testing
  # survival_type = "PFS" # for testing
  # survival_censor = "PFS_censorship" # for testing
  
  clinical.data_ = as.data.frame(clinical.data_)
  clinical.data_ = clinical.data_[which(!is.na(clinical.data_[,survival_type])),]
  clinical.data_ = clinical.data_[which(!is.na(clinical.data_[,survival_censor])),]
  all_ = permut_logrank(survival = clinical.data_[,survival_type],
                        tmb = clinical.data_$TMB,
                        censorship = clinical.data_[,survival_censor])
  real_all = logrank_percutoff(survival = clinical.data_[,survival_type],
                               tmb = clinical.data_$TMB,
                               censorship = clinical.data_[,survival_censor])
  
  clinical.data_1 = clinical.data_[which(clinical.data_$Stratif_type == stratification_1),]
  stratif_1 = permut_logrank(survival = clinical.data_1[,survival_type],
                             tmb = clinical.data_1$TMB,
                             censorship = clinical.data_1[,survival_censor])
  real_stratif_1 = logrank_percutoff(survival = clinical.data_1[,survival_type],
                                     tmb = clinical.data_1$TMB,
                                     censorship = clinical.data_1[,survival_censor])
  
  clinical.data_2 = clinical.data_[which(clinical.data_$Stratif_type == stratification_2),]
  stratif_2 = permut_logrank(survival = clinical.data_2[,survival_type],
                             tmb = clinical.data_2$TMB,
                             censorship = clinical.data_2[,survival_censor])
  real_stratif_2 = logrank_percutoff(survival = clinical.data_2[,survival_type],
                                     tmb = clinical.data_2$TMB,
                                     censorship = clinical.data_2[,survival_censor]) 
  
  
  to_plot = rbind(cbind(all_, "all"),
                  cbind(stratif_1, stratification_1),
                  cbind(stratif_2, stratification_2))
  to_plot = data.frame(to_plot, stringsAsFactors = F)
  colnames(to_plot) = c("p_value", "stratification")
  to_plot$p_value = as.numeric(as.character(to_plot$p_value))
  to_plot$stratification = (as.character(to_plot$stratification))
  
  df = data.frame(p_value = c(real_all, real_stratif_1, real_stratif_2), 
                  stratification = c("all", stratification_1, stratification_2))
  
  to_plot$stratification = factor(to_plot$stratification, levels=c(
                                                   'acral/ mucosal',
                                                   'skin/ occult',
                                                   'never',
                                                   'former/ current',
                                                   'all'))
  to_plot[which(to_plot$p_value < 10^-10),1] = 10^-10
  to_plot2 = to_plot
  p<-ggplot(to_plot2, aes(x=stratification, y=p_value, fill=stratification, levels = stratification)) +
    # scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + 
    geom_violin(trim=T, scale = "width") + yscale("log10", .format = TRUE) +
    geom_hline(aes(yintercept = 0.05), color = "chartreuse3",  linetype="dashed") +
    geom_point(data = df, col = 'darkred', shape = 18, size = 5) + 
    coord_flip(ylim = c(10^-10,1))
  print(p)
  
  pval_all = sum(to_plot[which(to_plot$stratification == "all"),"p_value"] <= real_all)/ length(to_plot[which(to_plot$stratification == "all"),"p_value"])
  pval_stratif_1 = sum(to_plot[which(to_plot$stratification == stratification_1),"p_value"] <= real_stratif_1)/ length(to_plot[which(to_plot$stratification == stratification_1),"p_value"])
  pval_stratif_2 = sum(to_plot[which(to_plot$stratification == stratification_2),"p_value"] <= real_stratif_2)/ length(to_plot[which(to_plot$stratification == stratification_2),"p_value"])
  
  
  return(list(c(pval_all, pval_stratif_1, pval_stratif_2),
              p, to_plot))
  
}
permut_logrank <-function(survival, tmb, censorship)
{
  shuffled_pvalues = c()
  for (i in 1:1000)
  {
    set.seed(i)
    shuffled_pvalues = c(shuffled_pvalues,
                         logrank_percutoff(survival, sample(tmb), censorship))
    
  }
  return(shuffled_pvalues)
}
logrank_percutoff <-function(survival, tmb, censorship)
{
  pvalues = c()
  data_ = data.frame(cbind(survival, tmb, censorship))
  for (i in sort(tmb))
  {
    if (!(i %in% c(max(tmb), min(tmb))))
    {
      data_tmp = data_
      data_tmp$tmb_categ = NA
      data_tmp$tmb_categ = ifelse(data_tmp$tmb < i ,"lTMB",
                                  ifelse(data_tmp$tmb > i,"hTMB",
                                         NA))
      res_surv <- survdiff(Surv(survival, censorship) ~ tmb_categ,  data = data_tmp,
                           rho = 0)
      res_p.value <- 1 - pchisq(res_surv$chisq, length(res_surv$n) - 1)
      pvalues = c(pvalues, res_p.value)
    }
  }
  return(min(pvalues))
}

# Figure 3A and 3B
plot_auc <- function(data_)
{
  data_ = data_tmb_wes # For testing
  data_auc = c()
  data_auc$tmb = as.numeric(as.character(data_$mutation_rate))
  data_auc$response = as.character(data_$response)
  data_auc$dataset = as.character(data_$dataset)
  data_auc = data.frame(data_auc)
  data_auc = data_auc[complete.cases(data_auc),]
  
  cp <- cutpointr(data_auc, tmb,response,  dataset)
  opt_cut <- cutpointr(data_auc, tmb, response, dataset, metric = youden)
  p = plot_roc(opt_cut) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  print(p)
  return(c(opt_cut$AUC, opt_cut$optimal_cutpoint))
}

# Figure 3C
FDA_youden_cutoffs <- function(youden_indexes)
{
  data_ = c()
  data_$TMB = c(as.numeric(as.character(mel1$mutation_rate)),
                as.numeric(as.character(mel2$mutation_rate)),
                as.numeric(as.character(lung1$mutation_rate)),
                as.numeric(as.character(lung2$mutation_rate)))
  data_$response = c((as.character(mel1$response)),
                     (as.character(mel2$response)),
                     (as.character(lung1$response)),
                     (as.character(lung2$response)))
  data_$dataset = c(rep( "mel1", length(as.character(mel1$response))),
                    rep( "mel2", length(as.character(mel2$response))),
                    rep( "lung1", length(as.character(lung1$response))),
                    rep( "lung2", length(as.character(lung2$response))))
  data_$youden = c(rep( youden_indexes[1], length(as.character(mel1$response))),
                   rep( youden_indexes[2], length(as.character(mel2$response))),
                   rep( youden_indexes[3], length(as.character(lung1$response))),
                   rep( youden_indexes[4], length(as.character(lung2$response))))
  data_ = data.frame(data_)
  data_ = data_[complete.cases(data_),]
  data_$dataset = factor(data_$dataset, levels=c('mel1',
                                                 'mel2',
                                                 'lung1',
                                                 'lung2'))
  
  p <- ggboxplot(data_, x = "response", y = "TMB",
                 color = "response", palette =  c( "darkgoldenrod2", "cadetblue3"), 
                 # add = "jitter",
                 facet.by = "dataset", short.panel.labs = T) 
  # Use only p.format as label. Remove method name.
  p = p + yscale("log10", .format = TRUE)  + coord_cartesian(ylim = c(1,10^4.5))
  p = p + stat_compare_means(comparisons = list(c("NR", "R")),
                             label = "p.format", method = "wilcox.test")
  p = p+geom_hline(aes(yintercept = 10), color = "black",  linetype="dashed")+
    facet_wrap(~dataset,  ncol=1)
  p = p + geom_hline(aes(yintercept = youden),  color = "black",  linetype="dashed") + coord_flip()
  print(p)
  
}

# Figure 3D
misclassified_pats <- function(youden_indexes, mel1, mel2, lung1, lung2)
{
  NR_getting_trt = c()
  R_not_getting_trt = c()
  
  j = 1
  for (i in list(mel1, mel2, lung1, lung2))
  {
    NR_getting_trt = rbind(NR_getting_trt,
                           data.frame(Proportion = length(i[which((i$mutation_rate>=10) & (i$response == "NR")),1])/length(i$response[which(i$response == "NR")]),
                                      row.names = paste("FDA", i$dataset[1])))
    
    NR_getting_trt = rbind(NR_getting_trt,
                           data.frame(Proportion = length(i[which((i$mutation_rate>=youden_indexes[j]) & (i$response == "NR")),1])/length(i$response[which(i$response == "NR")]),
                                      row.names = paste("Youden", i$dataset[1])))
    
    R_not_getting_trt = rbind(R_not_getting_trt,
                           data.frame(Proportion = length(i[which((i$mutation_rate<10) & (i$response == "R")),1])/length(i$response[which(i$response == "R")]),
                                      row.names = paste("FDA", i$dataset[1])))
    
    R_not_getting_trt = rbind(R_not_getting_trt,
                           data.frame(Proportion = length(i[which((i$mutation_rate<youden_indexes[j]) & (i$response == "R")),1])/length(i$response[which(i$response == "R")]),
                                      row.names = paste("Youden", i$dataset[1])))
    j = j+1
  }
  
  NR_getting_trt$x_label = rownames(NR_getting_trt)
  R_not_getting_trt$x_label = rownames(R_not_getting_trt)
  
  NR_getting_trt$Proportion =  round(NR_getting_trt$Proportion, 2)
  R_not_getting_trt$Proportion =  round( R_not_getting_trt$Proportion, 2)
  
  R_not_getting_trt$x_label = factor(R_not_getting_trt$x_label, 
                                     levels =c("Youden lung2",
                                               "FDA lung2",
                                               "Youden lung1",
                                               "FDA lung1",
                                               "Youden mel2",
                                               "FDA mel2",
                                               "Youden mel1",
                                               "FDA mel1"))
  
  
  NR_getting_trt$x_label = factor(NR_getting_trt$x_label, 
                                      levels =c("Youden lung2",
                                                "FDA lung2",
                                                "Youden lung1",
                                                "FDA lung1",
                                                "Youden mel2",
                                                "FDA mel2",
                                                "Youden mel1",
                                                "FDA mel1"))
  
  p1 = ggplot(data=R_not_getting_trt, aes(x=x_label, y=Proportion)) +
    geom_bar(stat="identity")+
    geom_text(aes(y=Proportion, label=Proportion), vjust=1.6, 
              size=3.5)+
    scale_fill_brewer(palette="Paired")+ coord_flip()+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  
    ggtitle("Responders not getting treatment")
  
  p2 = ggplot(data=NR_getting_trt, aes(x=x_label, y=Proportion)) +
    geom_bar(stat="identity")+
    geom_text(aes(y=Proportion, label=Proportion), vjust=1.6, 
              size=3.5)+
    scale_fill_brewer(palette="Paired")+ coord_flip()+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    ggtitle("Non responders getting treatment")
  

  p = gridExtra::grid.arrange(p1, p2, nrow = 1, ncol =2)
  return(p)
}

# Figure 4B
immunogenicity_model <- function()
{
  N=1:10000;
  p=0.22;
  kcut=1
  P = 1- ppois(kcut, N*p)
  
  cols = viridis(42) 
  plot(N, P, log = "x", xlab = "TMB", ylab = "Probability of being immunogenic", type = "l", ylim = c(10^-6,1), col = cols[1])
  i=1
  for (p_tmp in c(0.22, 0.64))
  {
    P = 1- ppois(kcut, N*p_tmp)
    lines(N, P, col = cols[i], lwd = 2)
    i = i+1
  }
  
  
  kcut=2
  cols = heat.colors(42) 
  i=1
  for (p_tmp in c(0.22, 0.64))
  {
    P = 1- ppois(kcut, N*p_tmp)
    lines(N, P, col = cols[i], lwd = 2)
    i = i+1
  }
}

# Figure S1
copd_tcga_analysis <- function(data_)
{
  my_comparisons <- list(c("0", "1"))
  ggboxplot(data_, x = "Fev1.fvc.ratio.postbroncholiator", y = "tmb",
            color = "Fev1.fvc.ratio.postbroncholiator", palette =  c( "darkgoldenrod2", "cadetblue3", "aquamarine4"), 
            add = "jitter") + yscale("log10", .format = TRUE) +   
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
}

# Figure S3
permutation_analysis_targeted <- function(data_)
{
  df = c()
  to_plot = c()
  
  for (i in as.character(unique(data_$Cancer.Type)))
  {
    data_tmp = data_[which(as.character(data_$Cancer.Type) == i),]
    all_ = permut_logrank(survival = as.numeric(as.character(data_tmp$SURVIVAL_MONTHS)),
                          tmb = as.numeric(as.character(data_tmp$TMB)),
                          censorship = as.numeric(as.character(data_tmp$SURVIVAL_EVENT)))
    real_all = logrank_percutoff(survival = as.numeric(as.character(data_tmp$SURVIVAL_MONTHS)),
                                 tmb = as.numeric(as.character(data_tmp$TMB)),
                                 censorship = as.numeric(as.character(data_tmp$SURVIVAL_EVENT)))
    
    to_plot_tmp = rbind(cbind(all_, i))
    to_plot_tmp = data.frame(to_plot_tmp, stringsAsFactors = F)
    colnames(to_plot_tmp) = c("p_value", "stratification")
    to_plot_tmp$p_value = as.numeric(as.character(to_plot_tmp$p_value))
    to_plot_tmp$stratification = (as.character(to_plot_tmp$stratification))
    
    df_tmp = data.frame(p_value = c(real_all), 
                        stratification = c(i))
    
    to_plot = rbind(to_plot, to_plot_tmp)
    df = rbind(df, df_tmp)
  }
  
  for (i in as.character(unique(data_$Cancer.Type)))
  {
    pval_all = sum(to_plot[which(to_plot$stratification == i),"p_value"] <= df[which(df$stratification == i), "p_value"])/ length(to_plot[which(to_plot$stratification == i),"p_value"])
    print(paste(i, pval_all))
  }
  
  to_plot[which(to_plot$p_value < 10^-10),1] = 10^-10
  to_plot2 = to_plot
  
  p<-ggplot(to_plot2, aes(x=stratification, y=p_value, fill=stratification, levels = stratification)) +
    # scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + 
    geom_violin(trim=FALSE) + yscale("log10", .format = TRUE) +
    geom_point(data = df, col = 'darkred', shape = 18, size = 5) + coord_flip(ylim = c(10^-20,1))
  print(p)
  
  return(pval_all)
}

# Figure S4
tmb_icb_19cancers <- function(data_)
{
  df_ = data.frame(data_)
  
  df_$General_Tumor_Type = gsub("Cancer", "", df_$General_Tumor_Type)
  df_$General_Tumor_Type = gsub("Carcinoma", "", df_$General_Tumor_Type)
  df_$General_Tumor_Type = gsub("  ", "", df_$General_Tumor_Type)
  
  rownames(df_) = df_$General_Tumor_Type
  
  mono_all = cor.test(df_$TMB_.median., df_$Single)
  print(paste("Monotherapy RR vs TMB for all cancers: p=", round(mono_all$p.value,4), ", r=", round(mono_all$estimate,4)))
  
  combi_all = cor.test(df_$TMB_.median., df_$Dual)
  print(paste("Combination therapy RR vs TMB for all cancers: p=", round(combi_all$p.value,4), ", r=", round(combi_all$estimate,4)))
  
  df_tmp = df_[-which(rownames(df_) %in% c("Melanoma", #"Ocular/Uveal melanoma",
                                                           "Colorectal– MSI", "Colorectal– MSS")), ]
  mono_all_no_outliers = cor.test(df_tmp$TMB_.median., df_tmp$Single)
  print(paste("Monotherapy RR vs TMB for all cancers except CRC and MELANOMA: p=", round(mono_all_no_outliers$p.value,4), ", r=", round(mono_all_no_outliers$estimate,4)))
  combi_all_no_outliers = cor.test(df_tmp$TMB_.median., df_tmp$Dual)
  print(paste("Combination therapy RR vs TMB for all cancers except CRC and MELANOMA: p=", round(combi_all_no_outliers$p.value,4), ", r=", round(combi_all_no_outliers$estimate,4)))
  
  p1 = ggplot(df_, aes(TMB_.median., Single)) + ylim(0,0.6) +
    geom_point(color='#2980B9', size = 4) + 
    theme_classic(base_size = 10) +
    scale_x_continuous(trans='log10', limits = c(1,50)) +
    scale_y_continuous(trans='log10', limits = c(0.005,0.6)) +
    geom_text_repel(aes(label = General_Tumor_Type, force = 1),
                    size = 3.5) + scale_x_continuous(trans='log10')
  
  p2 = ggplot(df_, aes(TMB_.median., Dual)) + ylim(0,0.6) +
    geom_point(color='#2980B9', size = 4) + 
    theme_classic(base_size = 10) +
    scale_x_continuous(trans='log10', limits = c(1,50)) +
    scale_y_continuous(trans='log10', limits = c(0.005,0.6)) +
    geom_text_repel(aes(label = General_Tumor_Type, force = 1),
                    size = 3.5) + scale_x_continuous(trans='log10')
  
 
  p = gridExtra::grid.arrange(p1, p2, nrow = 1, ncol =2)
  
  return(p)
}

# Figure S5
immunogenicity_model_params <- function()
{
    N=1:1000;
    
    par(mfrow = c(2,2))
    par(mar = c(2,2,1,1))
    
    # k variation weak
    cols = viridis(1000) 
    plot(1,1, log = "x", xlab = "TMB", ylab = "Probability of being immunogenic", type = "l", ylim = c(10^-6,1),
         xlim = c(1,10000), col = cols[1])
    p_tmp=0.22
    i=1
    for (kcut in seq(1, 100, 0.1))
    {
      P = 1- ppois(kcut, N*p_tmp)
      lines(N, P, col = cols[i], lwd = 2)
      i = i+1
    }
    
    # k variation strong
    cols = viridis(1000) 
    x = plot(1,1, log = "x", xlab = "TMB", ylab = "Probability of being immunogenic", type = "l", ylim = c(10^-6,1),
         xlim = c(1,10000), col = cols[1])
    p_tmp=0.64
    i=1
    for (kcut in seq(1, 100, 0.1))
    {
      P = 1- ppois(kcut, N*p_tmp)
      lines(N, P, col = cols[i], lwd = 2)
      i = i+1
    }
    
    # p variation kcut = 1
    plot(1,1, log = "x", xlab = "TMB", ylab = "Probability of being immunogenic", type = "l", ylim = c(10^-6,1),
         xlim = c(1,10000), col = cols[1])
    kcut=1
    cols = viridis(1000) 
    i=1
    for (p_tmp in seq(0.01,1,0.001))
    {
      P = 1- ppois(kcut, N*p_tmp)
      lines(N, P, col = cols[i], lwd = 2)
      i = i+1
    }
    
    # p variation kcut = 2
    plot(1,1, log = "x", xlab = "TMB", ylab = "Probability of being immunogenic", type = "l", ylim = c(10^-6,1),
         xlim = c(1,10000), col = cols[1])
    kcut=2
    cols = viridis(1000) 
    i=1
    for (p_tmp in seq(0.01,1,0.001))
    {
      P = 1- ppois(kcut, N*p_tmp)
      lines(N, P, col = cols[i], lwd = 2)
      i = i+1
    }
    par(mfrow = c(1,1))
}

