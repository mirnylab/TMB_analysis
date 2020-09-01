library("cutpointr")
library(data.table)
library(viridis)
library(ggplot2)
library(ggpubr)
library(maftools)
library(survival)
library(ggrepel)


source("/home/cgurjao/Dropbox (Partners HealthCare)/Mirnylab/Function_GURJAO_2020.R")

data_tmb19cancers = read.delim("/home/cgurjao/Dropbox (Partners HealthCare)/Mirnylab/Github_submission/data_tmb19cancers.csv", stringsAsFactors = F)
data_tmb_targeted = read.delim("/home/cgurjao/Dropbox (Partners HealthCare)/Mirnylab/Github_submission/data_tmb_targeted.csv", stringsAsFactors = F)
data_tmb_wes = read.delim("/home/cgurjao/Dropbox (Partners HealthCare)/Mirnylab/Github_submission/data_tmb_wes", stringsAsFactors = F)
data_copd = read.delim("/home/cgurjao/Dropbox (Partners HealthCare)/Mirnylab/Github_submission/data_copd", stringsAsFactors = F)

mel1 = data_tmb_wes[which(data_tmb_wes$dataset == "mel1"),]
mel2 = data_tmb_wes[which(data_tmb_wes$dataset == "mel2"),]
lung1 = data_tmb_wes[which(data_tmb_wes$dataset == "lung1"),]
lung2 = data_tmb_wes[which(data_tmb_wes$dataset == "lung2"),]

# Figure 1A
plot_all_cancers(data_tmb_wes)

# Figure 1B
mel1_stratif = plot_stratif(mel1, title = "mel1")
mel2_stratif = plot_stratif(mel2, title = "mel2")
p = gridExtra::grid.arrange(mel1_stratif, mel2_stratif, 
                            nrow = 1, ncol = 2)

# Figure 1C
lung1_stratif = plot_stratif(lung1, title = "lung1")
lung2_stratif = plot_stratif(lung2, title = "lung2")
p = gridExtra::grid.arrange(lung1_stratif, lung2_stratif, 
                            nrow = 1, ncol = 2)

# Figure 2A
mel1_raw = plot_raw_data(mel1, title = "mel1")
mel2_raw = plot_raw_data(mel2, title = "mel2")
lung1_raw = plot_raw_data(lung1, title = "lung1")
lung2_raw = plot_raw_data(lung2, title = "lung2")
gridExtra::grid.arrange(mel1_raw, mel2_raw,lung1_raw, lung2_raw, nrow = 2, ncol =2)

# Figure 2B
NA

# Figure 2C
mel1_permut = permutation_analysis(clinical.data_ = mel1,
                                   stratification_1 = ("skin/ occult"),
                                   stratification_2 = ("acral/ mucosal"),
                                   survival_type = "PFS", 
                                   survival_censor = "PFS_censorship")
mel2_permut = permutation_analysis(clinical.data_ = mel2,
                                   stratification_1 = ("skin/ occult"),
                                   stratification_2 = ("acral/ mucosal"),
                                   survival_type = "PFS", 
                                   survival_censor = "PFS_censorship")

# Figure 2D
lung1_permut = permutation_analysis(clinical.data_ = lung1,
                                    stratification_1 = ("former/ current"),
                                    stratification_2 = ("never"),
                                    survival_type = "PFS", 
                                    survival_censor = "PFS_censorship")
lung2_permut = permutation_analysis(clinical.data_ = lung2,
                                    stratification_1 = ("former/ current"),
                                    stratification_2 = ("never"),
                                    survival_type = "PFS", 
                                    survival_censor = "PFS_censorship")


# Figure 3A
youden_indexes = c()
youden_indexes[1:2] = plot_auc(rbind(mel1, mel2))[3:4]

# Figure 3B
youden_indexes[3:4] = plot_auc(rbind(lung1, lung2))[3:4]

# Figure 3C
FDA_youden_cutoffs(youden_indexes)

# Figure 3D
misclassified_pats(youden_indexes, mel1, mel2, lung1, lung2)

# Figure 4A
NA

# Figure 4B
immunogenicity_model()

# Figure S1
copd_tcga_analysis(data_copd)

# Figure S2A
permut_mel1_OS = permutation_analysis(clinical.data_ = mel1,
                                   stratification_1 = ("skin/ occult"),
                                   stratification_2 = ("acral/ mucosal"), 
                                   survival_type = "OS", 
                                   survival_censor = "OS_censorship")
permut_mel2_OS = permutation_analysis(clinical.data_ = mel2,
                                   stratification_1 = ("skin/ occult"),
                                   stratification_2 = ("acral/ mucosal"), 
                                   survival_type = "OS", 
                                   survival_censor = "OS_censorship")

# Figure S2B
permut_lung1_OS = permutation_analysis(clinical.data_ = lung1,
                                    stratification_1 = ("former/ current"),
                                    stratification_2 = ("never"), 
                                    survival_type = "OS", 
                                    survival_censor = "OS_censorship")

# Figure S3
permut_targeted = permutation_analysis_targeted(data_tmb_targeted)

# Figure S4
tmb_icb_19cancers(data_tmb19cancers)

# Figure S5
immunogenicity_model_params()
