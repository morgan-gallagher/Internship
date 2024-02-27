
############################# LOAD IN LIBRARIES ################################
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiome")

library(BiocManager)
library(microbiome)
library(dplyr)
library(ggplot2)
library(phyloseq)
install.packages(
  "microViz",
  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
)
library(microViz)

################################ GET DATA ######################################

#get data from microbiome
data(dietswap)
#print(dietswap)
#class(sample_data(dietswap))

#patient data frame
sam_data = sample_data(dietswap) |>
  data.frame()

#create count data frame
otu_df <- otu_table(dietswap) |>
  as.matrix() |>
  as.data.frame()

#taxonomy data frame 
taxonomy = tax_table(dietswap) |>
  as.data.frame()

################################ DATA CLEAN #####################################
#check if any na values
colSums(is.na(otu_df))


################################################################################


#abundance of microbes in data set
abundance_bacteria1 <- otu_df |>
  #sum across rows in data frame, sum across each bacteria for all samples
  rowSums(na.rm = TRUE, dims = 1) |>
  as.data.frame()


#rename column to counts  
abundance_bacteria1 <- abundance_bacteria1 |>
  rename(`counts` = `rowSums(otu_df, na.rm = TRUE, dims = 1)`)

#add column that contains names of bacteria
abundance_bacteria1$bacteria <- rownames(abundance_bacteria1)

##Bar graph that contains abundance of bacteria, bar not appropriate 
#p <- ggplot(abundance_bacteria1, aes(x = bacteria, y = counts)) +
#  geom_bar(stat = "identity", fill = "blue", alpha = 0.7) +
#  labs(title = "Bacteria Counts", x = "bacteria", y = "Counts") +
#  theme_minimal()
#p


#indices that are african american
african_amer_ind <- which(sam_data$nationality == "AAM")

#sample numbers that are african american
aam_sample_nums <- rownames(sam_data[african_amer_ind, ])


#tallys of bacterial in african americans
abundance_aam <- otu_df[, aam_sample_nums] |>
  rowSums(na.rm = TRUE, dims = 1) |>
  as.data.frame()


abundance_aam <- abundance_aam %>%
  rename(Counts = `rowSums(otu_df[, aam_sample_nums], na.rm = TRUE, dims = 1)`)

abundance_aam$Genus <- rownames(abundance_aam)

abundance_aam <- merge(abundance_aam,taxonomy,by="Genus")

aam_phylum_counts <- abundance_aam %>%
  group_by(Phylum) %>%
  summarize(total_counts = sum(Counts)) %>%
  arrange(desc(total_counts))

aam_family_counts <- abundance_aam %>%
  group_by(Family) %>%
  summarize(total_counts = sum(Counts)) %>%
  arrange(desc(total_counts))


#sample numbers for afr patients
afr_sample_nums <- rownames(data[which(data$nationality == "AFR" &
                                         (data$timepoint == 1 |
                                            data$timepoint == 2 |
                                            data$timepoint == 3))])

#tallys of bacterial in african americans
abundance_afr <- otu_df[, afr_sample_nums] |>
  rowSums(na.rm = TRUE, dims = 1) |>
  as.data.frame()


abundance_afr$Genus <- rownames(abundance_afr)

abundance_afr <- abundance_afr %>%
  rename(Counts = `rowSums(otu_df[, afr_sample_nums], na.rm = TRUE, dims = 1)`)

abundance_afr <- merge(abundance_afr,taxonomy,by="Genus")


afr_phylum_counts <- abundance_afr %>%
  group_by(Phylum) %>%
  summarize(total_counts = sum(Counts)) %>%
  arrange(desc(total_counts))

afr_family_counts <- abundance_afr %>%
  group_by(Family) %>%
  summarize(total_counts = sum(Counts)) %>%
  arrange(desc(total_counts))

#diet <- psmelt(dietswap)
#head(diet)


#plots proportion of each nationality that are lean, overweight, obese
plot_frequencies(x = data, 
                 Groups = "nationality", Factor = "bmi_group") +
  
  # below are ggplot2 functions for prettier visual
  labs(fill = "bmi_group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))


AAM_dataset <- diet |>
  filter(nationality == "AAM")

AFR_dataset <- diet |>
  filter(nationality == "AFR")

AAM_data <- subset_samples(dietswap, nationality == "AAM") |>
  core(detection = 50, prevalence = 50/100)

diet_relav <- transform(AAM_data, "compositional")


plot_composition(diet_relav,
                 taxonomic.level = "Genus",
                 average_by = "timepoint", # average by timepoint
                 otu.sort = "abundance",
                 x.label = "timepoint") +
  labs(x = "Time point",
       y = "Abundance",
       title = "African American") +
  theme(axis.text.x = element_text(angle = 0,
                                   hjust = 0.5
  ))

######THINGS FROM INTERNET######################################################

######MICROBIOME COMPOSITION####################################################

AAM_dataset <- subset_samples(dietswap, nationality == "AAM")
AFR_dataset <- subset_samples(dietswap, nationality == "AFR")

View(sam_data(AAM_dataset))

home_dataset <- subset_samples(dietswap, (timepoint == 1 | timepoint == 2 | timepoint == 3))


#proportions of bacteria in genuses
AAM_microbiome <- AAM_dataset %>%
  ps_filter(timepoint == 1) %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 15, other_name = "Other",
    taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
    palette = distinct_palette(n = 15, add = "grey90"),
    merge_other = TRUE, bar_outline_colour = "black"
  ) +
  coord_flip() +
  facet_wrap("bmi_group", nrow = 1, scales = "free") +
  ggtitle("Microbial Abundance by BMI Group - AAM")

AFR_microbiome <- AFR_dataset %>%
  ps_filter(timepoint == 1) %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 15, other_name = "Other",
    taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
    palette = distinct_palette(n = 15, add = "grey90"),
    merge_other = TRUE, bar_outline_colour = "black"
  ) +
  coord_flip() +
  facet_wrap("bmi_group", nrow = 1, scales = "free") +
  ggtitle("Microbial Abundance by BMI Group - AFR")

ALL_microbiome <- dietswap %>%
  ps_filter(timepoint == 1) %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 15, other_name = "Other",
    taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
    palette = distinct_palette(n = 15, add = "grey90"),
    merge_other = FALSE, bar_outline_colour = "darkgrey"
  ) +
  coord_flip() +
  facet_wrap("nationality", nrow = 1, scales = "free") +
  ggtitle("Microbial Abundance Across Nationality")

#library(ggpubr)
#ggarrange(AFR_microbiome, AAM_microbiome, nrow = 2, 
#          labels = c("  Native African","African American"), 
#          vjust = 1,
#          font.label = list(size = 12, face = "bold"))


p <- comp_barplot(dietswap,
                  n_taxa = 15, tax_level = "Genus",
                  bar_outline_colour = "black", merge_other = TRUE,
                  sample_order = "aitchison", group_by = "sex"
)

# plot them side by side with patchwork package
patch <- patchwork::wrap_plots(p, ncol = 2, guides = "collect")
patch & coord_flip() # make bars in all plots horizontal (note: use & instead of +)

############ ORDINATION PLOT ####################################################


############ ON ALL DATA ####################################################

core_diet <- core(x = dietswap, 
                  detection =  50, # treshold for presence for each OTU in a sample
                  prevalence = 50/100) # treshold for the fraction of OTUs that exceeds detection treshold
core_diet 

View(otu_table(core_diet))

diet_log <- transform_sample_counts(core_diet, function(x) log(1 + x))

ord <- ordinate(physeq = diet_log, 
                method = "MDS", 
                distance = "bray")
evals <- ord$values$Eigenvalues

pca_nationality <- plot_ordination(physeq = diet_log,
                ordination = ord, 
                color = "nationality",
                axes = c(4,5)) +
                geom_point(size = 2) +
                coord_fixed(sqrt(evals[2] / evals[1])) +
  ggtitle("PCA Plot by Nationality")


plot_scree(ordination = ord, axes = c(1:20))


pca_bmi <- plot_ordination(physeq = diet_log,
                                   ordination = ord, 
                                   color = "bmi_group") +
  geom_point(size = 2) +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  ggtitle("PCA Plot by BMI")


library(cluster)

sil_width <- silhouette(clustering_vector, dist = dist_matrix)

############ AGGREGATE FAMILY MICROVIZ   #######################################
unconstrained_aitchison_pca <- dietswap %>%
  tax_agg("Family") %>%
  tax_transform("clr") %>%
  ord_calc()

pca_plot <- unconstrained_aitchison_pca %>%
  ord_plot(
    plot_taxa = 1:6, colour = "bmi_group", size = 1.5,
    tax_vec_length = 0.325,
    tax_lab_style = tax_lab_style(max_angle = 90, aspect_ratio = 1),
    auto_caption = 8
  )

pca_plot

customised_plot <- pca_plot +
  stat_ellipse(aes(linetype = bmi_group, colour = bmi_group), linewidth = 0.3) +
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = "bottom") +
  coord_fixed(ratio = 1, clip = "off") # makes rotated labels align correctly

customised_plot

############ TEST VAR ASSOCIATION WITH OVERALL COMPOSITION   ###################

#permanova test
aitchison_dists <- dietswap %>%
  tax_transform("identity", rank = "Family") %>%
  dist_calc("aitchison")

# the more permutations you request, the longer it takes
# but also the more stable and precise your p-values become
aitchison_perm <- aitchison_dists %>%
  dist_permanova(
    seed = 1234, # for set.seed to ensure reproducibility of random process
    n_processes = 1, n_perms = 99, # you should use at least 999!
    variables = "bmi_group"
  )


# view the permanova results
perm_get(aitchison_perm) %>% as.data.frame()


# view the info stored about the distance calculation
info_get(aitchison_perm)


########################### ALPHA DIVERSITY   ##################################

tab <- dietswap |>
  ps_filter(timepoint == 1) |>
  alpha(index = "all")
  
  alpha(dietswap, index = "all")
head(tab)

p.shannon.nationality <- boxplot_alpha(ps_filter(dietswap, timepoint == 1), 
                                       index = "shannon", 
                                       x_var = "nationality", 
                                       fill.colors = c(AAM = "cyan4", AFR = 'deeppink4')) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +  # Remove both x-axis labels and title
  labs(y = "Shannon diversity")




p.shannon.sex <- boxplot_alpha(ps_filter(dietswap, timepoint == 1), 
                               index = "shannon", 
                               x_var = "sex", 
                               fill.colors = c(female = "cyan4", 
                                               male = 'deeppink4')) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = "nationality", y = "Shannon diversity")



p.shannon.bmi <- boxplot_alpha(ps_filter(dietswap, timepoint == 1), 
                               index = "shannon", 
                               x_var = "bmi_group", 
                               fill.colors = c(lean = "cyan4", 
                                               overweight = 'deeppink4',
                                               obese = "purple1")) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = "nationality", y = "Shannon diversity")

#questions to look into
#comparing counts of microbiome data by nationality
#stratifying by bmi group, microbiome 
#top phylum and family by nationality 
#most common phylym and family across the board


################################################################################
# TESTING BART: DATA PREPARATION

###Matrix Formation###

#filter data from home environment
he_data <- sam_data |>
  filter(timepoint == "1" |
           timepoint == "2" |
           timepoint == "3")

#filter data from swapped diet
swap_data <- sam_data |>
  filter(timepoint == "4" |
           timepoint == "5" |
           timepoint == "6")

#check how many samples collected per patient
frequency_df <- sam_data |>
  dplyr::group_by(subject)|>
  summarize(
    n = n()) |>
  as.data.frame()


#FUNCTION: summarize patient samples and nationalities
summarize_patient_samples <- function(data) {
  summarized_data <- data |>
    group_by(subject) |>
    summarize(samples = list(unique(sample)),
              nationality = first(nationality))
  
  return(summarized_data)
}

#FUNCTION: create otu average across each patient
create_otu_avg <- function(count_df, patient_samples, remove_byu = FALSE, remove_dwk = FALSE) {
  
  #create empty dataframe
  otu_df_avg <- data.frame(subject = character(), stringsAsFactors = FALSE)
  
  #remove byu from the home environ if using HE data from patient samples df
  if (remove_byu == TRUE) {
    patient_samples <- patient_samples |>
      filter(subject != 'byu')
  }
  
  if(remove_dwk == TRUE) {
    patient_samples <- patient_samples |>
      filter(subject != 'dwk')
  }
    
  #fill in the otu_avg df
  for (patient in 1:nrow(patient_samples)) {
    subject <- patient_samples$subject[patient]
    samples <- patient_samples$samples[[patient]]
    
    # subset otu_df for the patient samples, only bacteria columns
    subset_df <- count_df %>%
      select(one_of(samples), bacteria)
    
    # row-wise averages
    avg_counts <- rowMeans(subset_df[, -ncol(subset_df)], na.rm = TRUE)
    
    # dataframe for the patient with subject and average counts
    patient_df <- data.frame(subject = subject,t(avg_counts))
    
    # merge patient dataframe to otu_df_avg
    otu_df_avg <- bind_rows(otu_df_avg, patient_df)
  }
  
  nationality_data <- patient_samples[c("subject", "nationality")]
  
  #add nationality column
  otu_df_avg <- merge(otu_df_avg, nationality_data, by = "subject", all.x = TRUE)
  
  return(otu_df_avg)
  
}

#add bacteria column to otu_df
otu_df$bacteria <- rownames(otu_df)

#create summary of patients
patient_samples <- summarize_patient_samples(sam_data)

#create otu average df
otu_df_avg <- create_otu_avg(count_df = otu_df, patient_samples = patient_samples, remove_byu = FALSE)


################################################################################
# TESTING BART: SPLIT DATA

#Design matrix
X <- subset(otu_df_avg, select = -c(nationality, subject))

#response
y <- otu_df_avg$nationality

#shuffled response to break relationship
y_shuffle <- sample(y)

################################################################################
# TESTING BART: BUILD MODEL

options(java.parameters = "-Xmx5000m")
library(bartMachine)

#classical model
bart_machine <- bartMachine(X, y_shuffle)
bart_machine
bart_machine$misclassification_error
#out of sample model
cross_val <- k_fold_cv(X,y_shuffle, k_folds = 10)
cross_val$confusion_matrix
cross_val$misclassification_error

#using predict function
predict(bart_machine_cv, X[1:7, ])


#grid search over parameters
bart_machine_cv <- bartMachineCV(X, y_shuffle)
#win: k: 2, m:50
bart_machine_cv$misclassification_error


########################### VARIABLE SELECTION  ###############################

vs <- var_selection_by_permute(bart_machine, 
                               bottom_margin = 10, num_permute_samples = 10)
vs
vs$important_vars_local_col_nums

#variable selection for cross validated bart_machine
vs_cv <- var_selection_by_permute(bart_machine_cv,
                                  bottom_margin = 10, num_permute_samples = 10)

length(vs_cv$important_vars_local_col_nums)
vs_cv$important_vars_local_names

#bartmachine cross val only using selected variables
bart_machine_vs_cv <- bartMachineCV(X[ ,vs_cv$important_vars_local_col_nums], y)

#bartmachine vanilla using selected variables
bart_machine_vs <- bartMachine(X[ ,vs_cv$important_vars_local_col_nums], y)

#bartmachine out of sample using selected variables
bart_machine_kfold <- k_fold_cv(X[ ,vs_cv$important_vars_local_col_nums], y)

#bartmachine cross val only using global se var
bart_machine_global_se_cv <- bartMachineCV(X[ ,vs_cv$important_vars_global_se_col_nums], y)
bart_machine_global_se_cv$misclassification_error

#bartmachine vanilla only using global se var
bart_machine_global_se <- bartMachine(X[ ,vs_cv$important_vars_global_se_col_nums], y)
bart_machine_global_se$misclassification_error

#bartmachine kfold only using global se var
bart_machine_global_se_kfold <- k_fold_cv(X[ ,vs_cv$important_vars_global_se_col_nums], y)
bart_machine_global_se_kfold$misclassification_error

######################### VARIABLE SELECTION PLOTS #############################

#partial dependence plot
pd_plot(bart_machine_cv, j = "Clostridium.difficile.et.rel.")

#covariate importance
cov_importance_test(bart_machine_cv, covariates = "Tannerella.et.rel.")

#proportion of times each variable was used as splitting
investigate_var_importance(bart_machine_cv, num_replicates_for_avg = 20)


plot_convergence_diagnostics(bart_machine_cv)

#not working
#var_selection_by_permute_cv(bart_machine_cv)

######################### TESTING UNPROCESSED DATA #############################

#use otu data to make training data
X_all <- as.data.frame(t(otu_df[, -which(names(otu_df) == "bacteria")]))

#make target category
y_all <- sam_data$nationality

#shuffled data
y_all_shuffle <- sample(y_all)

#create model on unprocessed data
bart_machine <- bartMachine(X_all, y_all)
bart_machine$misclassification_error

#create model on unprocessed data out of sample
bart_machine_cv <- k_fold_cv(X_all,y_all_shuffle, k_folds = 10)
bart_machine_cv$misclassification_error

#################### CORRELATIONS BETWEEN BACTERIA #############################

#correlation table with bacteria
cc <- associate(t(otu_df[, -which(names(otu_df) == 'bacteria')]), 
                t(otu_df[, -which(names(otu_df) == 'bacteria')]), 
                method = 'spearman', 
                p.adj.threshold = 0.01, 
                n.signif = 5, 
                filter.self.correlations = TRUE)


#list the top 20 positively correlated bacteria
correlation_pos <- cc |>
  filter(Correlation > 0) |>
  head(20)

#list top 20 negatively correlated bacteria
correlation_neg <- cc|>
  filter(Correlation < 0) |>
  head(20)

set.seed(11)
n = 200
p = 10
X = data.frame(matrix(runif(n * p), ncol = p))
y = 10 * sin(pi* X[ ,1] * X[,2]) +20 * (X[,3] -.5)^2 + 10 * X[ ,4] + 5 * X[,5] + rnorm(n)
##build BART regression model
bart_machine = bartMachine(X, y, num_trees = 20)
#investigate interactions
interaction_investigator(bart_machine)


########################## TIME PERIODS ########################################

#transform raw data into proportions
proportional_table <- as.data.frame(prop.table(as.matrix(otu_df[, !colnames(otu_df) %in% "bacteria"]), margin = 2))

#add bacteria column
proportional_table$bacteria <- rownames(otu_df)

#average patient prop so one observation per patient
patient_samples <- summarize_patient_samples(he_data)
otu_df_avg <- create_otu_avg(count_df = proportional_table, patient_samples = patient_samples, remove_byu = TRUE, remove_dwk = TRUE)

#create swapped proportion dataframe
patient_samples2 <- summarize_patient_samples(swap_data)

View(swap_data)

target_otu_avg <- create_otu_avg(count_df = proportional_table, patient_samples = patient_samples2, remove_byu = TRUE, remove_dwk = TRUE)

  
#create empty dataframe
target_otu_avg2 <- data.frame(subject = character(), stringsAsFactors = FALSE)
  
patient_samples <- patient_samples[patient_samples$subject != 'dwk', ]
patient_samples <- patient_samples[patient_samples$subject != 'byu', ]
  

View(patient_samples)
patient_samples$samples


  #fill in the otu_avg df
  for (patient in 1:nrow(patient_samples)) {
    subject <- patient_samples$subject[patient]
    samples <- patient_samples$samples[[patient]]
    print(samples)
    
    # subset otu_df for the patient samples, only bacteria columns
    subset_df <- proportional_table %>%
      select(one_of(samples), bacteria)
    
    # row-wise averages
    avg_counts <- rowMeans(subset_df[, -ncol(subset_df)], na.rm = TRUE)
    
    # dataframe for the patient with subject and average counts
    patient_df <- data.frame(subject = subject,t(avg_counts))
    
    # merge patient dataframe to otu_df_avg
    target_otu_avg2 <- bind_rows(target_otu_avg2, patient_df)
  }
  

#check same subjects in both
target_otu_avg2$subject == target_otu_avg$subject


#before and after proportion matrices
target = target_otu_avg[, -1]
microb = target_otu_avg2[, -1]

target <- target[, !names(target) %in% "nationality"]

#substitute epsilon for zero values
epsilon <- 1e-10

#log transform data
microb_log <- round(log(microb + epsilon), digits = 5)
target_log <- round(log(target + epsilon), digits = 5)

#find difference of log
difference_log = target_log - microb_log
#percent_diff_log = difference_log / microb_log

#factorize response variable
row_sums <- rowSums(abs(difference_log))
intervals <- seq(0, 220, by = 10)
factor_intervals <- cut(row_sums, breaks = intervals, labels = paste0("(", intervals[-length(intervals)], ", ", intervals[-1], "]"))
factorized_row_sums <- factor(factor_intervals)

#turn factor levels into numeric values for response
factor_levels <- levels(factorized_row_sums)
numeric_values <- seq_along(factor_levels)
numeric_row_sums <- numeric_values[as.integer(factorized_row_sums)]



#BART
#feature matrix: original averaged microbiome, target_otu_avg2
#target: factorized_row_sums
bart_machine_perc_change <- mbart2(x.train = target_otu_avg2[, -1], y.train = numeric_row_sums, nskip = 1, ndpost = 1, type = "lbart")

View(bart_machine_perc_change)

bart_machine_perc_change$prob.train

predict(bart_machine_perc_change, target_otu_avg2[, -1])

dim(target_otu_avg2[, -1])
length(factorized_row_sums)



