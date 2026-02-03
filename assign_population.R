# assign population labels to study data
#-------------------------
# 3. assign_population.R
# Assign super-population labels to study samples using a
#	reference trained Random Forest classifier at a
#	selected probability threshold.
# Input:
#   - study PCs based on reference PC loadings (PLINK --score output)
#   - Reference PC eigenvalues (used for PC rescaling)
#   - Trained RF model (from train_rf_classifier.R)
# Output:
#   - Table with per-sample ancestry probabilities and assignments
#   - PCA plots with assigned population labels
#---------------------------------------------

library(dplyr)
library(data.table)
library(ggplot2)

setwd("/path/to/dir")

RF_MODEL_FILE   <- "/path/to/rf_model.RData"
STUDY_PCS_FILE  <- "/path/to/study_pca.sscore"
EIGENVAL_FILE   <- "/path/to/ref_pca.eigenval"

N_PCS <- 10
PROB_THRESHOLDS <- c(0.9, 0.8, 0.7, 0.5) # could use one if you like

#-------- load the trained rf model 
load(RF_MODEL_FILE)

#--------- load study PCs
# Expected PLINK --score format: FID IID <other columns> PC1 PC2 ...
study_raw  <- read.table(STUDY_PCS_FILE, header = TRUE)

# identify PC columns (PLINK output can vary slightly - i don't actually know why)
pc_cols <- grep("^PC", names(study_raw), value = TRUE)

study_pcs <- study_raw %>%
  select(FID, IID, all_of(pc_cols[1:N_PCS]))

# ------------- rescale study PCs to reference scale
# Required when PCs are projected using allele weights
# See: https://www.cog-genomics.org/plink/2.0/score#pca_project
# Also Kate Shim's comment here ==> https://groups.google.com/g/plink2-users/c/W6DL5-hs_Q4/m/pMwsSYxtAwAJ
# Christopher Chang suggests another way to do this but i haven't tried it - https://groups.google.com/g/plink2-users/c/W6DL5-hs_Q4/m/b_o3JMrxAwAJ

# Here, each PC is divided by:
#   -sqrt(lambda_k) / 2 where lambda_k is the k-th eigenvalue from the reference PCA
# when you visualise study + ref PCs on the same plot you can see it 
# ps, you don't need this if you generated the PCs with the data (ref+study data) combined 

eigen_vals <- read.table(EIGENVAL_FILE)

study_pcs_adj <- study_pcs %>%
  mutate(across(
    starts_with("PC"),
    ~ .x / (-sqrt(eigen_vals[as.numeric(sub("PC", "", cur_column())), 1]) / 2)
  ))

#--------- predict ancestry probabilities
# apply rf_model 
pred_probs <- predict(rf_model, study_pcs_adj, type="prob")

study_pred <- bind_cols(study_pcs_adj, pred_prob)

#-------------assign population using probability thresholds
# For each individual, the random forest outputs posterior probabilities
# of belonging to each ancestry group

# so here assign an ancestry label only if the highest posterior probability
# exceeds a specified confidence threshold

# if no ancestry meets the threshold, the individual is left unassigned (NA)

# it generates separate columns for each thresholds e.g pred_0.9, pred_0.8, ...
# ps, thresholds were selected at the top (line 32??) of this script

assign_by_threshold <- function(df, threshold) {
  pop_cols <- colnames(pred_probs)

  df %>%
    mutate(
      assigned_pop = apply(
        select(., all_of(pop_cols)),
        1,
        function(x) {
          if (max(x) >= threshold) {
            pop_cols[which.max(x)]
          } else {
            NA
          }
        }
      )
    ) %>%
    rename(!!paste0("pred_", threshold) := assigned_pop)
}

for (t in PROB_THRESHOLDS) {
  study_pred <- assign_by_threshold(study_pred, t)
}

#--------- i like to see how the different thresholds separate samples
# colour Palette
palette <- c(
  "AFR" = "#E41A1C",
  "EUR" = "#377EB8", 
  "AMR" = "#4DAF4A", 
  "EAS" = "#FF7F00",
  "SAS" = "#984EA3",
  "NA" = "grey70"
)

for (t in PROB_THRESHOLDS) {

  col <- paste0("pred_", t)

  # PCA plot
  p <- ggplot(study_pred, aes(PC1, PC2, color = .data[[col]])) +
    geom_point(size = 1) +
    theme_classic() +
    scale_color_manual(values = palette, name = "Population") +
    labs(title = paste0("P > ", t))

  ggsave(paste0("prob", t, ".png"), p, width = 6, height = 4, dpi = 300)
}

# ------ save table of outputs - includes columns for all thresholds used
write.table(study_pred,file = "ancestry_assignments.tsv",sep = "\t",row.names = FALSE,quote = FALSE)
