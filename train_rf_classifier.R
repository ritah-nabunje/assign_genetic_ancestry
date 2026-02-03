# Be ware of typos!! 

#-----------------------
# 2. train_rf_classifier.R
#Train a Random Forest classifier to predict ancestry super-populations
#   using principal components derived from a reference panel (e.g. 1000G)
# Input:
#   - PCA eigenvectors from the reference panel
#   - Reference sample metadata with population labels
# Output:
#   - Trained RF model (RData)
#   - PCA plot of reference samples
#   - Model accuracy plot
#   - Confusion matrix heatmap
#-------------------------------

library(caret)
library(randomForest)
library(pROC)
library(plyr)
library(ggplot2)

setwd("/path/to/dir")

REF_PCS_FILE   <- "/path/to/ref_pca.eigenvec"
REF_META_FILE  <- "/path/to/20130606_1000GP_sample_info.csv"

N_PCS <- 10
SEED  <- 123 # for reproducibility

#------- load reference PCs 
# expected PLINK format: FID IID PC1 PC2...PCN
ref_pcs<- read.table(REF_PCS_FILE, header=TRUE, comment.char="")
colnames(ref_pcs)[1] <- "FID"

#------- load reference metadata 
ref_meta_raw  <-read.csv(REF_META_FILE)

# keep required columns
# 1 = sample ID and 3 = Population codes (e.g. YRI, GBR)
ref_meta <- ref_meta_raw[,c(1,3)]
colnames(ref_meta) <- c("V2", "population")

# merge pcs with metadata 
ref_data  <- merge(ref_pcs, ref_meta, by="V2")

# map populations to super-populations
# based on 1000 Genomes population definitions:
# more info on the population codes: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/README_populations.md
pop_to_superpop <- list(
  AFR = c("ACB", "ASW", "GWD", "ESN", "LWK", "MSL", "YRI"),
  EUR = c("GBR", "CEU", "FIN", "TSI", "IBS"),
  EAS = c("CDX", "CHB", "CHS", "JPT", "KHV"),
  SAS = c("BEB", "GIH", "ITU", "PJL", "STU"),
  AMR = c("CLM", "MXL", "PEL", "PUR")
)

ref_data$superpop <- NA

for (sp in names(pop_to_superpop)) {
  ref_data$superpop[ref_data$population %in% pop_to_superpop[[sp]]] <- sp
}

ref_data$superpop <- as.factor(ref_data$superpop)

# visualise reference PCA
p <- ggplot(ref_data, aes(PC1, PC2, color = superpop)) +
  geom_point()
ggsave("ref_pca.png", p, width = 6, height = 4, dpi = 300)

#-------------Training/test data split
set.seed(SEED)  # for reproducibility

splitIndex <- createDataPartition(ref_data$superpop, p = 0.8, list = FALSE)

train_data <- ref_data[splitIndex, ]
test_data <- ref_data[-splitIndex, ]

#--------------- Random Forest training
# define the training control for cross validation
train_control <- trainControl(method = "repeatedcv",
                              number = 5,
                              repeats = 5,
                              classProbs = TRUE)

# parameter tuning - how many trees to try
tuning_parameters <- expand.grid(mtry = c(1:N_PCS))

# build/train RF model
set.seed(SEED)
#rf_model <- train(superpop ~ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
#                  data = train_data,
#                  method = "rf",
#                  trControl = train_control,
#                  verbose = TRUE,
#                  tuneGrid=tuning_parameters,
#                  metric="Accuracy")

rf_model <- train(
  superpop ~ .,
  data = train_data[, c(paste0("PC", 1:N_PCS), "superpop")],
  method = "rf",
  trControl = train_control,
  tuneGrid = tune_grid,
  metric = "Accuracy"
)

print(rf_model)

# ----------- model diagnostics
png(file="rf_model_accuracy.png",width=600, height=400)
plot(rf_model, metric = "Accuracy")
dev.off()

#-------------model performance on test data
test_preds  <- predict(rf_model, test_data)

# confusion matrix (comparing predictions with original labels)
conf_matrix <- confusionMatrix(test_preds, test_data$superpop)
print(conf_matrix)

conf_df <- as.data.frame(conf_matrix$table)
conf_df$Percentage <- conf_df$Freq / colSums(conf_table)[conf_df$Reference] * 100

# make heatmap of the confusion matrix
ggplot(conf_df, aes(x = Reference, y = Prediction, fill = Percentage)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high="#009194", limits = c(0, 100)) +
  geom_text(aes(label = sprintf("%.1f", Percentage)), vjust = 1) +
  theme_bw() +
  labs(title = "confusion matrix")
ggsave("confusion_matrix_rf.png", width = 6, height = 4, dpi = 300)

# save the trained model - if satsfied with the performace
save(rf_model, file="rf_model.RData")
