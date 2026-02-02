# Genetic ancestry assignment using PCA + Random Forest
Assign population labels to study samples using principal components derived from genotype data and a supervised Random Forest classifier trained on reference populations  (I used 1000 Genomes).

- Perform QC and SNP selection between reference and study data
- Compute reference PCs and project study samples in the same PC space
- Train a supervised population classifier
- Assign ancestry labels probabilistically with configurable confidence thresholds

## Structure  
1. qc_prune_pca.sh ==> QC, variant ID harmonisation and PCA
2. train_rf_classifier.R ==> Train Random Forest classifier
3. assign_populations.R ==> Assign population labels to study samples


> PS, this approach is intended for global ancestry stratification (AFR, EUR, EAS, SAS, AMR), not fine-scale population inference or admixture estimation.
