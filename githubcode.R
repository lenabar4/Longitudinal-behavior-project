#################################################################################
##########################Longitudinal behavior project##########################
#################################################################################

#Project: Early childhood multi-fluid metabolite signatures of emotional and behavioral 
#symptoms and their trajectories to adolescence: the Rhea Mother–Child Cohort 
#Author: Eleni Barmpa 
#### Libraries ####
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(gt)
library(writexl)
library(openxlsx)

library(stats)
library(glmnet)
library(lme4)
library(lmerTest)
library(broom.mixed)

library(mice)
library(miceadds)
library(randomForest)
library(missForest)

library(kSamples)
library(RNOmni)
library(corrplot)
library(pheatmap)
library(grid)

##########Setting working directory##########----
getwd()
setwd('set_path_to_data_folder')

##########Metabolomic dataset creation##########----
##########Serum metabolomic dataset preprocessing##########----

###Loading datasets (Serum metabolomic dataset from HELIX)
# Participants diagnosed as ADS-positive did not participate
# in the internalizing and externalizing assessments.
# These individuals are excluded from the analytic dataset.
# We remove these helixids 

RHEA_serum_data <- read_excel("data/datasets/Preprocessing/RHEA_serum_data.xlsx")

RHEA_serum_data <- RHEA_serum_data %>%
  filter(!helixid %in% c("RHExxxxxx", "RHExxxxxx", "RHExxxxxx"))  # Exclude specific helixid values
str(RHEA_serum_data)
summary(RHEA_serum_data)

# Calculate Median and format Lower Quartile (25%) and Upper Quartile (75%) in parentheses
# of each serum metabolite

summary_table_serum <- RHEA_serum_data %>%
  select(where(is.numeric)) %>%  # Ensure only numeric columns are summarized
  summarise(across(everything(),
                   list(Median = ~ round(median(., na.rm = TRUE), 3),
                        LowerQ = ~ round(quantile(., 0.25, na.rm = TRUE), 3),
                        UpperQ = ~ round(quantile(., 0.75, na.rm = TRUE), 3)),
                   .names = "{col}_{fn}")) %>%
  pivot_longer(everything(),
               names_to = "Metabolite_Stat",
               values_to = "Value") %>%
  separate(Metabolite_Stat, into = c("Metabolite", "Stat"), sep = "_") %>%
  pivot_wider(names_from = Stat, values_from = Value) %>%
  mutate(Summary = paste0(Median, " (", LowerQ, "-", UpperQ, ")")) %>%
  select(Metabolite, Summary)

View(summary_table_serum)
#Save the dataset with the summary statostics of the serum metabolites 
#write_xlsx(summary_table_serum, "d1_serum_summary_conc.xlsx")

#Categorise each serum metabolite based on the general category 

metabolite_columns <- colnames(RHEA_serum_data)[-1] 
metabolite_columns

metabolites <- as.data.frame(metabolite_columns)

# Update vectors for each category with the missing metabolites
amino_acids <- c("Ala", "Arg", "Asn", "Asp", "Cit", "Glu", "Gln", "Gly", "His", "Ile", "Leu", 
                 "Lys", "Met", "Orn", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")

biogenic_amines <- c("Ac-Orn", "alpha-AAA", "ADMA", "Carnosine", "Creatinine", "DOPA", "Dopamine", 
                     "Histamine", "Kynurenine", "Met-SO", "Nitro-Tyr", "PEA", "Putrescine", 
                     "Sarcosine", "Serotonin", "Spermidine", "Spermine", "c4-OH-Pro", "t4-OH-Pro", 
                     "SDMA", "Taurine", "total DMA")  # Added 'total DMA'

monosaccharides <- c("H1")

acylcarnitines <- c("C0", "C10:1", "C10:2", "C2", "C3", "C3:1", "C3-OH", "C3-DC (C4-OH)", 
                    "C4", "C4:1", "C4-OH (C3-DC)", "C5", "C5:1", "C5:1-DC", "C5-DC (C6-OH)", 
                    "C5-M-DC", "C5-OH (C3-DC-M)", "C6", "C6:1", "C6 (C4:1-DC)",  # Added 'C6 (C4:1-DC)'
                    "C7-DC", "C8", "C9", "C10", "C12", "C12-DC", "C12:1", "C14", "C14:1", 
                    "C14:1-OH", "C14:2", "C14:2-OH", "C16", "C16:1", "C16:1-OH", "C16:2", 
                    "C16:2-OH", "C16-OH", "C18", "C18:1", "C18:1-OH", "C18:2")

glycerophospholipids <- c("lysoPC a C14:0", "lysoPC a C16:0", "lysoPC a C16:1", "lysoPC a C17:0", 
                          "lysoPC a C18:0", "lysoPC a C18:1", "lysoPC a C18:2", "lysoPC a C20:3", 
                          "lysoPC a C20:4", "lysoPC a C24:0", "lysoPC a C26:0", "lysoPC a C26:1", 
                          "lysoPC a C28:0", "lysoPC a C28:1", "PC aa C24:0", "PC aa C26:0", "PC aa C28:1", 
                          "PC aa C30:0", "PC aa C30:2", "PC aa C32:0", "PC aa C32:1", "PC aa C32:2", 
                          "PC aa C32:3", "PC aa C34:1", "PC aa C34:2", "PC aa C34:3", "PC aa C34:4", 
                          "PC aa C36:0", "PC aa C36:1", "PC aa C36:2", "PC aa C36:3", "PC aa C36:4", 
                          "PC aa C36:5", "PC aa C36:6", "PC aa C38:0", "PC aa C38:1", "PC aa C38:3", 
                          "PC aa C38:4", "PC aa C38:5", "PC aa C38:6", "PC aa C40:1", "PC aa C40:2", 
                          "PC aa C40:3", "PC aa C40:4", "PC aa C40:5", "PC aa C40:6", "PC aa C42:0", 
                          "PC aa C42:1", "PC aa C42:2", "PC aa C42:4", "PC aa C42:5", "PC aa C42:6",
                          "PC ae C30:0", "PC ae C30:1", "PC ae C30:2", "PC ae C32:1", "PC ae C32:2", 
                          "PC ae C34:0", "PC ae C34:1", "PC ae C34:2", "PC ae C34:3", "PC ae C36:0", 
                          "PC ae C36:1", "PC ae C36:2", "PC ae C36:3", "PC ae C36:4", "PC ae C36:5", 
                          "PC ae C38:0", "PC ae C38:1", "PC ae C38:2", "PC ae C38:3", "PC ae C38:4", 
                          "PC ae C38:5", "PC ae C38:6", "PC ae C40:1", "PC ae C40:2", "PC ae C40:3", 
                          "PC ae C40:4", "PC ae C40:5", "PC ae C40:6", "PC ae C42:0", "PC ae C42:1", 
                          "PC ae C42:2", "PC ae C42:3", "PC ae C42:4", "PC ae C42:5", "PC ae C44:3", 
                          "PC ae C44:4", "PC ae C44:5", "PC ae C44:6")

sphingolipids <- c("SM (OH) C14:1", "SM (OH) C16:1", "SM (OH) C22:1", "SM (OH) C22:2", "SM (OH) C24:1", 
                   "SM C16:0", "SM C16:1", "SM C18:0", "SM C18:1", "SM C20:2", "SM C24:0", "SM C24:1", 
                   "SM C26:0", "SM C26:1")

# Assign categories based on metabolite names
metabolites$Category <- ifelse(metabolites$metabolite_columns %in% amino_acids, "Amino Acids",
                               ifelse(metabolites$metabolite_columns %in% biogenic_amines, "Biogenic Amines", 
                                      ifelse(metabolites$metabolite_columns %in% monosaccharides, "Monosaccharides", 
                                             ifelse(metabolites$metabolite_columns %in% acylcarnitines, "Acylcarnitines", 
                                                    ifelse(metabolites$metabolite_columns %in% glycerophospholipids, "Glycerophospholipids", 
                                                           ifelse(metabolites$metabolite_columns %in% sphingolipids, "Sphingolipids", NA))))))

# View the updated data frame with the new category column
print(metabolites)

# Merge categories back with the RHEA_serum_data
RHEA_serum_data_with_categories <- RHEA_serum_data %>%
  select(helixid, everything())

summary_table_serum_metabolites_category <- RHEA_serum_data_with_categories %>%
  select(-helixid) %>%  # Remove helixid column
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Value") %>%  # Reshape data
  left_join(metabolites, by = c("Metabolite" = "metabolite_columns")) %>%  # Add category information
  group_by(Metabolite, Category) %>%  # Group by Metabolite and Category
  summarise(
    Mean = round(mean(Value, na.rm = TRUE), 3),
    Median = round(median(Value, na.rm = TRUE), 3),
    LowerQ = round(quantile(Value, 0.25, na.rm = TRUE), 3),  # 25% IQR
    UpperQ = round(quantile(Value, 0.75, na.rm = TRUE), 3),  # 75% IQR
    SD = round(sd(Value, na.rm = TRUE), 2),  # Standard deviation
    n = sum(!is.na(Value)),  # Sample size
    CI_Lower = round(Mean - 1.96 * (SD / sqrt(n)), 3),  # Confidence Interval (Lower)
    CI_Upper = round(Mean + 1.96 * (SD / sqrt(n)), 3),  # Confidence Interval (Upper)
    .groups = 'drop'  # Remove grouping
  ) %>%
  mutate(Summary = paste0(Median, " (", LowerQ, "-", UpperQ, ")")) %>%  # Format IQR in the same column as Median
  select(Metabolite, Category, Mean, Summary, SD, n, CI_Lower, CI_Upper)  # Reorder columns

# Print the final summary table
print(summary_table_serum_metabolites_category)

# Save the summary table as an Excel file
#write_xlsx(summary_table_serum_metabolites_category, "d1_serum_summary_conc.xlsx")

# Count the number of metabolites in each category
category_counts <- metabolites %>%
  group_by(Category) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))  # Optional: Sort by count
# Print the summary
print(category_counts)


# Create a bar graph for metabolite counts per category - Not used 
ggplot(category_counts, aes(x = reorder(Category, -Count), y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Metabolite Counts per Category",
    x = "Category",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_brewer(palette = "Set3")

#####################################################################################################
#####################################################################################################
# Check the data  

#Check for duplicate ids 
duplicate_helixids_serum <- RHEA_serum_data %>%
  group_by(helixid) %>%
  filter(n() > 1)
# Print duplicate helixids, if any
if (nrow(duplicate_helixids_serum) > 0) {
  print("Duplicate helixid values found:")
  print(duplicate_helixids_serum)
} else {
  print("No duplicate helixid values found.")
}

#calculate the negative values if they exist 
negative_values_per_column_serum <- colSums(RHEA_serum_data[, -1] < 0, na.rm = TRUE)
# Print the result
print(negative_values_per_column_serum)
ids_with_negatives_serum <- RHEA_serum_data$helixid[apply(RHEA_serum_data[, -1] < 0, 1, any, na.rm = TRUE)]
# Print the IDs
cat("IDs with negative values:\n")
print(ids_with_negatives_serum)

# Calculate the total number of missing values
total_missing_values_serum <- sum(is.na(RHEA_serum_data))
cat("Total number of missing values in the dataset:", total_missing_values_serum, "\n")# Calculate the total number of zero values in the dataset, excluding the first column (ID)

# Calculate the total number of zero values 
total_zero_values_serum <- sum(RHEA_serum_data[, -1] == 0, na.rm = TRUE)
cat("Total number of zero values in the dataset (excluding ID column):", total_zero_values_serum, "\n")

#Which columns have zero values?
columns_with_zero_values <- names(RHEA_serum_data)[-1][colSums(RHEA_serum_data[, -1] == 0, na.rm = TRUE) > 0]
columns_with_zero_values ##only 2 columns 

#How many zero values there are in every column and where 
for (column in columns_with_zero_values) {
  cat("Column", column, "has", sum(RHEA_serum_data[[column]] == 0, na.rm = TRUE), "zero values\n")
}
#####Column PC aa C30:2 has 21 zero values
#####Column PC ae C38:1 has 2 zero values


#Find the exact locations of zero values 
zero_locations <- which(RHEA_serum_data == 0, arr.ind = TRUE)
zero_locations <- zero_locations[zero_locations[,2] != 1, ] # Exclude ID column
print(zero_locations)

# Create a new dataset only with the columns that have zero values
serum_columns_with_zero_values <- RHEA_serum_data %>%
  select(helixid, all_of(columns_with_zero_values))
View(serum_columns_with_zero_values)


# Filter out zero values for "PC aa C30:2"
filtered_data_PC30_2 <- serum_columns_with_zero_values %>%
  filter(`PC aa C30:2` != 0) %>%
  select(helixid, `PC aa C30:2`)

# Create a density plot for "PC aa C30:2" (excluding zero values)
median_before <- median(filtered_data_PC30_2$`PC aa C30:2`)

p1 <- ggplot(filtered_data_PC30_2, aes(x = `PC aa C30:2`)) +
  geom_density(fill = "red", alpha = 0.5) +
  labs(title = "Density plot for PC aa C30:2 (excluding zero values)",
       x = "PC aa C30:2",
       y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
p1

# Filter out zero values for "PC ae C38:1"
filtered_data_PC38_1 <- serum_columns_with_zero_values %>%
  filter(`PC ae C38:1` != 0) %>%
  select(helixid, `PC ae C38:1`)

# Create a density plot for "PC ae C38:1" (excluding zero values)
p2 <- ggplot(filtered_data_PC38_1, aes(x = `PC ae C38:1`)) +
  geom_density(fill = "red", alpha = 0.5) +
  labs(title = "Density plot for PC ae C38:1 (excluding zero values)",
       x = "PC ae C38:1",
       y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
# Print both plots
print(p1)
print(p2)


#Replace the zero with NAs
serum_NA_replaced <- RHEA_serum_data
serum_NA_replaced[serum_NA_replaced == 0] <- NA
View(serum_NA_replaced)

# Calculate the total number of missing values in the dataset, excluding the first column (ID)
total_NA <- sum(is.na(serum_NA_replaced[, -1]))
cat("Total number of missing values in the dataset (excluding ID column):", total_NA, "\n")

# Calculate the percentage of missing values for each metabolite
missing_percentage_serum <- colSums(is.na(serum_NA_replaced[, -1])) / nrow(serum_NA_replaced) * 100
missing_df <- data.frame(
  Metabolite = names(missing_percentage_serum),
  Missingness = missing_percentage_serum
)
missing_df
missing_df_filtered <- missing_df %>% filter(Missingness > 0)
missing_df_filtered

#Metabolite Missingness
#PC aa C30:2 PC aa C30:2        10.5
#PC ae C38:1 PC ae C38:1         1.0

# Plot the percentage of missing values using ggplot2
ggplot(missing_df_filtered, aes(x = reorder(Metabolite, -Missingness), y = Missingness)) +
  geom_bar(stat = "identity", fill = "blue") +
  geom_hline(yintercept = 20, color = "red", linetype = "dashed") + # Add the cutoff line
  annotate("text", x = Inf, y = 22, label = "Percentage 20", hjust = 1.1, color = "red") + # Add annotation
  theme_minimal() +
  labs(title = "Percentage of Missing Values for Each Metabolite (Missingness > 0%)",
       x = "Metabolite",
       y = "Percentage of Missingness") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


#############################################################################################
#############################################################################################
#################Imputation of missing values ###############################################

# Filter columns with NAs
missing_df_filtered <- missing_df %>% filter(Missingness > 0)
columns_with_NAs <- missing_df_filtered$Metabolite

# Create a new data frame with columns that have NAs
RHEA_columns_with_NAs <- serum_NA_replaced %>%
  select(helixid, all_of(columns_with_NAs))
View(RHEA_columns_with_NAs)
str(RHEA_columns_with_NAs)

# Metabolites with >0% and <20% missing values are retained
# and will be imputed using a Random Forest imputation approach.
# Metabolites with ≥20% missingness are excluded from analysis.

data_to_impute <- serum_NA_replaced[, -1]
data_to_impute
data_to_impute= data_to_impute %>% mutate_if(is.character, as.numeric)
data_to_impute = as.data.frame(sapply(data_to_impute, function(x) as.numeric(x)))
data_to_impute = data_to_impute %>% mutate_if(is.character, as.numeric)
View(data_to_impute)
str(data_to_impute)

# Reproducibility: Set Random Seed
set.seed(1)
imputed_data <- missForest::missForest(data_to_impute)
RHEA_imputed_data <-imputed_data$ximp
RHEA_imputed_data

# Add the helixid column to the imputed data
# Assuming you have a vector with helixid values
helixid <- serum_NA_replaced$helixid  # Extracting helixid from original dataset
RHEA_imputed_data <- as.data.frame(RHEA_imputed_data) 
RHEA_imputed_data <- cbind(helixid, RHEA_imputed_data)

#now the imputed data are ready
#save the imputed data set 
#write.xlsx(RHEA_imputed_data, file = "RHEA_imputed_serum_data.xlsx")

# Inspect Imputed Values
# ----------------------------------------------------------
# Extract and examine imputed values from the multiple
# imputation object to evaluate plausibility and range.
imputed_values <- data.frame(
  Row = zero_locations,
  Column = colnames(data_to_impute)[zero_locations[, 2]],
  Original_NA = NA,
  Imputed_Value = apply(zero_locations, 1, function(loc) RHEA_imputed_data[loc[1], loc[2]])
)
# Print the imputed values
print(imputed_values)

# Create the density plots one by one before and after imputation 
RHEA_imputed_data_columns <- RHEA_imputed_data %>%
  select(helixid, `PC aa C30:2`, `PC ae C38:1`)
View(RHEA_imputed_data_columns)
str(RHEA_imputed_data_columns)
median_after <- median(RHEA_imputed_data_columns$`PC aa C30:2`, na.rm = TRUE)

p3 <- ggplot(RHEA_imputed_data_columns, aes(x = `PC aa C30:2`)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Density Plot for PC aa C30:2 after imputation", 
       x = "PC aa C30:2", 
       y = "Density") +
  theme_minimal()

# Create density plot for 'PC ae C38:1'
p4 <- ggplot(RHEA_imputed_data_columns, aes(x = `PC ae C38:1`)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Density Plot for PC ae C38:1  after imputation", 
       x = "PC ae C38:1", 
       y = "Density") +
  theme_minimal()

print(p3)
print(p4)

# Create the density plots combined before and after imputation 
# Combine data for "PC aa C30:2"
combined_PC30_2 <- rbind(
  data.frame(Value = filtered_data_PC30_2$`PC aa C30:2`, Status = "Before Imputation"),
  data.frame(Value = RHEA_imputed_data_columns$`PC aa C30:2`, Status = "After Imputation")
)

# Plot for "PC aa C30:2"
p_combined_PC30_2 <- ggplot(combined_PC30_2, aes(x = Value, color = Status, fill = Status)) +
  geom_density(alpha = 0.4) +
  labs(title = "Density Plot for PC aa C30:2 Before and After Imputation",
       x = "PC aa C30:2",
       y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Combine data for "PC ae C38:1"
combined_PC38_1 <- rbind(
  data.frame(Value = filtered_data_PC38_1$`PC ae C38:1`, Status = "Before Imputation"),
  data.frame(Value = RHEA_imputed_data_columns$`PC ae C38:1`, Status = "After Imputation")
)

# Plot for "PC ae C38:1"
p_combined_PC38_1 <- ggplot(combined_PC38_1, aes(x = Value, color = Status, fill = Status)) +
  geom_density(alpha = 0.4) +
  labs(title = "Density Plot for PC ae C38:1 Before and After Imputation",
       x = "PC ae C38:1",
       y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Print the plots
print(p_combined_PC30_2)
print(p_combined_PC38_1)


# Perform the Anderson-Darling k-sample test
ad_test_result <- ad.test(
  filtered_data_PC30_2$`PC aa C30:2`, # Before imputation
  RHEA_imputed_data_columns$`PC aa C30:2` # After imputation
)

# Print the test results
print(ad_test_result)


ad_test_result_38_1 <- ad.test(
  filtered_data_PC38_1$`PC ae C38:1`,
  RHEA_imputed_data_columns$`PC ae C38:1`
)

print(ad_test_result_38_1)

# Inverse normal transformation 
# Apply the inverse normal transformation to each metabolite column

DB_norm = sapply(RHEA_imputed_data[, -1], function(x) RankNorm(x))
View(DB_norm)

# Convert DB_norm to a data frame and assign column names if necessary
DB_norm_df <- as.data.frame(DB_norm)
colnames(DB_norm_df) <- colnames(RHEA_imputed_data)[-1]

# Add the 'helixid' column back to the transformed data
RHEA_imputed_data_normalized <- cbind(helixid = RHEA_imputed_data_columns$helixid, DB_norm_df)

# View the transformed data
View(RHEA_imputed_data_normalized)
# Save the normalized data set 
#write.xlsx(RHEA_imputed_data_normalized, file = "RHEA_imputed_normalized_serum_data.xlsx")

# Check the distribution after normalization 
# Compute the Pearson correlation matrix for all metabolites
correlation_matrix <- cor(RHEA_imputed_data_normalized[, -1], method = "spearman")

# Compute the correlation matrix (assuming your data is in 'RHEA_imputed_data')
correlation_matrix <- cor(RHEA_imputed_data_normalized[, -1], method = "spearman")

# Step 1: Compute Spearman correlation matrix
correlation_matrix <- cor(RHEA_imputed_data_normalized[, -1], method = "spearman", use = "pairwise.complete.obs")


# Quality Control: Check for Negative Values in Metabolite Features
# ----------------------------------------------------------
# After preprocessing and imputation, we verify that no
# negative values remain in the metabolite feature matrix.
# The ID column (helixid) is excluded from this check.

correlation_matrix <- cor(RHEA_imputed_data_normalized[, -1], method = "spearman", use = "all.obs")
View(correlation_matrix)

# Abbreviate to ~6 chars (tweak minlength if needed)
lab_row <- abbreviate(rownames(correlation_matrix), minlength = 6)
lab_col <- abbreviate(colnames(correlation_matrix), minlength = 6)

# Abbreviate to ~6 chars (tweak minlength if needed)
lab_row <- abbreviate(rownames(correlation_matrix), minlength = 6)
lab_col <- abbreviate(colnames(correlation_matrix), minlength = 6)

p <- pheatmap(
  correlation_matrix,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  labels_row = lab_row,
  labels_col = lab_col,
  main = "Spearman Correlation Heatmap of Serum Metabolites",
  fontsize = 6,
  fontsize_row = 4,
  fontsize_col = 4,
  angle_col = 45
)

grid.force()
grid.gedit(gPath("row_names"), gp = gpar(fontface = "bold"))
grid.gedit(gPath("col_names"), gp = gpar(fontface = "bold"))

################################################################################
################################################################################
##########Urine metabolomic dataset preprocessing##########----

URINE <- read_excel("data/datasets/Preprocessing/HELIX_urine_metabol_data_IC_v2_CHL_2017_01_26_onlyRhea.xlsx", 
                    sheet = "creatinine norm")


# Extract relevant parts from Urine Sample ID
URINE_ <- URINE %>%
  mutate(
    Child_ID_from_ID = str_extract(`Urine Sample ID`, "(?<=RHE-)[0-9]+"),
    suffix = str_extract(`Urine Sample ID`, "1[a-zx]")
  )

filtered_URINE <- URINE_ %>%
  filter(
    suffix == "1x" |                                     # keep all 1x
      (suffix == "1a" & Child_ID_from_ID %in% 
         Child_ID_from_ID[suffix %in% "1b"])               # keep 1a if 1b exists (morning and night void samples)
  ) %>%
  distinct(Child_ID_from_ID, suffix, .keep_all = TRUE)   # remove 1b duplicates if both 1a and 1b exist


filtered_URINE <- filtered_URINE %>%
  mutate(helixid = paste0("RHE", Child_ID)) %>%
  select(helixid, everything(), -`Urine Sample ID`, -Centre, -Child_ID, -RunOrder, -suffix, -Child_ID_from_ID)

#setwd("")

##########Metabolomic dataset creation##########----
##########Urine metabolomic dataset preprocessing##########----

###Loading datasets (Urine metabolomic dataset from HELIX)
# Participants diagnosed as ADS-positive did not participate
# in the internalizing and externalizing assessments.
# These individuals are excluded from the analytic dataset.
# We remove these helixids 

filtered_URINE <- filtered_URINE %>%
  filter(!helixid %in% c("RHExxxxxx", "RHExxxxxx", "RHExxxxxx"))

#write_xlsx(filtered_URINE, "data/datasets/RHEA_urine_data.xlsx")

#Add the concentation of creatinine to have it in the summary table of metabolites
concentration <- read_excel("data/datasets/Preprocessing/HELIX_urine_metabol_data_IC_v2_CHL_2017_01_26_onlyRhea.xlsx", 
                            sheet = "concentration")
# Extract suffix and Child_ID
concentration <- concentration %>%
  mutate(
    suffix = str_extract(`Urine Sample ID`, "1[a-zx]"),
    Child_ID_from_ID = str_extract(`Urine Sample ID`, "(?<=RHE-)[0-9]+")
  )


# Apply filtering: keep all 1x, and 1a if there's also a 1b
filtered_concentration <- concentration %>%
  filter(
    suffix == "1x" |
      (suffix == "1a" & Child_ID_from_ID %in% Child_ID_from_ID[suffix == "1b"])
  ) %>%
  distinct(Child_ID_from_ID, suffix, .keep_all = TRUE)

# Create helixid and remove unwanted helixid rows
creatinine_concentration <- filtered_concentration %>%
  mutate(helixid = paste0("RHE", Child_ID)) %>%
  filter(!helixid %in% c("RHExxxxxx", "RHExxxxxx", "RHExxxxxx")) %>%
  select(helixid, Creatinine)

#write_xlsx(creatinine_concentration, "data/datasets/creatinine_concentration_urine.xlsx")

#MERGE THE CREATININE CONCETRATION WITH THE REST DATASET THAT IT IS NORMALIZED FOR THE CREATININE 
# Merge creatinine concentration into the filtered urine data

filtered_URINE_with_creatinine <- filtered_URINE %>%
  left_join(creatinine_concentration, by = "helixid")

#write_xlsx(filtered_URINE_with_creatinine, "data/datasets/RHEA_data_urine_withcreatinine.xlsx")

# Create summary with only median, Q1, and Q3
urine_summary_median_IQR <- filtered_URINE_with_creatinine %>%
  select(where(is.numeric)) %>%  # Keep only numeric variables
  pivot_longer(cols = everything(), names_to = "metabolite", values_to = "value") %>%
  group_by(metabolite) %>%
  summarise(
    median = median(value, na.rm = TRUE),
    Q1 = quantile(value, 0.25, na.rm = TRUE),
    Q3 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    summary = sprintf("%.3f (%.3f–%.3f)", median, Q1, Q3)
  ) %>%
  select(metabolite, summary)

# View the final table
View(urine_summary_median_IQR)
#write.csv(urine_summary_median_IQR, "urine_summary_median_IQR.csv", row.names = FALSE)

##########################################################################################
################### Read the prepared metabolomic urine dataset ##########################

RHEA_urine_data <- read_excel("data/datasets/RHEA_urine_data.xlsx")

##########################################################################################
##########################################################################################
# Check the data  
#Search for missing data 
missing_per_column <- colSums(is.na(RHEA_urine_data))
print(missing_per_column)

# Total number of missing values in the entire dataset
total_missing <- sum(is.na(RHEA_urine_data))
cat("Total missing values in the dataset:", total_missing, "\n")

# Optional: show only variables with any missing values
missing_per_column[missing_per_column > 0]

# Search for zero values 

RHEA_urine_data %>%
  select(where(is.numeric)) %>%
  summarise(total_zeros = sum(across(everything(), ~sum(. == 0, na.rm = TRUE))))

RHEA_urine_data%>%
  filter(if_any(where(is.numeric), ~ . == 0))

zero_summary <-RHEA_urine_data %>%
  select(where(is.numeric)) %>%
  pivot_longer(cols = everything(), names_to = "metabolite", values_to = "value") %>%
  group_by(metabolite) %>%
  summarise(
    zero_count = sum(value == 0, na.rm = TRUE),
    non_missing = sum(!is.na(value)),
    zero_percent = round(100 * zero_count / non_missing, 2)
  ) %>%
  filter(zero_count > 0) %>%
  arrange(desc(zero_percent))

helixid_zero_details <- RHEA_urine_data %>%
  select(helixid, where(is.numeric)) %>%
  pivot_longer(-helixid, names_to = "metabolite", values_to = "value") %>%
  filter(value == 0) %>%
  arrange(helixid)

View(helixid_zero_details)

# Get long format with helixid, metabolite, and value
zero_details <- RHEA_urine_data %>%
  select(helixid, where(is.numeric)) %>%
  pivot_longer(-helixid, names_to = "metabolite", values_to = "value") %>%
  filter(value == 0)

# Summarise number of zeros and metabolites per helixid
helixid_zero_summary <- zero_details %>%
  group_by(helixid) %>%
  summarise(
    zero_count = n(),
    metabolites_with_zero = paste(unique(metabolite), collapse = ", ")
  ) %>%
  arrange(desc(zero_count))

# View the result
View(helixid_zero_summary)

#Check for negative values  
RHEA_urine_data%>%
  select(where(is.numeric)) %>%
  summarise(across(everything(), ~sum(. < 0, na.rm = TRUE))) %>%
  pivot_longer(cols = everything(), names_to = "metabolite", values_to = "negative_count") %>%
  filter(negative_count > 0) %>%
  arrange(desc(negative_count))

# Create density plots of metabolites with zero values 
zero_metabolites <- zero_summary$metabolite

# Loop to create density plots with red fill + legend
for (met in zero_metabolites) {
  plot_data <- RHEA_urine_data %>%
    select(helixid, all_of(met)) %>%
    filter(.data[[met]] != 0) %>%
    mutate(group = "Before imputation (zeros excluded)")  # for legend
  
  if (nrow(plot_data) > 1) {
    p <- ggplot(plot_data, aes(x = .data[[met]], fill = group)) +
      geom_density(alpha = 0.5) +
      scale_fill_manual(values = c("Before imputation (zeros excluded)" = "red")) +
      labs(
        title = paste("Density Plot:", met),
        fill = "Legend",
        x = met,
        y = "Density"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "right"
      )
    
    print(p)
  }
}


# Replace zero values with missing values 
# Replace all 0 values with NA (excluding helixid)
urine_NA_replaced <- RHEA_urine_data
urine_NA_replaced[urine_NA_replaced == 0] <- NA

# Calculate total missing values (excluding helixid)
total_NA_urine <- sum(is.na(urine_NA_replaced[, -1]))  # assumes helixid is the first column
cat("Total number of missing values in the urine dataset (excluding helixid):", total_NA_urine, "\n")

# Calculate percentage of missing values per metabolite
missing_percentage_urine <- colSums(is.na(urine_NA_replaced[, -1])) / nrow(urine_NA_replaced) * 100

# Create a data frame
missing_df_urine <- data.frame(
  Metabolite = names(missing_percentage_urine),
  Missingness = missing_percentage_urine
)

# Filter to keep only metabolites with missingness > 0
missing_df_urine_filtered <- missing_df_urine %>%
  filter(Missingness > 0)

# View the result
View(missing_df_urine_filtered)


# Plot for the percentage of missigness of the metabolites 
ggplot(missing_df_urine_filtered, aes(x = reorder(Metabolite, -Missingness), y = Missingness)) +
  geom_bar(stat = "identity", fill = "blue") +
  geom_hline(yintercept = 20, color = "red", linetype = "dashed") +  # 20% cutoff line
  annotate("text", x = Inf, y = 22, label = "20% threshold", hjust = 1.1, color = "red") +
  theme_minimal() +
  labs(title = "Percentage of Missing Values per Metabolite (Urine Dataset)",
       x = "Metabolite",
       y = "Percentage of Missingness") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


#N-methylpicolinic acid has more than 20% of missigness so it should be removed 
urine_cleaned <- urine_NA_replaced %>%
  select(-`N-methylpicolinic acid`, -SamplingType)

#Random forest imputation of the missing metabolites  
# Separate ID and numeric data
data_to_impute_URINE <- urine_cleaned %>%
  select(-helixid)

# Ensure all columns are numeric (just in case)
data_to_impute_URINE <- data_to_impute_URINE %>%
  mutate_if(is.character, as.numeric) %>%
  mutate_if(is.factor, as.numeric) %>%
  as.data.frame()

# Set seed for reproducibility
set.seed(1)

# Apply Random Forest Imputation
imputed_result_urine <- missForest(data_to_impute_URINE)

# Extract imputed data
urine_imputed_data <- imputed_result_urine$ximp

# Add helixid back (and SamplingType if needed later)
urine_imputed_data <- cbind(helixid = urine_cleaned$helixid, urine_imputed_data)

# View the result
View(urine_imputed_data)

any_missing <- any(is.na(urine_imputed_data))
cat("Are there any missing values after imputation? ", any_missing, "\n")

# Optional: total count of NAs (should be 0 ideally)
total_missing <- sum(is.na(urine_imputed_data))
cat("Total number of missing values: ", total_missing, "\n")

# Count exact zeros
total_zeros <- sum(urine_imputed_data == 0, na.rm = TRUE)
cat("Total number of zeros: ", total_zeros, "\n")

#write.xlsx(urine_imputed_data, file = "RHEA_imputed_urine_data.xlsx")

RHEA_imputed_urine_data <- read_excel("data/datasets/Preprocessing/RHEA_imputed_urine_data.xlsx")

#I want to check what values the imputed values took 
# Extract imputed values from the imputed dataset

# Remove ID and SamplingType from both for comparison
original_values <- urine_NA_replaced %>%
  select(-helixid, -SamplingType)

imputed_values <- urine_imputed_data %>%
  select(-helixid)

# Create a logical mask of where imputation occurred
imputed_mask <- is.na(original_values)

# Extract imputed values only
imputed_only <- imputed_values[imputed_mask]

# Summary of imputed values
summary(imputed_values)

# Create the density plots before and after imputation 

# Prepare numeric-only parts (excluding helixid and SamplingType)
original_data <- urine_NA_replaced %>% select(-helixid)
imputed_data <- urine_imputed_data %>% select(-helixid)

# Identify where imputation happened
imputed_mask <- is.na(original_data)

# Get list of metabolites with at least one imputed value
imputed_metabolites <- colnames(original_data)[colSums(imputed_mask) > 0]

# Loop to create plots with legend
for (met in imputed_metabolites) {
  plot_data <- data.frame(value = imputed_data[[met]]) %>%
    mutate(group = "After imputation")
  
  if (sum(!is.na(plot_data$value)) > 1) {
    p <- ggplot(plot_data, aes(x = value, fill = group)) +
      geom_density(alpha = 0.5) +
      scale_fill_manual(values = c("After imputation" = "blue")) +
      labs(
        title = paste("Density Plot:", met),
        fill = "Legend",
        x = met,
        y = "Density"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right"
      )
    
    print(p)
  }
}

##################################################################################
# Compine the 2 plot before and after imputation together
# Extract the numeric part of both datasets (no helixid or SamplingType)
original_data <- urine_NA_replaced %>% select(-helixid, -SamplingType)
imputed_data <- urine_imputed_data %>% select(-helixid)

# Identify which metabolites had missing values (i.e., were imputed)
imputed_mask <- is.na(original_data)
imputed_metabolites <- colnames(original_data)[colSums(imputed_mask) > 0]

# Loop through each metabolite and plot before vs after imputation
for (met in imputed_metabolites) {
  before_values <- original_data[[met]][!is.na(original_data[[met]])]
  after_values <- imputed_data[[met]]
  
  # Skip if insufficient values
  if (length(before_values) < 2 || length(after_values) < 2) next
  
  # Combine into one data frame
  combined <- rbind(
    data.frame(Value = before_values, Status = "Before Imputation"),
    data.frame(Value = after_values, Status = "After Imputation")
  )
  
  # Plot
  p <- ggplot(combined, aes(x = Value, color = Status, fill = Status)) +
    geom_density(alpha = 0.4) +
    labs(
      title = paste("Density Plot Before vs After Imputation:", met),
      x = met,
      y = "Density"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(p)
}

# Anderson darling test before and after imputation 
# Prepare numeric-only data
original_data <- urine_NA_replaced %>% select(-helixid, -SamplingType)
imputed_data <- urine_imputed_data %>% select(-helixid)

# Identify imputed metabolites (those that had NA originally)
imputed_mask <- is.na(original_data)
imputed_metabolites <- colnames(original_data)[colSums(imputed_mask) > 0]



# Define clean datasets
original_data <- urine_NA_replaced %>% select(-helixid, -SamplingType)
imputed_data <- urine_imputed_data %>% select(-helixid)

# Identify metabolites that were imputed
imputed_mask <- is.na(original_data)
imputed_metabolites <- colnames(original_data)[colSums(imputed_mask) > 0]

manual_imputed_metabolites <- c(
  "3-hydroxybutyrate/3-aminoisobutyrate",
  "Proline betaine",
  "Succinate",
  "Scyllo-inositol",
  "N-methyl-2-pyridone-5-carboxamide",
  "Sucrose",
  "3-aminoisobutyrate",
  "Carnitine",
  "N-acetyl neuraminic acid",
  "Taurine",
  "Trimethylamine oxide",
  "3-Indoxylsulfate",
  "5-oxoproline",
  "Acetate",
  "Dimethylamine",
  "Glutamine",
  "N1-methyl nicotinamide"
)

# Loop through metabolites and print full test results with custom labels
for (met in manual_imputed_metabolites) {
  before_values <- urine_NA_replaced[[met]][!is.na(urine_NA_replaced[[met]])]
  after_values <- urine_imputed_data[[met]]
  
  if (length(before_values) < 5 || length(after_values) < 5) next
  
  cat("==================================================\n")
  cat("Anderson-Darling Test for:", met, "\n")
  cat("==================================================\n")
  
  result <- tryCatch({
    ad.test(before_values, after_values)
  }, error = function(e) return(NULL))
  
  if (!is.null(result)) {
    # Modify version labels
    rownames(result$ad) <- c("Before Imputation", "After Imputation")
    print(result)
  } else {
    cat("Test could not be performed for", met, "\n")
  }

# Inverse normal transformation normalization 
  # Load the required package
  # Step 1: Apply the inverse normal transformation (excluding helixid)
  DB_norm_urine <- sapply(urine_imputed_data[, -1], function(x) RankNorm(x))
  
  # Step 2: Convert to data frame
  DB_norm_urine_df <- as.data.frame(DB_norm_urine)
  
  # Step 3: Assign original column names back (excluding helixid)
  colnames(DB_norm_urine_df) <- colnames(urine_imputed_data)[-1]
  
  # Step 4: Add helixid back
  urine_imputed_data_normalized <- cbind(helixid = urine_imputed_data$helixid, DB_norm_urine_df)
  
  cat("\n\n")
}

total_zeros <- sum(urine_imputed_data_normalized == 0, na.rm = TRUE)
cat("Total number of zeros: ", total_zeros, "\n")


zeros_per_col <- colSums(X == 0, na.rm = TRUE)
metabolites_with_zeros <- zeros_per_col[zeros_per_col > 0]
print(metabolites_with_zeros)
names(metabolites_with_zeros)

#write.xlsx(urine_imputed_data_normalized, file = "RHEA_imputed_normalized_urine_data.xlsx")

# Correlation of the metabolites 
# Step 1: Compute the Spearman correlation matrix for normalized urine data
correlation_matrix_urine <- cor(
  urine_imputed_data_normalized[, -1],  # exclude helixid
  method = "spearman",
  use = "pairwise.complete.obs"
)

# View correlation matrix (optional)
View(correlation_matrix_urine)

# Plot using pheatmap
pheatmap(
  correlation_matrix_urine,
  color = colorRampPalette(c("blue", "white", "red"))(100),  # Blue to red gradient
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Spearman Correlation Heatmap of Urine Metabolites",
  fontsize = 5
)

################################################################################
################################################################################
############Preprocessing of outcomes Internalizing, Externalizing, ADHD########

getwd()
setwd('set_path_to_data_folder')

#Read data from internalizing, externalizing, adhd from RHEA study 
inter_externdata <- read_excel("data/datasets/intern_extern_data/primary dataset/inter_externdata.xlsx")

#Filter the outcome data for participants 6 years old that is the baseline age
inter_externdata_age6 <- inter_externdata %>%
  mutate(helixid = paste0("RHE", childid)) %>%
  select(helixid, sex, age6, cbcl_pro_ext6, cbcl_pro_int6, diagnosis_6to15, adhd_pro6, preterm) %>%
  rename(
  cbcl_ext6 = cbcl_pro_ext6)
#save the data set for 6 years 
#write.xlsx(inter_externdata_age6, file = "inter_externdata_age6.xlsx")

#Load the serum metabolomic data after the preprocessing (imputed and normalized)
serum_final_data <- read_excel("data/datasets/Preprocessing/RHEA_imputed_normalized_serum_data.xlsx")

#Load the dataset that we created filtered at the age 6 
inter_externdata_age6 <- read_excel("data/datasets/intern_extern_data/primary dataset/inter_externdata_age6.xlsx")

# General view of the dataset 
# Summary statistics for the outcome 
summary(inter_externdata_age6)
numeric_cols <- inter_externdata_age6 %>%
  select(where(is.numeric))
summary(numeric_cols)
descriptive_stats <- psych::describe(inter_externdata_age6 %>% select_if(is.numeric))
print(descriptive_stats)

# Count the number of the males and felames in the dataset 
sex_counts <- inter_externdata_age6 %>%
  count(sex)  # Count occurrences of each gender
# Print the count table
print(sex_counts)

# Select only numeric columns
numeric_data <- inter_externdata_age6 %>% select_if(is.numeric)

# Identify in which ids you have missing values in every test
# Identify rows where cbcl_adhd6 is NA
missing_helixids <- inter_externdata_age6 %>%
  filter(is.na(adhd_pro6)) %>%
  select(helixid)
# Print the helixids with missing cbcl_adhd6 values
print(missing_helixids)


#Identify rows where cbcl_ext6 is NA
missing_helixids_cbcl_ext6 <- inter_externdata_age6 %>%
  filter(is.na(cbcl_ext6)) %>%
  select(helixid)
# Print the helixids with missing cbcl_adhd6 values
print(missing_helixids_cbcl_ext6)


#Identify rows where cbcl_int6 is NA
missing_helixids_cbcl_int6 <- inter_externdata_age6 %>%
  filter(is.na(cbcl_pro_int6)) %>%
  select(helixid)
# Print the helixids with missing cbcl_adhd6 values
print(missing_helixids_cbcl_int6)

# Calculate the total number of missing values in numeric columns
total_missing_values <- sum(is.na(numeric_data))
cat("Total number of missing values in numeric columns:", total_missing_values, "\n")

# Calculate the total number of zero values in numeric columns
total_zero_values <- sum(numeric_data == 0, na.rm = TRUE)
cat("Total number of zero values in numeric columns:", total_zero_values, "\n")

# Identify numeric columns with zero values
columns_with_zero_values <- names(numeric_data)[
  colSums(numeric_data == 0, na.rm = TRUE) > 0
]
cat("Numeric columns with zero values:", paste(columns_with_zero_values, collapse = ", "), "\n")

# Identify numeric columns with missing values
columns_with_missing_values <- names(numeric_data)[
  colSums(is.na(numeric_data)) > 0
]
cat("Numeric columns with missing values:", paste(columns_with_missing_values, collapse = ", "), "\n")

################################################################################
####### Create the baseline dataset for the ourcomes Internalizing, Externalizing, ADHD#####

# Create data set for each outcome but without the missing values (Internalizing)
# Create a dataset with helixid and cbcl_pro_int6 without missing values
cbcl_pro_int6_clean <- inter_externdata_age6 %>%
  select(helixid, cbcl_pro_int6) %>%  # Select only helixid and cbcl_pro_int6
  filter(!is.na(cbcl_pro_int6))  # Remove rows with missing values in cbcl_pro_int6
summary(cbcl_pro_int6_clean)

#195 helixids
# Save the dataset to a file
#write.xlsx(cbcl_pro_int6_clean, file = "cbcl_pro_int6_clean.xlsx")

# Create a dataset with helixid and cbcl_pro_ext6 without missing values (Externalizing)
cbcl_pro_ext6_clean <- inter_externdata_age6 %>%
  select(helixid, cbcl_ext6) %>%  # Select only helixid and cbcl_pro_ext6
  filter(!is.na(cbcl_ext6))       # Remove rows with missing values in cbcl_pro_ext6
summary(cbcl_pro_ext6_clean)
#195 helixids
#write.xlsx(cbcl_pro_ext6_clean, file = "cbcl_pro_ext6_clean.xlsx")


# Create a dataset with helixid and ADHD6 without missing values (ADHD)
# 193 helixids
cbcl_pro_adhd6_clean <- inter_externdata_age6 %>%
  select(helixid, adhd_pro6) %>%  # Select only helixid and cbcl_pro_ext6
  filter(!is.na(adhd_pro6))       # Remove rows with missing values in cbcl_pro_ext6
summary(cbcl_pro_adhd6_clean)
# 193 helixids
#write.xlsx(cbcl_pro_adhd6_clean, file = "cbcl_pro_adhd6_clean.xlsx")



# Create density plots for the ourcomes after excluding the missing values 
# Select only numeric columns
numeric_cols <- inter_externdata_age6 %>%
  select(where(is.numeric)) %>%
  colnames()  # Ensure we are iterating over column names

# Loop through each numeric column to create a separate density plot
for (col in numeric_cols) {
  # Filter out rows with missing values in the current column
  filtered_data <- inter_externdata_age6 %>%
    filter(!is.na(.data[[col]]))  # Use .data[[col]] safely
  
  # Create the density plot
  p <- ggplot(filtered_data, aes(x = .data[[col]])) +  # Use .data[[col]] inside aes()
    geom_density(fill = "skyblue", color = "black", alpha = 0.5) +
    labs(title = paste("Density Plot of", col, "(excluding missing values)"), x = col, y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(p)  # Print each plot individually
}


# Apply Anderson-Darling test to each numeric column (excluding missing values)
ad_test_results <- sapply(numeric_cols, function(col) {
  data <- inter_externdata_age6[[col]]
  data <- data[!is.na(data)]  # Exclude missing values
  ad.test(data)$p.value       # Return only the p-value
})

# Create a summary table for better readability
ad_test_summary <- data.frame(
  Column = numeric_cols,
  P_Value = ad_test_results
)

# Print the summary
print(ad_test_summary)

# I want to see if they have common ids between outcome dataset and serum metabolomic dataset 
common_helixids <- intersect(inter_externdata_age6$helixid, serum_final_data$helixid)

# Display common helixid
summary(common_helixids)
# Exclude 2 ids because they have ASD
# Identify the helixid values that are not common (unique to each dataset)

unique_to_serum_final_data <- setdiff(serum_final_data$helixid, inter_externdata_age6$helixid)
# Display the unique helix IDs for each dataset
cat("Helix IDs unique to inter_externdata:\n")
print(unique_to_serum_final_data)

########################################################################################################

# Apply square root transformation and round to 4 decimal places for cbcl_pro_int6_clean
sqrt_transformed_cbcl_pro_int6 <- cbcl_pro_int6_clean %>%
  mutate(sqrt_cbcl_pro_int6 = round(sqrt(cbcl_pro_int6), 4))

#write_xlsx(cbcl_pro_int6_clean, "sqrt_transformed_cbcl_pro_int6.xlsx")

# Apply square root transformation and round to 4 decimal places for cbcl_pro_ext6_clean
sqrt_transformed_cbcl_pro_ext6 <- cbcl_pro_ext6_clean %>%
  mutate(sqrt_cbcl_ext6 = round(sqrt(cbcl_ext6), 4))

#write_xlsx(sqrt_transformed_cbcl_pro_ext6, "sqrt_transformed_cbcl_pro_ext6.xlsx")

# Apply square root transformation and round to 4 decimal places for cbcl_pro_adhd6_clean
sqrt_transformed_cbcl_pro_adhd6 <- cbcl_pro_adhd6_clean %>%
  mutate(sqrt_adhd_pro6 = round(sqrt(adhd_pro6), 4))
#write_xlsx(sqrt_transformed_cbcl_pro_adhd6, "sqrt_transformed_cbcl_pro_adhd6.xlsx")

# Density plot for square root transformed cbcl_pro_int6
ggplot(sqrt_transformed_cbcl_pro_int6, aes(x = sqrt_cbcl_pro_int6)) +
  geom_density(fill = "blue", color = "black", alpha = 0.5) +
  labs(title = "Density Plot of sqrt(cbcl_pro_int6)", x = "sqrt(cbcl_pro_int6)", y = "Density") +
  theme_minimal()

# Density plot for square root transformed cbcl_ext6
ggplot(sqrt_transformed_cbcl_pro_ext6, aes(x = sqrt_cbcl_ext6)) +
  geom_density(fill = "blue", color = "black", alpha = 0.5) +
  labs(title = "Density Plot of sqrt(cbcl_ext6)", x = "sqrt(cbcl_ext6)", y = "Density") +
  theme_minimal()

# Density plot for square root transformed adhd_pro6
ggplot(sqrt_transformed_cbcl_pro_adhd6, aes(x = sqrt_adhd_pro6)) +
  geom_density(fill = "blue", color = "black", alpha = 0.5) +
  labs(title = "Density Plot of sqrt(adhd_pro6)", x = "sqrt(adhd_pro6)", y = "Density") +
  theme_minimal()

#Combine density plots for cbcl_pro_int6
ggplot() +
  geom_density(data = cbcl_pro_int6_clean, aes(x = cbcl_pro_int6, color = "Before"), fill = "skyblue", alpha = 0.4) +
  geom_density(data = sqrt_transformed_cbcl_pro_int6, aes(x = sqrt_cbcl_pro_int6, color = "After"), fill = "blue", alpha = 0.4) +
  labs(title = "Density Plot: Before and After Square Root Transformation (cbcl_pro_int6)",
       x = "cbcl_pro_int6 / sqrt(cbcl_pro_int6)",
       y = "Density") +
  scale_color_manual(values = c("Before" = "skyblue", "After" = "blue")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank())

# Combine density plots for cbcl_ext6
ggplot() +
  geom_density(data = cbcl_pro_ext6_clean, aes(x = cbcl_ext6, color = "Before"), fill = "skyblue", alpha = 0.4) +
  geom_density(data = sqrt_transformed_cbcl_pro_ext6, aes(x = sqrt_cbcl_ext6, color = "After"), fill = "blue", alpha = 0.4) +
  labs(title = "Density Plot: Before and After Square Root Transformation (cbcl_ext6)",
       x = "cbcl_ext6 / sqrt(cbcl_ext6)",
       y = "Density") +
  scale_color_manual(values = c("Before" = "skyblue", "After" = "blue")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank())

# Combine density plots for adhd_pro6
ggplot() +
  geom_density(data = cbcl_pro_adhd6_clean, aes(x = adhd_pro6, color = "Before"), fill = "skyblue", alpha = 0.4) +
  geom_density(data = sqrt_transformed_cbcl_pro_adhd6, aes(x = sqrt_adhd_pro6, color = "After"), fill = "blue", alpha = 0.4) +
  labs(title = "Density Plot: Before and After Square Root Transformation (adhd_pro6)",
       x = "adhd_pro6 / sqrt(adhd_pro6)",
       y = "Density") +
  scale_color_manual(values = c("Before" = "skyblue", "After" = "blue")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank())

##################################################################################
# Perform Anderson-Darling test for normality on original and transformed data
ad_test_cbcl_pro_int6 <- ad.test(cbcl_pro_int6_clean$cbcl_pro_int6)
ad_test_cbcl_ext6 <- ad.test(cbcl_pro_ext6_clean$cbcl_ext6)
ad_test_adhd_pro6 <- ad.test(cbcl_pro_adhd6_clean$adhd_pro6)

cat("Anderson-Darling Test Results for Original Data:\n")
print(ad_test_cbcl_pro_int6)
print(ad_test_cbcl_ext6)
print(ad_test_adhd_pro6)

# Apply Anderson-Darling test to each transformed column
ad_test_transformed_cbcl_pro_int6 <- ad.test(sqrt_transformed_cbcl_pro_int6$sqrt_cbcl_pro_int6)
ad_test_transformed_cbcl_ext6 <- ad.test(sqrt_transformed_cbcl_pro_ext6$sqrt_cbcl_ext6)
ad_test_transformed_cbcl_pro_adhd6 <- ad.test(sqrt_transformed_cbcl_pro_adhd6$sqrt_adhd_pro6)

# Create a summary table for better readability
ad_test_summary_transformed <- data.frame(
  Column = c("sqrt_cbcl_pro_int6", "sqrt_cbcl_ext6", "sqrt_adhd_pro6"),
  P_Value = c(
    ad_test_transformed_cbcl_pro_int6$p.value,
    ad_test_transformed_cbcl_ext6$p.value,
    ad_test_transformed_cbcl_pro_adhd6$p.value
  )
)

# Print the summary
cat("Anderson-Darling Test Results for Transformed Data:\n")
print(ad_test_summary_transformed)

##############################################################################################
#################################################################################
##########Merge  metabolomic datasets with the outcomes at baseline at age 6  ########
############### Serum #################################################################

# Merge square root transformed cbcl_pro_int6_clean with serum_final_data
merged_sqrt_cbcl_pro_int6 <- sqrt_transformed_cbcl_pro_int6 %>%
  inner_join(serum_final_data, by = "helixid")

# Save the dataset to a file (optional)
#write.xlsx(merged_sqrt_cbcl_pro_int6, file = "merged_sqrt_cbcl_pro_int6.xlsx")

# Merge square root transformed cbcl_pro_ext6_clean with serum_final_data
merged_sqrt_cbcl_pro_ext6 <- sqrt_transformed_cbcl_pro_ext6 %>%
  inner_join(serum_final_data, by = "helixid")

# Save the dataset to a file (optional)
#write.xlsx(merged_sqrt_cbcl_pro_ext6, file = "merged_sqrt_cbcl_pro_ext6.xlsx")

# Merge square root transformed cbcl_pro_adhd6_clean with serum_final_data
merged_sqrt_cbcl_pro_adhd6 <- sqrt_transformed_cbcl_pro_adhd6 %>%
  inner_join(serum_final_data, by = "helixid")

# Save the dataset to a file (optional)
#write.xlsx(merged_sqrt_cbcl_pro_adhd6, file = "merged_sqrt_cbcl_pro_adhd6.xlsx")

# Print summaries of the merged datasets
cat("Merged Dataset for sqrt(cbcl_pro_int6):\n")
print(summary(merged_sqrt_cbcl_pro_int6))

cat("Merged Dataset for sqrt(cbcl_ext6):\n")
print(summary(merged_sqrt_cbcl_pro_ext6))

cat("Merged Dataset for sqrt(adhd_pro6):\n")
print(summary(merged_sqrt_cbcl_pro_adhd6))


################################# Urine ########################################
RHEA_imputed_normalized_urine_data <- read_excel("RHEA_imputed_normalized_urine_data.xlsx")
sqrt_cbcl_pro_ext6$helixid <- as.character(sqrt_cbcl_pro_ext6$helixid)
RHEA_imputed_normalized_urine_data$helixid <- as.character(RHEA_imputed_normalized_urine_data$helixid)
# Merge by helixid
merged_sqrt_urine_ext6 <- left_join(sqrt_cbcl_pro_ext6, RHEA_imputed_normalized_urine_data, by = "helixid")
#write_xlsx(merged_sqrt_urine_ext6, path = "merged_sqrt_urine_ext6.xlsx")

sqrt_cbcl_pro_int6$helixid <- as.character(sqrt_cbcl_pro_int6$helixid)
RHEA_imputed_normalized_urine_data$helixid <- as.character(RHEA_imputed_normalized_urine_data$helixid)
# Merge by helixid
merged_sqrt_urine_int6 <- left_join(sqrt_cbcl_pro_int6, RHEA_imputed_normalized_urine_data, by = "helixid")
#write_xlsx(merged_sqrt_urine_int6, path = "merged_sqrt_urine_int6.xlsx")

sqrt_pro_adhd6$helixid <- as.character(sqrt_pro_adhd6$helixid)
RHEA_imputed_normalized_urine_data$helixid <- as.character(RHEA_imputed_normalized_urine_data$helixid)
# Merge by helixid
merged_sqrt_urine_adhd6 <- left_join(sqrt_pro_adhd6, RHEA_imputed_normalized_urine_data, by = "helixid")
#write_xlsx(merged_sqrt_urine_adhd6, path = "merged_sqrt_urine_adhd6.xlsx")


#################### Lasso/Ridge Analysis- Linear regression with elastic net ######################
####################################################################################################

####################### Serum Analysis #############################################################
# Clean R Environment
############################################################
# Remove all objects from the current workspace to ensure
# a clean and reproducible execution of the script.
rm(list = ls())
gc()
getwd()
setwd('set_path_to_data_folder')

############################################################
# Elastic Net Regression (Alpha & Lambda Selection)
############################################################
# Objective:
# To identify optimal alpha (mixing parameter) and lambda
# (regularization strength) for predicting symptom scores
# using metabolite features.
#
# Models are fitted separately for each symptom domain:
#   - Internalizing
#   - Externalizing
#   - ADHD
# ------------------------------------------------------------------
# NOTE:
# The linear regression model is repeated separately for each symptom
# and for each biofluid.
#
# Before running a new model, the R environment is cleared to ensure:
# 1) No objects from previous analyses remain in memory
# 2) No parameter carryover occurs
# 3) Each model is run in a clean, independent workspace
#
# This guarantees full analytical independence between runs.
# ------------------------------------------------------------------
############################################################
# Select symptom score for analysis.
# Options:
#   "sqrt_cbcl_pro_int6"  -> Internalizing
#   "sqrt_cbcl_ext6"      -> Externalizing
#   "sqrt_adhd_pro6"      -> ADHD


# Create Predictor Matrix (Metabolites Only)
# Internalizing symptoms for serum metabolites 
# Load the merged datasets with the square root transformed symptpoms and the metabolomics (imputed/normalized)
merged_sqrt_cbcl_pro_int6 <- read_excel("data/datasets/intern_extern_data/merget outcome-metab-normal/merged_sqrt_cbcl_pro_int6.xlsx")
X <- as.matrix(merged_sqrt_cbcl_pro_int6 %>%
                 select(where(is.numeric)) %>%
                 select(-sqrt_cbcl_pro_int6, -cbcl_pro_int6))

symptom <- "sqrt_cbcl_pro_int6"
# Outcome vector
y <- merged_sqrt_cbcl_pro_int6[[symptom]]

#Initialize Storage Objects
############################################################
# These lists will store model outputs for each symptom:
#   - fit:      cross-validated glmnet models
#   - coef:     selected coefficients at optimal lambda
#   - alpha:    optimal alpha value (Elastic Net mixing parameter)
#   - lambda:   optimal lambda value (regularization strength)
#   - rmse:     cross-validated prediction error (RMSE)
############################################################

fit    <- list()
coef   <- list()
alpha  <- list()
lambda <- list()
rmse   <- list()

# 10-iteration loop for alpha and lambda tuning
set.seed(1)
for (i in 1:10) {
  cat("Iteration:", i, "Symptom:", symptom, "\n")
  
  # Split data into training (90%) and testing (10%) sets
  train_idx <- sample(1:nrow(X), size = 0.9 * nrow(X), replace = FALSE)
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  X_test <- X[-train_idx, ]
  y_test <- y[-train_idx]
  
  # Tune alpha and lambda using cross-validation
  cv_results <- lapply(seq(0.1, 1, by = 0.05), function(alpha_val) {
    cv_model <- cv.glmnet(X_train, y_train, alpha = alpha_val, nfolds = 10, family = "gaussian", type.measure = "mse")
    list(alpha = alpha_val, lambda = cv_model$lambda.min, cv_model = cv_model)
  })
  
  # Retrieve best alpha and lambda based on minimum cross-validation error
  best_result <- cv_results[[which.min(sapply(cv_results, function(res) min(res$cv_model$cvm)))]]
  alpha[[i]] <- best_result$alpha
  lambda[[i]] <- best_result$lambda
  
  # Fit the Elastic Net model using the best alpha and lambda
  fit[[i]] <- glmnet(X_train, y_train, alpha = alpha[[i]], lambda = lambda[[i]], family = "gaussian")
  
  # Extract non-zero coefficients
  coef[[i]] <- coef(fit[[i]])
  
  # Predict on the test set and calculate RMSE
  y_pred <- predict(fit[[i]], newx = X_test)
  rmse[[i]] <- sqrt(mean((y_test - y_pred)^2))
}

# Consolidate results into a data frame
results_inter <- data.frame(
  RMSE = unlist(rmse),
  Alpha = unlist(alpha),
  Lambda = unlist(lambda)
)

# Save results to a CSV file
#write.xlsx(results_inter, "serum_int_alpha_lambda.xlsx", rowNames = FALSE)

cat("Results saved for symptom:", symptom, "\n")
# Select best alpha based on lowest CV error
best_combination_intern <- results_inter[which.min(results_inter$RMSE), ]
print(best_combination_intern)

# Model Parameters (From Tuning)
# ------------------------------------------------------------
# Apply the best alpha and lambda values 
alpha_value  <- 0.4
lambda_value <- 0.1691768

# Prepare data 
X <- as.matrix(merged_sqrt_cbcl_pro_int6 %>%
                 select(where(is.numeric)) %>%
                 select(-sqrt_cbcl_pro_int6, -cbcl_pro_int6))
symptom <- "sqrt_cbcl_pro_int6"
y <- merged_sqrt_cbcl_pro_int6[[symptom]]

# Debugging: Check data
print(summary(X))
print(summary(y))
if (any(is.na(X)) || any(is.na(y))) {
  stop("Missing values detected in X or y.")
}

# Initialize a list to store coefficients from iterations
my_coef <- list()

# Set seed for reproducibility
set.seed(124876)

# Loop for coefficient selection
for (i in 1:10) {
  print(paste("Iteration:", i))
  
  # Split data into training (90%) and testing (10%) sets
  train_idx <- sample(1:nrow(X), size = 0.9 * nrow(X), replace = FALSE)
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  
  # Fit Elastic Net model using pre-selected alpha and lambda
  fit <- glmnet(X_train, y_train, family = "gaussian", alpha = alpha_value, lambda = lambda_value, intercept = FALSE)
  
  # Debugging: Check coefficients
  print(coef(fit))
  
  # Extract non-zero coefficients
  my_coef[[i]] <- as.matrix(coef(fit))
  my_coef[[i]] <- my_coef[[i]][my_coef[[i]] != 0, ]
}

# Combine all coefficients across iterations
final_coef <- c()
for (i in 1:10) {
  final_coef <- c(final_coef, my_coef[[i]])
}

# Debugging: Check if coefficients are retained
print("Final coefficients:")
print(final_coef)
# Combine all coefficients across iterations
final_coef <- c()
for (i in 1:length(my_coef)) {
  final_coef <- c(final_coef, my_coef[[i]])
}

# Define a dynamic threshold for filtering over 6 iteration 
final_coef <- final_coef[names(final_coef) %in%
                           names(table(names(final_coef))[table(names(final_coef)) > 6])]

# Debug: Check filtered coefficients
if (length(final_coef) == 0) {
  stop("No coefficients met the frequency threshold. Adjust filtering criteria.")
}
print("Filtered coefficients:")
print(final_coef)

# Create a matrix for coefficients
coef <- matrix(NA, nrow = 10, ncol = length(table(names(final_coef))))
for (i in seq_along(table(names(final_coef)))) {
  c <- final_coef[names(final_coef) %in% names(table(names(final_coef)))[i]]
  coef[1:length(c), i] <- c
}
colnames(coef) <- names(table(names(final_coef)))

# Calculate mean and standard deviation of selected coefficients
est <- apply(coef, 2, mean, na.rm = TRUE)
est_sd <- apply(coef, 2, sd, na.rm = TRUE)

# Create a data frame for the coefficients
DF_COEF_int6 <- data.frame(est, est_sd)
# Calculate 95% Confidence Intervals
lower <- est - 1.96 * est_sd
upper <- est + 1.96 * est_sd

# Create a data frame for the coefficients including 95% CIs
DF_COEF_int6 <- data.frame(
  Mean = est,
  SD = est_sd,
  Lower95CI = lower,
  Upper95CI = upper
)
# Print the resulting data frame
print(DF_COEF_int6)

DF_COEF_int6$Metabolite <- rownames(DF_COEF_int6)
#write.xlsx(DF_COEF_int6, file = "Serum_int_coef.xlsx", rowNames = TRUE)

# Reorder for plotting (optional: by mean effect size)
DF_COEF_int6 <- DF_COEF_int6 %>%
  arrange(desc(Mean))

# Create a label column with Mean and 95% CI
DF_COEF_int6$Label <- sprintf("%.2f [%.2f, %.2f]", 
                              DF_COEF_int6$Mean, 
                              DF_COEF_int6$Lower95CI, 
                              DF_COEF_int6$Upper95CI)

# Forest plot with labels to the right of the point
p <-ggplot(DF_COEF_int6, aes(x = reorder(Metabolite, Mean), y = Mean)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Lower95CI, ymax = Upper95CI), width = 0.5) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(
    title = "Forest Plot of Serum Metabolites Associated with Internalizing Symptoms",
    x = "Metabolite",
    y = "Mean Coefficient (95% CI)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0),
    axis.text.y = element_text(color = "black")  # y-axis labels in black
  )

# ============================================================
# RMSE Evaluation Using Pre-Estimated Coefficients
# 10 Iterations (90% Train / 10% Test)
# ============================================================
# Coefficients and their names
metabolite_names <- names(est)
mean_values <- as.numeric(est)

# Ensure the column names in X match the names in the coefficients
if (!all(names(est) %in% colnames(X))) {
  stop("Mismatch between coefficient names and column names in X.")
}

# Initialize a list to store RMSE values
rmse <- list()

# Set seed for reproducibility
set.seed(123)

# Loop for 10 iterations to calculate RMSE
for (i in 1:10) {
  # Create train-test split (90% train, 10% test)
  ind <- sample(1:nrow(X), size = nrow(X) * 0.9, replace = FALSE)
  train <- X[ind, ]
  test <- X[-ind, ]
  y_train <- y[ind]
  y_test <- y[-ind]
  
  # Calculate predicted scores for the test set
  pred <- rowSums(sapply(names(est), function(var) {
    if (var %in% colnames(test)) {
      test[, var] * est[var]
    } else {
      0  # If variable is missing, contribute 0
    }
  }))
  
  # Calculate RMSE for the test set
  rmse[[i]] <- sqrt(mean((y_test - pred)^2, na.rm = TRUE))
}

# Calculate the mean and 95% confidence interval of RMSE
mean_rmse <- mean(unlist(rmse))
sd_rmse <- sd(unlist(rmse))
sqrt_n <- sqrt(10)

# 95% Confidence intervals
upper <- mean_rmse + 1.96 * (sd_rmse / sqrt_n)
lower <- mean_rmse - 1.96 * (sd_rmse / sqrt_n)

# Output results
cat("Mean RMSE:", mean_rmse, "\n")
cat("95% CI of RMSE:", lower, "to", upper, "\n")


# Round to 3 decimals
mean_rmse_round <- round(mean_rmse, 3)
lower_round <- round(lower, 3)
upper_round <- round(upper, 3)

# Print to console
cat("Mean RMSE:", mean_rmse_round, "\n")
cat("95% CI of RMSE:", lower_round, "to", upper_round, "\n")

# Save to Excel

results <- data.frame(
  Mean_RMSE = mean_rmse_round,
  Lower95CI = lower_round,
  Upper95CI = upper_round
)

#write_xlsx(results, "all figures tables/serum_int_RMSE_results.xlsx")

# ============================================================
# Internalizing Signature Score Calculation
# ============================================================
X <- as.data.frame(X)

required_columns <- c("C16:1-OH", "C3", "C3-DC (C4-OH)", "Gln", "lysoPC a C17:0", "Glu")

# Check if all required columns exist in X
missing_columns <- setdiff(required_columns, colnames(X))
if (length(missing_columns) > 0) {
  stop(paste("ERROR: The following required columns are missing in X:", 
             paste(missing_columns, collapse = ", ")))
}

# Compute Internalizing Signature Score
X$score_internalizing <- (-0.1262098 * X[["C16:1-OH"]]) +
  (-0.08253165 * X[["C3"]]) +
  (-0.07723611 * X[["C3-DC (C4-OH)"]]) +
  (-0.06750363 * X[["Gln"]]) +
  (-0.06060119 * X[["lysoPC a C17:0"]]) +
  ( 0.07189446 * X[["Glu"]])

# Display the calculated Internalizing signature scores
print(X$score_internalizing)
summary(X$score_internalizing)

# Compute z-score normalization for Internalizing Signature Score
X$zscore_internalizing <- (X$score_internalizing - mean(X$score_internalizing, na.rm = TRUE)) / 
  sd(X$score_internalizing, na.rm = TRUE)
# Show summary statistics of the z-score
summary(X$zscore_internalizing)

# ============================================================
# Pearson Correlation Analysis
# ============================================================
correlation_result_intern <- cor(X$score_internalizing, y, method = "pearson", use = "complete.obs")

# ============================================================
# Pearson Correlation Stability (1000 Iterations)
# 90% Subsampling Procedure
# ============================================================
set.seed(238)

# Create a dataframe to store Pearson correlation results
pear <- data.frame(matrix(nrow = 1000, ncol = 2)) 
colnames(pear) <- c("estimate", "p.value")  

# Run 1000 iterations for cross-validation
for (i in 1:1000) {
  
  # Sample 90% of the data for training
  ind <- sample(seq_len(nrow(X)), size = round(0.9 * nrow(X)), replace = FALSE)
  subset90 <- X[ind, ]
  y_subset <- y[ind]
  
  # Ensure 'score_internalizing' and 'y_subset' are numeric
  subset90$score_internalizing <- as.numeric(subset90$score_internalizing)
  y_subset <- as.numeric(y_subset)
  
  # Compute Pearson correlation
  cor_test <- cor.test(subset90$score_internalizing, y_subset, method = "pearson")
  
  # Store correlation coefficient and p-value
  pear[i,1] <- cor_test$estimate  
  pear[i,2] <- cor_test$p.value    
}

# Compute summary statistics
df_pearson <- data.frame(
  Mean = mean(pear$estimate, na.rm = TRUE),
  SD = sd(pear$estimate, na.rm = TRUE)
)

# Compute 95% Confidence Interval
df_pearson$Lower95CI <- df_pearson$Mean - 1.96 * (df_pearson$SD / sqrt(1000))
df_pearson$Upper95CI <- df_pearson$Mean + 1.96 * (df_pearson$SD / sqrt(1000))
df_pearson$p_value <- mean(pear$p.value, na.rm = TRUE) 

# Print final Pearson correlation results
print(df_pearson) 

# ============================================================
# R² Stability Analysis (1000 Iterations)
# ============================================================
set.seed(238)

# Create a dataframe to store R² results
squared <- data.frame(matrix(nrow = 1000, ncol = 2))  
colnames(squared) <- c("r.squared", "p.value")  

# Run 1000 iterations for cross-validation
for (i in 1:1000) {
  
  # Sample 90% of the data for training
  ind <- sample(seq_len(nrow(X)), size = round(0.9 * nrow(X)), replace = FALSE)
  subset90 <- X[ind, ]
  y_subset <- y[ind]
  
  # Ensure 'score_internalizing' and 'y_subset' are numeric
  subset90$score_internalizing <- as.numeric(subset90$score_internalizing)
  y_subset <- as.numeric(y_subset)
  
  # Fit a linear regression model
  model <- lm(y_subset ~ subset90$score_internalizing)
  
  # Store R-squared and p-value
  squared[i,1] <- summary(model)$r.squared  
  squared[i,2] <- summary(model)$coefficients[2,4]  # P-value of the coefficient
}

# Compute summary statistics
df_r2 <- data.frame(
  Mean = mean(squared$r.squared, na.rm = TRUE),
  SD = sd(squared$r.squared, na.rm = TRUE)
)

# Compute 95% Confidence Interval
df_r2$Lower95CI <- df_r2$Mean - 1.96 * (df_r2$SD / sqrt(1000))
df_r2$Upper95CI <- df_r2$Mean + 1.96 * (df_r2$SD / sqrt(1000))
df_r2$p_value <- mean(squared$p.value, na.rm = TRUE) 

# Print final R-squared results
print(df_r2)

# ============================================================
# Final Performance Metrics Export
# ============================================================

# --- Round helper ---
round_df <- function(df, digits = 3) {
  as.data.frame(lapply(df, function(x) if (is.numeric(x)) round(x, digits) else x),
                stringsAsFactors = FALSE)
}

# --- Pearson full ---
pearson_full_df <- data.frame(
  Metric    = "Pearson_r_full",
  Mean      = round(as.numeric(correlation_result_intern), 3),
  SD        = NA,
  Lower95CI = NA,
  Upper95CI = NA,
  p_value   = NA
)

# --- Pearson CV summary ---
pearson_cv_df <- cbind(Metric = "Pearson_CV", round_df(df_pearson, 3))

# --- R² CV summary ---
r2_cv_df <- cbind(Metric = "R2_CV", round_df(df_r2, 3))

# --- RMSE summary (you already computed mean_rmse, lower, upper) ---
rmse_df <- data.frame(
  Metric    = "RMSE",
  Mean      = round(mean_rmse, 3),
  SD        = NA,
  Lower95CI = round(lower, 3),
  Upper95CI = round(upper, 3),
  p_value   = NA
)

# --- Combine all metrics ---
all_metrics <- rbind(pearson_full_df, pearson_cv_df, r2_cv_df, rmse_df)

# --- Save to Excel ---
#write_xlsx(all_metrics, "all figures tables/serum_int_pred_metrics.xlsx")

# ============================================================
# Spearman Correlation Matrix
# Selected Metabolites
# ============================================================
selected_metabs <- c("C16:1-OH", "C3", "C3-DC (C4-OH)", "Gln", "lysoPC a C17:0", "Glu")
data_selected <- merged_sqrt_cbcl_pro_int6[, selected_metabs]

cor_matrix <- cor(data_selected, use = "pairwise.complete.obs", method = "spearman")
library(corrplot)
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, addCoef.col = "black")

###############################################################################################
# Externalizing symptoms for serum metabolites 
rm(list = ls())
# Load the merged datasets with the square root transformed symptpoms and the metabolomics (imputed/normalized)
merged_sqrt_cbcl_pro_ext6 <- read_excel("data/datasets/intern_extern_data/merget outcome-metab-normal/merged_sqrt_cbcl_pro_ext6.xlsx")

# Create Predictor Matrix (Metabolites Only)
X <- as.matrix(merged_sqrt_cbcl_pro_ext6 %>%
                 select(where(is.numeric)) %>%
                 select(-sqrt_cbcl_ext6, -cbcl_ext6))
symptom <- "sqrt_cbcl_ext6"
y <- merged_sqrt_cbcl_pro_ext6[[symptom]]  # Extract response variable for the selected symptom

# Initialize lists to store results
fit <- list()
coef <- list()
alpha <- list()
lambda <- list()
rmse <- list()

# 10-iteration loop for alpha and lambda tuning
set.seed(1)
for (i in 1:10) {
  cat("Iteration:", i, "Symptom:", symptom, "\n")
  
  # Split data into training (90%) and testing (10%) sets
  train_idx <- sample(1:nrow(X), size = 0.9 * nrow(X), replace = FALSE)
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  X_test <- X[-train_idx, ]
  y_test <- y[-train_idx]
  
  # Tune alpha and lambda using cross-validation
  cv_results <- lapply(seq(0.1, 1, by = 0.05), function(alpha_val) {
    cv_model <- cv.glmnet(X_train, y_train, alpha = alpha_val, nfolds = 10, family = "gaussian", type.measure = "mse")
    list(alpha = alpha_val, lambda = cv_model$lambda.min, cv_model = cv_model)
  })
  
  # Retrieve best alpha and lambda based on minimum cross-validation error
  best_result <- cv_results[[which.min(sapply(cv_results, function(res) min(res$cv_model$cvm)))]]
  alpha[[i]] <- best_result$alpha
  lambda[[i]] <- best_result$lambda
  
  # Fit the Elastic Net model using the best alpha and lambda
  fit[[i]] <- glmnet(X_train, y_train, alpha = alpha[[i]], lambda = lambda[[i]], family = "gaussian")
  
  # Extract non-zero coefficients
  coef[[i]] <- coef(fit[[i]])
  
  # Predict on the test set and calculate RMSE
  y_pred <- predict(fit[[i]], newx = X_test)
  rmse[[i]] <- sqrt(mean((y_test - y_pred)^2))
}

# Consolidate results into a data frame
results_extern<- data.frame(
  RMSE = unlist(rmse),
  Alpha = unlist(alpha),
  Lambda = unlist(lambda)
)

# Save results to a CSV file
#output_path_exter <- paste0("C:/Users/lenab/OneDrive/Υπολογιστής/thesis/data/datasets/alpha_lambda_cbcbl_pro_ext6(10iter).csv")
#write.csv(results_extern, file = output_path_exter, row.names = FALSE)

cat("Results saved for symptom:", symptom, "\n")

best_combination_extern <- results_extern[which.min(results_extern$RMSE), ]

print(best_combination_extern)


# Define alpha and lambda based on the tuning results

alpha <- 0.75 
lambda <- 0.1249258


# Prepare data for Elastic Net
X <- as.matrix(merged_sqrt_cbcl_pro_ext6 %>%
                 select(where(is.numeric)) %>%
                 select(-sqrt_cbcl_ext6, -cbcl_ext6))
symptom <- "sqrt_cbcl_ext6"
y <- merged_sqrt_cbcl_pro_ext6[[symptom]]  # Extract response variable for the selected symptom

# Debugging: Check data
print(summary(X))
print(summary(y))
if (any(is.na(X)) || any(is.na(y))) {
  stop("Missing values detected in X or y.")
}

# Initialize a list to store coefficients from iterations
my_coef <- list()

# Set seed for reproducibility
set.seed(124876)

# Loop for coefficient selection
for (i in 1:10) {
  print(paste("Iteration:", i))
  
  # Split data into training (90%) and testing (10%) sets
  train_idx <- sample(1:nrow(X), size = 0.9 * nrow(X), replace = FALSE)
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  
  # Fit Elastic Net model using pre-selected alpha and lambda
  fit <- glmnet(X_train, y_train, family = "gaussian", alpha = alpha, lambda = lambda, intercept = FALSE)
  
  # Extract non-zero coefficients
  my_coef[[i]] <- as.matrix(coef(fit))
  my_coef[[i]] <- my_coef[[i]][my_coef[[i]] != 0, ]
}

# Combine all coefficients across iterations
final_coef <- unlist(my_coef)

# Filter coefficients that appear in more than 5 iterations
final_coef_table <- table(names(final_coef))
filtered_coef <- names(final_coef_table[final_coef_table > 6])

# Debug: Check filtered coefficients
if (length(filtered_coef) == 0) {
  stop("No coefficients met the frequency threshold (more than 5 iterations).")
}
print("Filtered coefficients:")
print(filtered_coef)

# Create a matrix for coefficients
coef_matrix <- matrix(NA, nrow = 10, ncol = length(filtered_coef))
colnames(coef_matrix) <- filtered_coef

# Populate the coefficient matrix for selected coefficients
for (i in 1:10) {
  for (j in filtered_coef) {
    if (j %in% names(my_coef[[i]])) {
      coef_matrix[i, j] <- my_coef[[i]][j]
    }
  }
}

# Calculate mean and standard deviation of selected coefficients
est <- apply(coef_matrix, 2, mean, na.rm = TRUE)
est_sd <- apply(coef_matrix, 2, sd, na.rm = TRUE)
metabolite_names <- names(est)
mean_values <- as.numeric(est)

# Calculate 95% Confidence Intervals
lower <- est - 1.96 * est_sd
upper <- est + 1.96 * est_sd

# Create a data frame for the coefficients including 95% CIs
DF_COEF_ext6 <- data.frame(
  Mean = est,
  SD = est_sd,
  Lower95CI = lower,
  Upper95CI = upper
)

# Print the resulting data frame
print(DF_COEF_ext6)

#write.xlsx(DF_COEF_ext6, file = "Serum_extern_COEF_ext6.xlsx", rowNames = TRUE)


#Create forest plot of the selected metabolite 
DF_COEF_ext6$Metabolite <- rownames(DF_COEF_ext6)

# Reorder for plotting (optional: by mean effect size)
DF_COEF_ext6 <- DF_COEF_ext6 %>%
  arrange(desc(Mean))

# Create a label column with Mean and 95% CI
DF_COEF_ext6$Label <- sprintf("%.2f [%.2f, %.2f]", 
                              DF_COEF_ext6$Mean, 
                              DF_COEF_ext6$Lower95CI, 
                              DF_COEF_ext6$Upper95CI)

ggplot(DF_COEF_ext6, aes(x = reorder(Metabolite, Mean), y = Mean)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Lower95CI, ymax = Upper95CI), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(
    title = "Forest Plot of Serum Metabolites Associated with Externalizing Symptoms",
    x = "Metabolite",
    y = "Mean Coefficient (95% CI)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0),
    axis.text.y = element_text(color = "black") # <-- y-axis labels (metabolites)
  )

# Calculate the prediction RMSE 

metabolite_names <- names(est)
mean_values <- as.numeric(est)

# Ensure the column names in X match the names in the coefficients
if (!all(names(est) %in% colnames(X))) {
  stop("Mismatch between coefficient names and column names in X.")
}

# Initialize a list to store RMSE values
rmse <- list()

# Set seed for reproducibility
set.seed(123)

# Loop for 10 iterations to calculate RMSE
for (i in 1:10) {
  # Create train-test split (90% train, 10% test)
  ind <- sample(1:nrow(X), size = nrow(X) * 0.9, replace = FALSE)
  train <- X[ind, ]
  test <- X[-ind, ]
  y_train <- y[ind]
  y_test <- y[-ind]
  
  # Calculate predicted scores for the test set
  pred <- rowSums(sapply(names(est), function(var) {
    if (var %in% colnames(test)) {
      test[, var] * est[var]
    } else {
      0  # If variable is missing, contribute 0
    }
  }))
  
  # Calculate RMSE for the test set
  rmse[[i]] <- sqrt(mean((y_test - pred)^2, na.rm = TRUE))
}

# Convert list to numeric vector
rmse_values <- unlist(rmse)

# Calculate mean and SD
mean_rmse <- mean(rmse_values)
sd_rmse   <- sd(rmse_values)

# 95% CI
ci_margin <- 1.96 * (sd_rmse / sqrt(length(rmse_values)))
lower     <- mean_rmse - ci_margin
upper     <- mean_rmse + ci_margin



# Round to 3 decimals
mean_rmse_round <- round(mean_rmse, 3)
lower_round <- round(lower, 3)
upper_round <- round(upper, 3)

# Print to console
cat("Mean RMSE:", mean_rmse_round, "\n")
cat("95% CI of RMSE:", lower_round, "to", upper_round, "\n")

results <- data.frame(
  Mean_RMSE = mean_rmse_round,
  Lower95CI = lower_round,
  Upper95CI = upper_round
)

#write_xlsx(results, "RMSE_results.xlsx")


###################################################################################

# Calculate score and z score 
X <- as.data.frame(X)

# Define required metabolites
required_columns <- c(
  "ADMA", "C10:2", "C16:1-OH", "Glu", 
  "lysoPC a C17:0", "PC aa C38:4", "SDMA"
)

# Check if all required columns exist in X
missing_columns <- setdiff(required_columns, colnames(X))
if (length(missing_columns) > 0) {
  stop(paste("ERROR: The following required columns are missing in X:", 
             paste(missing_columns, collapse = ", ")))
}

# Compute Externalizing Signature Score
X$score_externalizing <- (0.11755724  * X[["ADMA"]]) +
  (0.07197259  * X[["C10:2"]]) +
  (-0.18074912 * X[["C16:1-OH"]]) + 
  (0.13879549  * X[["Glu"]]) + 
  (-0.07380767 * X[["lysoPC a C17:0"]]) + 
  (0.20066734  * X[["PC aa C38:4"]]) +   
  (-0.09243611 * X[["SDMA"]])

# Display the calculated Externalizing signature scores
print(X$score_externalizing)
summary(X$score_externalizing)

# Compute Z-score normalization for Externalizing Signature Score
X$zscore_externalizing <- (X$score_externalizing - mean(X$score_externalizing, na.rm = TRUE)) / 
  sd(X$score_externalizing, na.rm = TRUE)

# Show summary statistics of the z-score
summary(X$zscore_externalizing)

####################################################################################
# Pearson colleration 
correlation_result_extern <- cor(X$score_externalizing, y, method = "pearson", use = "complete.obs")

######################################################################################
set.seed(238)
# Create a dataframe to store Pearson correlation results for Externalizing
pear_externalizing <- data.frame(matrix(nrow = 1000, ncol = 2)) 
colnames(pear_externalizing) <- c("estimate", "p.value")  

# Run 1000 iterations for cross-validation
for (i in 1:1000) {
  
  # Sample 90% of the data for training
  ind <- sample(seq_len(nrow(X)), size = round(0.9 * nrow(X)), replace = FALSE)
  subset90 <- X[ind, ]
  y_subset <- y[ind]  # Ensure the correct response variable is used
  
  # Ensure 'score_externalizing' and 'y_subset' are numeric
  subset90$score_externalizing <- as.numeric(subset90$score_externalizing)
  y_subset <- as.numeric(y_subset)
  
  # Compute Pearson correlation
  cor_test <- cor.test(subset90$score_externalizing, y_subset, method = "pearson")
  
  # Store correlation coefficient and p-value
  pear_externalizing[i,1] <- cor_test$estimate  
  pear_externalizing[i,2] <- cor_test$p.value    
}

# Compute summary statistics
df_pearson_externalizing <- data.frame(
  Mean = mean(pear_externalizing$estimate, na.rm = TRUE),
  SD = sd(pear_externalizing$estimate, na.rm = TRUE)
)

# Compute 95% Confidence Interval
df_pearson_externalizing$Lower95CI <- df_pearson_externalizing$Mean - 1.96 * (df_pearson_externalizing$SD / sqrt(1000))
df_pearson_externalizing$Upper95CI <- df_pearson_externalizing$Mean + 1.96 * (df_pearson_externalizing$SD / sqrt(1000))
df_pearson_externalizing$p_value <- mean(pear_externalizing$p.value, na.rm = TRUE) 

# Print final Pearson correlation results
print(df_pearson_externalizing)

############## R² for Externalizing Signature Score ####################################
# Set seed for reproducibility
set.seed(238)

# Create a dataframe to store R² results for Externalizing
squared_externalizing <- data.frame(matrix(nrow = 1000, ncol = 2))  
colnames(squared_externalizing) <- c("r.squared", "p.value")  

# Run 1000 iterations for cross-validation
for (i in 1:1000) {
  
  # Sample 90% of the data for training
  ind <- sample(seq_len(nrow(X)), size = round(0.9 * nrow(X)), replace = FALSE)
  subset90 <- X[ind, ]
  y_subset <- y[ind]
  
  # Ensure 'score_externalizing' and 'y_subset' are numeric
  subset90$score_externalizing <- as.numeric(subset90$score_externalizing)
  y_subset <- as.numeric(y_subset)
  
  # Fit a linear regression model
  model <- lm(y_subset ~ subset90$score_externalizing)
  
  # Store R-squared and p-value
  squared_externalizing[i,1] <- summary(model)$r.squared  
  squared_externalizing[i,2] <- summary(model)$coefficients[2,4]  # P-value of the coefficient
}

# Compute summary statistics
df_r2_externalizing <- data.frame(
  Mean = mean(squared_externalizing$r.squared, na.rm = TRUE),
  SD = sd(squared_externalizing$r.squared, na.rm = TRUE)
)

# Compute 95% Confidence Interval
df_r2_externalizing$Lower95CI <- df_r2_externalizing$Mean - 1.96 * (df_r2_externalizing$SD / sqrt(1000))
df_r2_externalizing$Upper95CI <- df_r2_externalizing$Mean + 1.96 * (df_r2_externalizing$SD / sqrt(1000))
df_r2_externalizing$p_value <- mean(squared_externalizing$p.value, na.rm = TRUE) 

# Print final R-squared results
print(df_r2_externalizing)

###################################################
# --- Round helper ---
round_df <- function(df, digits = 3) {
  as.data.frame(lapply(df, function(x) if (is.numeric(x)) round(x, digits) else x),
                stringsAsFactors = FALSE)
}

# --- Pearson full ---
pearson_full_df <- data.frame(
  Metric    = "Pearson_r_full",
  Mean      = round(as.numeric(correlation_result_extern), 3),
  SD        = NA,
  Lower95CI = NA,
  Upper95CI = NA,
  p_value   = NA
)

# --- Pearson CV summary ---
pearson_cv_df <- cbind(Metric = "Pearson_CV", round_df(df_pearson_externalizing, 3))

# --- R² CV summary ---
r2_cv_df <- cbind(Metric = "R2_CV", round_df(df_r2_externalizing, 3))

# --- RMSE summary (you already computed mean_rmse, lower, upper) ---
rmse_df <- data.frame(
  Metric    = "RMSE",
  Mean      = round(mean_rmse, 3),
  SD        = NA,
  Lower95CI = round(lower, 3),
  Upper95CI = round(upper, 3),
  p_value   = NA
)

# --- Combine all metrics ---
all_metrics <- rbind(pearson_full_df, pearson_cv_df, r2_cv_df, rmse_df)

# --- Save to Excel ---
library(writexl)
write_xlsx(all_metrics, "externalizing_metrics_with_RMSE.xlsx")

# Check save location
getwd()


#CORRELATION MATRIX

selected_metabs <- c("ADMA", "C10:2", "C16:1-OH", "Glu", "lysoPC a C17:0", "PC aa C38:4", "SDMA")
data_selected <- merged_sqrt_cbcl_pro_ext6[, selected_metabs]

cor_matrix <- cor(data_selected, use = "pairwise.complete.obs", method = "spearman")
library(corrplot)
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, addCoef.col = "black")

###################################################################################################
###############################################################################################
# ADHD symptoms for serum metabolites 
# Clean the environment 
rm(list = ls())
# Load the merged datasets with the square root transformed symptpoms and the metabolomics (imputed/normalized)
merged_sqrt_cbcl_pro_adhd6 <- read_excel("data/datasets/intern_extern_data/merget outcome-metab-normal/merged_sqrt_cbcl_pro_adhd6.xlsx")

View(merged_sqrt_cbcl_pro_adhd6)


X <- as.matrix(merged_sqrt_cbcl_pro_adhd6 %>%
                 select(where(is.numeric)) %>%
                 select(-adhd_pro6, -sqrt_adhd_pro6))

symptom <- "sqrt_adhd_pro6"
y <- merged_sqrt_cbcl_pro_adhd6[[symptom]]  # Extract response variable for the selected symptom

# Initialize lists to store results
fit <- list()
coef <- list()
alpha <- list()
lambda <- list()
rmse <- list()

# 10-iteration loop for alpha and lambda tuning
set.seed(1)
for (i in 1:10) {
  cat("Iteration:", i, "Symptom:", symptom, "\n")
  
  # Split data into training (90%) and testing (10%) sets
  train_idx <- sample(1:nrow(X), size = 0.9 * nrow(X), replace = FALSE)
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  X_test <- X[-train_idx, ]
  y_test <- y[-train_idx]
  
  # Tune alpha and lambda using cross-validation
  cv_results <- lapply(seq(0.1, 1, by = 0.05), function(alpha_val) {
    cv_model <- cv.glmnet(X_train, y_train, alpha = alpha_val, nfolds = 10, family = "gaussian", type.measure = "mse")
    list(alpha = alpha_val, lambda = cv_model$lambda.min, cv_model = cv_model)
  })
  
  # Retrieve best alpha and lambda based on minimum cross-validation error
  best_result <- cv_results[[which.min(sapply(cv_results, function(res) min(res$cv_model$cvm)))]]
  alpha[[i]] <- best_result$alpha
  lambda[[i]] <- best_result$lambda
  
  # Fit the Elastic Net model using the best alpha and lambda
  fit[[i]] <- glmnet(X_train, y_train, alpha = alpha[[i]], lambda = lambda[[i]], family = "gaussian")
  
  # Extract non-zero coefficients
  coef[[i]] <- coef(fit[[i]])
  
  # Predict on the test set and calculate RMSE
  y_pred <- predict(fit[[i]], newx = X_test)
  rmse[[i]] <- sqrt(mean((y_test - y_pred)^2))
}

# Consolidate results into a data frame
results_adhd <- data.frame(
  RMSE = unlist(rmse),
  Alpha = unlist(alpha),
  Lambda = unlist(lambda)
)

# Save results to a CSV file
#write.xlsx(results_adhd, "all figures tables/serum_adhd_alpha_lambda.xlsx", rowNames = FALSE)


cat("Results saved for symptom:", symptom, "\n")

best_combination_adhd <- results_adhd[which.min(results_adhd$RMSE), ]
print(best_combination_adhd)

######################################################################################
# FOR ADHD

alpha <- 0.15  
lambda <- 1.179094 

# Prepare data for Elastic Net
X <- as.matrix(merged_sqrt_cbcl_pro_adhd6 %>%
                 select(where(is.numeric)) %>%
                 select(-adhd_pro6, -sqrt_adhd_pro6))
symptom <- "sqrt_adhd_pro6"
y <- merged_sqrt_cbcl_pro_adhd6[[symptom]]

# Debugging: Check data
print(summary(X))
print(summary(y))
if (any(is.na(X)) || any(is.na(y))) {
  stop("Missing values detected in X or y.")
}

# Initialize a list to store coefficients from iterations
my_coef <- list()

# Set seed for reproducibility
set.seed(124876)

# Loop for coefficient selection
for (i in 1:10) {
  print(paste("Iteration:", i))
  
  # Split data into training (90%) and testing (10%) sets
  train_idx <- sample(1:nrow(X), size = 0.9 * nrow(X), replace = FALSE)
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  
  # Fit Elastic Net model using pre-selected alpha and lambda
  fit <- glmnet(X_train, y_train, family = "gaussian", alpha = alpha, lambda = lambda, intercept = FALSE)
  
  # Extract non-zero coefficients
  my_coef[[i]] <- as.matrix(coef(fit))
  my_coef[[i]] <- my_coef[[i]][my_coef[[i]] != 0, ]
}

# Combine all coefficients across iterations
final_coef <- unlist(my_coef)

# Filter coefficients that appear in more than 6 iterations
final_coef_table <- table(names(final_coef))
filtered_coef <- names(final_coef_table[final_coef_table > 6])

# Debug: Check filtered coefficients
if (length(filtered_coef) == 0) {
  stop("No coefficients met the frequency threshold (more than x iterations).")
}
print("Filtered coefficients:")
print(filtered_coef)

# Create a matrix for coefficients
coef_matrix <- matrix(NA, nrow = 10, ncol = length(filtered_coef))
colnames(coef_matrix) <- filtered_coef

# Populate the coefficient matrix for selected coefficients
for (i in 1:10) {
  for (j in filtered_coef) {
    if (j %in% names(my_coef[[i]])) {
      coef_matrix[i, j] <- my_coef[[i]][j]
    }
  }
}

# Calculate mean and standard deviation of selected coefficients
est <- apply(coef_matrix, 2, mean, na.rm = TRUE)
est_sd <- apply(coef_matrix, 2, sd, na.rm = TRUE)
metabolite_names <- names(est)
mean_values <- as.numeric(est)


# Calculate 95% Confidence Intervals
lower <- est - 1.96 * est_sd
upper <- est + 1.96 * est_sd

# Create a data frame for the coefficients including 95% CIs
DF_COEF_adhd6 <- data.frame(
  Mean = est,
  SD = est_sd,
  Lower95CI = lower,
  Upper95CI = upper
)

# Print the resulting data frame
print(DF_COEF_adhd6)
#write.xlsx(DF_COEF_adhd6, file = "serum_adhd_coef.xlsx", rowNames = TRUE)

library(ggplot2)

DF_COEF_adhd6$Metabolite <- rownames(DF_COEF_adhd6)


# Reorder for plotting (optional: by mean effect size)
DF_COEF_adhd6 <- DF_COEF_adhd6 %>%
  arrange(desc(Mean))

# Create a label column with Mean and 95% CI
DF_COEF_adhd6$Label <- sprintf("%.2f [%.2f, %.2f]", 
                               DF_COEF_adhd6$Mean, 
                               DF_COEF_adhd6$Lower95CI, 
                               DF_COEF_adhd6$Upper95CI)

# Forest plot with labels to the right of the point
ggplot(DF_COEF_adhd6, aes(x = reorder(Metabolite, Mean), y = Mean)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Lower95CI, ymax = Upper95CI), width = 0.2) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(
    title = "Forest Plot of Serum Metabolites Associated with ADHD",
    x = "Metabolite",
    y = "Mean Coefficient (95% CI)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0),
    axis.text.y = element_text(color = "black")  # y-axis labels in black
  )

################################################################################
#CALCULATE THE PREDICTION RMSE 

# Ensure the column names in X match the names in the coefficients
if (!all(names(est) %in% colnames(X))) {
  stop("Mismatch between coefficient names and column names in X.")
}

# Initialize a list to store RMSE values
rmse <- list()

# Set seed for reproducibility
set.seed(123)

# Loop for 10 iterations to calculate RMSE
for (i in 1:10) {
  # Create train-test split (90% train, 10% test)
  ind <- sample(1:nrow(X), size = nrow(X) * 0.9, replace = FALSE)
  train <- X[ind, ]
  test <- X[-ind, ]
  y_train <- y[ind]
  y_test <- y[-ind]
  
  # Calculate predicted scores for the test set
  pred <- rowSums(sapply(names(est), function(var) {
    if (var %in% colnames(test)) {
      test[, var] * est[var]
    } else {
      0  # If variable is missing, contribute 0
    }
  }))
  
  # Calculate RMSE for the test set
  rmse[[i]] <- sqrt(mean((y_test - pred)^2, na.rm = TRUE))
}

# Calculate the mean and 95% confidence interval of RMSE
mean_rmse <- mean(unlist(rmse))
sd_rmse <- sd(unlist(rmse))
sqrt_n <- sqrt(10)  # Number of iterations

# 95% Confidence intervals
upper <- mean_rmse + 1.96 * (sd_rmse / sqrt_n)
lower <- mean_rmse - 1.96 * (sd_rmse / sqrt_n)

# Output results
cat("Mean RMSE (ADHD):", mean_rmse, "\n")
cat("95% CI of RMSE:", lower, "to", upper, "\n")

# Save to Excel
# Round to 3 decimals
mean_rmse_round <- round(mean_rmse, 3)
lower_round <- round(lower, 3)
upper_round <- round(upper, 3)

# Print to console
cat("Mean RMSE:", mean_rmse_round, "\n")
cat("95% CI of RMSE:", lower_round, "to", upper_round, "\n")
results <- data.frame(
  Mean_RMSE = mean_rmse_round,
  Lower95CI = lower_round,
  Upper95CI = upper_round
)

#write_xlsx(results, "all figures tables/serum_adhd_RMSE_results.xlsx")

################################################################################

#CALCULATE SCORE AND ZSCORE 

X <- as.data.frame(X)

# Define required metabolite for ADHD signature
required_columns <- c("C16:1-OH")

# Check if all required columns exist in X
missing_columns <- setdiff(required_columns, colnames(X))
if (length(missing_columns) > 0) {
  stop(paste("ERROR: The following required columns are missing in X:", 
             paste(missing_columns, collapse = ", ")))
}

# Compute ADHD Signature Score (Only using C16:1-OH)
X$score_adhd <- (-0.04002156 * X[["C16:1-OH"]])

# Display the calculated ADHD signature scores
print(X$score_adhd)
summary(X$score_adhd)

# Compute Z-score normalization for ADHD Signature Score
X$zscore_adhd <- (X$score_adhd - mean(X$score_adhd, na.rm = TRUE)) / 
  sd(X$score_adhd, na.rm = TRUE)

# Show summary statistics of the z-score
summary(X$zscore_adhd)

##############################################################################
#SIMPLE COLLERATION 
correlation_result_adhd <- cor(X$score_adhd, y, method = "pearson", use = "complete.obs")

############################################################################

set.seed(238)
# Create a dataframe to store Pearson correlation results for ADHD
pear_adhd <- data.frame(matrix(nrow = 1000, ncol = 2)) 
colnames(pear_adhd) <- c("estimate", "p.value")  

# Run 1000 iterations for cross-validation
for (i in 1:1000) {
  
  # Sample 90% of the data for training
  ind <- sample(seq_len(nrow(X)), size = round(0.9 * nrow(X)), replace = FALSE)
  subset90 <- X[ind, ]
  y_subset <- y[ind]  # Ensure the correct response variable is used
  
  # Ensure 'score_adhd' and 'y_subset' are numeric
  subset90$score_adhd <- as.numeric(subset90$score_adhd)
  y_subset <- as.numeric(y_subset)
  
  # Compute Pearson correlation
  cor_test <- cor.test(subset90$score_adhd, y_subset, method = "pearson")
  
  # Store correlation coefficient and p-value
  pear_adhd[i,1] <- cor_test$estimate  
  pear_adhd[i,2] <- cor_test$p.value    
}

# Compute summary statistics
df_pearson_adhd <- data.frame(
  Mean = mean(pear_adhd$estimate, na.rm = TRUE),
  SD = sd(pear_adhd$estimate, na.rm = TRUE)
)

# Compute 95% Confidence Interval
df_pearson_adhd$Lower95CI <- df_pearson_adhd$Mean - 1.96 * (df_pearson_adhd$SD / sqrt(1000))
df_pearson_adhd$Upper95CI <- df_pearson_adhd$Mean + 1.96 * (df_pearson_adhd$SD / sqrt(1000))
df_pearson_adhd$p_value <- mean(pear_adhd$p.value, na.rm = TRUE) 

# Print final Pearson correlation results
print(df_pearson_adhd)

###############################################################################
# R² for ADHD Signature Score ####################################
# Set seed for reproducibility

set.seed(238)

# Create a dataframe to store R² results for ADHD
squared_adhd <- data.frame(matrix(nrow = 1000, ncol = 2))  
colnames(squared_adhd) <- c("r.squared", "p.value")  

# Run 1000 iterations for cross-validation
for (i in 1:1000) {
  
  # Sample 90% of the data for training
  ind <- sample(seq_len(nrow(X)), size = round(0.9 * nrow(X)), replace = FALSE)
  subset90 <- X[ind, ]
  y_subset <- y[ind]
  
  # Ensure 'score_adhd' and 'y_subset' are numeric
  subset90$score_adhd <- as.numeric(subset90$score_adhd)
  y_subset <- as.numeric(y_subset)
  
  # Fit a linear regression model
  model <- lm(y_subset ~ subset90$score_adhd)
  
  # Store R-squared and p-value
  squared_adhd[i,1] <- summary(model)$r.squared  
  squared_adhd[i,2] <- summary(model)$coefficients[2,4]  # P-value of the coefficient
}

# Compute summary statistics
df_r2_adhd <- data.frame(
  Mean = mean(squared_adhd$r.squared, na.rm = TRUE),
  SD = sd(squared_adhd$r.squared, na.rm = TRUE)
)

# Compute 95% Confidence Interval
df_r2_adhd$Lower95CI <- df_r2_adhd$Mean - 1.96 * (df_r2_adhd$SD / sqrt(1000))
df_r2_adhd$Upper95CI <- df_r2_adhd$Mean + 1.96 * (df_r2_adhd$SD / sqrt(1000))
df_r2_adhd$p_value <- mean(squared_adhd$p.value, na.rm = TRUE) 

# Print final R² results
print(df_r2_adhd)

################################################################################

# --- Round helper ---
round_df <- function(df, digits = 3) {
  as.data.frame(lapply(df, function(x) if (is.numeric(x)) round(x, digits) else x),
                stringsAsFactors = FALSE)
}

# --- Pearson full ---
pearson_full_df <- data.frame(
  Metric    = "Pearson_r_full",
  Mean      = round(as.numeric(correlation_result_adhd), 3),
  SD        = NA,
  Lower95CI = NA,
  Upper95CI = NA,
  p_value   = NA
)

# --- Pearson CV summary ---
pearson_cv_df <- cbind(Metric = "Pearson_CV", round_df(df_pearson_adhd, 3))

# --- R² CV summary ---
r2_cv_df <- cbind(Metric = "R2_CV", round_df(df_r2_adhd, 3))

# --- RMSE summary (you already computed mean_rmse, lower, upper) ---
rmse_df <- data.frame(
  Metric    = "RMSE",
  Mean      = round(mean_rmse, 3),
  SD        = NA,
  Lower95CI = round(lower, 3),
  Upper95CI = round(upper, 3),
  p_value   = NA
)

# --- Combine all metrics ---
all_metrics <- rbind(pearson_full_df, pearson_cv_df, r2_cv_df, rmse_df)

# --- Save to Excel ---
library(writexl)
#write_xlsx(all_metrics, "all figures tables/serum_adhd_pred_metrics.xlsx")

# Check save location
getwd()

#################### Lasso/Ridge Analysis- Linear regression with elastic net ######################
####################################################################################################
####################### Urine Analysis #############################################################
#We are following the same logic for the urine metabolomic data separate internalizing, externalizing, ADHD 
# Internalizing symptoms for urine metabolites 
# Clean the environment 
rm(list = ls())

merged_sqrt_urine_int6 <- read_excel("merged_sqrt_urine_int6.xlsx")

#LASSO ELASTIC NET REGRESSION TO FIND THE ALPHA AND LAMDA FOR INTERNALIZING SYMPTOMS
X <- as.matrix(
  merged_sqrt_urine_int6 %>%
    select(where(is.numeric)) %>%
    select(-sqrt_cbcl_pro_int6, -cbcl_pro_int6)  # exclude both transformed and raw outcome
)

# Step 2: Define the outcome variable
symptom <- "sqrt_cbcl_pro_int6"
y <- merged_sqrt_urine_int6[[symptom]]

# Step 3: Initialize storage lists
fit <- list()
coef <- list()
alpha <- list()
lambda <- list()
rmse <- list()

# Step 4: Run 10 iterations for alpha/lambda tuning
set.seed(1)
for (i in 1:10) {
  cat("Iteration:", i, "Symptom:", symptom, "\n")
  
  # Split into training (90%) and testing (10%)
  train_idx <- sample(1:nrow(X), size = 0.9 * nrow(X), replace = FALSE)
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  X_test <- X[-train_idx, ]
  y_test <- y[-train_idx]
  
  # Tune alpha and lambda
  cv_results <- lapply(seq(0.1, 1, by = 0.05), function(alpha_val) {
    cv_model <- cv.glmnet(X_train, y_train, alpha = alpha_val, nfolds = 10,
                          family = "gaussian", type.measure = "mse")
    list(alpha = alpha_val, lambda = cv_model$lambda.min, cv_model = cv_model)
  })
  
  # Find best result based on lowest CV error
  best_result <- cv_results[[which.min(sapply(cv_results, function(res) min(res$cv_model$cvm)))]]
  alpha[[i]] <- best_result$alpha
  lambda[[i]] <- best_result$lambda
  
  # Fit Elastic Net model
  fit[[i]] <- glmnet(X_train, y_train, alpha = alpha[[i]], lambda = lambda[[i]], family = "gaussian")
  
  # Extract non-zero coefficients
  coef[[i]] <- coef(fit[[i]])
  
  # Predict and calculate RMSE
  y_pred <- predict(fit[[i]], newx = X_test)
  rmse[[i]] <- sqrt(mean((y_test - y_pred)^2))
}

# Step 5: Create results summary table
results_urine_int <- data.frame(
  RMSE = unlist(rmse),
  Alpha = unlist(alpha),
  Lambda = unlist(lambda)
)


# Save results to a CSV file
#write.xlsx(results_urine_int, "all figures tables/urine_int_alpha_lambda.xlsx", rowNames = FALSE)

cat("Results saved for symptom:", symptom, "\n")

best_combination_intern <- results_urine_int[which.min(results_urine_int$RMSE), ]
print(best_combination_intern)


cat("Results saved for symptom:", symptom, "\n")

################################################################################
#COEFFICIENTS DETECTION 
#INTERNALIZING 

# Define alpha and lambda based on the tuning results
alpha <- 0.85  # Adjust alpha if needed
lambda <- 0.08963741
# Adjust lambda if needed


X <- as.matrix(merged_sqrt_urine_int6 %>%
                 select(where(is.numeric)) %>%
                 select(-sqrt_cbcl_pro_int6, -cbcl_pro_int6))
symptom <- "sqrt_cbcl_pro_int6"
y <- merged_sqrt_urine_int6[[symptom]]

# Debugging: Check data
print(summary(X))
print(summary(y))
if (any(is.na(X)) || any(is.na(y))) {
  stop("Missing values detected in X or y.")
}

# Initialize a list to store coefficients from iterations
my_coef <- list()

# Set seed for reproducibility
set.seed(124876)

# Loop for coefficient selection
for (i in 1:10) {
  print(paste("Iteration:", i))
  
  # Split data into training (90%) and testing (10%) sets
  train_idx <- sample(1:nrow(X), size = 0.9 * nrow(X), replace = FALSE)
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  
  # Fit Elastic Net model using pre-selected alpha and lambda
  fit <- glmnet(X_train, y_train, family = "gaussian", alpha = alpha, lambda = lambda, intercept = FALSE)
  
  # Debugging: Check coefficients
  print(coef(fit))
  
  # Extract non-zero coefficients
  my_coef[[i]] <- as.matrix(coef(fit))
  my_coef[[i]] <- my_coef[[i]][my_coef[[i]] != 0, ]
}

# Combine all coefficients across iterations
final_coef <- c()
for (i in 1:10) {
  final_coef <- c(final_coef, my_coef[[i]])
}

# Debugging: Check if coefficients are retained
print("Final coefficients:")
print(final_coef)
# Combine all coefficients across iterations
final_coef <- c()
for (i in 1:length(my_coef)) {
  final_coef <- c(final_coef, my_coef[[i]])
}

# Define a dynamic threshold for filtering


final_coef <- final_coef[names(final_coef) %in%
                           names(table(names(final_coef))[table(names(final_coef)) > 6])]

# Filter coefficients based on frequency
#final_coef <- final_coef[names(final_coef) %in% names(table(names(final_coef))[table(names(final_coef)) > threshold])]

# Debug: Check filtered coefficients
if (length(final_coef) == 0) {
  stop("No coefficients met the frequency threshold. Adjust filtering criteria.")
}
print("Filtered coefficients:")
print(final_coef)

# Create a matrix for coefficients
coef <- matrix(NA, nrow = 10, ncol = length(table(names(final_coef))))
for (i in seq_along(table(names(final_coef)))) {
  c <- final_coef[names(final_coef) %in% names(table(names(final_coef)))[i]]
  coef[1:length(c), i] <- c
}
colnames(coef) <- names(table(names(final_coef)))

# Calculate mean and standard deviation of selected coefficients
est <- apply(coef, 2, mean, na.rm = TRUE)
est_sd <- apply(coef, 2, sd, na.rm = TRUE)

# Create a data frame for the coefficients
DF_COEF_int6 <- data.frame(est, est_sd)
# Calculate 95% Confidence Intervals
lower <- est - 1.96 * est_sd
upper <- est + 1.96 * est_sd

# Create a data frame for the coefficients including 95% CIs
DF_COEF_int6 <- data.frame(
  Mean = est,
  SD = est_sd,
  Lower95CI = lower,
  Upper95CI = upper
)
# Print the resulting data frame
print(DF_COEF_int6)

DF_COEF_int6$Metabolite <- rownames(DF_COEF_int6)
#write.xlsx(DF_COEF_int6, file = "urine_int_coef.xlsx", rowNames = TRUE)


#Forest plot of the selected metabolites 
# Create a column for predictor names
DF_COEF_int6$Metabolite <- rownames(DF_COEF_int6)

# Reorder metabolites by effect size (optional)
DF_COEF_int6 <- DF_COEF_int6 %>%
  arrange(desc(Mean))

# Plot
ggplot(DF_COEF_int6, aes(x = reorder(Metabolite, Mean), y = Mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = Lower95CI, ymax = Upper95CI), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(
    title = "Forest Plot of Urine Metabolites Assosiated with Internalizing symptoms",
    x = "Metabolite",
    y = "Mean Coefficient (95% CI)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = -4),
    axis.text.y = element_text(color = "black")  # y-axis labels in black
  )


#CALCULATE PREDICTION RMSE 

# Coefficients and their names
metabolite_names <- names(est)
mean_values <- as.numeric(est)

# Ensure the column names in X match the names in the coefficients
if (!all(names(est) %in% colnames(X))) {
  stop("Mismatch between coefficient names and column names in X.")
}

# Initialize a list to store RMSE values
rmse <- list()

# Set seed for reproducibility
set.seed(123)

# Loop for 10 iterations to calculate RMSE
for (i in 1:10) {
  # Create train-test split (90% train, 10% test)
  ind <- sample(1:nrow(X), size = nrow(X) * 0.9, replace = FALSE)
  train <- X[ind, ]
  test <- X[-ind, ]
  y_train <- y[ind]
  y_test <- y[-ind]
  
  # Calculate predicted scores for the test set
  pred <- rowSums(sapply(names(est), function(var) {
    if (var %in% colnames(test)) {
      test[, var] * est[var]
    } else {
      0  # If variable is missing, contribute 0
    }
  }))
  
  # Calculate RMSE for the test set
  rmse[[i]] <- sqrt(mean((y_test - pred)^2, na.rm = TRUE))
}

# Calculate the mean and 95% confidence interval of RMSE
mean_rmse <- mean(unlist(rmse))
sd_rmse <- sd(unlist(rmse))
sqrt_n <- sqrt(10)

# 95% Confidence intervals
upper <- mean_rmse + 1.96 * (sd_rmse / sqrt_n)
lower <- mean_rmse - 1.96 * (sd_rmse / sqrt_n)

# Output results
cat("Mean RMSE:", mean_rmse, "\n")
cat("95% CI of RMSE:", lower, "to", upper, "\n")

#########################################################################################
#CALCULATE THE SCORE AND THE Z SCORE 

X <- as.data.frame(X)

required_columns <- c(
  "3-hydroxyisovalerate",
  "5-oxoproline",
  "Dimethylamine",
  "Glutamine",
  "p-cresol sulfate",
  "Tyrosine"
)
# Check if all required columns exist in X
missing_columns <- setdiff(required_columns, colnames(X))
if (length(missing_columns) > 0) {
  stop(paste("ERROR: The following required columns are missing in X:", 
             paste(missing_columns, collapse = ", ")))
}

# Compute Internalizing Signature Score
X$score_internalizing <- (
  -0.03230146 * X[["3-hydroxyisovalerate"]] +
    -0.04665104 * X[["5-oxoproline"]] +
    -0.04922661 * X[["Dimethylamine"]] +
    0.06439756 * X[["Glutamine"]] +
    -0.06812796 * X[["p-cresol sulfate"]] +
    0.12402503 * X[["Tyrosine"]]
)

# Display the calculated Internalizing signature scores
print(X$score_internalizing)
summary(X$score_internalizing)

# Compute z-score normalization for Internalizing Signature Score
X$zscore_internalizing <- (X$score_internalizing - mean(X$score_internalizing, na.rm = TRUE)) / 
  sd(X$score_internalizing, na.rm = TRUE)
# Show summary statistics of the z-score
summary(X$zscore_internalizing)

#########################################################################################

#SIMPLE PEARSON COLLERATION 
correlation_result_intern <- cor(X$score_internalizing, y, method = "pearson", use = "complete.obs")

####################################################################################
#PEARSON COLLERATION WITH 1000 ITERATIONS 

set.seed(238)

# Create a dataframe to store Pearson correlation results
pear <- data.frame(matrix(nrow = 1000, ncol = 2)) 
colnames(pear) <- c("estimate", "p.value")  

# Run 1000 iterations for cross-validation
for (i in 1:1000) {
  
  # Sample 90% of the data for training
  ind <- sample(seq_len(nrow(X)), size = round(0.9 * nrow(X)), replace = FALSE)
  subset90 <- X[ind, ]
  y_subset <- y[ind]
  
  # Ensure 'score_internalizing' and 'y_subset' are numeric
  subset90$score_internalizing <- as.numeric(subset90$score_internalizing)
  y_subset <- as.numeric(y_subset)
  
  # Compute Pearson correlation
  cor_test <- cor.test(subset90$score_internalizing, y_subset, method = "pearson")
  
  # Store correlation coefficient and p-value
  pear[i,1] <- cor_test$estimate  
  pear[i,2] <- cor_test$p.value    
}

# Compute summary statistics
df_pearson <- data.frame(
  Mean = mean(pear$estimate, na.rm = TRUE),
  SD = sd(pear$estimate, na.rm = TRUE)
)

# Compute 95% Confidence Interval
df_pearson$Lower95CI <- df_pearson$Mean - 1.96 * (df_pearson$SD / sqrt(1000))
df_pearson$Upper95CI <- df_pearson$Mean + 1.96 * (df_pearson$SD / sqrt(1000))
df_pearson$p_value <- mean(pear$p.value, na.rm = TRUE) 

# Print final Pearson correlation results
print(df_pearson) 

########################################################################################
############## R2 ################################################################
# Set seed for reproducibility
set.seed(238)

# Create a dataframe to store R² results
squared <- data.frame(matrix(nrow = 1000, ncol = 2))  
colnames(squared) <- c("r.squared", "p.value")  

# Run 1000 iterations for cross-validation
for (i in 1:1000) {
  
  # Sample 90% of the data for training
  ind <- sample(seq_len(nrow(X)), size = round(0.9 * nrow(X)), replace = FALSE)
  subset90 <- X[ind, ]
  y_subset <- y[ind]
  
  # Ensure 'score_internalizing' and 'y_subset' are numeric
  subset90$score_internalizing <- as.numeric(subset90$score_internalizing)
  y_subset <- as.numeric(y_subset)
  
  # Fit a linear regression model
  model <- lm(y_subset ~ subset90$score_internalizing)
  
  # Store R-squared and p-value
  squared[i,1] <- summary(model)$r.squared  
  squared[i,2] <- summary(model)$coefficients[2,4]  # P-value of the coefficient
}

# Compute summary statistics
df_r2 <- data.frame(
  Mean = mean(squared$r.squared, na.rm = TRUE),
  SD = sd(squared$r.squared, na.rm = TRUE)
)

# Compute 95% Confidence Interval
df_r2$Lower95CI <- df_r2$Mean - 1.96 * (df_r2$SD / sqrt(1000))
df_r2$Upper95CI <- df_r2$Mean + 1.96 * (df_r2$SD / sqrt(1000))
df_r2$p_value <- mean(squared$p.value, na.rm = TRUE) 

# Print final R-squared results
print(df_r2)


###############################################################################

# --- Round helper ---
round_df <- function(df, digits = 3) {
  as.data.frame(lapply(df, function(x) if (is.numeric(x)) round(x, digits) else x),
                stringsAsFactors = FALSE)
}

# --- Pearson full ---
pearson_full_df <- data.frame(
  Metric    = "Pearson_r_full",
  Mean      = round(as.numeric(correlation_result_intern), 3),
  SD        = NA,
  Lower95CI = NA,
  Upper95CI = NA,
  p_value   = NA
)

# --- Pearson CV summary ---
pearson_cv_df <- cbind(Metric = "Pearson_CV", round_df(df_pearson, 3))

# --- R² CV summary ---
r2_cv_df <- cbind(Metric = "R2_CV", round_df(df_r2, 3))

# --- RMSE summary (you already computed mean_rmse, lower, upper) ---
rmse_df <- data.frame(
  Metric    = "RMSE",
  Mean      = round(mean_rmse, 3),
  SD        = NA,
  Lower95CI = round(lower, 3),
  Upper95CI = round(upper, 3),
  p_value   = NA
)

# --- Combine all metrics ---
all_metrics <- rbind(pearson_full_df, pearson_cv_df, r2_cv_df, rmse_df)

# --- Save to Excel ---
library(writexl)
#write_xlsx(all_metrics, "all figures tables/urine_int_pred_metrics.xlsx")

#CORRELATION 

selected_metabs <- c("3-hydroxyisovalerate", 
                     "5-oxoproline", 
                     "Dimethylamine", 
                     "Glutamine", 
                     "p-cresol sulfate", 
                     "Tyrosine")

# Subset the original metabolite data
data_selected <- merged_sqrt_urine_int6[, selected_metabs]

# Compute Spearman correlation matrix
cor_matrix <- cor(data_selected, use = "pairwise.complete.obs", method = "spearman")

# Visualize the matrix
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, addCoef.col = "black")


# Externalizing symptoms for urine metabolites 
# Clean the environment 
rm(list = ls())
merged_sqrt_urine_ext6 <- read_excel("merged_sqrt_urine_ext6.xlsx")

X <- as.matrix(merged_sqrt_urine_ext6 %>%
                 select(where(is.numeric)) %>%
                 select(-sqrt_cbcl_ext6, -cbcl_ext6))
symptom <- "sqrt_cbcl_ext6"
y <- merged_sqrt_urine_ext6[[symptom]]  # Extract response variable for the selected symptom


# Initialize lists to store results
fit <- list()
coef <- list()
alpha <- list()
lambda <- list()
rmse <- list()

# 10-iteration loop for alpha and lambda tuning
set.seed(1)
for (i in 1:10) {
  cat("Iteration:", i, "Symptom:", symptom, "\n")
  
  # Split data into training (90%) and testing (10%) sets
  train_idx <- sample(1:nrow(X), size = 0.9 * nrow(X), replace = FALSE)
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  X_test <- X[-train_idx, ]
  y_test <- y[-train_idx]
  
  # Tune alpha and lambda using cross-validation
  cv_results <- lapply(seq(0.1, 1, by = 0.05), function(alpha_val) {
    cv_model <- cv.glmnet(X_train, y_train, alpha = alpha_val, nfolds = 10, family = "gaussian", type.measure = "mse")
    list(alpha = alpha_val, lambda = cv_model$lambda.min, cv_model = cv_model)
  })
  
  # Retrieve best alpha and lambda based on minimum cross-validation error
  best_result <- cv_results[[which.min(sapply(cv_results, function(res) min(res$cv_model$cvm)))]]
  alpha[[i]] <- best_result$alpha
  lambda[[i]] <- best_result$lambda
  
  # Fit the Elastic Net model using the best alpha and lambda
  fit[[i]] <- glmnet(X_train, y_train, alpha = alpha[[i]], lambda = lambda[[i]], family = "gaussian")
  
  # Extract non-zero coefficients
  coef[[i]] <- coef(fit[[i]])
  
  # Predict on the test set and calculate RMSE
  y_pred <- predict(fit[[i]], newx = X_test)
  rmse[[i]] <- sqrt(mean((y_test - y_pred)^2))
}

# Consolidate results into a data frame
results_extern<- data.frame(
  RMSE = unlist(rmse),
  Alpha = unlist(alpha),
  Lambda = unlist(lambda)
)


# Save results to a CSV file
#write.xlsx(results_extern, "all figures tables/urine_ext_alpha_lambda.xlsx", rowNames = FALSE)

cat("Results saved for symptom:", symptom, "\n")


best_combination_extern <- results_extern[which.min(results_extern$RMSE), ]

print(best_combination_extern)

########################################################################################
#FOR EXTERNALIZING 
# Define alpha and lambda based on the tuning results


alpha <- 0.95 
lambda <- 0.08163184



# Prepare data for Elastic Net
X <- as.matrix(merged_sqrt_urine_ext6 %>%
                 select(where(is.numeric)) %>%
                 select(-sqrt_cbcl_ext6, -cbcl_ext6))
symptom <- "sqrt_cbcl_ext6"
y <- merged_sqrt_urine_ext6[[symptom]]  # Extract response variable for the selected symptom

# Debugging: Check data
print(summary(X))
print(summary(y))
if (any(is.na(X)) || any(is.na(y))) {
  stop("Missing values detected in X or y.")
}

# Initialize a list to store coefficients from iterations
my_coef <- list()

# Set seed for reproducibility
set.seed(124876)

# Loop for coefficient selection
for (i in 1:10) {
  print(paste("Iteration:", i))
  
  # Split data into training (90%) and testing (10%) sets
  train_idx <- sample(1:nrow(X), size = 0.9 * nrow(X), replace = FALSE)
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  
  # Fit Elastic Net model using pre-selected alpha and lambda
  fit <- glmnet(X_train, y_train, family = "gaussian", alpha = alpha, lambda = lambda, intercept = FALSE)
  
  # Extract non-zero coefficients
  my_coef[[i]] <- as.matrix(coef(fit))
  my_coef[[i]] <- my_coef[[i]][my_coef[[i]] != 0, ]
}

# Combine all coefficients across iterations
final_coef <- unlist(my_coef)

# Filter coefficients that appear in more than 5 iterations
final_coef_table <- table(names(final_coef))
filtered_coef <- names(final_coef_table[final_coef_table > 6])

# Debug: Check filtered coefficients
if (length(filtered_coef) == 0) {
  stop("No coefficients met the frequency threshold (more than 5 iterations).")
}
print("Filtered coefficients:")
print(filtered_coef)

# Create a matrix for coefficients
coef_matrix <- matrix(NA, nrow = 10, ncol = length(filtered_coef))
colnames(coef_matrix) <- filtered_coef

# Populate the coefficient matrix for selected coefficients
for (i in 1:10) {
  for (j in filtered_coef) {
    if (j %in% names(my_coef[[i]])) {
      coef_matrix[i, j] <- my_coef[[i]][j]
    }
  }
}

# Calculate mean and standard deviation of selected coefficients
est <- apply(coef_matrix, 2, mean, na.rm = TRUE)
est_sd <- apply(coef_matrix, 2, sd, na.rm = TRUE)
metabolite_names <- names(est)
mean_values <- as.numeric(est)

# Calculate 95% Confidence Intervals
lower <- est - 1.96 * est_sd
upper <- est + 1.96 * est_sd

# Create a data frame for the coefficients including 95% CIs
DF_COEF_ext6 <- data.frame(
  Mean = est,
  SD = est_sd,
  Lower95CI = lower,
  Upper95CI = upper
)

# Print the resulting data frame
print(DF_COEF_ext6)
DF_COEF_ext6$Metabolite <- rownames(DF_COEF_ext6)
#write.xlsx(DF_COEF_ext6, file = "all figures tables/urine_ext_coef.xlsx", rowNames = TRUE)


#Forest plot of selected metabolites 
# Add metabolite names to the data frame
DF_COEF_ext6$Metabolite <- rownames(DF_COEF_ext6)

# Reorder metabolites by effect size (optional)
DF_COEF_ext6 <- DF_COEF_ext6 %>%
  arrange(desc(Mean))

# Create forest plot
ggplot(DF_COEF_ext6, aes(x = reorder(Metabolite, Mean), y = Mean)) +
  geom_point(color = "black") +
  geom_errorbar(aes(ymin = Lower95CI, ymax = Upper95CI), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(
    title = "Forest Plot of Urine Metabolites Associated with Externalizing Symptoms",
    x = "Metabolite",
    y = "Mean Coefficient (95% CI)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 4),
    axis.text.y = element_text(color = "black")  # y-axis labels in black
  )


#CALCULATE THE PREDICTION RMSE 
metabolite_names <- names(est)
mean_values <- as.numeric(est)

# Ensure the column names in X match the names in the coefficients
if (!all(names(est) %in% colnames(X))) {
  stop("Mismatch between coefficient names and column names in X.")
}

# Initialize a list to store RMSE values
rmse <- list()

# Set seed for reproducibility
set.seed(123)

# Loop for 10 iterations to calculate RMSE
for (i in 1:10) {
  # Create train-test split (90% train, 10% test)
  ind <- sample(1:nrow(X), size = nrow(X) * 0.9, replace = FALSE)
  train <- X[ind, ]
  test <- X[-ind, ]
  y_train <- y[ind]
  y_test <- y[-ind]
  
  # Calculate predicted scores for the test set
  pred <- rowSums(sapply(names(est), function(var) {
    if (var %in% colnames(test)) {
      test[, var] * est[var]
    } else {
      0  # If variable is missing, contribute 0
    }
  }))
  
  # Calculate RMSE for the test set
  rmse[[i]] <- sqrt(mean((y_test - pred)^2, na.rm = TRUE))
}

# Calculate the mean and 95% confidence interval of RMSE
mean_rmse <- mean(unlist(rmse))
sd_rmse <- sd(unlist(rmse))
sqrt_n <- sqrt(10)  # Number of iterations

# 95% Confidence intervals
upper <- mean_rmse + 1.96 * (sd_rmse / sqrt_n)
lower <- mean_rmse - 1.96 * (sd_rmse / sqrt_n)

# Output results
cat("Mean RMSE:", mean_rmse, "\n")
cat("95% CI of RMSE:", lower, "to", upper, "\n")

###################################################################################

# CALCULATE SCORE AND ZSCORE 
X <- as.data.frame(X)

# Define required metabolites
required_columns <-  c(
  "4-deoxythreonic acid",
  "Glycine",
  "Lysine",
  "N-acetyl neuraminic acid",
  "p-cresol sulfate",
  "Pantothenic acid",
  "Taurine",
  "Trimethylamine oxide",
  "Tyrosine"
)

# Check if all required columns exist in X
missing_columns <- setdiff(required_columns, colnames(X))
if (length(missing_columns) > 0) {
  stop(paste("ERROR: The following required columns are missing in X:", 
             paste(missing_columns, collapse = ", ")))
}

# Compute Externalizing Signature Score
X$score_externalizing <- (
  -0.1578368  * X[["4-deoxythreonic acid"]] +
    0.13536015 * X[["Glycine"]] +
    -0.07365586 * X[["Lysine"]] +
    -0.19192026 * X[["N-acetyl neuraminic acid"]] +
    -0.08695367 * X[["p-cresol sulfate"]] +
    -0.15753733 * X[["Pantothenic acid"]] +
    0.07641308 * X[["Taurine"]] +
    0.16974074 * X[["Trimethylamine oxide"]] +
    0.11098109 * X[["Tyrosine"]]
)
# Display the calculated Externalizing signature scores
print(X$score_externalizing)
summary(X$score_externalizing)

# Compute Z-score normalization for Externalizing Signature Score
X$zscore_externalizing <- (X$score_externalizing - mean(X$score_externalizing, na.rm = TRUE)) / 
  sd(X$score_externalizing, na.rm = TRUE)


# Show summary statistics of the z-score
summary(X$zscore_externalizing)

####################################################################################
#Simple pearson colleration 
correlation_result_extern <- cor(X$score_externalizing, y, method = "pearson", use = "complete.obs")

######################################################################################
set.seed(238)

# Create a dataframe to store Pearson correlation results for Externalizing
pear_externalizing <- data.frame(matrix(nrow = 1000, ncol = 2)) 
colnames(pear_externalizing) <- c("estimate", "p.value")  

# Run 1000 iterations for cross-validation
for (i in 1:1000) {
  
  # Sample 90% of the data for training
  ind <- sample(seq_len(nrow(X)), size = round(0.9 * nrow(X)), replace = FALSE)
  subset90 <- X[ind, ]
  y_subset <- y[ind]  # Ensure the correct response variable is used
  
  # Ensure 'score_externalizing' and 'y_subset' are numeric
  subset90$score_externalizing <- as.numeric(subset90$score_externalizing)
  y_subset <- as.numeric(y_subset)
  
  # Compute Pearson correlation
  cor_test <- cor.test(subset90$score_externalizing, y_subset, method = "pearson")
  
  # Store correlation coefficient and p-value
  pear_externalizing[i,1] <- cor_test$estimate  
  pear_externalizing[i,2] <- cor_test$p.value    
}

# Compute summary statistics
df_pearson_externalizing <- data.frame(
  Mean = mean(pear_externalizing$estimate, na.rm = TRUE),
  SD = sd(pear_externalizing$estimate, na.rm = TRUE)
)

# Compute 95% Confidence Interval
df_pearson_externalizing$Lower95CI <- df_pearson_externalizing$Mean - 1.96 * (df_pearson_externalizing$SD / sqrt(1000))
df_pearson_externalizing$Upper95CI <- df_pearson_externalizing$Mean + 1.96 * (df_pearson_externalizing$SD / sqrt(1000))
df_pearson_externalizing$p_value <- mean(pear_externalizing$p.value, na.rm = TRUE) 

# Print final Pearson correlation results
print(df_pearson_externalizing)

########################################################################################
############## R² for Externalizing Signature Score ####################################
# Set seed for reproducibility
set.seed(238)

# Create a dataframe to store R² results for Externalizing
squared_externalizing <- data.frame(matrix(nrow = 1000, ncol = 2))  
colnames(squared_externalizing) <- c("r.squared", "p.value")  

# Run 1000 iterations for cross-validation
for (i in 1:1000) {
  
  # Sample 90% of the data for training
  ind <- sample(seq_len(nrow(X)), size = round(0.9 * nrow(X)), replace = FALSE)
  subset90 <- X[ind, ]
  y_subset <- y[ind]
  
  # Ensure 'score_externalizing' and 'y_subset' are numeric
  subset90$score_externalizing <- as.numeric(subset90$score_externalizing)
  y_subset <- as.numeric(y_subset)
  
  # Fit a linear regression model
  model <- lm(y_subset ~ subset90$score_externalizing)
  
  # Store R-squared and p-value
  squared_externalizing[i,1] <- summary(model)$r.squared  
  squared_externalizing[i,2] <- summary(model)$coefficients[2,4]  # P-value of the coefficient
}

# Compute summary statistics
df_r2_externalizing <- data.frame(
  Mean = mean(squared_externalizing$r.squared, na.rm = TRUE),
  SD = sd(squared_externalizing$r.squared, na.rm = TRUE)
)


# Compute 95% Confidence Interval
df_r2_externalizing$Lower95CI <- df_r2_externalizing$Mean - 1.96 * (df_r2_externalizing$SD / sqrt(1000))
df_r2_externalizing$Upper95CI <- df_r2_externalizing$Mean + 1.96 * (df_r2_externalizing$SD / sqrt(1000))
df_r2_externalizing$p_value <- mean(squared_externalizing$p.value, na.rm = TRUE) 

# Print final R-squared results
print(df_r2_externalizing)

###############################################################################

# --- Round helper ---
round_df <- function(df, digits = 3) {
  as.data.frame(lapply(df, function(x) if (is.numeric(x)) round(x, digits) else x),
                stringsAsFactors = FALSE)
}

# --- Pearson full ---
pearson_full_df <- data.frame(
  Metric    = "Pearson_r_full",
  Mean      = round(as.numeric(correlation_result_extern), 3),
  SD        = NA,
  Lower95CI = NA,
  Upper95CI = NA,
  p_value   = NA
)

# --- Pearson CV summary ---
pearson_cv_df <- cbind(Metric = "Pearson_CV", round_df(df_pearson_externalizing, 3))

# --- R² CV summary ---
r2_cv_df <- cbind(Metric = "R2_CV", round_df(df_r2_externalizing, 3))

# --- RMSE summary (you already computed mean_rmse, lower, upper) ---
rmse_df <- data.frame(
  Metric    = "RMSE",
  Mean      = round(mean_rmse, 3),
  SD        = NA,
  Lower95CI = round(lower, 3),
  Upper95CI = round(upper, 3),
  p_value   = NA
)

# --- Combine all metrics ---
all_metrics <- rbind(pearson_full_df, pearson_cv_df, r2_cv_df, rmse_df)

# --- Save to Excel ---
#write_xlsx(all_metrics, "all figures tables/urine_ext_pred_metrics.xlsx")

#CORRELATION 

selected_metabs <- c("4-deoxythreonic acid",
                     "Glycine",
                     "Lysine",
                     "N-acetyl neuraminic acid",
                     "p-cresol sulfate",
                     "Pantothenic acid",
                     "Taurine",
                     "Trimethylamine oxide",
                     "Tyrosine")

# Subset your full metabolomics data
data_selected <- merged_sqrt_urine_ext6[, selected_metabs]

# Compute Spearman correlation matrix
cor_matrix <- cor(data_selected, use = "pairwise.complete.obs", method = "spearman")

# Visualize
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 70, addCoef.col = "black")


# ADHD symptoms for urine metabolites 
# Clean the environment 
rm(list = ls())

merged_sqrt_urine_adhd6 <- read_excel("merged_sqrt_urine_adhd6.xlsx")

X <- as.matrix(merged_sqrt_urine_adhd6 %>%
                 select(where(is.numeric)) %>%
                 select(-adhd_pro6, -sqrt_adhd_pro6))

symptom <- "sqrt_adhd_pro6"
y <- merged_sqrt_urine_adhd6[[symptom]]  # Extract response variable for the selected symptom

# Initialize lists to store results
fit <- list()
coef <- list()
alpha <- list()
lambda <- list()
rmse <- list()

# 10-iteration loop for alpha and lambda tuning
set.seed(1)
for (i in 1:10) {
  cat("Iteration:", i, "Symptom:", symptom, "\n")
  
  # Split data into training (90%) and testing (10%) sets
  train_idx <- sample(1:nrow(X), size = 0.9 * nrow(X), replace = FALSE)
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  X_test <- X[-train_idx, ]
  y_test <- y[-train_idx]
  
  # Tune alpha and lambda using cross-validation
  cv_results <- lapply(seq(0.1, 1, by = 0.05), function(alpha_val) {
    cv_model <- cv.glmnet(X_train, y_train, alpha = alpha_val, nfolds = 10, family = "gaussian", type.measure = "mse")
    list(alpha = alpha_val, lambda = cv_model$lambda.min, cv_model = cv_model)
  })
  
  # Retrieve best alpha and lambda based on minimum cross-validation error
  best_result <- cv_results[[which.min(sapply(cv_results, function(res) min(res$cv_model$cvm)))]]
  alpha[[i]] <- best_result$alpha
  lambda[[i]] <- best_result$lambda
  
  # Fit the Elastic Net model using the best alpha and lambda
  fit[[i]] <- glmnet(X_train, y_train, alpha = alpha[[i]], lambda = lambda[[i]], family = "gaussian")
  
  # Extract non-zero coefficients
  coef[[i]] <- coef(fit[[i]])
  
  # Predict on the test set and calculate RMSE
  y_pred <- predict(fit[[i]], newx = X_test)
  rmse[[i]] <- sqrt(mean((y_test - y_pred)^2))
}

# Consolidate results into a data frame
results_adhd <- data.frame(
  RMSE = unlist(rmse),
  Alpha = unlist(alpha),
  Lambda = unlist(lambda)
)


# Save results to a CSV file
#write.xlsx(results_adhd, "all figures tables/urine_adhd_alpha_lambda.xlsx", rowNames = FALSE)

cat("Results saved for symptom:", symptom, "\n")

best_combination_adhd <- results_adhd[which.min(results_adhd$RMSE), ]
print(best_combination_adhd)

######################################################################################
# FOR ADHD

alpha <- 0.3 
lambda <- 0.3932095


# Prepare data for Elastic Net
X <- as.matrix(merged_sqrt_urine_adhd6 %>%
                 select(where(is.numeric)) %>%
                 select(-adhd_pro6, -sqrt_adhd_pro6))
symptom <- "sqrt_adhd_pro6"
y <- merged_sqrt_urine_adhd6[[symptom]]

# Debugging: Check data
print(summary(X))
print(summary(y))
if (any(is.na(X)) || any(is.na(y))) {
  stop("Missing values detected in X or y.")
}

# Initialize a list to store coefficients from iterations
my_coef <- list()

# Set seed for reproducibility
set.seed(124876)

# Loop for coefficient selection
for (i in 1:10) {
  print(paste("Iteration:", i))
  
  # Split data into training (90%) and testing (10%) sets
  train_idx <- sample(1:nrow(X), size = 0.9 * nrow(X), replace = FALSE)
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  
  # Fit Elastic Net model using pre-selected alpha and lambda
  fit <- glmnet(X_train, y_train, family = "gaussian", alpha = alpha, lambda = lambda, intercept = FALSE)
  
  # Extract non-zero coefficients
  my_coef[[i]] <- as.matrix(coef(fit))
  my_coef[[i]] <- my_coef[[i]][my_coef[[i]] != 0, ]
}

# Combine all coefficients across iterations
final_coef <- unlist(my_coef)

# Filter coefficients that appear in more than 5 iterations
final_coef_table <- table(names(final_coef))
filtered_coef <- names(final_coef_table[final_coef_table > 6])

# Debug: Check filtered coefficients
if (length(filtered_coef) == 0) {
  stop("No coefficients met the frequency threshold (more than x iterations).")
}
print("Filtered coefficients:")
print(filtered_coef)

# Create a matrix for coefficients
coef_matrix <- matrix(NA, nrow = 10, ncol = length(filtered_coef))
colnames(coef_matrix) <- filtered_coef

# Populate the coefficient matrix for selected coefficients
for (i in 1:10) {
  for (j in filtered_coef) {
    if (j %in% names(my_coef[[i]])) {
      coef_matrix[i, j] <- my_coef[[i]][j]
    }
  }
}

# Calculate mean and standard deviation of selected coefficients
est <- apply(coef_matrix, 2, mean, na.rm = TRUE)
est_sd <- apply(coef_matrix, 2, sd, na.rm = TRUE)
metabolite_names <- names(est)
mean_values <- as.numeric(est)


# Calculate 95% Confidence Intervals
lower <- est - 1.96 * est_sd
upper <- est + 1.96 * est_sd

# Create a data frame for the coefficients including 95% CIs
DF_COEF_adhd6 <- data.frame(
  Mean = est,
  SD = est_sd,
  Lower95CI = lower,
  Upper95CI = upper
)

# Print the resulting data frame
print(DF_COEF_adhd6)

DF_COEF_adhd6$Metabolite <- rownames(DF_COEF_adhd6)
#write.xlsx(DF_COEF_adhd6, file = "all figures tables/urine_adhd_coef.xlsx", rowNames = TRUE)


# Forest plot of selected metabolites 
# Add metabolite names as a column
DF_COEF_adhd6$Metabolite <- rownames(DF_COEF_adhd6)

# Reorder for plotting (optional: by mean effect size)
DF_COEF_adhd6 <- DF_COEF_adhd6 %>%
  arrange(desc(Mean))

# Create forest plot
ggplot(DF_COEF_adhd6, aes(x = reorder(Metabolite, Mean), y = Mean)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Lower95CI, ymax = Upper95CI), width = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(
    title = "Forest Plot of Urine Metabolites Associated with ADHD",
    x = "Metabolite",
    y = "Mean Coefficient (95% CI)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0),
    axis.text.y = element_text(color = "black")  # y-axis labels in black
  )

################################################################################
#CALCULATE THE PREDICTION RMSE 

# Ensure the column names in X match the names in the coefficients
if (!all(names(est) %in% colnames(X))) {
  stop("Mismatch between coefficient names and column names in X.")
}

# Initialize a list to store RMSE values
rmse <- list()

# Set seed for reproducibility
set.seed(123)

# Loop for 10 iterations to calculate RMSE
for (i in 1:10) {
  # Create train-test split (90% train, 10% test)
  ind <- sample(1:nrow(X), size = nrow(X) * 0.9, replace = FALSE)
  train <- X[ind, ]
  test <- X[-ind, ]
  y_train <- y[ind]
  y_test <- y[-ind]
  
  # Calculate predicted scores for the test set
  pred <- rowSums(sapply(names(est), function(var) {
    if (var %in% colnames(test)) {
      test[, var] * est[var]
    } else {
      0  # If variable is missing, contribute 0
    }
  }))
  
  # Calculate RMSE for the test set
  rmse[[i]] <- sqrt(mean((y_test - pred)^2, na.rm = TRUE))
}

# Calculate the mean and 95% confidence interval of RMSE
mean_rmse <- mean(unlist(rmse))
sd_rmse <- sd(unlist(rmse))
sqrt_n <- sqrt(10)  # Number of iterations

# 95% Confidence intervals
upper <- mean_rmse + 1.96 * (sd_rmse / sqrt_n)
lower <- mean_rmse - 1.96 * (sd_rmse / sqrt_n)

# Output results
cat("Mean RMSE (ADHD):", mean_rmse, "\n")
cat("95% CI of RMSE:", lower, "to", upper, "\n")

################################################################################

#CALCULATE SCORE AND ZSCORE 

X <- as.data.frame(X)

# Define required metabolite for ADHD signature
required_columns <- c(
  "4-deoxythreonic acid",
  "Citrate",
  "Pantothenic acid"
)

# Check if all required columns exist in X
missing_columns <- setdiff(required_columns, colnames(X))
if (length(missing_columns) > 0) {
  stop(paste("ERROR: The following required columns are missing in X:", 
             paste(missing_columns, collapse = ", ")))
}

# Compute ADHD Signature Score (Only using C16:1-OH)
X$score_adhd <- (
  -0.05031965 * X[["4-deoxythreonic acid"]] +
    0.07722535 * X[["Citrate"]] +
    -0.14404295 * X[["Pantothenic acid"]]
)

# Display the calculated ADHD signature scores
print(X$score_adhd)
summary(X$score_adhd)

# Compute Z-score normalization for ADHD Signature Score
X$zscore_adhd <- (X$score_adhd - mean(X$score_adhd, na.rm = TRUE)) / 
  sd(X$score_adhd, na.rm = TRUE)

# Show summary statistics of the z-score
summary(X$zscore_adhd)

##############################################################################
#SIMPLE COLLERATION 
correlation_result_adhd <- cor(X$score_adhd, y, method = "pearson", use = "complete.obs")
############################################################################

############################################################################

set.seed(238)

# Create a dataframe to store Pearson correlation results for ADHD
pear_adhd <- data.frame(matrix(nrow = 1000, ncol = 2)) 
colnames(pear_adhd) <- c("estimate", "p.value")  

# Run 1000 iterations for cross-validation
for (i in 1:1000) {
  
  # Sample 90% of the data for training
  ind <- sample(seq_len(nrow(X)), size = round(0.9 * nrow(X)), replace = FALSE)
  subset90 <- X[ind, ]
  y_subset <- y[ind]  # Ensure the correct response variable is used
  
  # Ensure 'score_adhd' and 'y_subset' are numeric
  subset90$score_adhd <- as.numeric(subset90$score_adhd)
  y_subset <- as.numeric(y_subset)
  
  # Compute Pearson correlation
  cor_test <- cor.test(subset90$score_adhd, y_subset, method = "pearson")
  
  # Store correlation coefficient and p-value
  pear_adhd[i,1] <- cor_test$estimate  
  pear_adhd[i,2] <- cor_test$p.value    
}

# Compute summary statistics
df_pearson_adhd <- data.frame(
  Mean = mean(pear_adhd$estimate, na.rm = TRUE),
  SD = sd(pear_adhd$estimate, na.rm = TRUE)
)

# Compute 95% Confidence Interval
df_pearson_adhd$Lower95CI <- df_pearson_adhd$Mean - 1.96 * (df_pearson_adhd$SD / sqrt(1000))
df_pearson_adhd$Upper95CI <- df_pearson_adhd$Mean + 1.96 * (df_pearson_adhd$SD / sqrt(1000))
df_pearson_adhd$p_value <- mean(pear_adhd$p.value, na.rm = TRUE) 

# Print final Pearson correlation results
print(df_pearson_adhd)

###############################################################################
# R² for ADHD Signature Score ####################################
# Set seed for reproducibility

set.seed(238)

# Create a dataframe to store R² results for ADHD
squared_adhd <- data.frame(matrix(nrow = 1000, ncol = 2))  
colnames(squared_adhd) <- c("r.squared", "p.value")  

# Run 1000 iterations for cross-validation
for (i in 1:1000) {
  
  # Sample 90% of the data for training
  ind <- sample(seq_len(nrow(X)), size = round(0.9 * nrow(X)), replace = FALSE)
  subset90 <- X[ind, ]
  y_subset <- y[ind]
  
  # Ensure 'score_adhd' and 'y_subset' are numeric
  subset90$score_adhd <- as.numeric(subset90$score_adhd)
  y_subset <- as.numeric(y_subset)
  
  # Fit a linear regression model
  model <- lm(y_subset ~ subset90$score_adhd)
  
  # Store R-squared and p-value
  squared_adhd[i,1] <- summary(model)$r.squared  
  squared_adhd[i,2] <- summary(model)$coefficients[2,4]  # P-value of the coefficient
}

# Compute summary statistics
df_r2_adhd <- data.frame(
  Mean = mean(squared_adhd$r.squared, na.rm = TRUE),
  SD = sd(squared_adhd$r.squared, na.rm = TRUE)
)

# Compute 95% Confidence Interval
df_r2_adhd$Lower95CI <- df_r2_adhd$Mean - 1.96 * (df_r2_adhd$SD / sqrt(1000))
df_r2_adhd$Upper95CI <- df_r2_adhd$Mean + 1.96 * (df_r2_adhd$SD / sqrt(1000))
df_r2_adhd$p_value <- mean(squared_adhd$p.value, na.rm = TRUE) 

# Print final R² results
print(df_r2_adhd)

###############################################################################

# --- Round helper ---
round_df <- function(df, digits = 3) {
  as.data.frame(lapply(df, function(x) if (is.numeric(x)) round(x, digits) else x),
                stringsAsFactors = FALSE)
}

# --- Pearson full ---
pearson_full_df <- data.frame(
  Metric    = "Pearson_r_full",
  Mean      = round(as.numeric(correlation_result_adhd), 3),
  SD        = NA,
  Lower95CI = NA,
  Upper95CI = NA,
  p_value   = NA
)

# --- Pearson CV summary ---
pearson_cv_df <- cbind(Metric = "Pearson_CV", round_df(df_pearson_adhd, 3))

# --- R² CV summary ---
r2_cv_df <- cbind(Metric = "R2_CV", round_df(df_r2_adhd, 3))

# --- RMSE summary (you already computed mean_rmse, lower, upper) ---
rmse_df <- data.frame(
  Metric    = "RMSE",
  Mean      = round(mean_rmse, 3),
  SD        = NA,
  Lower95CI = round(lower, 3),
  Upper95CI = round(upper, 3),
  p_value   = NA
)

# --- Combine all metrics ---
all_metrics <- rbind(pearson_full_df, pearson_cv_df, r2_cv_df, rmse_df)

# --- Save to Excel ---
#write_xlsx(all_metrics, "all figures tables/urine_adhd_pred_metrics.xlsx")


#CORRELATION 


selected_metabs <- c("4-deoxythreonic acid", 
                     "Citrate", 
                     "Pantothenic acid")

# Subset the original data
data_selected <- merged_sqrt_urine_adhd6[, selected_metabs]

# Compute Spearman correlation matrix
cor_matrix <- cor(data_selected, use = "pairwise.complete.obs", method = "spearman")

# Visualize the matrix
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, addCoef.col = "black")


####################################################################################################
################################# Table 1 - Sample Characteristics #################################
####################################################################################################

# Urine datasets
merged_sqrt_urine_ext6  <- read_excel("merged_sqrt_urine_ext6.xlsx")
merged_sqrt_urine_int6  <- read_excel("merged_sqrt_urine_int6.xlsx")
merged_sqrt_urine_adhd6 <- read_excel("merged_sqrt_urine_adhd6.xlsx")

# Serum datasets
merged_sqrt_serum_ext6  <- read_excel("data/datasets/intern_extern_data/merget outcome-metab-normal/merged_sqrt_cbcl_pro_ext6.xlsx")
merged_sqrt_serum_int6  <- read_excel("data/datasets/intern_extern_data/merget outcome-metab-normal/merged_sqrt_cbcl_pro_int6.xlsx")
merged_sqrt_serum_adhd6 <- read_excel("data/datasets/intern_extern_data/merget outcome-metab-normal/merged_sqrt_cbcl_pro_adhd6.xlsx")

# Full covariates dataset
full_covariates_2 <- read_excel("Covariates dataset/primary datasets/full_covariates_2.xlsx")
#Rename childid column as helixid 
full_covariates_2 <- full_covariates_2 %>%
  rename(helixid = childid)

# ============================================================
# EXTRACT HELIX IDs FROM ALL METABOLOMIC RELATED DATASETS
# ============================================================

helix_ids_all <- c(
  merged_sqrt_urine_ext6$helixid,
  merged_sqrt_urine_int6$helixid,
  merged_sqrt_urine_adhd6$helixid,
  merged_sqrt_serum_ext6$helixid,
  merged_sqrt_serum_int6$helixid,
  merged_sqrt_serum_adhd6$helixid
)

# Keep unique IDs only
helix_ids_unique <- unique(helix_ids_all)

cat("Total unique HELIX IDs from metabolite datasets:",
    length(helix_ids_unique), "\n")


# ============================================================
# FILTER FULL COVARIATES DATASET
# ============================================================

full_covariates_filtered <- full_covariates_2 %>%
  filter(helixid %in% helix_ids_unique)

cat("Children retained in covariates dataset:",
    nrow(full_covariates_filtered), "\n")

# ============================================================
# VARIABLE SELECTION
# ============================================================
#
# From the filtered dataset, we retain:
#
# 1) Baseline characteristics at age 6
#    - Sex
#    - Baseline age
#    - Baseline symptom scores
#    - Maternal and anthropometric covariates
#
# 2) Follow-up symptom measures at ages 11 and 15

vars_keep <- c(
  # ID
  "helixid",
  
  # Baseline (age 6)
  "sex",
  "age6",
  "cbcl_pro_int6",
  "cbcl_pro_ext6",
  "adhd_pro6",
  
  # Follow-up symptoms
  "age11",
  "cbcl_pro_int11",
  "cbcl_pro_ext11",
  "adhd_pro11",
  
  "age15",
  "cbcl_pro_ext15",
  "adhd_pro15",
  "cbcl_pro_int15",
  
  # Covariates
  "mum_educ",
  "z_who_bmi_7y",
  "bmi_iotf_7y",
  "m_greek",
  "mage",
  "bmi7y"
)

# Create final analytic dataset
analytic_dataset <- full_covariates_filtered %>%
  select(all_of(vars_keep))

summary(analytic_dataset)

# One child had a missing value for age at the 6-year assessment,
# despite having ADHD symptom data recorded at that wave.
#
# After direct communication with the RHEA study team,
# the correct age at assessment was confirmed as 6.6 years.
#
# The missing age value was therefore manually completed
# using the verified value (6.6).

analytic_dataset <- analytic_dataset %>%
  mutate(age6 = ifelse(is.na(age6), 6.6, age6))

# Confirm no missing values remain
sum(is.na(analytic_dataset$age6))

#Merge bmi cuttoff categories to be 2 categories 
d1 <- analytic_dataset %>%
  mutate(
    bmi_cutoff_merge = case_when(
      bmi_iotf_7y %in% c("overweight", "obese") ~ "overweight/obese",
      bmi_iotf_7y == "underweight/normal" ~ "underweight/normal",
      TRUE ~ NA_character_  # For missing or unknown values
    )
  )

# CBCL symptom scores were categorized into Normal,
# Borderline, and Clinical ranges according to
# age- and sex-specific thresholds as defined in
# the CBCL manual.
#
# Sex-specific cutoffs were applied separately
# at ages 6, 11, and 15 years.
d2 <- d1 %>%
  mutate(
    ## --------------------
    ## INTERNALIZING SCORES
    ## --------------------
    
    # Age 6
    Diagnosis_internalizing_6 = case_when(
      sex == "Female" & cbcl_pro_int6 >= 0  & cbcl_pro_int6 <= 10 ~ "Normal",
      sex == "Female" & cbcl_pro_int6 >= 11 & cbcl_pro_int6 <= 13 ~ "Borderline",
      sex == "Female" & cbcl_pro_int6 >= 14 & cbcl_pro_int6 <= 64 ~ "Clinical",
      sex == "Male"   & cbcl_pro_int6 >= 0  & cbcl_pro_int6 <= 8  ~ "Normal",
      sex == "Male"   & cbcl_pro_int6 >= 9  & cbcl_pro_int6 <= 11 ~ "Borderline",
      sex == "Male"   & cbcl_pro_int6 >= 12 & cbcl_pro_int6 <= 64 ~ "Clinical",
      TRUE ~ NA_character_
    ),
    
    # Age 11
    Diagnosis_internalizing_11 = case_when(
      sex == "Female" & cbcl_pro_int11 >= 0  & cbcl_pro_int11 <= 10 ~ "Normal",
      sex == "Female" & cbcl_pro_int11 >= 11 & cbcl_pro_int11 <= 13 ~ "Borderline",
      sex == "Female" & cbcl_pro_int11 >= 14 & cbcl_pro_int11 <= 64 ~ "Clinical",
      sex == "Male"   & cbcl_pro_int11 >= 0  & cbcl_pro_int11 <= 8  ~ "Normal",
      sex == "Male"   & cbcl_pro_int11 >= 9  & cbcl_pro_int11 <= 11 ~ "Borderline",
      sex == "Male"   & cbcl_pro_int11 >= 12 & cbcl_pro_int11 <= 64 ~ "Clinical",
      TRUE ~ NA_character_
    ),
    
    # Age 15
    Diagnosis_internalizing_15 = case_when(
      sex == "Female" & cbcl_pro_int15 >= 0  & cbcl_pro_int15 <= 11 ~ "Normal",
      sex == "Female" & cbcl_pro_int15 >= 12 & cbcl_pro_int15 <= 14 ~ "Borderline",
      sex == "Female" & cbcl_pro_int15 >= 15 & cbcl_pro_int15 <= 64 ~ "Clinical",
      sex == "Male"   & cbcl_pro_int15 >= 0  & cbcl_pro_int15 <= 10 ~ "Normal",
      sex == "Male"   & cbcl_pro_int15 >= 11 & cbcl_pro_int15 <= 13 ~ "Borderline",
      sex == "Male"   & cbcl_pro_int15 >= 14 & cbcl_pro_int15 <= 64 ~ "Clinical",
      TRUE ~ NA_character_
    ),
    
    ## --------------------
    ## EXTERNALIZING SCORES
    ## --------------------
    
    # Age 6
    Diagnosis_externalizing_6 = case_when(
      sex == "Female" & cbcl_pro_ext6 >= 0  & cbcl_pro_ext6 <= 11 ~ "Normal",
      sex == "Female" & cbcl_pro_ext6 >= 12 & cbcl_pro_ext6 <= 14 ~ "Borderline",
      sex == "Female" & cbcl_pro_ext6 >= 15 & cbcl_pro_ext6 <= 70 ~ "Clinical",
      sex == "Male"   & cbcl_pro_ext6 >= 0  & cbcl_pro_ext6 <= 11 ~ "Normal",
      sex == "Male"   & cbcl_pro_ext6 >= 12 & cbcl_pro_ext6 <= 15 ~ "Borderline",
      sex == "Male"   & cbcl_pro_ext6 >= 16 & cbcl_pro_ext6 <= 70 ~ "Clinical",
      TRUE ~ NA_character_
    ),
    
    # Age 11
    Diagnosis_externalizing_11 = case_when(
      sex == "Female" & cbcl_pro_ext11 >= 0  & cbcl_pro_ext11 <= 11 ~ "Normal",
      sex == "Female" & cbcl_pro_ext11 >= 12 & cbcl_pro_ext11 <= 14 ~ "Borderline",
      sex == "Female" & cbcl_pro_ext11 >= 15 & cbcl_pro_ext11 <= 70 ~ "Clinical",
      sex == "Male"   & cbcl_pro_ext11 >= 0  & cbcl_pro_ext11 <= 11 ~ "Normal",
      sex == "Male"   & cbcl_pro_ext11 >= 12 & cbcl_pro_ext11 <= 15 ~ "Borderline",
      sex == "Male"   & cbcl_pro_ext11 >= 16 & cbcl_pro_ext11 <= 70 ~ "Clinical",
      TRUE ~ NA_character_
    ),
    
    # Age 15
    Diagnosis_externalizing_15 = case_when(
      sex == "Female" & cbcl_pro_ext15 >= 0  & cbcl_pro_ext15 <= 11 ~ "Normal",
      sex == "Female" & cbcl_pro_ext15 >= 12 & cbcl_pro_ext15 <= 15 ~ "Borderline",
      sex == "Female" & cbcl_pro_ext15 >= 16 & cbcl_pro_ext15 <= 70 ~ "Clinical",
      sex == "Male"   & cbcl_pro_ext15 >= 0  & cbcl_pro_ext15 <= 13 ~ "Normal",
      sex == "Male"   & cbcl_pro_ext15 >= 14 & cbcl_pro_ext15 <= 18 ~ "Borderline",
      sex == "Male"   & cbcl_pro_ext15 >= 19 & cbcl_pro_ext15 <= 70 ~ "Clinical",
      TRUE ~ NA_character_
    )
  )


#Combine the categories cuttoff 
d3 <- d2 %>%
  mutate(
    Diagnosis_internalizing_6_combined = case_when(
      Diagnosis_internalizing_6 == "Normal" ~ "Normal",
      Diagnosis_internalizing_6 %in% c("Borderline", "Clinical") ~ "Borderline/Clinical",
      TRUE ~ NA_character_
    ),
    
    Diagnosis_internalizing_11_combined = case_when(
      Diagnosis_internalizing_11 == "Normal" ~ "Normal",
      Diagnosis_internalizing_11 %in% c("Borderline", "Clinical") ~ "Borderline/Clinical",
      TRUE ~ NA_character_
    ),
    
    Diagnosis_internalizing_15_combined = case_when(
      Diagnosis_internalizing_15 == "Normal" ~ "Normal",
      Diagnosis_internalizing_15 %in% c("Borderline", "Clinical") ~ "Borderline/Clinical",
      TRUE ~ NA_character_
    ),
    
    Diagnosis_externalizing_6_combined = case_when(
      Diagnosis_externalizing_6 == "Normal" ~ "Normal",
      Diagnosis_externalizing_6 %in% c("Borderline", "Clinical") ~ "Borderline/Clinical",
      TRUE ~ NA_character_
    ),
    
    Diagnosis_externalizing_11_combined = case_when(
      Diagnosis_externalizing_11 == "Normal" ~ "Normal",
      Diagnosis_externalizing_11 %in% c("Borderline", "Clinical") ~ "Borderline/Clinical",
      TRUE ~ NA_character_
    ),
    
    Diagnosis_externalizing_15_combined = case_when(
      Diagnosis_externalizing_15 == "Normal" ~ "Normal",
      Diagnosis_externalizing_15 %in% c("Borderline", "Clinical") ~ "Borderline/Clinical",
      TRUE ~ NA_character_
    )
  )

# ============================================================
# CENTER MATERNAL AGE
# ============================================================
# Centering is performed using the sample mean.
d4 <- d3  %>% mutate(mother_age_centered = mage - mean(mage, na.rm = TRUE))

#write_xlsx(d4, "filtered_covariates_for_table1_version2.xlsx")

d_table <- read_excel("filtered_covariates_for_table1_version2.xlsx")

# Maternal characteristics
mother_cont <- c("mage")
mother_cat <- c("m_greek", "Parity", "mum_educ")

# Child age 6
child6_mean_sd <- c("age6", "z_who_bmi_7y")
child6_median <- c("cbcl_pro_int6", "cbcl_pro_ext6", "adhd_pro6", "bmi7y")
child6_cat <- c("Diagnosis_internalizing_6_combined", "Diagnosis_externalizing_6_combined", 
                "passive_smk_6y", "sex", "preterm", "bmi_cutoff_merge")

# Child age 11
child11_mean_sd <- c("age11")
child11_median <- c("cbcl_pro_int11", "cbcl_pro_ext11", "adhd_pro11")
child11_cat <- c("Diagnosis_internalizing_11_combined", "Diagnosis_externalizing_11_combined")

# Child age 15
child15_mean_sd <- c("age15")
child15_median <- c("cbcl_pro_int15", "cbcl_pro_ext15", "adhd_pro15")
child15_cat <- c("Diagnosis_internalizing_15_combined", "Diagnosis_externalizing_15_combined")

# === Helper summarization function ===

describe_variable <- function(data, var, type = c("mean_sd", "median_iqr", "categorical"), filter_var = NULL) {
  type <- match.arg(type)
  
  if (!is.null(filter_var)) {
    data <- data[!is.na(data[[filter_var]]), ]
  }
  
  x <- data[[var]]
  n_non_missing <- sum(!is.na(x))
  n_missing <- sum(is.na(x))
  
  if (type == "mean_sd") {
    mean_val <- round(mean(x, na.rm = TRUE), 2)
    sd_val <- round(sd(x, na.rm = TRUE), 2)
    desc <- paste0(mean_val, " (", sd_val, ")")
  } else if (type == "median_iqr") {
    med <- round(median(x, na.rm = TRUE), 2)
    q1 <- round(quantile(x, 0.25, na.rm = TRUE), 2)
    q3 <- round(quantile(x, 0.75, na.rm = TRUE), 2)
    desc <- paste0(med, " [", q1, ", ", q3, "]")
  } else if (type == "categorical") {
    counts <- table(x)
    percents <- prop.table(counts) * 100
    desc <- paste0(names(counts), ": ", as.integer(counts), " (", round(percents, 1), "%)", collapse = "; ")
  }
  
  tibble(
    Variable = var,
    Description = desc,
    Non_Missing = n_non_missing,
    Missing = n_missing
  )
}

# === Summarize ===

mother_summary <- bind_rows(
  lapply(mother_cont, \(v) describe_variable(d_table, v, "mean_sd")),
  lapply(mother_cat, \(v) describe_variable(d_table, v, "categorical"))
)

child6_summary <- bind_rows(
  lapply(child6_mean_sd, \(v) describe_variable(d_table, v, "mean_sd")),
  lapply(child6_median, \(v) describe_variable(d_table, v, "median_iqr")),
  lapply(child6_cat, \(v) describe_variable(d_table, v, "categorical"))
)

child11_summary <- bind_rows(
  lapply(child11_mean_sd, \(v) describe_variable(d_table, v, "mean_sd", filter_var = "age6")),
  lapply(child11_median, \(v) describe_variable(d_table, v, "median_iqr", filter_var = gsub("11", "6", v))),
  lapply(child11_cat, \(v) describe_variable(d_table, v, "categorical", filter_var = gsub("11", "6", v)))
)

child15_summary <- bind_rows(
  lapply(child15_mean_sd, \(v) describe_variable(d_table, v, "mean_sd", filter_var = "age6")),
  lapply(child15_median, \(v) describe_variable(d_table, v, "median_iqr", filter_var = gsub("15", "6", v))),
  lapply(child15_cat, \(v) describe_variable(d_table, v, "categorical", filter_var = gsub("15", "6", v)))
)

# === Final table ===

final_descriptive_table <- bind_rows(
  mutate(mother_summary, Section = "Maternal Characteristics"),
  mutate(child6_summary, Section = "Child at Age 6"),
  mutate(child11_summary, Section = "Child at Age 11"),
  mutate(child15_summary, Section = "Child at Age 15")
) %>%
  relocate(Section)

# View the final descriptive table
print(final_descriptive_table)

#write_xlsx(final_descriptive_table, path = "final_descriptive_table.xlsx")


##############################################################################
# ============================================================
# LONGITUDINAL DATA AVAILABILITY SUMMARY
# ============================================================
#
# This table summarizes the number of participants with
# available symptom data across assessment waves (ages 6, 11, 15)
# for each symptom domain (Internalizing, Externalizing, ADHD).
#
# Categories indicate whether participants contributed data:
# - Only at age 6
# - At ages 6 & 11
# - At ages 6 & 15
# - At ages 6, 11 & 15

# Flag internalizing scores for 6, 11, and 15
d2_internalizing_flags <- d_table %>%
  mutate(
    int6 = !is.na(cbcl_pro_int6),
    int11 = !is.na(cbcl_pro_int11),
    int15 = !is.na(cbcl_pro_int15)
  )

# === Internalizing Flags ===
int_flags <- d_table %>%
  mutate(
    int6 = !is.na(cbcl_pro_int6),
    int11 = !is.na(cbcl_pro_int11),
    int15 = !is.na(cbcl_pro_int15)
  ) %>%
  summarise(
    Domain = "Internalizing",
    `Only 6` = sum(int6 & !int11 & !int15),
    `6 & 11` = sum(int6 & int11 & !int15),
    `6 & 15` = sum(int6 & !int11 & int15),
    `6, 11 & 15` = sum(int6 & int11 & int15),
    `6 (any)` = sum(int6)
  )

# === Externalizing Flags ===
ext_flags <- d_table %>%
  mutate(
    ext6 = !is.na(cbcl_pro_ext6),
    ext11 = !is.na(cbcl_pro_ext11),
    ext15 = !is.na(cbcl_pro_ext15)
  ) %>%
  summarise(
    Domain = "Externalizing",
    `Only 6` = sum(ext6 & !ext11 & !ext15),
    `6 & 11` = sum(ext6 & ext11 & !ext15),
    `6 & 15` = sum(ext6 & !ext11 & ext15),
    `6, 11 & 15` = sum(ext6 & ext11 & ext15),
    `6 (any)` = sum(ext6)
  )

# === ADHD Flags ===
adhd_flags <- d_table %>%
  mutate(
    adhd6 = !is.na(adhd_pro6),
    adhd11 = !is.na(adhd_pro11),
    adhd15 = !is.na(adhd_pro15)
  ) %>%
  summarise(
    Domain = "ADHD",
    `Only 6` = sum(adhd6 & !adhd11 & !adhd15),
    `6 & 11` = sum(adhd6 & adhd11 & !adhd15),
    `6 & 15` = sum(adhd6 & !adhd11 & adhd15),
    `6, 11 & 15` = sum(adhd6 & adhd11 & adhd15),
    `6 (any)` = sum(adhd6)
  )

# === Combine all summaries ===
summary_all <- bind_rows(int_flags, ext_flags, adhd_flags)

# View the summary
print(summary_all)


# ============================================================
# BASELINE CHARACTERISTICS (AGE 6)
# ============================================================
#
# This table summarizes baseline characteristics (age 6)

child6_vars <- c("age6", "bmi7y", "sex", "bmi_cutoff_merge", "passive_smk_6y", 
                 "Diagnosis_internalizing_6_combined", "Diagnosis_externalizing_6_combined")

# Helper function to determine variable type
get_variable_type <- function(var) {
  if (var %in% c("age6", "z_who_bmi_7y")) {
    return("mean_sd")
  } else {
    return("categorical")
  }
}


# Summary for children with non-missing internalizing symptoms at age 6
internalizing6_summary <- bind_rows(
  lapply(child6_vars, function(v) {
    describe_variable(
      data = d_table,
      var = v,
      type = get_variable_type(v),
      filter_var = "cbcl_pro_int6"
    )
  })
) %>%
  mutate(Section = "Age 6 - Non-missing Internalizing")

# Summary for children with non-missing externalizing symptoms at age 6
externalizing6_summary <- bind_rows(
  lapply(child6_vars, function(v) {
    describe_variable(
      data = d_table,
      var = v,
      type = get_variable_type(v),
      filter_var = "cbcl_pro_ext6"
    )
  })
) %>%
  mutate(Section = "Age 6 - Non-missing Externalizing")

# Combine summaries
filtered_summary <- bind_rows(internalizing6_summary, externalizing6_summary) %>%
  relocate(Section)

# View the final filtered descriptive summary
print(filtered_summary)


# ============================================================
# LONGITUDINAL INTERNALIZING SAMPLE – DESCRIPTIVE TABLE
# ============================================================
#
# This table summarizes maternal and child characteristics
# for children with:
#
# 1) Non-missing internalizing symptoms at age 6, AND
# 2) At least one follow-up internalizing assessment
#    at age 11 or age 15.

# ============================================================
# CLEAN R ENVIRONMENT
# ============================================================

rm(list = ls())
gc()

# === INTERNALIZING-RELATED VARIABLE GROUPS ===
d1 <- read_excel("filtered_covariates_for_table1_version2.xlsx")

d1_int <- d1 %>%
  filter(
    !is.na(cbcl_pro_int6) &
      (!is.na(cbcl_pro_int11) | !is.na(cbcl_pro_int15))
  )

# Maternal characteristics (related to internalizing risk)
mother_cont <- c("mage")
mother_cat  <- c("mum_educ", "m_greek")

# Child age 6 (internalizing-focused)
child6_mean_sd <- c("age6", "z_who_bmi_7y")
child6_median  <- c("cbcl_pro_int6", "bmi7y")
child6_cat     <- c("Diagnosis_internalizing_6_combined",
                    "sex",
                    "bmi_cutoff_merge",
                    "passive_smk_6y")

# Child age 11 (only internalizing variables)
child11_mean_sd <- c("age11")
child11_median  <- c("cbcl_pro_int11")
child11_cat     <- c("Diagnosis_internalizing_11_combined")

# Child age 15 (only internalizing variables)
child15_mean_sd <- c("age15")
child15_median  <- c("cbcl_pro_int15")
child15_cat     <- c("Diagnosis_internalizing_15_combined")

# === Helper summarization function ===

describe_variable <- function(data, var, type = c("mean_sd", "median_iqr", "categorical"), filter_var = NULL) {
  type <- match.arg(type)
  
  if (!is.null(filter_var)) {
    data <- data[!is.na(data[[filter_var]]), ]
  }
  
  x <- data[[var]]
  n_non_missing <- sum(!is.na(x))
  n_missing     <- sum(is.na(x))
  
  if (type == "mean_sd") {
    desc <- sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
  } else if (type == "median_iqr") {
    q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
    desc <- sprintf("%.2f [%.2f, %.2f]", q[2], q[1], q[3])
  } else if (type == "categorical") {
    counts <- table(x)
    perc   <- round(prop.table(counts) * 100, 1)
    desc <- paste0(names(counts), ": ", counts, " (", perc, "%)", collapse = "; ")
  }
  
  tibble(
    Variable = var,
    Description = desc,
    Non_Missing = n_non_missing,
    Missing = n_missing
  )
}

# === GENERATE SUMMARY TABLES ===

mother_summary <- bind_rows(
  lapply(mother_cont, \(v) describe_variable(d1_int, v, "mean_sd")),
  lapply(mother_cat,  \(v) describe_variable(d1_int, v, "categorical"))
)

child6_summary <- bind_rows(
  lapply(child6_mean_sd, \(v) describe_variable(d1_int, v, "mean_sd")),
  lapply(child6_median,  \(v) describe_variable(d1_int, v, "median_iqr")),
  lapply(child6_cat,     \(v) describe_variable(d1_int, v, "categorical"))
)

child11_summary <- bind_rows(
  lapply(child11_mean_sd, \(v) describe_variable(d1_int, v, "mean_sd")),
  lapply(child11_median,  \(v) describe_variable(d1_int, v, "median_iqr")),
  lapply(child11_cat,     \(v) describe_variable(d1_int, v, "categorical"))
)

child15_summary <- bind_rows(
  lapply(child15_mean_sd, \(v) describe_variable(d1_int, v, "mean_sd")),
  lapply(child15_median,  \(v) describe_variable(d1_int, v, "median_iqr")),
  lapply(child15_cat,     \(v) describe_variable(d1_int, v, "categorical"))
)

# === FINAL COMBINED TABLE ===

final_descriptive_table <- bind_rows(
  mutate(mother_summary, Section = "Maternal Characteristics"),
  mutate(child6_summary, Section = "Child at Age 6"),
  mutate(child11_summary, Section = "Child at Age 11"),
  mutate(child15_summary, Section = "Child at Age 15")
) %>%
  relocate(Section)

# === Output the table ===
print(final_descriptive_table)

#write_xlsx(final_descriptive_table, path = "baseline_characteristics_internal_163.xlsx")

# ============================================================
# EXCLUDED LONGITUDINAL SAMPLE – DESCRIPTIVE TABLE
# ============================================================
#
# This table summarizes baseline (age 6) maternal and child
# characteristics for participants who:
#
# - Had internalizing symptom data at age 6, BUT
# - Had no internalizing follow-up data at age 11 or 15.

d1_excluded <- d1 %>%
  filter(!is.na(cbcl_pro_int6) & is.na(cbcl_pro_int11) & is.na(cbcl_pro_int15))

# === Descriptive summaries for EXCLUDED group ===
mother_summary_excl <- bind_rows(
  lapply(mother_cont, \(v) describe_variable(d1_excluded, v, "mean_sd")),
  lapply(mother_cat,  \(v) describe_variable(d1_excluded, v, "categorical"))
)

child6_summary_excl <- bind_rows(
  lapply(child6_mean_sd, \(v) describe_variable(d1_excluded, v, "mean_sd")),
  lapply(child6_median,  \(v) describe_variable(d1_excluded, v, "median_iqr")),
  lapply(child6_cat,     \(v) describe_variable(d1_excluded, v, "categorical"))
)

# === Combine into one table for excluded group ===
excluded_descriptive_table <- bind_rows(
  mutate(mother_summary_excl, Section = "Maternal Characteristics"),
  mutate(child6_summary_excl, Section = "Child at Age 6")
) %>%
  relocate(Section)

# === View the excluded table ===
library(writexl)
print(excluded_descriptive_table)
#write_xlsx(excluded_descriptive_table, path = "baseline_characteristics_internal_excluded32.xlsx")


# Follow same procedure and for the other symptoms externalizing and ADHD 
#Externalizing follow up data and excluded 

d1 <- read_excel("filtered_covariates_for_table1_version2.xlsx")

# === Define externalizing-related variables ===

# Maternal
mother_cont <- c("mage")
mother_cat  <- c("mum_educ", "m_greek")

# Child at age 6 (focus on externalizing-relevant covariates)
child6_mean_sd <- c("age6", "z_who_bmi_7y")
child6_median  <- c("cbcl_pro_ext6", "bmi7y")
child6_cat     <- c("Diagnosis_externalizing_6_combined", "sex", "bmi_cutoff_merge", "passive_smk_6y")

# Age 11
child11_mean_sd <- c("age11")
child11_median  <- c("cbcl_pro_ext11")
child11_cat     <- c("Diagnosis_externalizing_11_combined")

# Age 15
child15_mean_sd <- c("age15")
child15_median  <- c("cbcl_pro_ext15")
child15_cat     <- c("Diagnosis_externalizing_15_combined")

# === Helper summarization function ===

describe_variable <- function(data, var, type = c("mean_sd", "median_iqr", "categorical"), filter_var = NULL) {
  type <- match.arg(type)
  
  if (!is.null(filter_var)) {
    data <- data[!is.na(data[[filter_var]]), ]
  }
  
  x <- data[[var]]
  n_non_missing <- sum(!is.na(x))
  n_missing     <- sum(is.na(x))
  
  if (type == "mean_sd") {
    desc <- sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
  } else if (type == "median_iqr") {
    q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
    desc <- sprintf("%.2f [%.2f, %.2f]", q[2], q[1], q[3])
  } else if (type == "categorical") {
    counts <- table(x)
    perc   <- round(prop.table(counts) * 100, 1)
    desc <- paste0(names(counts), ": ", counts, " (", perc, "%)", collapse = "; ")
  }
  
  tibble(
    Variable = var,
    Description = desc,
    Non_Missing = n_non_missing,
    Missing = n_missing
  )
}

# === INCLUDED: Children with externalizing at 6 and at 11 or 15 ===

d1_ext_included <- d1 %>%
  filter(!is.na(cbcl_pro_ext6) & (!is.na(cbcl_pro_ext11) | !is.na(cbcl_pro_ext15)))

# Summarize for included
mother_summary_inc <- bind_rows(
  lapply(mother_cont, \(v) describe_variable(d1_ext_included, v, "mean_sd")),
  lapply(mother_cat,  \(v) describe_variable(d1_ext_included, v, "categorical"))
)

child6_summary_inc <- bind_rows(
  lapply(child6_mean_sd, \(v) describe_variable(d1_ext_included, v, "mean_sd")),
  lapply(child6_median,  \(v) describe_variable(d1_ext_included, v, "median_iqr")),
  lapply(child6_cat,     \(v) describe_variable(d1_ext_included, v, "categorical"))
)

child11_summary_inc <- bind_rows(
  lapply(child11_mean_sd, \(v) describe_variable(d1_ext_included, v, "mean_sd")),
  lapply(child11_median,  \(v) describe_variable(d1_ext_included, v, "median_iqr")),
  lapply(child11_cat,     \(v) describe_variable(d1_ext_included, v, "categorical"))
)

child15_summary_inc <- bind_rows(
  lapply(child15_mean_sd, \(v) describe_variable(d1_ext_included, v, "mean_sd")),
  lapply(child15_median,  \(v) describe_variable(d1_ext_included, v, "median_iqr")),
  lapply(child15_cat,     \(v) describe_variable(d1_ext_included, v, "categorical"))
)

included_ext_summary <- bind_rows(
  mutate(mother_summary_inc, Section = "Maternal Characteristics"),
  mutate(child6_summary_inc, Section = "Child at Age 6"),
  mutate(child11_summary_inc, Section = "Child at Age 11"),
  mutate(child15_summary_inc, Section = "Child at Age 15")
) %>%
  relocate(Section)


# === EXCLUDED: Had externalizing at 6, but missing at both 11 and 15 ===

d1_ext_excluded <- d1 %>%
  filter(!is.na(cbcl_pro_ext6) & is.na(cbcl_pro_ext11) & is.na(cbcl_pro_ext15))

mother_summary_exc <- bind_rows(
  lapply(mother_cont, \(v) describe_variable(d1_ext_excluded, v, "mean_sd")),
  lapply(mother_cat,  \(v) describe_variable(d1_ext_excluded, v, "categorical"))
)

child6_summary_exc <- bind_rows(
  lapply(child6_mean_sd, \(v) describe_variable(d1_ext_excluded, v, "mean_sd")),
  lapply(child6_median,  \(v) describe_variable(d1_ext_excluded, v, "median_iqr")),
  lapply(child6_cat,     \(v) describe_variable(d1_ext_excluded, v, "categorical"))
)

excluded_ext_summary <- bind_rows(
  mutate(mother_summary_exc, Section = "Maternal Characteristics"),
  mutate(child6_summary_exc, Section = "Child at Age 6")
) %>%
  relocate(Section)

# === View summaries ===
print(included_ext_summary)
print(excluded_ext_summary)

#write_xlsx(excluded_ext_summary, path = "baseline_characteristics_external_excluded32.xlsx")
#write_xlsx(included_ext_summary, path = "baseline_characteristics_external_included163.xlsx")



#ADHD follow up data and excluded 
###########################################################################################

d1 <- read_excel("filtered_covariates_for_table1_version2.xlsx")

# === Define ADHD-related variables ===

# Maternal
mother_cont <- c("mage")
mother_cat  <- c("mum_educ", "m_greek")

# Age 6 (ADHD relevant)
child6_mean_sd <- c("age6", "z_who_bmi_7y")
child6_median  <- c("adhd_pro6", "bmi7y")
child6_cat     <- c("sex", "passive_smk_6y", "bmi_cutoff_merge")

# Age 11
child11_mean_sd <- c("age11")
child11_median  <- c("adhd_pro11")

# Age 15
child15_mean_sd <- c("age15")
child15_median  <- c("adhd_pro15")

# === Helper summarization function ===

describe_variable <- function(data, var, type = c("mean_sd", "median_iqr", "categorical"), filter_var = NULL) {
  type <- match.arg(type)
  
  if (!is.null(filter_var)) {
    data <- data[!is.na(data[[filter_var]]), ]
  }
  
  x <- data[[var]]
  n_non_missing <- sum(!is.na(x))
  n_missing     <- sum(is.na(x))
  
  if (type == "mean_sd") {
    desc <- sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
  } else if (type == "median_iqr") {
    q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
    desc <- sprintf("%.2f [%.2f, %.2f]", q[2], q[1], q[3])
  } else if (type == "categorical") {
    counts <- table(x)
    perc   <- round(prop.table(counts) * 100, 1)
    desc <- paste0(names(counts), ": ", counts, " (", perc, "%)", collapse = "; ")
  }
  
  tibble(
    Variable = var,
    Description = desc,
    Non_Missing = n_non_missing,
    Missing = n_missing
  )
}

# === INCLUDED: ADHD at 6 and at 11 or 15 ===
d1_adhd_included <- d1 %>%
  filter(!is.na(adhd_pro6) & (!is.na(adhd_pro11) | !is.na(adhd_pro15)))

mother_summary_inc <- bind_rows(
  lapply(mother_cont, \(v) describe_variable(d1_adhd_included, v, "mean_sd")),
  lapply(mother_cat,  \(v) describe_variable(d1_adhd_included, v, "categorical"))
)

child6_summary_inc <- bind_rows(
  lapply(child6_mean_sd, \(v) describe_variable(d1_adhd_included, v, "mean_sd")),
  lapply(child6_median,  \(v) describe_variable(d1_adhd_included, v, "median_iqr")),
  lapply(child6_cat,     \(v) describe_variable(d1_adhd_included, v, "categorical"))
)

child11_summary_inc <- bind_rows(
  lapply(child11_mean_sd, \(v) describe_variable(d1_adhd_included, v, "mean_sd")),
  lapply(child11_median,  \(v) describe_variable(d1_adhd_included, v, "median_iqr"))
)

child15_summary_inc <- bind_rows(
  lapply(child15_mean_sd, \(v) describe_variable(d1_adhd_included, v, "mean_sd")),
  lapply(child15_median,  \(v) describe_variable(d1_adhd_included, v, "median_iqr"))
)

included_adhd_summary <- bind_rows(
  mutate(mother_summary_inc, Section = "Maternal Characteristics"),
  mutate(child6_summary_inc, Section = "Child at Age 6"),
  mutate(child11_summary_inc, Section = "Child at Age 11"),
  mutate(child15_summary_inc, Section = "Child at Age 15")
) %>%
  relocate(Section)


# === EXCLUDED: ADHD at 6, but missing at both 11 and 15 ===
d1_adhd_excluded <- d1 %>%
  filter(!is.na(adhd_pro6) & is.na(adhd_pro11) & is.na(adhd_pro15))

mother_summary_exc <- bind_rows(
  lapply(mother_cont, \(v) describe_variable(d1_adhd_excluded, v, "mean_sd")),
  lapply(mother_cat,  \(v) describe_variable(d1_adhd_excluded, v, "categorical"))
)

child6_summary_exc <- bind_rows(
  lapply(child6_mean_sd, \(v) describe_variable(d1_adhd_excluded, v, "mean_sd")),
  lapply(child6_median,  \(v) describe_variable(d1_adhd_excluded, v, "median_iqr")),
  lapply(child6_cat,     \(v) describe_variable(d1_adhd_excluded, v, "categorical"))
)

excluded_adhd_summary <- bind_rows(
  mutate(mother_summary_exc, Section = "Maternal Characteristics"),
  mutate(child6_summary_exc, Section = "Child at Age 6")
) %>%
  relocate(Section)

# === View summaries ===
print(included_adhd_summary)
print(excluded_adhd_summary)

#write_xlsx(excluded_adhd_summary, path = "baseline_characteristics_adhd_excluded31.xlsx")
#write_xlsx(included_adhd_summary, path = "baseline_characteristics_adhd_included162.xlsx")

# ============================================================
# SELECTION BIAS / ATTRITION ANALYSIS SUPPLEMENTARY TABLE S2
# ============================================================
#
# This table compares baseline characteristics between:
#
# 1) Participants included in the longitudinal  analysis
#    (baseline symptom data + at least one follow-up at 11 or 15), and
#
# 2) Participants excluded from the longitudinal analysis
#    (baseline symptom data only, no follow-up).
#
# Continuous variables are compared using Wilcoxon tests.
# Categorical variables are compared using Chi-square or
# Fisher’s exact tests when expected cell counts < 5.


# ADHD 
#######################################################################################
#This part of code is to assess potentioanl selection bias tables between participants included in the 
#lognitudinal analysis and with those not included

d1_adhd_included <- d1 %>%
  filter(!is.na(adhd_pro6) & (!is.na(adhd_pro11) | !is.na(adhd_pro15))) %>%
  mutate(group = "Included")

# ADHD excluded: baseline but no follow-up
d1_adhd_excluded <- d1 %>%
  filter(!is.na(adhd_pro6) & is.na(adhd_pro11) & is.na(adhd_pro15)) %>%
  mutate(group = "Excluded")

# Combine into one dataset for comparison
adhd_compare <- bind_rows(d1_adhd_included, d1_adhd_excluded)




mean_sd_vars <- c("age6", "z_who_bmi_7y", "mage")
median_iqr_vars <- c("adhd_pro6", "bmi7y")

# ---- MEAN (SD) ----
compare_mean <- lapply(mean_sd_vars, function(v) {
  means <- adhd_compare %>%
    group_by(group) %>%
    summarise(
      Mean = round(mean(.data[[v]], na.rm = TRUE), 2),
      SD = round(sd(.data[[v]], na.rm = TRUE), 2)
    )
  
  test <- wilcox.test(as.formula(paste(v, "~ group")), data = adhd_compare)
  
  tibble(
    Variable = v,
    Included = paste0(means$Mean[means$group == "Included"], " (", means$SD[means$group == "Included"], ")"),
    Excluded = paste0(means$Mean[means$group == "Excluded"], " (", means$SD[means$group == "Excluded"], ")"),
    p_value = signif(test$p.value, 3)
  )
})

# ---- MEDIAN [IQR] ----
compare_median <- lapply(median_iqr_vars, function(v) {
  medians <- adhd_compare %>%
    group_by(group) %>%
    summarise(
      Median = round(median(.data[[v]], na.rm = TRUE), 2),
      IQR = paste0("[", round(quantile(.data[[v]], 0.25, na.rm = TRUE), 2), ", ",
                   round(quantile(.data[[v]], 0.75, na.rm = TRUE), 2), "]")
    )
  
  test <- wilcox.test(as.formula(paste(v, "~ group")), data = adhd_compare)
  
  tibble(
    Variable = v,
    Included = paste0(medians$Median[medians$group == "Included"], " ", medians$IQR[medians$group == "Included"]),
    Excluded = paste0(medians$Median[medians$group == "Excluded"], " ", medians$IQR[medians$group == "Excluded"]),
    p_value = signif(test$p.value, 3)
  )
})

# Combine both
compare_cont <- bind_rows(compare_mean, compare_median)



compare_cat <- lapply(cat_vars, function(v) {
  
  tab <- table(adhd_compare[[v]], adhd_compare$group)
  
  # Compute expected values manually
  row_totals <- rowSums(tab)
  col_totals <- colSums(tab)
  total <- sum(tab)
  expected <- outer(row_totals, col_totals) / total
  
  # Choose test WITHOUT triggering the warning
  if (any(expected < 5)) {
    test <- fisher.test(tab)
  } else {
    test <- chisq.test(tab)
  }
  
  # Proportion table
  prop_tab <- prop.table(tab, 2) * 100
  
  # Loop over levels
  levels_v <- rownames(tab)
  summary_list <- sapply(levels_v, function(lvl) {
    inc_n <- if ("Included" %in% colnames(tab)) tab[lvl, "Included"] else 0
    exc_n <- if ("Excluded" %in% colnames(tab)) tab[lvl, "Excluded"] else 0
    inc_p <- if ("Included" %in% colnames(tab)) round(prop_tab[lvl, "Included"], 1) else 0
    exc_p <- if ("Excluded" %in% colnames(tab)) round(prop_tab[lvl, "Excluded"], 1) else 0
    
    paste0(lvl, ": ", inc_n, " (", inc_p, "%) / ", exc_n, " (", exc_p, "%)")
  })
  
  tibble(
    Variable = v,
    Summary = paste(summary_list, collapse = " | "),
    p_value = signif(test$p.value, 3)
  )
  
}) %>% bind_rows()


##################################################################################
#This part of code is to assess potentioanl selection bias tables between participants included in the 
#lognitudinal analysis and with those not included

#### Internalizing ########
d1_internalizing_included <- d1 %>%
  filter(!is.na(cbcl_pro_int6) & (!is.na(cbcl_pro_int11) | !is.na(cbcl_pro_int15))) %>%
  mutate(group = "Included")

d1_internalizing_excluded <- d1 %>%
  filter(!is.na(cbcl_pro_int6) & is.na(cbcl_pro_int11) & is.na(cbcl_pro_int15)) %>%
  mutate(group = "Excluded")

internalizing_compare <- bind_rows(d1_internalizing_included, d1_internalizing_excluded)

# Variables to compare
mean_sd_vars <- c("age6", "z_who_bmi_7y", "mage")
median_iqr_vars <- c("cbcl_pro_int6", "bmi7y")
cat_vars <- c("sex", "bmi_cutoff_merge", "mum_educ", "passive_smk_6y", "m_greek", "Diagnosis_internalizing_6_combined", "Diagnosis_externalizing_6_combined")

# Mean (SD)
compare_mean_int <- lapply(mean_sd_vars, function(v) {
  means <- internalizing_compare %>%
    group_by(group) %>%
    summarise(
      Mean = round(mean(.data[[v]], na.rm = TRUE), 2),
      SD = round(sd(.data[[v]], na.rm = TRUE), 2)
    )
  
  test <- wilcox.test(as.formula(paste(v, "~ group")), data = internalizing_compare)
  
  tibble(
    Variable = v,
    Included = paste0(means$Mean[means$group == "Included"], " (", means$SD[means$group == "Included"], ")"),
    Excluded = paste0(means$Mean[means$group == "Excluded"], " (", means$SD[means$group == "Excluded"], ")"),
    p_value = signif(test$p.value, 3)
  )
})

# Median [IQR]
compare_median_int <- lapply(median_iqr_vars, function(v) {
  medians <- internalizing_compare %>%
    group_by(group) %>%
    summarise(
      Median = round(median(.data[[v]], na.rm = TRUE), 2),
      IQR = paste0("[", round(quantile(.data[[v]], 0.25, na.rm = TRUE), 2), ", ",
                   round(quantile(.data[[v]], 0.75, na.rm = TRUE), 2), "]")
    )
  
  test <- wilcox.test(as.formula(paste(v, "~ group")), data = internalizing_compare)
  
  tibble(
    Variable = v,
    Included = paste0(medians$Median[medians$group == "Included"], " ", medians$IQR[medians$group == "Included"]),
    Excluded = paste0(medians$Median[medians$group == "Excluded"], " ", medians$IQR[medians$group == "Excluded"]),
    p_value = signif(test$p.value, 3)
  )
})

compare_cont_int <- bind_rows(compare_mean_int, compare_median_int)

# Categorical variables
compare_cat_int <- lapply(cat_vars, function(v) {
  tab <- table(internalizing_compare[[v]], internalizing_compare$group)
  
  expected <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  test <- if (any(expected < 5)) fisher.test(tab) else chisq.test(tab)
  
  prop_tab <- prop.table(tab, 2) * 100
  levels_v <- rownames(tab)
  
  summary_list <- sapply(levels_v, function(lvl) {
    inc_n <- if ("Included" %in% colnames(tab)) tab[lvl, "Included"] else 0
    exc_n <- if ("Excluded" %in% colnames(tab)) tab[lvl, "Excluded"] else 0
    inc_p <- if ("Included" %in% colnames(tab)) round(prop_tab[lvl, "Included"], 1) else 0
    exc_p <- if ("Excluded" %in% colnames(tab)) round(prop_tab[lvl, "Excluded"], 1) else 0
    
    paste0(lvl, ": ", inc_n, " (", inc_p, "%) / ", exc_n, " (", exc_p, "%)")
  })
  
  tibble(
    Variable = v,
    Summary = paste(summary_list, collapse = " | "),
    p_value = signif(test$p.value, 3)
  )
}) %>% bind_rows()

##############################################################################################
#This part of code is to assess potentioanl selection bias tables between participants included in the 
#lognitudinal analysis and with those not included
#Externalizing 
d1_externalizing_included <- d1 %>%
  filter(!is.na(cbcl_pro_ext6) & (!is.na(cbcl_pro_ext11) | !is.na(cbcl_pro_ext15))) %>%
  mutate(group = "Included")

d1_externalizing_excluded <- d1 %>%
  filter(!is.na(cbcl_pro_ext6) & is.na(cbcl_pro_ext11) & is.na(cbcl_pro_ext15)) %>%
  mutate(group = "Excluded")

externalizing_compare <- bind_rows(d1_externalizing_included, d1_externalizing_excluded)

# Median var changes
median_iqr_vars <- c("cbcl_pro_ext6", "bmi7y")

# Mean (SD)
compare_mean_ext <- lapply(mean_sd_vars, function(v) {
  means <- externalizing_compare %>%
    group_by(group) %>%
    summarise(
      Mean = round(mean(.data[[v]], na.rm = TRUE), 2),
      SD = round(sd(.data[[v]], na.rm = TRUE), 2)
    )
  
  test <- wilcox.test(as.formula(paste(v, "~ group")), data = externalizing_compare)
  
  tibble(
    Variable = v,
    Included = paste0(means$Mean[means$group == "Included"], " (", means$SD[means$group == "Included"], ")"),
    Excluded = paste0(means$Mean[means$group == "Excluded"], " (", means$SD[means$group == "Excluded"], ")"),
    p_value = signif(test$p.value, 3)
  )
})

# Median [IQR]
compare_median_ext <- lapply(median_iqr_vars, function(v) {
  medians <- externalizing_compare %>%
    group_by(group) %>%
    summarise(
      Median = round(median(.data[[v]], na.rm = TRUE), 2),
      IQR = paste0("[", round(quantile(.data[[v]], 0.25, na.rm = TRUE), 2), ", ",
                   round(quantile(.data[[v]], 0.75, na.rm = TRUE), 2), "]")
    )
  
  test <- wilcox.test(as.formula(paste(v, "~ group")), data = externalizing_compare)
  
  tibble(
    Variable = v,
    Included = paste0(medians$Median[medians$group == "Included"], " ", medians$IQR[medians$group == "Included"]),
    Excluded = paste0(medians$Median[medians$group == "Excluded"], " ", medians$IQR[medians$group == "Excluded"]),
    p_value = signif(test$p.value, 3)
  )
})

compare_cont_ext <- bind_rows(compare_mean_ext, compare_median_ext)

# Categorical variables
compare_cat_ext <- lapply(cat_vars, function(v) {
  tab <- table(externalizing_compare[[v]], externalizing_compare$group)
  
  expected <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  test <- if (any(expected < 5)) fisher.test(tab) else chisq.test(tab)
  
  prop_tab <- prop.table(tab, 2) * 100
  levels_v <- rownames(tab)
  
  summary_list <- sapply(levels_v, function(lvl) {
    inc_n <- if ("Included" %in% colnames(tab)) tab[lvl, "Included"] else 0
    exc_n <- if ("Excluded" %in% colnames(tab)) tab[lvl, "Excluded"] else 0
    inc_p <- if ("Included" %in% colnames(tab)) round(prop_tab[lvl, "Included"], 1) else 0
    exc_p <- if ("Excluded" %in% colnames(tab)) round(prop_tab[lvl, "Excluded"], 1) else 0
    
    paste0(lvl, ": ", inc_n, " (", inc_p, "%) / ", exc_n, " (", exc_p, "%)")
  })
  
  tibble(
    Variable = v,
    Summary = paste(summary_list, collapse = " | "),
    p_value = signif(test$p.value, 3)
  )
}) %>% bind_rows()


# ============================================================
# SUPPLEMENTARY TABLE S5
# Interpretation of Square Root–Transformed CBCL Outcomes
#
# This table presents raw means and standard deviations of CBCL
# symptom scores across assessment waves (ages 6, 11, 15).
#
# To aid interpretation of regression coefficients estimated
# on the square root–transformed scale, we back-transform
# +0.5 and +1.0 √-unit increases to the raw score metric.
############################################################################]

cbcl_vars <- c(
  "cbcl_pro_int6", "cbcl_pro_ext6",
  "cbcl_pro_int11", "cbcl_pro_ext11",
  "cbcl_pro_int15", "cbcl_pro_ext15",
  "adhd_pro6", "adhd_pro11", "adhd_pro15"
)

# 4. Define outcome labels and time points
outcomes <- c(
  "Internalizing", "Externalizing",
  "Internalizing", "Externalizing",
  "Internalizing", "Externalizing",
  "ADHD", "ADHD", "ADHD"
)

timepoints <- c(
  "6", "6",
  "11", "11",
  "15", "15",
  "6", "11", "15"
)

# 5. Calculate raw means and SDs
raw_means <- sapply(d1[cbcl_vars], mean, na.rm = TRUE)
raw_sds <- sapply(d1[cbcl_vars], sd, na.rm = TRUE)

# 6. Create base table
table_S5 <- data.frame(
  Outcome = outcomes,
  Age = timepoints,
  Raw_Mean = round(raw_means, 2),
  Raw_SD = round(raw_sds, 2),
  stringsAsFactors = FALSE
)

# 7. Calculate square root of mean
table_S5$Sqrt_Mean <- round(sqrt(table_S5$Raw_Mean), 2)

# 8. Back-transform to raw scale for +0.5 and +1.0 √ units
table_S5$Raw_at_plus_0.5 <- round((table_S5$Sqrt_Mean + 0.5)^2, 2)
table_S5$Delta_Raw_0.5 <- round(table_S5$Raw_at_plus_0.5 - table_S5$Raw_Mean, 2)

table_S5$Raw_at_plus_1.0 <- round((table_S5$Sqrt_Mean + 1.0)^2, 2)
table_S5$Delta_Raw_1.0 <- round(table_S5$Raw_at_plus_1.0 - table_S5$Raw_Mean, 2)

# 9. Create combined Mean (SD) column
table_S5$Raw_Mean_SD <- paste0(table_S5$Raw_Mean, " (", table_S5$Raw_SD, ")")

# 10. Select and reorder columns
final_table <- table_S5 %>%
  mutate(`Outcome (Domain)` = paste(Outcome)) %>%
  select(`Outcome (Domain)`, `Age (Years)` = Age,
         `Mean (SD)` = Raw_Mean_SD,
         `√-Transformed Mean` = Sqrt_Mean,
         `Raw Score Change for +0.5 √ Unit` = Delta_Raw_0.5,
         `Raw Score Change for +1.0 √ Unit` = Delta_Raw_1.0)

# 11. Reorder rows
desired_order <- c(
  "Internalizing", "Internalizing", "Internalizing",
  "Externalizing", "Externalizing", "Externalizing",
  "ADHD", "ADHD", "ADHD"
)
ages_order <- c("6", "11", "15")
final_table <- final_table %>%
  arrange(factor(`Outcome (Domain)`, levels = unique(desired_order)),
          factor(`Age (Years)`, levels = ages_order))

# 12. Print the final table
print(final_table)

# 13. Export to Excel
#write_xlsx(final_table, "Table_S5_CBCL_Sqrt_Interpretation_with_SqrtMean.xlsx")


# CREATE Z-STANDARDIZED METABOLITE SIGNATURE SCORES
# ============================================================
#
# For each biofluid (serum, urine) and each symptom domain
# (internalizing, externalizing, ADHD), we constructed
# metabolite signature scores based on the metabolites
# selected from the Elastic Net regression models.
#
# Metabolites were:
# 1) Extracted from the corresponding dataset
# 2) Standardized (z-scored)
# 3) Combined using their Elastic Net coefficients

#Internalizing serum 

# ------------------------------------------------------------
# Elastic Net Coefficients (Serum – Internalizing)
# ------------------------------------------------------------


mean_results_elastic_net_internalizing_serum <- data.frame(
  Metabolite = c("C16:1-OH", "C3", "C3-DC (C4-OH)", "Gln", "lysoPC a C17:0", "Glu"),
  Mean = c(-0.1262098, -0.08253165, -0.07723611, -0.06750363, -0.06060119, 0.07189446)
)
# ------------------------------------------------------------
# Load Serum Dataset (Imputed & Normalized)
# ------------------------------------------------------------
metabolite_int_serum  <- read_excel("data/datasets/intern_extern_data/merget outcome-metab-normal/merged_sqrt_cbcl_pro_int6.xlsx")

# Each metabolite is assigned its corresponding coefficient
# from the Elastic Net model
metabolite_coefficients_int_serum <- setNames(
  mean_results_elastic_net_internalizing_serum$Mean,
  mean_results_elastic_net_internalizing_serum$Metabolite
)

selected_metabolites_int_serum <- names(metabolite_coefficients_int_serum)[
  names(metabolite_coefficients_int_serum) %in% colnames(metabolite_int_serum)
]

if (length(selected_metabolites_int_serum) == 0) {
  stop("No matching internalizing serum metabolites found in the dataset!")
}

# Calculate the metabolomic score
# For each participant:
#   metabolite_score = Σ (metabolite_value × Elastic Net coefficient)
#
# rowwise() ensures that the weighted sum is computed per individual.

metabolite_int_serum_scored <- metabolite_int_serum %>%
  rowwise() %>%
  mutate(metabolite_score_int6 = sum(
    c_across(all_of(selected_metabolites_int_serum)) *
      metabolite_coefficients_int_serum[selected_metabolites_int_serum],
    na.rm = TRUE
  )) %>%
  ungroup() %>%
  mutate(time = 6)

# Select relevant columns
dataset_with_score_intern_serum <- metabolite_int_serum_scored %>%
  select(helixid, time, cbcl_pro_int6, sqrt_cbcl_pro_int6, metabolite_score_int6)

# Z-standardize the metabolite score
dataset_with_score_intern_serum <- dataset_with_score_intern_serum %>%
  mutate(z_metabolite_score_int6 = as.numeric(scale(metabolite_score_int6)))

# Save to Excel
#write_xlsx(dataset_with_score_intern_serum, "dataset_with_score_time6_intern_serum.xlsx")

# Remove all objects from the workspace
rm(list = ls())

# ============================================================
# EXTERNALIZING – SERUM METABOLITE SIGNATURE (AGE 6)
# ============================================================
# Define Elastic Net coefficients
mean_results_elastic_net_externalizing_serum <- data.frame(
  Metabolite = c("ADMA", "C10:2", "C16:1-OH", "Glu", "lysoPC a C17:0", "PC aa C38:4", "SDMA"),
  Mean = c(0.11755724, 0.07197259, -0.18074912, 0.13879549, -0.07380767, 0.20066734, -0.09243611)
)

# Load serum dataset (externalizing age 6)
metabolite_ext_serum <- read_excel("data/datasets/intern_extern_data/merget outcome-metab-normal/merged_sqrt_cbcl_pro_ext6.xlsx")

# Load serum dataset (externalizing age 6)
metabolite_coefficients_ext_serum <- setNames(
  mean_results_elastic_net_externalizing_serum$Mean,
  mean_results_elastic_net_externalizing_serum$Metabolite
)

# Identify metabolites present in dataset
selected_metabolites_ext_serum <- names(metabolite_coefficients_ext_serum)[
  names(metabolite_coefficients_ext_serum) %in% colnames(metabolite_ext_serum)
]

if (length(selected_metabolites_ext_serum) == 0) {
  stop("No matching externalizing serum metabolites found in the dataset!")
}

# Compute weighted metabolite score
metabolite_ext_serum_scored <- metabolite_ext_serum %>%
  rowwise() %>%
  mutate(metabolite_score_ext6 = sum(
    c_across(all_of(selected_metabolites_ext_serum)) *
      metabolite_coefficients_ext_serum[selected_metabolites_ext_serum],
    na.rm = TRUE
  )) %>%
  ungroup() %>%
  mutate(time = 6)

# Select output columns
dataset_with_score_extern_serum <- metabolite_ext_serum_scored %>%
  select(helixid, time, cbcl_ext6, sqrt_cbcl_ext6, metabolite_score_ext6)

# Standardize metabolite score (z-score)
dataset_with_score_extern_serum <- dataset_with_score_extern_serum %>%
  mutate(z_metabolite_score_ext6 = as.numeric(scale(metabolite_score_ext6)))

#write_xlsx(dataset_with_score_extern_serum, "dataset_with_score_time6_extern_serum.xlsx")

# ============================================================
# ADHD – SERUM METABOLITE SIGNATURE (AGE 6)
# ============================================================
## Remove all objects from the workspace
rm(list = ls())

# Elastic Net coefficient(s) for ADHD (serum)
mean_results_elastic_net_adhd_serum <- data.frame(
  Metabolite = c("C16:1-OH"),
  Mean = c(-0.04002156)
)

# Load serum dataset (ADHD age 6)
metabolite_adhd_serum <- read_excel("data/datasets/intern_extern_data/merget outcome-metab-normal/merged_sqrt_cbcl_pro_adhd6.xlsx")

# Create named coefficient vector
metabolite_coefficients_adhd_serum <- setNames(
  mean_results_elastic_net_adhd_serum$Mean,
  mean_results_elastic_net_adhd_serum$Metabolite
)

# Identify metabolites available in dataset
selected_metabolites_adhd_serum <- names(metabolite_coefficients_adhd_serum)[
  names(metabolite_coefficients_adhd_serum) %in% colnames(metabolite_adhd_serum)
]

if (length(selected_metabolites_adhd_serum) == 0) {
  stop("No matching ADHD serum metabolites found in the dataset!")
}

# Compute weighted metabolite score
metabolite_adhd_serum_scored <- metabolite_adhd_serum %>%
  rowwise() %>%
  mutate(metabolite_score_adhd6 = sum(
    c_across(all_of(selected_metabolites_adhd_serum)) *
      metabolite_coefficients_adhd_serum[selected_metabolites_adhd_serum],
    na.rm = TRUE
  )) %>%
  ungroup() %>%
  mutate(time = 6)

# Retain relevant variables
dataset_with_score_adhd_serum <- metabolite_adhd_serum_scored %>%
  select(helixid, time, adhd_pro6, sqrt_adhd_pro6, metabolite_score_adhd6)

# Standardize metabolite score
dataset_with_score_adhd_serum <- dataset_with_score_adhd_serum %>%
  mutate(z_metabolite_score_adhd6 = as.numeric(scale(metabolite_score_adhd6)))


# Optional: Save to Excel
#write_xlsx(dataset_with_score_adhd_serum, "dataset_with_score_time6_adhd_serum.xlsx")

# ============================================================
# URINE METABOLITE SIGNATURES
# ============================================================
## Remove all objects from the workspace
rm(list = ls())

# Elastic Net coefficients for internalizing (urine)
mean_results_elastic_net_internalizing <- data.frame(
  Metabolite = c("3-hydroxyisovalerate", "5-oxoproline", "Dimethylamine", 
                 "Glutamine", "p-cresol sulfate", "Tyrosine"),
  Mean = c(-0.03230146, -0.04665104, -0.04922661, 
           0.06439756, -0.06812796, 0.12402503)
)

# Load urine dataset (internalizing age 6)
# Includes raw and square-root transformed outcomes
# and inverse-normal transformed metabolite concentrations
metabolite_int_urine<- read_excel("merged_sqrt_urine_int6.xlsx")


# Create named coefficient vector
metabolite_coefficients_int <- setNames(mean_results_elastic_net_internalizing$Mean, 
                                        mean_results_elastic_net_internalizing$Metabolite)

# Identify metabolites present in dataset
selected_metabolites_int <- names(metabolite_coefficients_int)[names(metabolite_coefficients_int) %in% colnames(metabolite_int_urine)]

if (length(selected_metabolites_int) == 0) {
  stop("No matching internalizing metabolites found in the dataset!")
}


# Compute weighted metabolite score
metabolite_int_urine_scored <- metabolite_int_urine %>%
  rowwise() %>%
  mutate(metabolite_score_int6 = sum(c_across(all_of(selected_metabolites_int)) * 
                                       metabolite_coefficients_int[selected_metabolites_int], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(time = 6)

# Retain relevant variables
dataset_with_score_intern_urine <- metabolite_int_urine_scored %>%
  select(helixid, time, cbcl_pro_int6, sqrt_cbcl_pro_int6, metabolite_score_int6)

# Standardize metabolite score
dataset_with_score_intern_urine <- dataset_with_score_intern_urine %>%
  mutate(z_metabolite_score_int6 = as.numeric(scale(metabolite_score_int6)))

#write_xlsx(dataset_with_score_intern_urine, "dataset_with_score_time6_intern_urine.xlsx")

## Remove all objects from the workspace
rm(list = ls())
# ============================================================
# EXTERNALIZING – URINE METABOLITE SIGNATURE (AGE 6)
# ============================================================

# Elastic Net coefficients for externalizing (urine)
mean_results_elastic_net_externalizing <- data.frame(
  Metabolite = c("4-deoxythreonic acid", "Glycine", "Lysine", 
                 "N-acetyl neuraminic acid", "p-cresol sulfate", 
                 "Pantothenic acid", "Taurine", "Trimethylamine oxide", "Tyrosine"),
  Mean = c(-0.1578368, 0.13536015, -0.07365586, 
           -0.19192026, -0.08695367, -0.15753733, 
           0.07641308, 0.16974074, 0.11098109)
)

# Load urine dataset (externalizing age 6)
# Includes raw and square-root transformed outcomes
# and inverse-normal transformed metabolite concentrations
metabolite_ext_urine <- read_excel("merged_sqrt_urine_ext6.xlsx")

# Create named coefficient vector
metabolite_coefficients_ext <- setNames(mean_results_elastic_net_externalizing$Mean, 
                                        mean_results_elastic_net_externalizing$Metabolite)

# Identify metabolites present in dataset
selected_metabolites_ext <- names(metabolite_coefficients_ext)[names(metabolite_coefficients_ext) %in% colnames(metabolite_ext_urine)]

if (length(selected_metabolites_ext) == 0) {
  stop("No matching externalizing metabolites found in the dataset!")
}

# Compute weighted metabolite score
metabolite_ext_urine_scored <- metabolite_ext_urine %>%
  rowwise() %>%
  mutate(metabolite_score_ext6 = sum(
    c_across(all_of(selected_metabolites_ext)) * 
      metabolite_coefficients_ext[selected_metabolites_ext], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(time = 6)

# Retain relevant variables
dataset_with_score_extern_urine <- metabolite_ext_urine_scored %>%
  select(helixid, time, cbcl_ext6, sqrt_cbcl_ext6, metabolite_score_ext6)

# Standardize metabolite score
dataset_with_score_extern_urine <- dataset_with_score_extern_urine %>%
  mutate(z_metabolite_score_ext6 = as.numeric(scale(metabolite_score_ext6)))

#write_xlsx(dataset_with_score_extern_urine, "dataset_with_score_time6_extern_urine.xlsx")

## Remove all objects from the workspace
rm(list = ls())

# ============================================================
# ADHD – URINE METABOLITE SIGNATURE (AGE 6)
# ============================================================

# Elastic Net coefficients for ADHD (urine)
mean_results_elastic_net_adhd <- data.frame(
  Metabolite = c("4-deoxythreonic acid", "Citrate", "Pantothenic acid"),
  Mean = c(-0.05031965, 0.07722535, -0.14404295)
)

# Load urine dataset (ADHD age 6)
# Includes raw and square-root transformed outcomes
# and inverse-normal transformed metabolite concentrations
metabolite_adhd_urine <- read_excel("merged_sqrt_urine_adhd6.xlsx")

# Create named coefficient vector
metabolite_coefficients_adhd <- setNames(mean_results_elastic_net_adhd$Mean, 
                                         mean_results_elastic_net_adhd$Metabolite)

# Identify metabolites present in dataset
selected_metabolites_adhd <- names(metabolite_coefficients_adhd)[names(metabolite_coefficients_adhd) %in% colnames(metabolite_adhd_urine)]

if (length(selected_metabolites_adhd) == 0) {
  stop("No matching ADHD metabolites found in the dataset!")
}

# Compute weighted metabolite score
metabolite_adhd_urine_scored <- metabolite_adhd_urine %>%
  rowwise() %>%
  mutate(metabolite_score_adhd6 = sum(
    c_across(all_of(selected_metabolites_adhd)) * 
      metabolite_coefficients_adhd[selected_metabolites_adhd],
    na.rm = TRUE
  )) %>%
  ungroup() %>%
  mutate(time = 6)

# Retain relevant variables
dataset_with_score_adhd_urine <- metabolite_adhd_urine_scored %>%
  select(helixid, time, adhd_pro6, sqrt_adhd_pro6, metabolite_score_adhd6)

# Standardize metabolite score
dataset_with_score_adhd_urine <- dataset_with_score_adhd_urine %>%
  mutate(z_metabolite_score_adhd6 = as.numeric(scale(metabolite_score_adhd6)))

#write_xlsx(dataset_with_score_adhd_urine, "dataset_with_score_time6_adhd_urine.xlsx")

# ============================================================
# DATA PREPARATION FOR MULTIPLE IMPUTATION (MICE)
# ============================================================
# Steps performed:
# - Load baseline covariate dataset
# - Merge serum and urine metabolite signature scores
# - Apply square-root transformation to follow-up outcomes
# - Recode categorical variables as factors

# Remove all objects from the workspace
rm(list = ls())

# ------------------------------------------------------------
# Load baseline covariate dataset
# ------------------------------------------------------------
full <- read_excel("filtered_covariates_for_table1_version2.xlsx")

#Read the datasets with the z_metabolite score for each symtpom and each biofluid 
#Read the dataset with metabolomic score from serum 
serum_adhd <- read_excel("dataset_with_score_time6_adhd_serum.xlsx")
serum_extern <- read_excel("dataset_with_score_time6_extern_serum.xlsx")
serum_intern <- read_excel("dataset_with_score_time6_intern_serum.xlsx")

#Rename the column of z_metabolite_score to represent the serum biofluid 
serum_adhd <- serum_adhd %>% rename(serum_z_metabolite_score_adhd6 = z_metabolite_score_adhd6)
serum_extern <- serum_extern %>% rename(serum_z_metabolite_score_ext6 = z_metabolite_score_ext6)
serum_intern <- serum_intern %>%rename(serum_z_metabolite_score_int6 = z_metabolite_score_int6)


#Read the dataset with metabolomic score from urine 
urine_adhd <- read_excel("dataset_with_score_time6_adhd_urine.xlsx")
urine_extern <- read_excel("dataset_with_score_time6_extern_urine.xlsx")
urine_intern <- read_excel("dataset_with_score_time6_intern_urine.xlsx")

#Rename the column of z_metabolite_score to represent the urine biofluid 
urine_adhd <- urine_adhd %>% rename(urine_z_metabolite_score_adhd6 = z_metabolite_score_adhd6)
urine_extern <- urine_extern %>% rename(urine_z_metabolite_score_ext6 = z_metabolite_score_ext6)
urine_intern <- urine_intern %>%rename(urine_z_metabolite_score_int6 = z_metabolite_score_int6)

#Merge the dataset with the covariates with the z_metabolite scores 
full_new <- full %>%
  left_join(urine_adhd %>% select(helixid, urine_z_metabolite_score_adhd6, sqrt_adhd_pro6),
            by = "helixid") %>%
  left_join(urine_extern %>% select(helixid, urine_z_metabolite_score_ext6, sqrt_cbcl_ext6),
            by = "helixid") %>%
  left_join(urine_intern %>% select(helixid, urine_z_metabolite_score_int6, sqrt_cbcl_pro_int6),
            by = "helixid") %>%
  left_join(serum_adhd %>% select(helixid, serum_z_metabolite_score_adhd6),
            by = "helixid") %>%
  left_join(serum_extern %>% select(helixid, serum_z_metabolite_score_ext6),
            by = "helixid") %>%
  left_join(serum_intern %>% select(helixid, serum_z_metabolite_score_int6),
            by = "helixid")


#We need to square root transform the outcome variables
full_new_ <- full_new %>%
  mutate(
    sqrt_cbcl_ext11 = sqrt(cbcl_pro_ext11),
    sqrt_cbcl_ext15 = sqrt(cbcl_pro_ext15),
    sqrt_cbcl_int11 = sqrt(cbcl_pro_int11),
    sqrt_cbcl_int15 = sqrt(cbcl_pro_int15),
    sqrt_cbcl_adhd11 = sqrt(adhd_pro11),
    sqrt_cbcl_adhd15 = sqrt(adhd_pro15)
  ) %>%
  select(-gender)

#Create factors 
full_new_ <- full_new_ %>%
  mutate(
    sex = factor(sex),
    m_greek = factor(m_greek),
    passive_smk_6y = factor(passive_smk_6y),
    bmi_cutoff_merge = factor(bmi_cutoff_merge)
  )

full_new_$mum_educ <- ordered(full_new_$mum_educ,
                              levels = c("Low", "Medium", "High"))


str(full_new_[, c("sex","mum_educ","m_greek","passive_smk_6y","bmi_cutoff_merge")])

#write_xlsx(full_new_, "filtered_covariates_plus_z_metabolites_scores.xlsx")

# ============================================================
# MULTIPLE IMPUTATION USING MICE (WIDE FORMAT)
# ============================================================
# ============================================================
# LOAD DATA FOR IMPUTATION
# ============================================================
d1 <- read_excel("filtered_covariates_plus_z_metabolites_scores.xlsx")
summary(d1)

# Create wide-format dataset for imputation
wideformat <- d1 %>%
  select(
    helixid,
    serum_z_metabolite_score_ext6,
    serum_z_metabolite_score_int6,
    serum_z_metabolite_score_adhd6,
    urine_z_metabolite_score_ext6,
    urine_z_metabolite_score_int6,
    urine_z_metabolite_score_adhd6,
    sex, 
    age6,
    sqrt_cbcl_ext6,
    sqrt_adhd_pro6,
    sqrt_cbcl_pro_int6,
    age11,
    sqrt_cbcl_ext11,
    sqrt_cbcl_int11,
    sqrt_cbcl_adhd11,
    age15,
    sqrt_cbcl_ext15,
    sqrt_cbcl_int15,
    sqrt_cbcl_adhd15,
    passive_smk_6y, z_who_bmi_7y, mum_educ, mother_age_centered, m_greek, bmi_cutoff_merge,
    Diagnosis_internalizing_6_combined, Diagnosis_internalizing_11_combined, Diagnosis_internalizing_15_combined,
    Diagnosis_externalizing_6_combined, Diagnosis_externalizing_11_combined, Diagnosis_externalizing_15_combined
  )

wideformat <- wideformat %>%
  mutate(
    sex = factor(sex),
    passive_smk_6y = factor(passive_smk_6y),
    mum_educ = factor(mum_educ), 
    bmi_cutoff_merge = factor(bmi_cutoff_merge)
  )


#Create long format 
#Long format
long_age6 <- wideformat %>%
  transmute(
    helixid,
    time = 6,
    age = age6,
    sqrt_cbcl_ext = sqrt_cbcl_ext6,
    sqrt_cbcl_int = sqrt_cbcl_pro_int6,
    sqrt_cbcl_adhd = sqrt_adhd_pro6,
    serum_z_metabolite_score_ext6,
    serum_z_metabolite_score_adhd6,
    serum_z_metabolite_score_int6,
    urine_z_metabolite_score_ext6,
    urine_z_metabolite_score_int6,
    urine_z_metabolite_score_adhd6,
    sex, m_greek, passive_smk_6y, 
    z_who_bmi_7y, 
    mum_educ, mother_age_centered, bmi_cutoff_merge
  )

# Age 11
long_age11<- wideformat %>%
  transmute(
    helixid,
    time = 11,
    age = age11,
    sqrt_cbcl_ext = sqrt_cbcl_ext11,
    sqrt_cbcl_int = sqrt_cbcl_int11,
    sqrt_cbcl_adhd = sqrt_cbcl_adhd11,
    serum_z_metabolite_score_ext6,
    serum_z_metabolite_score_adhd6,
    serum_z_metabolite_score_int6,
    urine_z_metabolite_score_ext6,
    urine_z_metabolite_score_int6,
    urine_z_metabolite_score_adhd6,
    sex, m_greek, passive_smk_6y, 
    z_who_bmi_7y,
    mum_educ, mother_age_centered, bmi_cutoff_merge
  )

# Age 15
long_age15 <- wideformat %>%
  transmute(
    helixid,
    time = 15,
    age = age15,
    sqrt_cbcl_ext = sqrt_cbcl_ext15,
    sqrt_cbcl_int = sqrt_cbcl_int15,
    sqrt_cbcl_adhd = sqrt_cbcl_adhd15,
    serum_z_metabolite_score_ext6,
    serum_z_metabolite_score_adhd6,
    serum_z_metabolite_score_int6,
    urine_z_metabolite_score_ext6,
    urine_z_metabolite_score_int6,
    urine_z_metabolite_score_adhd6,
    sex, m_greek, passive_smk_6y, 
    z_who_bmi_7y,
    mum_educ, mother_age_centered, bmi_cutoff_merge
  )

# Combine all three timepoints
long_all <- bind_rows(long_age6, long_age11, long_age15) %>%
  arrange(helixid, time)

long_all <- long_all %>%
  mutate(
    time = factor(time),
    sex = factor(sex),
    passive_smk_6y = factor(passive_smk_6y,
                            levels = c("No", "Yes")),
    mum_educ = factor(mum_educ,
                      levels = c("Low", "Medium", "High"),
                      ordered = TRUE),
    bmi_cutoff_merge = factor(
      case_when(
        bmi_cutoff_merge %in% c("underweight/normal") ~ "underweight/normal",
        bmi_cutoff_merge %in% c("overweight/obese") ~ "overweight/obese",
        TRUE ~ NA_character_
      ),
      levels = c("underweight/normal", "overweight/obese")
    )
  )

#write_xlsx(long_all, "filtered_covariates_all_long_format.xlsx")



# Indicators for missing baseline metabolite scores

wideformat$miss_serum_ext6  <- is.na(wideformat$serum_z_metabolite_score_ext6)
wideformat$miss_serum_int6  <- is.na(wideformat$serum_z_metabolite_score_int6)
wideformat$miss_serum_adhd6 <- is.na(wideformat$serum_z_metabolite_score_adhd6)

wideformat$miss_urine_ext6  <- is.na(wideformat$urine_z_metabolite_score_ext6)
wideformat$miss_urine_int6  <- is.na(wideformat$urine_z_metabolite_score_int6)
wideformat$miss_urine_adhd6 <- is.na(wideformat$urine_z_metabolite_score_adhd6)


#Start setup imputation run
imp0 <- mice(wideformat, maxit = 0, defaultMethod = c("norm", "logreg", "polyreg", "polr"))
imp0$loggedEvents
meth <- imp0$method
pred <- imp0$predictorMatrix



#Don't impute the helixid 
meth["helixid"] <- ""
pred[, "helixid"] <- 0
pred["helixid", ] <- 0

#Don't impute the indicator 
indicators <- c("miss_serum_ext6",
                "miss_serum_int6",
                "miss_serum_adhd6",
                "miss_urine_ext6",
                "miss_urine_int6",
                "miss_urine_adhd6",
                "Diagnosis_internalizing_6_combined",
                "Diagnosis_internalizing_11_combined",
                "Diagnosis_internalizing_15_combined",
                "Diagnosis_externalizing_6_combined",
                "Diagnosis_externalizing_11_combined",
                "Diagnosis_externalizing_15_combined")

meth[indicators] <- ""
pred[, indicators] <- 0
pred[indicators, ] <- 0


vars_to_impute <- c(
  "sqrt_cbcl_ext6",
  "sqrt_adhd_pro6",
  "sqrt_cbcl_pro_int6",
  "age11",
  "sqrt_cbcl_ext11",
  "sqrt_cbcl_int11",
  "sqrt_cbcl_adhd11",
  "age15",
  "sqrt_cbcl_ext15",
  "sqrt_cbcl_int15",
  "sqrt_cbcl_adhd15",
  "serum_z_metabolite_score_ext6",
  "serum_z_metabolite_score_int6",
  "serum_z_metabolite_score_adhd6",
  "urine_z_metabolite_score_ext6",
  "urine_z_metabolite_score_int6",
  "urine_z_metabolite_score_adhd6"
)

meth[!names(meth) %in% vars_to_impute] <- ""

# Constrain imputed ages to plausible ranges
meth["age11"] <- "~ I(runif(sum(!ry), 10, 11.5))"
meth["age15"] <- "~ I(runif(sum(!ry), 14, 15.5))"


predictors_to_use <- c(
  "serum_z_metabolite_score_ext6",
  "serum_z_metabolite_score_int6",
  "serum_z_metabolite_score_adhd6",
  "urine_z_metabolite_score_ext6",
  "urine_z_metabolite_score_int6",
  "urine_z_metabolite_score_adhd6",
  "mother_age_centered", "z_who_bmi_7y", "sex", "passive_smk_6y", "mum_educ", "m_greek")

pred[vars_to_impute, predictors_to_use] <- 1

#Visit sequence
visSeq <- imp0$visitSequence

# ============================================================
# RUN MULTIPLE IMPUTATION
# ============================================================
imp.test <- mice(wideformat, method = meth, predictorMatrix = pred, visitSequence = visSeq, maxit = 30, m = 30, printFlag = TRUE, seed = 2022) #imputing long df

imp.test$loggedEvents

# Check one of the datasets 
d5 <- complete(imp.test, 5)

#Force symptoms to imputed from values 0-100 so they can be positive. 
post <- imp.test$post

symptoms <- c("sqrt_cbcl_ext6",
              "sqrt_cbcl_ext11",
              "sqrt_cbcl_ext15",
              "sqrt_cbcl_pro_int6",
              "sqrt_cbcl_int11",
              "sqrt_cbcl_int15",
              "sqrt_adhd_pro6",
              "sqrt_cbcl_adhd11",
              "sqrt_cbcl_adhd15")

post[symptoms] <- "imp[[j]][, i] <- pmax(imp[[j]][, i], 0)"

imp2 <- update(imp.test, post = post, maxit = 30, seed = 123)


#Use the indicator to set again missing the baseline metabolite score and the symptoms 
#that are aligned to this metabolite score 
completed_list <- complete(imp2, "all")

for (i in 1:length(completed_list)) {
  
  d <- completed_list[[i]]
  
  # =========================
  # Reset metabolites
  # =========================
  
  # Serum
  d$serum_z_metabolite_score_ext6[d$miss_serum_ext6 == TRUE] <- NA
  d$serum_z_metabolite_score_int6[d$miss_serum_int6 == TRUE] <- NA
  d$serum_z_metabolite_score_adhd6[d$miss_serum_adhd6 == TRUE] <- NA
  
  # Urine
  d$urine_z_metabolite_score_ext6[d$miss_urine_ext6 == TRUE] <- NA
  d$urine_z_metabolite_score_int6[d$miss_urine_int6 == TRUE] <- NA
  d$urine_z_metabolite_score_adhd6[d$miss_urine_adhd6 == TRUE] <- NA
  
  
  # =========================
  # Align symptoms to exposure
  # =========================
  
  # Externalizing
  d[d$miss_serum_ext6 == TRUE | d$miss_urine_ext6 == TRUE,
    c("sqrt_cbcl_ext6",
      "sqrt_cbcl_ext11",
      "sqrt_cbcl_ext15"
    )] <- NA
  
  # Internalizing
  d[d$miss_serum_int6 == TRUE | d$miss_urine_int6 == TRUE,
    c("sqrt_cbcl_pro_int6",
      "sqrt_cbcl_int11",
      "sqrt_cbcl_int15"
    )] <- NA
  
  # ADHD
  d[d$miss_serum_adhd6 == TRUE | d$miss_urine_adhd6 == TRUE,
    c("sqrt_adhd_pro6",
      "sqrt_cbcl_adhd11",
      "sqrt_cbcl_adhd15"
    )] <- NA
  
  
  completed_list[[i]] <- d
}


#Veryfy the number of the full data for each symtpom
colSums(!is.na(d[, c(
  "sqrt_cbcl_ext6",
  "sqrt_cbcl_ext11",
  "sqrt_cbcl_ext15",
  "sqrt_cbcl_pro_int6",
  "sqrt_cbcl_int11",
  "sqrt_cbcl_int15",
  "sqrt_adhd_pro6",
  "sqrt_cbcl_adhd11",
  "sqrt_cbcl_adhd15"
)]))

summary(d)



#Plot diagnostics-check convergence 
plot(imp2)
densityplot(
  imp2,
  ~ sqrt_cbcl_int11 +
    sqrt_cbcl_int15 +
    sqrt_cbcl_ext11 +
    sqrt_cbcl_ext15 +
    sqrt_cbcl_adhd11+
    sqrt_cbcl_adhd15+
    age11 +
    age15,
  layout = c(3,2)
)
#savehistory("mice_imputation_code.txt")
load("imputed_datasets_wide.RData")

#Convert from wide format to long format the imputed datasets
wide_to_long <- function(completed_list) {
  
  long_age6 <- wideformat %>%
    transmute(
      helixid,
      time = 6,
      age = age6,
      sqrt_cbcl_ext  = sqrt_cbcl_ext6,
      sqrt_cbcl_int  = sqrt_cbcl_pro_int6,
      sqrt_cbcl_adhd = sqrt_adhd_pro6,
      serum_z_metabolite_score_ext6,
      serum_z_metabolite_score_int6,
      serum_z_metabolite_score_adhd6,
      urine_z_metabolite_score_ext6,
      urine_z_metabolite_score_int6,
      urine_z_metabolite_score_adhd6,
      sex, m_greek, passive_smk_6y,
      z_who_bmi_7y, bmi_cutoff_merge,
      mum_educ, mother_age_centered
    )
  
  long_age11 <- wideformat %>%
    transmute(
      helixid,
      time = 11,
      age = age11,
      sqrt_cbcl_ext  = sqrt_cbcl_ext11,
      sqrt_cbcl_int  = sqrt_cbcl_int11,
      sqrt_cbcl_adhd = sqrt_cbcl_adhd11,
      serum_z_metabolite_score_ext6,
      serum_z_metabolite_score_int6,
      serum_z_metabolite_score_adhd6,
      urine_z_metabolite_score_ext6,
      urine_z_metabolite_score_int6,
      urine_z_metabolite_score_adhd6,
      sex, m_greek, passive_smk_6y,
      z_who_bmi_7y, bmi_cutoff_merge,
      mum_educ, mother_age_centered
    )
  
  long_age15 <- wideformat %>%
    transmute(
      helixid,
      time = 15,
      age = age15,
      sqrt_cbcl_ext  = sqrt_cbcl_ext15,
      sqrt_cbcl_int  = sqrt_cbcl_int15,
      sqrt_cbcl_adhd = sqrt_cbcl_adhd15,
      serum_z_metabolite_score_ext6,
      serum_z_metabolite_score_int6,
      serum_z_metabolite_score_adhd6,
      urine_z_metabolite_score_ext6,
      urine_z_metabolite_score_int6,
      urine_z_metabolite_score_adhd6,
      sex, m_greek, passive_smk_6y,
      z_who_bmi_7y, bmi_cutoff_merge,
      mum_educ, mother_age_centered
    )
  
  bind_rows(long_age6, long_age11, long_age15) %>%
    arrange(helixid, time) %>%
    mutate(
      
      # Time factor
      time = factor(time),
      
      # Safe factor conversion
      sex = factor(sex, levels = c("Male", "Female")),
      
      m_greek = factor(m_greek, levels = c("Greek", "Other")),
      
      passive_smk_6y = factor(
        passive_smk_6y,
        levels = c("No", "Yes")
      ),
      
      mum_educ = factor(
        mum_educ,
        levels = c("Low", "Medium", "High"),
        ordered = TRUE
      ),
      # Binary BMI grouping
      bmi_cutoff_merge = factor(
        case_when(
          bmi_cutoff_merge %in% c("underweight/normal") ~ "underweight/normal",
          bmi_cutoff_merge %in% c("overweight/obese") ~ "overweight/obese",
          TRUE ~ NA_character_
        ),
        levels = c("underweight/normal", "overweight/obese")
      )
    )
}

long_list <- lapply(completed_list, wide_to_long)


#Check
long5 <- long_list[[5]]
head(long5)
summary(long5)
str(long5)


#save(long_list, file = "imputed_datasets_long.RData")
#save(completed_list, file = "imputed_datasets_wide.RData")


# ==============================================================
# Linear Mixed Models with Imputed Data
# ==============================================================
# Each entry links a behavioral outcome with its corresponding
# metabolite exposure variable (serum or urine)
model_specs <- list(
  serum_ext  = list(outcome = "sqrt_cbcl_ext",  exposure = "serum_z_metabolite_score_ext6"),
  serum_int  = list(outcome = "sqrt_cbcl_int",  exposure = "serum_z_metabolite_score_int6"),
  serum_adhd = list(outcome = "sqrt_cbcl_adhd", exposure = "serum_z_metabolite_score_adhd6"),
  urine_ext  = list(outcome = "sqrt_cbcl_ext",  exposure = "urine_z_metabolite_score_ext6"),
  urine_int  = list(outcome = "sqrt_cbcl_int",  exposure = "urine_z_metabolite_score_int6"),
  urine_adhd = list(outcome = "sqrt_cbcl_adhd", exposure = "urine_z_metabolite_score_adhd6")
)

# ==============================================================
# Fit Linear Mixed-Effects Models Across Imputed Datasets
# ==============================================================
mi_results_list <- lapply(names(model_specs), function(name) {
  
  spec <- model_specs[[name]]
  
  d_first <- long_list[[1]]
  
  n_subjects <- d_first %>%
    filter(
      !is.na(.data[[spec$outcome]]),
      !is.na(.data[[spec$exposure]])
    ) %>%
    distinct(helixid) %>%
    nrow()
  
  n_observations <- d_first %>%
    filter(
      !is.na(.data[[spec$outcome]]),
      !is.na(.data[[spec$exposure]])
    ) %>%
    nrow()
  
  
  models <- lapply(long_list, function(d) {
    
    formula <- as.formula(
      paste0(spec$outcome, " ~ ",
             spec$exposure, " * time + ",
             "sex + mother_age_centered + z_who_bmi_7y + ",
             "passive_smk_6y + mum_educ + m_greek + ",
             "(1 | helixid)")
    )
    
    lmer(formula, data = d, REML = FALSE)
  })
  
  pooled <- pool(models)
  
  summary(pooled, conf.int = TRUE) %>%
    mutate(
      model = name,
      n_subjects = n_subjects,
      n_observations = n_observations,
      .before = 1
    )
})

mi_results <- bind_rows(mi_results_list) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

mi_results


#Validate with the non imputed linear mixed model 
d_long_raw <- read_excel("filtered_covariates_all_long_format.xlsx")

sapply(
  d_long_raw[, c("sex", "m_greek", "mum_educ", "passive_smk_6y", "bmi_cutoff_merge")],
  class
)

d_long_raw <- read_excel("filtered_covariates_all_long_format.xlsx") %>%
  mutate(
    sex = factor(sex, levels = c("Male", "Female")),
    m_greek = factor(m_greek, levels = c("Greek", "Other")),
    mum_educ = factor(mum_educ, levels = c("Low", "Medium", "High")),
    passive_smk_6y = factor(passive_smk_6y, levels = c("No", "Yes")),
    bmi_cutoff_merge = factor(
      bmi_cutoff_merge,
      levels = c("underweight/normal", "overweight/obese")
    )
  )

#Define Complete Case at Subject Level

cc_results_list <- lapply(names(model_specs), function(name) {
  
  spec <- model_specs[[name]]
  
  complete_ids <- d_long_raw %>%
    group_by(helixid) %>%
    filter(n_distinct(time) == 3) %>%
    filter(
      !is.na(.data[[spec$outcome]]),
      !is.na(.data[[spec$exposure]]),
      !is.na(sex),
      !is.na(mother_age_centered),
      !is.na(z_who_bmi_7y),
      !is.na(passive_smk_6y),
      !is.na(mum_educ),
      !is.na(m_greek)
    ) %>%
    summarise(n_rows = n()) %>%
    filter(n_rows == 3) %>%
    pull(helixid)
  
  n_subjects <- length(unique(complete_ids))
  d_complete <- d_long_raw %>%
    filter(helixid %in% complete_ids)
  n_obs <- nrow(d_complete)
  formula <- as.formula(
    paste0(spec$outcome, " ~ ",
           spec$exposure, " * time + ",
           "sex + mother_age_centered + z_who_bmi_7y + ",
           "passive_smk_6y + mum_educ + m_greek + ",
           "(1 | helixid)")
  )
  
  model <- lmer(formula, data = d_complete, REML = FALSE)
  
  tidy(model, effects = "fixed", conf.int = TRUE) %>%
    mutate(
      model = name,
      n_subjects = n_subjects,
      n_observations = n_obs,
      .before = 1
    )
})

cc_results <- bind_rows(cc_results_list) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

cc_results

# Standardize MI results
mi_plot <- mi_results %>%
  rename(
    conf.low = `2.5 %`,
    conf.high = `97.5 %`
  ) %>%
  mutate(type = "MI")

# Standardize CC results
cc_plot <- cc_results %>%
  mutate(type = "Complete Case")

# Combine
plot_data <- bind_rows(mi_plot, cc_plot)

#Keep only the main exposure 
plot_data <- plot_data %>%
  filter(!grepl(":", term)) %>%                 # remove interactions
  filter(grepl("metabolite_score", term)) %>%   # keep exposure only
  mutate(
    exposure_type = ifelse(grepl("serum", model), "Serum", "Urine"),
    outcome = case_when(
      grepl("ext", model)  ~ "Externalizing",
      grepl("int", model)  ~ "Internalizing",
      grepl("adhd", model) ~ "ADHD"
    )
  )


plot_data$outcome <- factor(
  plot_data$outcome,
  levels = c("Internalizing", "Externalizing", "ADHD")
)

ggplot(plot_data,
       aes(x = type,
           y = estimate,
           color = type)) +
  
  geom_point(size = 3,
             position = position_dodge(width = 0.4)) +
  
  geom_errorbar(
    aes(ymin = conf.low,
        ymax = conf.high),
    width = 0.15,
    linewidth = 1,
    position = position_dodge(width = 0.4)
  ) +
  
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "black",
             linewidth = 0.8) +
  
  facet_grid(rows = vars(exposure_type),
             cols = vars(outcome)) +
  
  scale_color_manual(
    values = c(
      "MI" = "#0072B2",            # blue
      "Complete Case" = "#D55E00"  # orange/red
    )
  ) +
  
  labs(
    title = expression(beta ~ "Coefficient (95% CI) for Metabolite Score"),
    x = "",
    y = expression(beta ~ "Coefficient (95% CI)"),
    color = "Analysis"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    strip.background = element_rect(fill = "gray90", color = "gray90"),
    strip.text = element_text(size = 18, face = "bold", color = "black"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 14)
  )

#########################################################################################
############Continue with Linear Mixed  Model Analysis with the imputed outcomes####
# ==============================================================
# Clean R Environment
# ==============================================================
rm(list = ls())
load("imputed_datasets_long.RData")
sapply(
  long_list[[1]][, c("sex", "m_greek", "mum_educ", "passive_smk_6y", "bmi_cutoff_merge", "time")],
  is.factor
)

# Each entry links a behavioral outcome with its corresponding
# metabolite exposure variable (serum or urine)
model_specs <- list(
  serum_ext  = list(outcome = "sqrt_cbcl_ext",  exposure = "serum_z_metabolite_score_ext6"),
  serum_int  = list(outcome = "sqrt_cbcl_int",  exposure = "serum_z_metabolite_score_int6"),
  serum_adhd = list(outcome = "sqrt_cbcl_adhd", exposure = "serum_z_metabolite_score_adhd6"),
  urine_ext  = list(outcome = "sqrt_cbcl_ext",  exposure = "urine_z_metabolite_score_ext6"),
  urine_int  = list(outcome = "sqrt_cbcl_int",  exposure = "urine_z_metabolite_score_int6"),
  urine_adhd = list(outcome = "sqrt_cbcl_adhd", exposure = "urine_z_metabolite_score_adhd6")
)


# ==============================================================
# Run Simple Linear Mixed Models
# ==============================================================

mi_simple_results_list <- lapply(names(model_specs), function(name) {
  
  spec <- model_specs[[name]]
  
  models <- lapply(long_list, function(d) {
    
    formula <- as.formula(
      paste0(spec$outcome, " ~ ",
             spec$exposure, " * factor(time) + ",
             "sex + (1 | helixid)")
    )
    
    lmer(formula, data = d, REML = FALSE)
  })
  
  pooled <- pool(models)
  
  summary(pooled, conf.int = TRUE) %>%
    mutate(
      model_type = "simple_model",
      model = name,
      .before = 1
    )
})

# ==============================================================
# Combine Results
# ==============================================================

mi_simple_results <- bind_rows(mi_simple_results_list) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

mi_simple_results

# ==============================================================
# Fully Adjusted Linear Mixed-Effects Models Across Imputed Datasets
# ==============================================================

mi_full_results_list <- lapply(names(model_specs), function(name) {
  
  spec <- model_specs[[name]]
  
  d_first <- long_list[[1]]
  
  n_subjects <- d_first %>%
    filter(
      !is.na(.data[[spec$outcome]]),
      !is.na(.data[[spec$exposure]])
    ) %>%
    distinct(helixid) %>%
    nrow()
  
  n_observations <- d_first %>%
    filter(
      !is.na(.data[[spec$outcome]]),
      !is.na(.data[[spec$exposure]])
    ) %>%
    nrow()
  
  models <- lapply(long_list, function(d) {
    
    formula <- as.formula(
      paste0(spec$outcome, " ~ ",
             spec$exposure, " * factor(time) + ",
             "sex + ",
             "m_greek + ",
             "mother_age_centered + ",
             "z_who_bmi_7y + ",
             "mum_educ + ",
             "passive_smk_6y + ",
             "(1 | helixid)")
    )
    
    lmer(formula, data = d, REML = FALSE)
  })
  
  pooled <- pool(models)
  
  summary(pooled, conf.int = TRUE) %>%
    mutate(
      model_type = "full_model",
      model = name,
      n_subjects = n_subjects,
      n_observations = n_observations,
      .before = 1
    )
})

mi_full_results <- bind_rows(mi_full_results_list) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

mi_full_results

################################################################################
# ==============================================================
# Test Exposure × Sex Interaction Across Imputed Datasets
# ==============================================================

mi_sex_interaction_results_list <- lapply(names(model_specs), function(name) {
  
  spec <- model_specs[[name]]
  
  models <- lapply(long_list, function(d) {
    
    formula <- as.formula(
      paste0(spec$outcome, " ~ ",
             spec$exposure, " * sex + ",
             "m_greek + ",
             "mother_age_centered + ",
             "z_who_bmi_7y + ",
             "mum_educ + ",
             "passive_smk_6y + ",
             "(1 | helixid)")
    )
    
    lmer(formula, data = d, REML = FALSE)
  })
  
  pooled <- pool(models)
  
  summary(pooled, conf.int = TRUE) %>%
    mutate(
      model_type = "sex_interaction_model",
      model = name,
      .before = 1
    )
})

mi_sex_interaction_results <- bind_rows(mi_sex_interaction_results_list) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

mi_sex_interaction_results

# ==============================================================
# Sex-Stratified Models Across Imputed Datasets
# ==============================================================
# Check the levels of sex variable
levels(long_list[[1]]$sex)
sex_levels <- c("Male", "Female")

mi_sex_stratified_results_list <- lapply(names(model_specs), function(name) {
  
  spec <- model_specs[[name]]
  
  lapply(sex_levels, function(sex_group) {
    
    subject_counts <- c()
    
    models <- lapply(long_list, function(d) {
      
      d_strat <- d %>% 
        filter(sex == sex_group) %>%
        filter(
          !is.na(.data[[spec$outcome]]),
          !is.na(.data[[spec$exposure]]),
          !is.na(m_greek),
          !is.na(mother_age_centered),
          !is.na(z_who_bmi_7y),
          !is.na(mum_educ),
          !is.na(passive_smk_6y),
          !is.na(helixid)
        )
      
      subject_counts <<- c(subject_counts, dplyr::n_distinct(d_strat$helixid))
      
      formula <- as.formula(
        paste0(spec$outcome, " ~ ",
               spec$exposure, " + ",
               "m_greek + ",
               "mother_age_centered + ",
               "z_who_bmi_7y + ",
               "mum_educ + ",
               "passive_smk_6y + ",
               "(1 | helixid)")
      )
      
      lmer(formula, data = d_strat, REML = FALSE)
    })
    
    pooled <- pool(models)
    
    summary(pooled, conf.int = TRUE) %>%
      mutate(
        model_type = "sex_stratified_model",
        model = name,
        sex = sex_group,
        N_subjects = round(mean(subject_counts), 1),
        N_subjects_min = min(subject_counts),
        N_subjects_max = max(subject_counts),
        .before = 1
      )
    
  }) %>% bind_rows()
})

mi_sex_stratified_results <- bind_rows(mi_sex_stratified_results_list) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

mi_sex_stratified_results

#==============================================================
# Test Exposure × BMI Interaction Across Imputed Datasets
# ==============================================================

mi_bmi_interaction_results_list <- lapply(names(model_specs), function(name) {
  
  spec <- model_specs[[name]]
  
  models <- lapply(long_list, function(d) {
    
    formula <- as.formula(
      paste0(spec$outcome, " ~ ",
             spec$exposure, " * bmi_cutoff_merge + ",
             "sex + ",
             "m_greek + ",
             "mother_age_centered + ",
             "mum_educ + ",
             "passive_smk_6y + ",
             "(1 | helixid)")
    )
    
    lmer(formula, data = d, REML = FALSE)
  })
  
  pooled <- pool(models)
  
  summary(pooled, conf.int = TRUE) %>%
    mutate(
      model_type = "bmi_interaction_model",
      model = name,
      .before = 1
    )
})

mi_bmi_interaction_results <- bind_rows(mi_bmi_interaction_results_list) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

mi_bmi_interaction_results


#==============================================================
# Stratified Analysis by BMI Category Across Imputed Datasets
# ==============================================================

bmi_levels <- c("underweight/normal", "overweight/obese")

mi_bmi_stratified_results_list <- lapply(names(model_specs), function(name) {
  
  spec <- model_specs[[name]]
  
  lapply(bmi_levels, function(bmi_group) {
    
    n_subjects_vec <- c()
    
    models <- lapply(long_list, function(d) {
      
      d_strat <- d %>% 
        filter(bmi_cutoff_merge == bmi_group) %>%
        filter(
          !is.na(.data[[spec$outcome]]),
          !is.na(.data[[spec$exposure]]),
          !is.na(sex),
          !is.na(mother_age_centered),
          !is.na(mum_educ),
          !is.na(passive_smk_6y),
          !is.na(helixid)
        )
      
      n_subjects_vec <<- c(n_subjects_vec, dplyr::n_distinct(d_strat$helixid))
      
      formula <- as.formula(
        paste0(spec$outcome, " ~ ",
               spec$exposure, " + ",
               "sex + ",
               "mother_age_centered + ",
               "mum_educ + ",
               "passive_smk_6y + ",
               "(1 | helixid)")
      )
      
      lmer(formula, data = d_strat, REML = FALSE)
    })
    
    pooled <- pool(models)
    
    summary(pooled, conf.int = TRUE) %>%
      mutate(
        model_type = "bmi_stratified_model",
        model = name,
        bmi_group = bmi_group,
        N_subjects = round(mean(n_subjects_vec), 1),
        .before = 1
      )
    
  }) %>% bind_rows()
})

mi_bmi_stratified_results <- bind_rows(mi_bmi_stratified_results_list) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

mi_bmi_stratified_results

#



# ==============================================================
# Test Exposure × Age Interaction Across Imputed Datasets
# ==============================================================

mi_age_interaction_results_list <- lapply(names(model_specs), function(name) {
  
  spec <- model_specs[[name]]
  
  models <- lapply(long_list, function(d) {
    
    formula <- as.formula(
      paste0(spec$outcome, " ~ ",
             spec$exposure, " * age + ",
             "sex + ",
             "m_greek + ",
             "mother_age_centered + ",
             "z_who_bmi_7y + ",
             "mum_educ + ",
             "passive_smk_6y + ",
             "(1 | helixid)")
    )
    lmer(formula, data = d, REML = FALSE)
    
  })
  
  pooled <- pool(models)
  
  summary(pooled, conf.int = TRUE) %>%
    mutate(
      model_type = "age_interaction_model",
      model = name,
      .before = 1
    )
})

mi_age_interaction_results <- bind_rows(mi_age_interaction_results_list) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

mi_age_interaction_results


# ==============================================================
# Age-Specific Linear Regression Models (Wide Imputed Datasets)
# ==============================================================
load("imputed_datasets_wide.RData")

#Transfrom the variables in factors 
completed_list <- lapply(completed_list, function(d) {
  d %>%
    mutate(
      sex = factor(sex, levels = c("Male", "Female")),
      m_greek = factor(m_greek, levels = c("Greek", "Other")),
      mum_educ = factor(mum_educ, levels = c("Low", "Medium", "High")),
      passive_smk_6y = factor(passive_smk_6y, levels = c("No", "Yes")),
      bmi_cutoff_merge = factor(
        bmi_cutoff_merge,
        levels = c("underweight/normal", "overweight/obese")
      )
    )
})
#Check if they have transformed correctly
sapply(
  completed_list[[1]][, c("sex", "m_greek", "mum_educ", "passive_smk_6y", "bmi_cutoff_merge")],
  class
)

wide_list <- completed_list
age_outcomes <- list(
  "6" = c(ext = "sqrt_cbcl_ext6",
          int = "sqrt_cbcl_pro_int6",
          adhd = "sqrt_adhd_pro6"),
  
  "11" = c(ext = "sqrt_cbcl_ext11",
           int = "sqrt_cbcl_int11",
           adhd = "sqrt_cbcl_adhd11"),
  
  "15" = c(ext = "sqrt_cbcl_ext15",
           int = "sqrt_cbcl_int15",
           adhd = "sqrt_cbcl_adhd15")
)

ages <- c("6","11","15")

mi_age_specific_results_list <- lapply(names(model_specs), function(name) {
  
  spec <- model_specs[[name]]
  
  outcome_type <- sub(".*_", "", name)   # ext / int / adhd
  
  lapply(ages, function(age_val) {
    
    outcome_var <- age_outcomes[[age_val]][outcome_type]
    
    models <- lapply(wide_list, function(d) {
      
      formula <- as.formula(
        paste0(outcome_var, " ~ ",
               spec$exposure, " + ",
               "sex + ",
               "m_greek + ",
               "mother_age_centered + ",
               "mum_educ + ",
               "passive_smk_6y")
      )
      
      lm(formula, data = d)
      
    })
    
    pooled <- pool(models)
    
    summary(pooled, conf.int = TRUE) %>%
      mutate(
        model_type = "age_specific_linear_model",
        model = name,
        age = age_val,
        .before = 1
      )
    
  }) %>% bind_rows()
})

mi_age_specific_results <- bind_rows(mi_age_specific_results_list) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

mi_age_specific_results


# ==============================================================
# Prepare Data for Plot
# ==============================================================

plot_data <- mi_age_specific_results %>%
  filter(grepl("metabolite_score", term)) %>%
  mutate(
    Age = factor(age,
                 levels = c("6","11","15"),
                 labels = c("Age 6","Age 11","Age 15"))
  )
# ==============================================================
# Plot
# ==============================================================

plots <- lapply(unique(plot_data$model), function(m) {
  
  df <- plot_data %>% filter(model == m)
  
  ggplot(df, aes(x = Age, y = estimate, group = 1)) +
    
    geom_line(color = "steelblue", linewidth = 1) +
    
    geom_point(size = 4, color = "steelblue") +
    
    geom_errorbar(
      aes(ymin = `2.5 %`, ymax = `97.5 %`),
      width = 0.2,
      color = "darkblue"
    ) +
    
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    
    labs(
      title = paste("Metabolite Score Association with", m),
      x = "Age at Assessment",
      y = "Beta Coefficient (95% CI)"
    ) +
    
    theme_minimal(base_size = 14)
  
})

names(plots) <- unique(plot_data$model)

plots$serum_ext
plots$serum_int
plots$serum_adhd
plots$urine_ext
plots$serum_int
plots$urine_adhd


#Binary analysis 
d1_binary <- read_excel("filtered_covariates_plus_z_metabolites_scores.xlsx")

d1_binary <- d1_binary %>%
  mutate(
    
    # ==========================================================
    # Externalizing Diagnosis (Binary)
    # ==========================================================
    
    diagnosis_binary_ext6 = case_when(
      Diagnosis_externalizing_6_combined == "Normal" ~ 0,
      Diagnosis_externalizing_6_combined == "Borderline/Clinical" ~ 1,
      TRUE ~ NA_real_
    ),
    
    diagnosis_binary_ext11 = case_when(
      Diagnosis_externalizing_11_combined == "Normal" ~ 0,
      Diagnosis_externalizing_11_combined == "Borderline/Clinical" ~ 1,
      TRUE ~ NA_real_
    ),
    
    diagnosis_binary_ext15 = case_when(
      Diagnosis_externalizing_15_combined == "Normal" ~ 0,
      Diagnosis_externalizing_15_combined == "Borderline/Clinical" ~ 1,
      TRUE ~ NA_real_
    ),
    
    
    # ==========================================================
    # Internalizing Diagnosis (Binary)
    # ==========================================================
    
    diagnosis_binary_int6 = case_when(
      Diagnosis_internalizing_6_combined == "Normal" ~ 0,
      Diagnosis_internalizing_6_combined == "Borderline/Clinical" ~ 1,
      TRUE ~ NA_real_
    ),
    
    diagnosis_binary_int11 = case_when(
      Diagnosis_internalizing_11_combined == "Normal" ~ 0,
      Diagnosis_internalizing_11_combined == "Borderline/Clinical" ~ 1,
      TRUE ~ NA_real_
    ),
    
    diagnosis_binary_int15 = case_when(
      Diagnosis_internalizing_15_combined == "Normal" ~ 0,
      Diagnosis_internalizing_15_combined == "Borderline/Clinical" ~ 1,
      TRUE ~ NA_real_
    )
  )

#Check if the variables are factors 
sapply(
  d1_binary[, c("sex", "m_greek", "mum_educ", "passive_smk_6y", "bmi_cutoff_merge")],
  is.factor
)

#Transform the variables as factors 
d1_binary <- d1_binary %>%
  mutate(
    sex = factor(sex, levels = c("Male", "Female")),
    m_greek = factor(m_greek, levels = c("Greek", "Other")),
    mum_educ = factor(mum_educ, levels = c("Low", "Medium", "High")),
    passive_smk_6y = factor(passive_smk_6y, levels = c("No", "Yes")),
    bmi_cutoff_merge = factor(
      bmi_cutoff_merge,
      levels = c("underweight/normal", "overweight/obese")
    ),
    
    Diagnosis_internalizing_6_combined = factor(
      Diagnosis_internalizing_6_combined,
      levels = c("Normal", "Borderline/Clinical")
    ),
    Diagnosis_internalizing_11_combined = factor(
      Diagnosis_internalizing_11_combined,
      levels = c("Normal", "Borderline/Clinical")
    ),
    Diagnosis_internalizing_15_combined = factor(
      Diagnosis_internalizing_15_combined,
      levels = c("Normal", "Borderline/Clinical")
    ),
    
    Diagnosis_externalizing_6_combined = factor(
      Diagnosis_externalizing_6_combined,
      levels = c("Normal", "Borderline/Clinical")
    ),
    Diagnosis_externalizing_11_combined = factor(
      Diagnosis_externalizing_11_combined,
      levels = c("Normal", "Borderline/Clinical")
    ),
    Diagnosis_externalizing_15_combined = factor(
      Diagnosis_externalizing_15_combined,
      levels = c("Normal", "Borderline/Clinical")
    )
  )

d1_binary <- d1_binary %>%
  mutate(
    diagnosis_binary_int6 = factor(diagnosis_binary_int6, levels = c(0, 1), labels = c("Normal", "Borderline/Clinical")),
    diagnosis_binary_int11 = factor(diagnosis_binary_int11, levels = c(0, 1), labels = c("Normal", "Borderline/Clinical")),
    diagnosis_binary_int15 = factor(diagnosis_binary_int15, levels = c(0, 1), labels = c("Normal", "Borderline/Clinical")),
    
    diagnosis_binary_ext6 = factor(diagnosis_binary_ext6, levels = c(0, 1), labels = c("Normal", "Borderline/Clinical")),
    diagnosis_binary_ext11 = factor(diagnosis_binary_ext11, levels = c(0, 1), labels = c("Normal", "Borderline/Clinical")),
    diagnosis_binary_ext15 = factor(diagnosis_binary_ext15, levels = c(0, 1), labels = c("Normal", "Borderline/Clinical"))
  )

#Confirm if they have transformed correct 
sapply(
  d1_binary[, c(
    "Diagnosis_internalizing_6_combined",
    "Diagnosis_internalizing_11_combined",
    "Diagnosis_internalizing_15_combined",
    "Diagnosis_externalizing_6_combined",
    "Diagnosis_externalizing_11_combined",
    "Diagnosis_externalizing_15_combined"
  )],
  is.factor
)
sapply(
  d1_binary[, c("sex", "m_greek", "mum_educ", "passive_smk_6y", "bmi_cutoff_merge")],
  is.factor
)
# ==============================================================
# Define Logistic Regression Specifications
# ==============================================================

diagnosis_specs <- list(
  serum_ext = list(
    outcome_prefix = "diagnosis_binary_ext",
    exposure = "serum_z_metabolite_score_ext6"
  ),
  urine_ext = list(
    outcome_prefix = "diagnosis_binary_ext",
    exposure = "urine_z_metabolite_score_ext6"
  ),
  serum_int = list(
    outcome_prefix = "diagnosis_binary_int",
    exposure = "serum_z_metabolite_score_int6"
  ),
  urine_int = list(
    outcome_prefix = "diagnosis_binary_int",
    exposure = "urine_z_metabolite_score_int6"
  )
)

ages <- c(6, 11, 15)

logistic_results_list <- lapply(names(diagnosis_specs), function(name) {
  
  spec <- diagnosis_specs[[name]]
  
  lapply(ages, function(age_val) {
    
    outcome_var <- paste0(spec$outcome_prefix, age_val)
    
    vars_needed <- c(
      outcome_var,
      spec$exposure,
      "sex",
      "mother_age_centered",
      "passive_smk_6y",
      "z_who_bmi_7y"
    )
    
    d_model <- d1_binary %>%
      filter(if_all(all_of(vars_needed[-1]), ~ !is.na(.)))
    
    n_subjects <- d_model %>%
      distinct(helixid) %>%
      nrow()
    
    formula <- as.formula(
      paste0(outcome_var, " ~ ",
             spec$exposure, " + ",
             "sex + ",
             "mother_age_centered + ",
             "z_who_bmi_7y + ",
             "passive_smk_6y")
    )
    
    model <- glm(
      formula,
      data = d_model,
      family = binomial(link = "logit")
    )
    
    tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
      mutate(
        model_type = "diagnosis_logistic_model",
        model = name,
        age = age_val,
        n_subjects = n_subjects,
        .before = 1
      )
    
  }) %>% bind_rows()
})

logistic_results <- bind_rows(logistic_results_list) %>%
  mutate(
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    statistic = round(statistic, 3),
    p.value = round(p.value, 3),
    conf.low = round(conf.low, 3),
    conf.high = round(conf.high, 3)
  )

logistic_results

#Plot of odds ratio 

logistic_results_raw <- bind_rows(logistic_results_list)

results <- logistic_results_raw %>%
  filter(term %in% c(
    "serum_z_metabolite_score_ext6",
    "urine_z_metabolite_score_ext6",
    "serum_z_metabolite_score_int6",
    "urine_z_metabolite_score_int6"
  )) %>%
  mutate(
    Age = age,
    OR = estimate,
    Lower_CI = conf.low,
    Upper_CI = conf.high,
    Outcome = case_when(
      grepl("_ext$", model) ~ "Externalizing",
      grepl("_int$", model) ~ "Internalizing"
    ),
    Source = case_when(
      grepl("^serum", model) ~ "Serum",
      grepl("^urine", model) ~ "Urine"
    ),
    Age = factor(Age, levels = c(6, 11, 15)),
    Outcome = factor(Outcome, levels = c("Internalizing", "Externalizing")),
    Source = factor(Source, levels = c("Serum", "Urine"))
  )

ggplot(results, aes(x = Age, y = OR, color = Outcome, group = interaction(Outcome, Source))) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI),
                width = 0.2, color = "darkblue", linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
  facet_grid(rows = vars(Source), cols = vars(Outcome)) +
  scale_color_manual(values = c("Internalizing" = "steelblue",
                                "Externalizing" = "steelblue")) +
  labs(
    title = "Odds Ratios for Borderline/Clinical Symptoms by Age, Outcome, and Biospecimen",
    x = "Age at assessment",
    y = "Odds Ratio (95% CI)"
  ) +
  ylim(0, 7) +
  theme_minimal(base_size = 16) +
  theme(
    strip.background = element_rect(fill = "gray90", color = "gray90"),
    strip.text = element_text(size = 22, face = "bold", color = "black"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 22),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "none"
  )

#Effect modification for sex and bmi 
run_stratified_linear_models_mi <- function(complete_list, outcomes, metabolite_var) {
  
  fit_lm_stratum <- function(data_subset, outcome, metabolite_var, strat_type = c("sex", "bmi")) {
    strat_type <- match.arg(strat_type)
    
    if (strat_type == "sex") {
      vars_needed <- c(
        outcome, metabolite_var,
        "mother_age_centered", "mum_educ",
        "passive_smk_6y", "z_who_bmi_7y"
      )
      
      d_model <- data_subset %>%
        filter(if_all(all_of(vars_needed), ~ !is.na(.)))
      
      formula <- as.formula(
        paste0(
          outcome, " ~ ", metabolite_var,
          " + mother_age_centered + mum_educ + passive_smk_6y + z_who_bmi_7y"
        )
      )
    }
    
    if (strat_type == "bmi") {
      vars_needed <- c(
        outcome, metabolite_var,
        "sex", "mother_age_centered",
        "mum_educ", "passive_smk_6y"
      )
      
      d_model <- data_subset %>%
        filter(if_all(all_of(vars_needed), ~ !is.na(.)))
      
      formula <- as.formula(
        paste0(
          outcome, " ~ ", metabolite_var,
          " + sex + mother_age_centered + mum_educ + passive_smk_6y"
        )
      )
    }
    
    lm(formula, data = d_model)
  }
  
  run_one_stratum <- function(stratum_name, filter_expr, outcome, strat_type) {
    models <- lapply(complete_list, function(d) {
      d_strat <- d %>% filter(!!rlang::parse_expr(filter_expr))
      fit_lm_stratum(d_strat, outcome, metabolite_var, strat_type)
    })
    
    pooled <- pool(as.mira(models))
    
    summary(pooled, conf.int = TRUE) %>%
      filter(term == metabolite_var) %>%
      mutate(
        Model = stratum_name,
        N_mean = round(mean(sapply(models, nobs)), 1),
        .before = 1
      )
  }
  
  bind_rows(
    run_one_stratum("sex_male_6",    'sex == "Male"', outcomes[1], "sex"),
    run_one_stratum("sex_female_6",  'sex == "Female"', outcomes[1], "sex"),
    run_one_stratum("sex_male_11",   'sex == "Male"', outcomes[2], "sex"),
    run_one_stratum("sex_female_11", 'sex == "Female"', outcomes[2], "sex"),
    run_one_stratum("sex_male_15",   'sex == "Male"', outcomes[3], "sex"),
    run_one_stratum("sex_female_15", 'sex == "Female"', outcomes[3], "sex"),
    
    run_one_stratum("bmi_normal_6",  'bmi_cutoff_merge == "underweight/normal"', outcomes[1], "bmi"),
    run_one_stratum("bmi_obese_6",   'bmi_cutoff_merge == "overweight/obese"', outcomes[1], "bmi"),
    run_one_stratum("bmi_normal_11", 'bmi_cutoff_merge == "underweight/normal"', outcomes[2], "bmi"),
    run_one_stratum("bmi_obese_11",  'bmi_cutoff_merge == "overweight/obese"', outcomes[2], "bmi"),
    run_one_stratum("bmi_normal_15", 'bmi_cutoff_merge == "underweight/normal"', outcomes[3], "bmi"),
    run_one_stratum("bmi_obese_15",  'bmi_cutoff_merge == "overweight/obese"', outcomes[3], "bmi")
  )
}


extract_pooled_linear_results <- function(results_df) {
  results_df %>%
    transmute(
      Model,
      Beta = round(estimate, 3),
      CI_Lower = round(`2.5 %`, 3),
      CI_Upper = round(`97.5 %`, 3),
      `Beta (95% CI)` = paste0(
        round(estimate, 3), " (",
        round(`2.5 %`, 3), ", ",
        round(`97.5 %`, 3), ")"
      ),
      p_value = round(p.value, 4),
      N = N_mean
    )
}

outcomes_ext  <- c("sqrt_cbcl_ext6",     "sqrt_cbcl_ext11",   "sqrt_cbcl_ext15")
outcomes_int  <- c("sqrt_cbcl_pro_int6", "sqrt_cbcl_int11",   "sqrt_cbcl_int15")
outcomes_adhd <- c("sqrt_adhd_pro6",     "sqrt_cbcl_adhd11",  "sqrt_cbcl_adhd15")

mi_serum_ext_raw <- run_stratified_linear_models_mi(
  complete_list = completed_list,
  outcomes = outcomes_ext,
  metabolite_var = "serum_z_metabolite_score_ext6"
)

mi_serum_int_raw <- run_stratified_linear_models_mi(
  complete_list = completed_list,
  outcomes = outcomes_int,
  metabolite_var = "serum_z_metabolite_score_int6"
)

mi_serum_adhd_raw <- run_stratified_linear_models_mi(
  complete_list = completed_list,
  outcomes = outcomes_adhd,
  metabolite_var = "serum_z_metabolite_score_adhd6"
)

mi_urine_ext_raw <- run_stratified_linear_models_mi(
  complete_list = completed_list,
  outcomes = outcomes_ext,
  metabolite_var = "urine_z_metabolite_score_ext6"
)

mi_urine_int_raw <- run_stratified_linear_models_mi(
  complete_list = completed_list,
  outcomes = outcomes_int,
  metabolite_var = "urine_z_metabolite_score_int6"
)

mi_urine_adhd_raw <- run_stratified_linear_models_mi(
  complete_list = completed_list,
  outcomes = outcomes_adhd,
  metabolite_var = "urine_z_metabolite_score_adhd6"
)

summary_serum_ext <- extract_pooled_linear_results(mi_serum_ext_raw) %>%
  mutate(Outcome = "Externalizing", Source = "Serum")

summary_serum_int <- extract_pooled_linear_results(mi_serum_int_raw) %>%
  mutate(Outcome = "Internalizing", Source = "Serum")

summary_serum_adhd <- extract_pooled_linear_results(mi_serum_adhd_raw) %>%
  mutate(Outcome = "ADHD", Source = "Serum")

summary_urine_ext <- extract_pooled_linear_results(mi_urine_ext_raw) %>%
  mutate(Outcome = "Externalizing", Source = "Urine")

summary_urine_int <- extract_pooled_linear_results(mi_urine_int_raw) %>%
  mutate(Outcome = "Internalizing", Source = "Urine")

summary_urine_adhd <- extract_pooled_linear_results(mi_urine_adhd_raw) %>%
  mutate(Outcome = "ADHD", Source = "Urine")

summary_all_stratified_sex_bmi_age <- bind_rows(
  summary_serum_ext, summary_serum_int, summary_serum_adhd,
  summary_urine_ext, summary_urine_int, summary_urine_adhd
) %>%
  mutate(
    Stratification = ifelse(grepl("sex", Model), "Sex", "BMI"),
    Subgroup = case_when(
      grepl("female", Model) ~ "Female",
      grepl("male", Model)   ~ "Male",
      grepl("normal", Model) ~ "Underweight/Normal",
      grepl("obese", Model)  ~ "Overweight/Obese"
    ),
    Age = case_when(
      grepl("_6$", Model)  ~ "6",
      grepl("_11$", Model) ~ "11",
      grepl("_15$", Model) ~ "15"
    )
  ) %>%
  select(Source, Outcome, Stratification, Subgroup, Age, N, `Beta (95% CI)`, p_value)

summary_all_stratified_sex_bmi_age


# ==============================================================
# Export All Model Results to a Multi-Sheet Excel Workbook
# ==============================================================
results_list <- list(
  "MI_Simple_Models" = mi_simple_results,
  "MI_Full_Models" = mi_full_results,
  "MI_Sex_Interaction" = mi_sex_interaction_results,
  "MI_Sex_Stratified" = mi_sex_stratified_results,
  "MI_BMI_Interaction" = mi_bmi_interaction_results,
  "MI_BMI_Stratified" = mi_bmi_stratified_results,
  "MI_Age_Interaction" = mi_age_interaction_results,
  "MI_Age_Specific" = mi_age_specific_results,
  "Logistic_Results" = logistic_results,
  "Stratified_Summary" = summary_all_stratified_sex_bmi_age
)

#write_xlsx(results_list, path = "all_model_results_with_imputed_scores.xlsx")


#\ END OF SCRIPT 


