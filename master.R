### Association of GI-related antibiotics and deprivation by LSOA in England, 2013-2022 ###
### Leah Linton 2024 ###


## Prelims -----------

# Clear environment
rm(list=ls())

# List of required packages
packages <- c("dplyr", "tidyr", "spData", "spdep", "sp", "sf", "geepack", "data.table", "future.apply", "car", "sandwich", "lmtest", "MASS", "ggplot2")

# Function to check and install missing packages
install_if_needed <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Install packages if needed
sapply(packages, install_if_needed)

# Load packages
library(dplyr)
library(tidyr)
library(spData)
library(spdep)
library(sp)
library(sf)
library(geepack)
library(data.table)
library(future.apply)
library(car)
library(sandwich)  
library(lmtest)  
library(MASS)
library(ggplot2)


## Read in and link data -----------

folder.path = "all GI antibiotics" # Change as required

# IF statement to aggregate datasets if all GI antibiotics selected
if(folder.path == "all GI antibiotics"){
  
  # Function to extract year and month from file name and read the CSV file
  read_and_add_date_info <- function(file) {
    # Extract year and month from the file name (assuming the format is consistent)
    # For example, p_1_07_202304.csv would give year = 2023, month = 04
    date_info <- strsplit(gsub("p_1_07_|\\.csv", "", file), "_")[[1]]
    year <- date_info[[5]]
    month <- date_info[[6]]
    
    # Read the CSV file
    data <- read.csv(file)
    
    # Add year and month columns to the data frame
    data$year <- as.integer(year)
    data$month <- as.character(month)
    
    return(data)
  }
  
  # Get list of all CSV files ending with "LSOA" in the data directory
  csv_files <- list.files(path = "clarithromycin", pattern = "LSOA\\.csv$", full.names = TRUE)
  
  # Read in all CSV files and row bind them
  combined_data1 <- bind_rows(lapply(csv_files, read_and_add_date_info))
  
  # Get list of all CSV files ending with "LSOA" in the data directory
  csv_files <- list.files(path = "azithromycin", pattern = "LSOA\\.csv$", full.names = TRUE)
  
  # Read in all CSV files and row bind them
  combined_data2 <- bind_rows(lapply(csv_files, read_and_add_date_info))
  
  # Get list of all CSV files ending with "LSOA" in the data directory
  csv_files <- list.files(path = "ciprofloxacin", pattern = "LSOA\\.csv$", full.names = TRUE)
  
  # Read in all CSV files and row bind them
  combined_data3 <- bind_rows(lapply(csv_files, read_and_add_date_info))
  
  # Keep required columns and calculate prescribing rate
  
  combined_data <- combined_data3 %>%
    left_join(combined_data2, by = c("lsoa11","all_ages","year","month")) %>%
    left_join(combined_data1, by = c("lsoa11","all_ages","year","month"))
  
  data.all <- combined_data %>%
    dplyr::select(lsoa11, all_ages, items.x, items.y, items, year, month) %>%
    mutate(items = ifelse(is.na(items), 0, items)) %>%
    mutate(items = items + items.x + items.y) %>%
    dplyr::select(-c(items.x, items.y)) %>%
    rename("lsoa" = "lsoa11",
           "pop.all" = "all_ages",
           "n.items" = "items") %>%
    filter(!grepl("^W", lsoa)) %>%
    mutate(r.items = 1000 * n.items/pop.all) %>%
    arrange(year, month)
  
} else {


# Get list of all CSV files ending with "LSOA" in the data directory
csv_files <- list.files(path = folder.path, pattern = "LSOA\\.csv$", full.names = TRUE)


# Function to extract year and month from file name and read the CSV file
read_and_add_date_info <- function(file) {
  # Extract year and month from the file name (assuming the format is consistent)
  # For example, p_1_07_202304.csv would give year = 2023, month = 04
  date_info <- strsplit(gsub("p_1_07_|\\.csv", "", file), "_")[[1]]
  year <- date_info[[5]]
  month <- date_info[[6]]
  
  # Read the CSV file
  data <- read.csv(file)
  
  # Add year and month columns to the data frame
  data$year <- as.integer(year)
  data$month <- as.character(month)
  
  return(data)
}

# Read in all CSV files and row bind them
combined_data <- bind_rows(lapply(csv_files, read_and_add_date_info))

# Keep required columns and calculate prescribing rate

data.all <- combined_data %>%
  dplyr::select(lsoa11, all_ages, items, year, month) %>%
  rename("lsoa" = "lsoa11",
         "pop.all" = "all_ages",
         "n.items" = "items") %>%
  
  filter(!grepl("^W", lsoa)) %>%
  mutate(r.items = 1000 * n.items/pop.all) %>%
  arrange(year, month)

}


# Read in geocoding data for LSOA and link

geocodes <- read.csv("lsoa_latlong.csv") %>%
  rename("lsoa" = "lsoa11cd") %>%
  dplyr::select(-LSOACD)

data.all <- data.all %>%
  left_join(geocodes, by = "lsoa")

# From https://opendatacommunities.org/resource?uri=http%3A%2F%2Fopendatacommunities.org%2Fdata%2Fsocietal-wellbeing%2Fimd2019%2Findices
# Read in IMD scores and link / calculate quintiles (1 is least deprived)

imd <- read.csv("imd2019lsoa.csv") %>%
  filter(Measurement == "Decile ") %>%
  filter(Indices.of.Deprivation == "a. Index of Multiple Deprivation (IMD)") %>%
  rename("lsoa" = "FeatureCode",
         "imd.decile" = "Value") %>%
  mutate(imd.quintile = case_when(
    imd.decile %in% c(1, 2) ~ 5,
    imd.decile %in% c(3, 4) ~ 4,
    imd.decile %in% c(5, 6) ~ 3,
    imd.decile %in% c(7, 8) ~ 2,
    imd.decile %in% c(9, 10) ~ 1
  )) %>%
  dplyr::select(lsoa, imd.quintile)

data.all <- data.all %>%
  left_join(imd, by = "lsoa")

# Add data on population density by LSOA from https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/lowersuperoutputareapopulationdensity
# mid 2017 estimates
# Density = people per square km

density <- read.csv("lsoapopdensity.csv") %>%
  dplyr::select(lsoa, lsoa.name, pop.density)

data.all <- data.all %>%
  left_join(density, by = "lsoa")

# Number of over 60s and femals data link (data from 2017)
# From https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland

females <- read.csv("popfemales.csv") %>%
  dplyr::select(lsoa, pop.females)

over60 <- read.csv("popage.csv") %>%
  rename("lsoa" = "Area.Codes")

# Get the column names that exist in the dataframe from X0 to X59
existing_columns <- intersect(colnames(over60), paste0("X", 0:59))

# Convert columns to numeric, applying as.numeric to each column
over60[existing_columns] <- lapply(over60[existing_columns], function(x) as.numeric(as.character(x)))

# Remove commas from All.Ages and convert to numeric for the specific column
over60$All.Ages <- as.numeric(gsub(",", "", as.character(over60$All.Ages)))


# Create a new column 'total_sum' that is the sum of columns X0 to X59, and subtract from total to get number of 60+
over60$pop.0to59 <- rowSums(over60[, paste0("X", 0:59)], na.rm = TRUE)
over60$pop.over60 <- over60$All.Ages - over60$pop.0to59

over60 <- over60 %>%
  rename("pop.total" = "All.Ages") %>%
  dplyr::select(lsoa, pop.over60, pop.total)

# link
data.all <- data.all %>%
  left_join(females, by = "lsoa") %>%
  left_join(over60, by = "lsoa")

# Remove commas from pop.females and convert to numeric for the specific column
data.all$pop.females <- as.numeric(gsub(",", "", as.character(data.all$pop.females)))


## Process data -----------------

# Calculate %females and %over60s as per 2017 populations

data.all <- data.all %>%
  mutate(pop.females = 100* pop.females/pop.total) %>%
  mutate(pop.over60 = 100* pop.over60/pop.total) %>%
  dplyr::select(-pop.total)

# Ensure each LSOA has data for every year and month

all_combinations <- data.all %>%
  dplyr::select(lsoa) %>% # Select LSOA column
  distinct() %>%   # Get unique LSOAs
  tidyr::expand(lsoa, year = 2013:2022, month = 1:12) %>% # Create all year and month combinations
  filter(grepl("^E", lsoa))

data.all <- data.all %>%
  mutate(month = as.numeric(month))

data.all <- all_combinations %>%
  left_join(data.all, by = c("lsoa", "year", "month"))

# Check for any NA values
any(is.na(data.all)) # If FALSE, no NA values

# Add linear time trend variable
data.all <- data.all %>%
  group_by(lsoa) %>%                         # Group by LSOA
  arrange(lsoa, year, month) %>%             # Ensure data is sorted by LSOA, year, and month
  mutate(time_trend = row_number()) %>%       # Create the time trend variable
  ungroup()                                  # Ungroup the data

# Add sine and cosine terms to account for seasonality
data.all$sin_month <- sin(2 * pi * data.all$month / 12)  # Sin component
data.all$cos_month <- cos(2 * pi * data.all$month / 12)  # Cosine component


## Regression ------

# Convert imd.quintile to an ordered factor (ordinal)
data.all$imd_quintile <- factor(data.all$imd.quintile, 
                                levels = c(1, 2, 3, 4, 5), 
                                ordered = TRUE)

# Convert to data table
setDT(data.all)

# Subset data to keep only relevant columns
data.subset <- data.all[, .(r.items, imd.quintile, pop.density, pop.females, pop.over60, time_trend, sin_month, cos_month, lsoa, year, month)]

# Remove commas from All.Ages and convert to numeric for the specific column
data.subset$pop.density <- as.numeric(gsub(",", "", as.character(data.subset$pop.density)))

# Aggregate data by lsoa and year to reduce  size
data.agg <- data.subset[, .(
  r.items = sum(r.items),  
  imd.quintile = mean(imd.quintile), # Sum r.items
  pop.density = mean(pop.density),  # Average pop.density
  pop.females = mean(pop.females),  # Average pop.females
  pop.over60 = mean(pop.over60),    # Average pop.over60
  time_trend = ((mean(time_trend) -6.5) /12) + 1 # Time trend
), by = .(lsoa, year)]  # Group by lsoa and year

# Ensure IMD is a factor
data.agg$imd.quintile <- factor(data.agg$imd.quintile, levels = c(1, 2, 3, 4, 5))


# Check for missing data
anyNA(data.agg)

# Check for multicollinearity (if VIF <5 no sigfnificant multicollinearity)
vif(lm(r.items ~ imd.quintile + pop.density + pop.females + pop.over60 + time_trend, data = data.agg))

# Scaling of variables
data.agg$pop.over60_scaled <- data.agg$pop.over60 / 10 # Scale so we look at the PRR of increasing over60s by 10%
data.agg$pop.density_scaled <- data.agg$pop.density / 1000 # Scale to to get PRR increase of increasing pop density by 1000
data.agg$pop.females_scaled <- data.agg$pop.females / 10 # Scale so we look at the PRR of increasing females by 10%

# Fit the negative binomial using glm (not quasipoisson as evidence of heteroscedasticity in resudual vs fitted plot)
# Confirmed no overdispersion (variance of deviance residuals is 0.94, so close to 1)
model.nb <- glm.nb(r.items ~ year + 
                  imd.quintile + 
                  pop.density_scaled + 
                  pop.females_scaled + 
                  pop.over60_scaled +
                  time_trend, 
                data = data.agg)

# Adjust standard errord for clustering by LSOA
clustered_se <- vcovCL(model.nb, cluster = ~ lsoa)

# Display results with cluster-robust standard errors
model_summary <- coeftest(model.nb, vcov = clustered_se)

# Extract coefficients and standard errors
coef <- model_summary[, 1]  # Coefficients (log-IRRs)
se <- model_summary[, 2]    # Robust standard errors
p_value <- model_summary[, 4]  # p-values

# Compute Prescribing Rate Ratios (PRR) by exponentiating the coefficients
PRR <- exp(coef)

# Compute the 95% confidence intervals
lower_ci <- exp(coef - 1.96 * se)  # Lower bound of the CI
upper_ci <- exp(coef + 1.96 * se)  # Upper bound of the CI

# Summary table of the results
irr_results <- data.frame(
  Coefficient = coef,
  PRR = PRR,
  Lower_95_CI = lower_ci,
  Upper_95_CI = upper_ci,
  p_value = p_value
)

# View the IRR and confidence intervals
print(irr_results)


# Save IRR results
irr_results_rounded <- irr_results
irr_results_rounded[, -1] <- round(irr_results_rounded[, -1], 3)  # Round results
irr_results_rounded <- irr_results_rounded %>%
  dplyr::select(-1) %>% # Remove coefficient column
  slice(-1) %>% # Remove intercept results
  mutate(p_value = ifelse(p_value == 0, "<0.001", as.character(p_value)))

write.csv(irr_results_rounded, file = paste0(folder.path," regression results.csv"))

## Plot -----------

# Create the coefficient table with IRRs and confidence intervals
coef_df <- data.frame(
  Predictor = rownames(irr_results),
  Coefficient = irr_results$Coefficient,
  IRR = irr_results$PRR,
  Lower_95_CI = irr_results$Lower_95_CI,
  Upper_95_CI = irr_results$Upper_95_CI
)

# Remove the intercept from the dataframe
coef_df <- coef_df[coef_df$Predictor != "(Intercept)", ]

# Create a mapping of variable names to meaningful labels
predictor_labels <- c(
  "year" = "Year\n(per 1 year increase)",
  "pop.density_scaled" = "Population density\n(per 1000 increase)",
  "pop.females_scaled" = "Percentage of females\n(per 10% increase)",
  "pop.over60_scaled" = "Percentage of population over 60yr\n(per 10% increase)",
  "imd.quintile2"= "IMD quintile 2\n(vs. quintile 1)",
  "imd.quintile3"= "IMD quintile 3\n(vs. quintile 1)",
  "imd.quintile4"= "IMD quintile 4\n(vs. quintile 1)",
  "imd.quintile5"= "IMD quintile 5\n(vs. quintile 1)"
  )

# Replace predictor names with meaningful labels
coef_df$Predictor <- predictor_labels[coef_df$Predictor]

# Reorder the levels of 'Predictor' so that Quintile 2 appears first, then 3, 4, 5, and others
coef_df$Predictor <- factor(coef_df$Predictor, 
                            levels = c("IMD quintile 5\n(vs. quintile 1)",
                                       "IMD quintile 4\n(vs. quintile 1)",
                                       "IMD quintile 3\n(vs. quintile 1)",
                                       "IMD quintile 2\n(vs. quintile 1)",
                                       "Percentage of population over 60yr\n(per 10% increase)",
                                       "Percentage of females\n(per 10% increase)", 
                                       "Population density\n(per 1000 increase)",
                                       "Year\n(per 1 year increase)"
                                       ))


# Plot the coefficients (IRRs) with their confidence intervals
ci.plot <- ggplot(coef_df, aes(x = Predictor, y = IRR)) +
  geom_point(colour = "#090088", size = 1) +
  geom_errorbar(aes(ymin = Lower_95_CI, ymax = Upper_95_CI), width = 0.2, colour = "#090088", linewidth = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black") +
  coord_flip() +  # Flip coordinates to make it horizontal
  labs(title = paste0("Predictor prescribing rate ratios for ", folder.path), y = "\nPRR", x = "Predictors\n") +
  theme_minimal()

# Save plot
ggsave(ci.plot, filename = paste0(folder.path, " PRR plot.png"), width = 8, height = 6, dpi = 600)
