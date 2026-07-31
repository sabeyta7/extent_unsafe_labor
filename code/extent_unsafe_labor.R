##=============================================================================##
##    Geographic Variation in Labor vs. Wage/Hour Violations - Hurdle Model    ##
##=============================================================================##

library(dplyr)
library(ggplot2)
library(sf)
library(INLA)
library(tidyr)
library(spdep)
library(gridExtra)
library(patchwork)
library(rmapshaper)

# Set up directories
current_dir <- getwd()
data_folder <- file.path(current_dir, "data")
code_folder <- file.path(current_dir, "code")
output_folder <- file.path(current_dir, "output")


# Load shapefile
spatial_data <- st_read(file.path(data_folder, "cross_sectional_data_multi_censored_26_update.shp"))

# Clean geometries
is_empty <- st_is_empty(spatial_data)
spatial_data_clean <- spatial_data[!is_empty, ]
spatial_data_clean <- st_make_valid(spatial_data_clean)

# Functions from code folder
source(file.path(code_folder, "functions.R"))

##=======================================##
##   1. PREPARE DATA FOR ALL MODELING    ##
##=======================================##

# Note that industries 99 and 81 are not included in the model due to collinearity

# Basic set up for the data and analysis:

# 1. Subsetting to the relevant variables. Quick couple notes on the variable names:
#    - `census_pop` and `census_wor` are the population and workforce counts from the census.
#    - `vio_overal`, `wh_overa_5`, `acc_overal`, and `insp_overa` are raw counts of violations, wage/hour, accidents, and inspections.
#       -Sub note: `wh_overa_5` is the wage/hour violations capped at the 99th percentile (max 85 for a given zip year). Other versions include the raw count which is just `wh_overal`, 
#       and `wh_overa_3` which is the Jenks natural breaks that caps at 51 for given zip year. Results with all three are nearly identical throughout with
#       slight differences in the coefficients.
#       -The two outocmes of wh_rate_tw (censoring based on employee and violations count) and wh_rate_ra (censoring based on employee and violation rate) are the two sensitivity versions
#    - `vio_rate`, `wh_rate_pe`, and `acc_rate` are the rates of violations, wage/hour, and accidents per establishment.
#    - The industry variables are industries per 1000 workers. I know the naming looks strange and inconsistent, but you are free to 
#      check. R cuts variable names from python if they are too long, hence why also wh_3, and 4 indicators
#    - The geographic identifiers include state and region.
#    - The inspection rate is also by establishment.

# 2. Making the spatial neighborhood structure for the eventual spatial modeling. As we have a binary component and a continuous components,
#    we make subsets of the data for each continuous version. In other words, each subset are the zip codes in which each outcome has more than 0,
#    hence the different subsets. The neighborhoods are made via a queen contiguity structure, and then we use the `knearneigh` function to fill in any islands.
#    Here we chose 3 nearest neighbors and made the relationship symetric. The neighborhood structure is then saved as an INLA graph file for use in the models at the end.
#       -Side tangent, it is only because of INLA that we had to use R in the first place. Woulda perferred just outright pythons, but what are you going to do besides complain
#        in the notes of your code?

# 3. Scaling the inspection rate variable in the overall data and all positive datasets. This is done to standardize the inspection rate across the datasets.


# Prepare data for modeling
model_data <- spatial_data_clean %>%
  dplyr::select(

    # Population vairables
    census_pop, census_wor,

    # Raw counts
    vio_overal, wh_overa_5, acc_overal, insp_overa, 

    # Response variables
    vio_rate, wh_rate_pe, acc_rate,
    
    # Industry composition
    ESTAB_11_P, ESTAB_21_P, ESTAB_22_P, ESTAB_23_P, ESTAB_31_1, 
    ESTAB_42_P, ESTAB_44_1, ESTAB_48_1, ESTAB_51_P, ESTAB_52_P, 
    ESTAB_53_P, ESTAB_54_P, ESTAB_55_P, ESTAB_56_P, ESTAB_61_P, 
    ESTAB_62_P, ESTAB_71_P, ESTAB_72_P, ESTAB_81_P, ESTAB_TOTA,

    # Geographic identifiers
    state, region_8,

    # Inspection rate
    insp_rate,

    # Geometry
    geometry

    # Sensitivity variables

    # # Natural breaks
    # wh_overa_4, wh_rate_na,
    # # Two way cap
    # wh_overa_6, wh_rate__1
  ) %>%
  na.omit()

# Create binary outcome variables (zero vs. non-zero) for hurdle models 
model_data$labor_binary <- as.numeric(model_data$vio_rate > 0)
model_data$wh_binary <- as.numeric(model_data$wh_rate_pe > 0)
model_data$acc_binary <- as.numeric(model_data$acc_rate > 0)

# Create positive-only datasets for the continuous component of the hurdle model
labor_positive <- model_data %>% filter(vio_rate > 0)
wh_positive <- model_data %>% filter(wh_rate_pe > 0)
acc_positive <- model_data %>% filter(acc_rate > 0)

# Zip code identifiers
model_data$zip_id <- 1:nrow(model_data)
labor_positive$zip_id_pos <- 1:nrow(labor_positive)
wh_positive$zip_id_pos <- 1:nrow(wh_positive)
acc_positive$zip_id_pos <- 1:nrow(acc_positive)

# Factor states and regions
model_data$state <- factor(model_data$state)
model_data$region_8 <- factor(model_data$region_8)
labor_positive$state <- factor(labor_positive$state)
wh_positive$state <- factor(wh_positive$state)
acc_positive$state <- factor(acc_positive$state)

# Create neighborhood structure for binary part
nb_queen <- poly2nb(model_data, queen = TRUE)
island_indices <- which(card(nb_queen) == 0)
if(length(island_indices) > 0) {
  coords <- st_coordinates(st_centroid(st_geometry(model_data)))
  knn <- knearneigh(coords, k = 3)
  knn_nb <- knn2nb(knn)
  
  for(i in island_indices) {
    nb_queen[[i]] <- knn_nb[[i]]
  }
  nb_queen <- make.sym.nb(nb_queen)
}

# Create INLA graph for binary part
nb2INLA(file.path(data_folder, "model_graph_binary"), nb_queen)
model_adj_binary <- paste(file.path(data_folder, "model_graph_binary"))

# For labor positive dataset (continuous part)
nb_queen_labor <- poly2nb(labor_positive, queen = TRUE)
island_indices_labor <- which(card(nb_queen_labor) == 0)
if(length(island_indices_labor) > 0) {
  coords_labor <- st_coordinates(st_centroid(st_geometry(labor_positive)))
  knn_labor <- knearneigh(coords_labor, k = 3)
  knn_nb_labor <- knn2nb(knn_labor)
  
  for(i in island_indices_labor) {
    nb_queen_labor[[i]] <- knn_nb_labor[[i]]
  }
  nb_queen_labor <- make.sym.nb(nb_queen_labor)
}

# Create INLA graph for labor positive
nb2INLA(file.path(data_folder, "model_graph_labor_pos"), nb_queen_labor)
model_adj_labor_pos <- paste(file.path(data_folder, "model_graph_labor_pos"))

# For wage/hour positive dataset
nb_queen_wh <- poly2nb(wh_positive, queen = TRUE)
island_indices_wh <- which(card(nb_queen_wh) == 0)
if(length(island_indices_wh) > 0) {
  coords_wh <- st_coordinates(st_centroid(st_geometry(wh_positive)))
  knn_wh <- knearneigh(coords_wh, k = 3)
  knn_nb_wh <- knn2nb(knn_wh)
  
  for(i in island_indices_wh) {
    nb_queen_wh[[i]] <- knn_nb_wh[[i]]
  }
  nb_queen_wh <- make.sym.nb(nb_queen_wh)
}

# Create INLA graph for wh positive
nb2INLA(file.path(data_folder, "model_graph_wh_pos"), nb_queen_wh)
model_adj_wh_pos <- paste(file.path(data_folder, "model_graph_wh_pos"))

# For accident rate positive dataset
nb_queen_acc <- poly2nb(acc_positive, queen = TRUE)
island_indices_acc <- which(card(nb_queen_acc) == 0)
if(length(island_indices_acc) > 0) {
  coords_acc <- st_coordinates(st_centroid(st_geometry(acc_positive)))
  knn_acc <- knearneigh(coords_acc, k = 3)
  knn_nb_acc <- knn2nb(knn_acc)
  
  for(i in island_indices_acc) {
    nb_queen_acc[[i]] <- knn_nb_acc[[i]]
  }
  nb_queen_acc <- make.sym.nb(nb_queen_acc)
}

# Create INLA graph for accident rate positive
nb2INLA(file.path(data_folder, "model_graph_acc_pos"), nb_queen_acc)
model_adj_acc_pos <- paste(file.path(data_folder, "model_graph_acc_pos"))

# Scale inspection rate in the overall data and all positive data sets
if("insp_rate" %in% colnames(model_data)) {
  model_data$insp_rate <- scale(model_data$insp_rate)
  labor_positive$insp_rate <- scale(labor_positive$insp_rate)
  wh_positive$insp_rate <- scale(wh_positive$insp_rate)
  acc_positive$insp_rate <- scale(acc_positive$insp_rate)
}

##=====================================##
##           2.Descvriptives           ##
##=====================================##
# Some packages screw with the creation of the neighborhood structure, so we need to load them after
library(viridis)
library(RColorBrewer)
library(factoextra) 
library(corrplot)   
library(textshape)
library(ggrepel)
library(tidytable)
library(ggridges)
library(knitr)
library(kableExtra)
library(pROC)
library(spatstat)
library(DCluster)
library(survey)
library(cowplot)
library(tidyverse)

# Basic descriptives. Pulling out summaries for the industry variable, plotting histograms, and creating a density plot for the most common industries

industry_vars <- grep("^ESTAB_(.*_P|.*_1)$", names(model_data), value = TRUE)

industry_data_long <- model_data %>%
  st_drop_geometry() %>%
  dplyr::select(all_of(industry_vars)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "industry_var",
    values_to = "percentage"
  ) %>%
  mutate(
    industry_name = case_when(
      industry_var == "ESTAB_11_P" ~ "Agriculture, Forestry, Fishing",
      industry_var == "ESTAB_21_P" ~ "Mining & Extraction",
      industry_var == "ESTAB_22_P" ~ "Utilities",
      industry_var == "ESTAB_23_P" ~ "Construction",
      industry_var == "ESTAB_31_1" ~ "Manufacturing",
      industry_var == "ESTAB_42_P" ~ "Wholesale Trade",
      industry_var == "ESTAB_44_1" ~ "Retail Trade",
      industry_var == "ESTAB_48_1" ~ "Transportation & Warehousing",
      industry_var == "ESTAB_51_P" ~ "Information",
      industry_var == "ESTAB_52_P" ~ "Finance & Insurance",
      industry_var == "ESTAB_53_P" ~ "Real Estate",
      industry_var == "ESTAB_54_P" ~ "Professional Services",
      industry_var == "ESTAB_55_P" ~ "Management",
      industry_var == "ESTAB_56_P" ~ "Administrative Services",
      industry_var == "ESTAB_61_P" ~ "Education",
      industry_var == "ESTAB_62_P" ~ "Healthcare",
      industry_var == "ESTAB_71_P" ~ "Arts & Entertainment",
      industry_var == "ESTAB_72_P" ~ "Accommodation & Food",
      industry_var == "ESTAB_81_P" ~ "Other Services",
      TRUE ~ industry_var
    )
  )

regional_descriptives <- model_data %>%
  st_drop_geometry() %>%
  group_by(region_8) %>%
  summarise(
    n_zip_codes = n(),
    
    total_violations = sum(vio_overal, na.rm = TRUE),
    total_wh_violations = sum(wh_overa_5, na.rm = TRUE),
    total_accidents = sum(acc_overal, na.rm = TRUE),
    total_inspections = sum(insp_overa, na.rm = TRUE),
    
    total_establishments = sum(ESTAB_TOTA, na.rm = TRUE),
    total_population = sum(census_pop, na.rm = TRUE),
    total_workers = sum(census_wor, na.rm = TRUE),
    
    agriculture_per_1k = weighted.mean(ESTAB_11_P, census_wor, na.rm = TRUE),
    mining_per_1k = weighted.mean(ESTAB_21_P, census_wor, na.rm = TRUE),
    utilities_per_1k = weighted.mean(ESTAB_22_P, census_wor, na.rm = TRUE),
    construction_per_1k = weighted.mean(ESTAB_23_P, census_wor, na.rm = TRUE),
    manufacturing_per_1k = weighted.mean(ESTAB_31_1, census_wor, na.rm = TRUE),
    wholesale_per_1k = weighted.mean(ESTAB_42_P, census_wor, na.rm = TRUE),
    retail_per_1k = weighted.mean(ESTAB_44_1, census_wor, na.rm = TRUE),
    transportation_per_1k = weighted.mean(ESTAB_48_1, census_wor, na.rm = TRUE),
    information_per_1k = weighted.mean(ESTAB_51_P, census_wor, na.rm = TRUE),
    finance_per_1k = weighted.mean(ESTAB_52_P, census_wor, na.rm = TRUE),
    real_estate_per_1k = weighted.mean(ESTAB_53_P, census_wor, na.rm = TRUE),
    professional_per_1k = weighted.mean(ESTAB_54_P, census_wor, na.rm = TRUE),
    management_per_1k = weighted.mean(ESTAB_55_P, census_wor, na.rm = TRUE),
    administrative_per_1k = weighted.mean(ESTAB_56_P, census_wor, na.rm = TRUE),
    education_per_1k = weighted.mean(ESTAB_61_P, census_wor, na.rm = TRUE),
    healthcare_per_1k = weighted.mean(ESTAB_62_P, census_wor, na.rm = TRUE),
    arts_per_1k = weighted.mean(ESTAB_71_P, census_wor, na.rm = TRUE),
    accommodation_per_1k = weighted.mean(ESTAB_72_P, census_wor, na.rm = TRUE),
    other_services_per_1k = weighted.mean(ESTAB_81_P, census_wor, na.rm = TRUE),
    avg_violation_rate = mean(vio_rate, na.rm = TRUE),
    avg_wh_rate = mean(wh_rate_pe, na.rm = TRUE),
    avg_accident_rate = mean(acc_rate, na.rm = TRUE),
    avg_inspection_rate = mean(insp_rate, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    regional_violation_rate = total_violations / total_establishments * 100,
    regional_wh_rate = total_wh_violations / total_establishments * 100,
    regional_accident_rate = total_accidents / total_establishments * 100,
    regional_inspection_rate = total_inspections / total_establishments * 100
  ) %>%
  arrange(desc(total_violations))


us_total <- model_data %>%
  st_drop_geometry() %>%
  summarise(
    region_8 = "US_TOTAL",
    n_zip_codes = n(),
    
    total_violations = sum(vio_overal, na.rm = TRUE),
    total_wh_violations = sum(wh_overa_5, na.rm = TRUE),
    total_accidents = sum(acc_overal, na.rm = TRUE),
    total_inspections = sum(insp_overa, na.rm = TRUE),
    
    total_establishments = sum(ESTAB_TOTA, na.rm = TRUE),
    total_population = sum(census_pop, na.rm = TRUE),
    total_workers = sum(census_wor, na.rm = TRUE),
    
    # Calculate industry percentages directly from raw data
    agriculture_per_1k = weighted.mean(ESTAB_11_P, census_wor, na.rm = TRUE),
    mining_per_1k = weighted.mean(ESTAB_21_P, census_wor, na.rm = TRUE),
    utilities_per_1k = weighted.mean(ESTAB_22_P, census_wor, na.rm = TRUE),
    construction_per_1k = weighted.mean(ESTAB_23_P, census_wor, na.rm = TRUE),
    manufacturing_per_1k = weighted.mean(ESTAB_31_1, census_wor, na.rm = TRUE),
    wholesale_per_1k = weighted.mean(ESTAB_42_P, census_wor, na.rm = TRUE),
    retail_per_1k = weighted.mean(ESTAB_44_1, census_wor, na.rm = TRUE),
    transportation_per_1k = weighted.mean(ESTAB_48_1, census_wor, na.rm = TRUE),
    information_per_1k = weighted.mean(ESTAB_51_P, census_wor, na.rm = TRUE),
    finance_per_1k = weighted.mean(ESTAB_52_P, census_wor, na.rm = TRUE),
    real_estate_per_1k = weighted.mean(ESTAB_53_P, census_wor, na.rm = TRUE),
    professional_per_1k = weighted.mean(ESTAB_54_P, census_wor, na.rm = TRUE),
    management_per_1k = weighted.mean(ESTAB_55_P, census_wor, na.rm = TRUE),
    administrative_per_1k = weighted.mean(ESTAB_56_P, census_wor, na.rm = TRUE),
    education_per_1k = weighted.mean(ESTAB_61_P, census_wor, na.rm = TRUE),
    healthcare_per_1k = weighted.mean(ESTAB_62_P, census_wor, na.rm = TRUE),
    arts_per_1k = weighted.mean(ESTAB_71_P, census_wor, na.rm = TRUE),
    accommodation_per_1k = weighted.mean(ESTAB_72_P, census_wor, na.rm = TRUE),
    other_services_per_1k = weighted.mean(ESTAB_81_P, census_wor, na.rm = TRUE)

  ) %>%
  mutate(
    regional_violation_rate = total_violations / total_establishments * 100,
    regional_wh_rate = total_wh_violations / total_establishments * 100,
    regional_accident_rate = total_accidents / total_establishments * 100,
    regional_inspection_rate = total_inspections / total_establishments * 100
  )

regional_with_us <- bind_rows(regional_descriptives, us_total)

print(regional_with_us, n = Inf, width = Inf)

reduced_descriptives <- regional_with_us %>%
  dplyr::select(
    region_8,
    n_zip_codes,
    agriculture_per_1k, mining_per_1k, utilities_per_1k, construction_per_1k,
    manufacturing_per_1k, wholesale_per_1k, retail_per_1k, transportation_per_1k,
    information_per_1k, finance_per_1k, real_estate_per_1k, professional_per_1k,
    management_per_1k, administrative_per_1k, education_per_1k, healthcare_per_1k,
    arts_per_1k, accommodation_per_1k, other_services_per_1k,
    regional_violation_rate, regional_wh_rate, regional_accident_rate, regional_inspection_rate
  ) %>%
  mutate(across(where(is.numeric), ~ round(., 2))) %>%
  rename(
    Region = region_8,
    `ZIP Codes` = n_zip_codes,
    Agriculture = agriculture_per_1k,
    Mining = mining_per_1k,
    Utilities = utilities_per_1k,
    Construction = construction_per_1k,
    Manufacturing = manufacturing_per_1k,
    Wholesale = wholesale_per_1k,
    Retail = retail_per_1k,
    Transportation = transportation_per_1k,
    Information = information_per_1k,
    Finance = finance_per_1k,
    `Real Estate` = real_estate_per_1k,
    Professional = professional_per_1k,
    Management = management_per_1k,
    Administrative = administrative_per_1k,
    Education = education_per_1k,
    Healthcare = healthcare_per_1k,
    Arts = arts_per_1k,
    Accommodation = accommodation_per_1k,
    `Other Services` = other_services_per_1k,
    `Violation Rate` = regional_violation_rate,
    `WH Rate` = regional_wh_rate,
    `Accident Rate` = regional_accident_rate,
    `Inspection Rate` = regional_inspection_rate
  ) %>%
  # Reorder so US_TOTAL comes first, then regions by your preferred order
  mutate(Region = factor(Region, levels = c(
    "US_TOTAL", "Pacific", "Southeast", "Great Lakes",
    "Mid-Atlantic", "Southwest", "Plains", "Mountain West", "New England"
  ))) %>%
  arrange(Region) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Metric")

# Fix column names to use region names instead of index numbers
colnames(reduced_descriptives) <- c("Metric", as.character(reduced_descriptives[1, -1]))
reduced_descriptives <- reduced_descriptives[-1, ]  # remove the Region row

print(reduced_descriptives)

write.csv(reduced_descriptives, file.path(output_folder, "regional_descriptives.csv"), row.names = FALSE)

##==================================================##
##            4. RELATIVE RISK ANALYSIS             ##
##==================================================##


# This section calculates relative risk measures for industry-region combinations. In other words, this section does several things:
# 1. For each region and each industry, it calculates the overall rates of OSHA violations, wage/hour violations, and accident rates.
# 2. It computes correlations between industry concentration and the binary outcomes (presence/absence of violations) as well as the continuous outcomes (severity rates) for each region-industry combination.
# 3. It also calculates the average severity rates for ZIP codes with above-median presence of each industry within that region.
# 4. Finally, it compiles all these metrics into a comprehensive data frame for further analysis or visualization.

# This seciton is ultimatly aiming to destill a set of largely desciriptive statistics that plainly indicate regional distinctions. To this end regional profiles
# largely encompass three main ideas:
# 1. Overall risk: How does the region compare to the national average in terms of likelihood and severity of violations/accidents?
# 2. Industry-specific risk: Within each region, how do zip codes with high concentrations of specific industries compare to the regional average?
# 3. Correlation analysis: What is the strength and direction of the relationship between industry concentration and both the likelihood and severity of violations/accidents?

# We can then also indicate that in a given region there is either consistency. For instnace, we can see if in a given region, one industry is consistently high risk across all three outcomes (for both presences and volume).



regions <- unique(model_data[["region_8"]])
results <- list()

industry_cols <- grep("^ESTAB_.*_P$|^ESTAB_.*_1$", names(model_data), value = TRUE)

for(region in regions) {
  region_data <- model_data[model_data[["region_8"]] == region, ]
  labor_positive_region <- labor_positive[labor_positive[["region_8"]] == region, ]
  wh_positive_region <- wh_positive[wh_positive[["region_8"]] == region, ]
  acc_positive_region <- acc_positive[acc_positive[["region_8"]] == region, ]
  
  national_labor_rate <- mean(model_data$labor_binary, na.rm = TRUE)
  national_wh_rate <- mean(model_data$wh_binary, na.rm = TRUE)
  national_acc_rate <- mean(model_data$acc_binary, na.rm = TRUE)
  
  national_labor_volume <- mean(labor_positive$vio_rate, na.rm = TRUE)
  national_wh_volume <- mean(wh_positive$wh_rate_pe, na.rm = TRUE)
  national_acc_volume <- mean(acc_positive$acc_rate, na.rm = TRUE)

  # Overall regional rates
  regional_labor_rate_overall <- mean(region_data$labor_binary, na.rm = TRUE)
  regional_wh_rate_overall <- mean(region_data$wh_binary, na.rm = TRUE)
  regional_acc_rate_overall <- mean(region_data$acc_binary, na.rm = TRUE)

  regional_labor_volume_overall <- mean(labor_positive_region$vio_rate, na.rm = TRUE)
  regional_wh_volume_overall <- mean(wh_positive_region$wh_rate_pe, na.rm = TRUE)
  regional_acc_volume_overall <- mean(acc_positive_region$acc_rate, na.rm = TRUE)
  
  for(industry in industry_cols) {
    ind_name <- gsub("ESTAB_|_P|_1", "", industry)
    
    # Binary correlations
    labor_binary_cor <- cor(region_data[[industry]], region_data$labor_binary, use = "pairwise.complete.obs")
    wh_binary_cor <- cor(region_data[[industry]], region_data$wh_binary, use = "pairwise.complete.obs")
    acc_binary_cor <- cor(region_data[[industry]], region_data$acc_binary, use = "pairwise.complete.obs")
    
    # Continuous correlations
    labor_cont_cor <- NA
    wh_cont_cor <- NA
    acc_cont_cor <- NA
    
    if(nrow(labor_positive_region) > 5) {
      labor_cont_cor <- cor(labor_positive_region[[industry]], labor_positive_region$vio_rate, use = "pairwise.complete.obs")
    }
    
    if(nrow(wh_positive_region) > 5) {
      wh_cont_cor <- cor(wh_positive_region[[industry]], wh_positive_region$wh_rate_pe, use = "pairwise.complete.obs")
    }
    
    if(nrow(acc_positive_region) > 5) {
      acc_cont_cor <- cor(acc_positive_region[[industry]], acc_positive_region$acc_rate, use = "pairwise.complete.obs")
    }
    
    results[[paste(region, ind_name, sep = "_")]] <- data.frame(
      region = region,
      industry = ind_name,
      industry_col = industry,
      
      # Overall regional rates
      regional_labor_rate_overall = regional_labor_rate_overall,
      regional_wh_rate_overall = regional_wh_rate_overall,
      regional_acc_rate_overall = regional_acc_rate_overall,
      
      # Overall regional volumes
      regional_labor_volume_overall = regional_labor_volume_overall,
      regional_wh_volume_overall = regional_wh_volume_overall,
      regional_acc_volume_overall = regional_acc_volume_overall,
      
      # Binary correlations
      labor_binary_correlation = labor_binary_cor,
      wh_binary_correlation = wh_binary_cor,
      acc_binary_correlation = acc_binary_cor,
      
      # Continuous correlations
      labor_cont_correlation = labor_cont_cor,
      wh_cont_correlation = wh_cont_cor,
      acc_cont_correlation = acc_cont_cor,
      
      # Relative risks compared to national rates
      relative_labor_risk = regional_labor_rate_overall / national_labor_rate,
      relative_wh_risk = regional_wh_rate_overall / national_wh_rate,
      relative_acc_risk = regional_acc_rate_overall / national_acc_rate,

      # Relative volumes compared to national volumes
      relative_labor_volume = regional_labor_volume_overall / national_labor_volume,
      relative_wh_volume = regional_wh_volume_overall / national_wh_volume,
      relative_acc_volume = regional_acc_volume_overall / national_acc_volume
    )
  }
}

relative_risk <- bind_rows(results)

regional_profiles <- relative_risk %>%
  group_by(region) %>%
  mutate(
    region_overall_binary_risk = mean(c(relative_labor_risk, relative_wh_risk, relative_acc_risk), na.rm = TRUE),
    region_overall_cont_risk = mean(c(relative_labor_volume, relative_wh_volume, relative_acc_volume), na.rm = TRUE)
  ) %>%
  arrange(desc(labor_binary_correlation)) %>%
  slice_head(n = 5) %>%
  dplyr::select(region, 
         region_overall_binary_risk, 
         region_overall_cont_risk,
         industry, 
         # Binary correlations
         labor_binary_correlation, 
         wh_binary_correlation, 
         acc_binary_correlation,
         # Continuous correlations
         labor_cont_correlation, 
         wh_cont_correlation, 
         acc_cont_correlation)

# Comprehensive region typology
region_typology <- data.frame()
for (reg in unique(relative_risk$region)) {
  region_data <- filter(relative_risk, region == reg)
  
  if (nrow(region_data) > 0) {
    # Binary correlation indices - finding the top industry for each outcome
    labor_binary_idx <- which.max(region_data$labor_binary_correlation)
    wh_binary_idx <- which.max(region_data$wh_binary_correlation)
    acc_binary_idx <- which.max(region_data$acc_binary_correlation)
    
    # Continuous correlation indices - finding the top industry for each outcome
    labor_cont_idx <- which.max(region_data$labor_cont_correlation)
    wh_cont_idx <- which.max(region_data$wh_cont_correlation)
    acc_cont_idx <- which.max(region_data$acc_cont_correlation)
    
    # Top industries for binary correlations
    top_labor_binary <- if (length(labor_binary_idx) > 0 && !is.na(labor_binary_idx)) region_data$industry[labor_binary_idx] else NA
    top_wh_binary <- if (length(wh_binary_idx) > 0 && !is.na(wh_binary_idx)) region_data$industry[wh_binary_idx] else NA
    top_acc_binary <- if (length(acc_binary_idx) > 0 && !is.na(acc_binary_idx)) region_data$industry[acc_binary_idx] else NA
    
    # Top industries for continuous correlations
    top_labor_cont <- if (length(labor_cont_idx) > 0 && !is.na(labor_cont_idx)) region_data$industry[labor_cont_idx] else NA
    top_wh_cont <- if (length(wh_cont_idx) > 0 && !is.na(wh_cont_idx)) region_data$industry[wh_cont_idx] else NA
    top_acc_cont <- if (length(acc_cont_idx) > 0 && !is.na(acc_cont_idx)) region_data$industry[acc_cont_idx] else NA
    
    # Binary overlap checks
    labor_wh_binary_same <- !is.na(top_labor_binary) && !is.na(top_wh_binary) && top_labor_binary == top_wh_binary
    labor_acc_binary_same <- !is.na(top_labor_binary) && !is.na(top_acc_binary) && top_labor_binary == top_acc_binary
    wh_acc_binary_same <- !is.na(top_wh_binary) && !is.na(top_acc_binary) && top_wh_binary == top_acc_binary
    all_binary_same <- labor_wh_binary_same && labor_acc_binary_same
    
    # Continuous overlap checks
    labor_wh_cont_same <- !is.na(top_labor_cont) && !is.na(top_wh_cont) && top_labor_cont == top_wh_cont
    labor_acc_cont_same <- !is.na(top_labor_cont) && !is.na(top_acc_cont) && top_labor_cont == top_acc_cont
    wh_acc_cont_same <- !is.na(top_wh_cont) && !is.na(top_acc_cont) && top_wh_cont == top_acc_cont
    all_cont_same <- labor_wh_cont_same && labor_acc_cont_same
    
    # Check overlap between binary and continuous
    labor_binary_cont_same <- !is.na(top_labor_binary) && !is.na(top_labor_cont) && top_labor_binary == top_labor_cont
    wh_binary_cont_same <- !is.na(top_wh_binary) && !is.na(top_wh_cont) && top_wh_binary == top_wh_cont
    acc_binary_cont_same <- !is.na(top_acc_binary) && !is.na(top_acc_cont) && top_acc_binary == top_acc_cont
        
    # Binary risk statistics
    mean_labor_binary_risk <- mean(region_data$relative_labor_risk, na.rm = TRUE)
    mean_wh_binary_risk <- mean(region_data$relative_wh_risk, na.rm = TRUE)
    mean_acc_binary_risk <- mean(region_data$relative_acc_risk, na.rm = TRUE)
    
    # Continuous risk statistics
    mean_labor_cont_risk <- mean(region_data$relative_labor_volume, na.rm = TRUE)
    mean_wh_cont_risk <- mean(region_data$relative_wh_volume, na.rm = TRUE)
    mean_acc_cont_risk <- mean(region_data$relative_acc_volume, na.rm = TRUE)
    
    # Overall risk measures
    overall_binary_risk <- mean(c(mean_labor_binary_risk, mean_wh_binary_risk, mean_acc_binary_risk), na.rm = TRUE)
    overall_cont_risk <- mean(c(mean_labor_cont_risk, mean_wh_cont_risk, mean_acc_cont_risk), na.rm = TRUE)
    
    # Create the result row
    result_row <- data.frame(
      region = reg,
      
      # Top industries for binary correlations
      top_labor_binary_industry = top_labor_binary,
      top_wh_binary_industry = top_wh_binary,
      top_acc_binary_industry = top_acc_binary,
      
      # Top industries for continuous correlations
      top_labor_cont_industry = top_labor_cont,
      top_wh_cont_industry = top_wh_cont,
      top_acc_cont_industry = top_acc_cont,
      
      # Binary overlap checks
      labor_wh_binary_same = labor_wh_binary_same,
      labor_acc_binary_same = labor_acc_binary_same,
      wh_acc_binary_same = wh_acc_binary_same,
      all_binary_same = all_binary_same,
      
      # Continuous overlap checks
      labor_wh_cont_same = labor_wh_cont_same,
      labor_acc_cont_same = labor_acc_cont_same,
      wh_acc_cont_same = wh_acc_cont_same,
      all_cont_same = all_cont_same,
      
      # Binary-continuous overlap
      labor_binary_cont_same = labor_binary_cont_same,
      wh_binary_cont_same = wh_binary_cont_same,
      acc_binary_cont_same = acc_binary_cont_same,
      
      # Overall risk measures
      overall_binary_risk = overall_binary_risk,
      overall_cont_risk = overall_cont_risk
    )
    
    # Append to results
    region_typology <- rbind(region_typology, result_row)
  }
}

# Create a summary table of the typology results
region_summary <- region_typology %>%
  dplyr::select(region, 
         top_labor_binary_industry, top_wh_binary_industry, top_acc_binary_industry,
         top_labor_cont_industry, top_wh_cont_industry, top_acc_cont_industry,
         labor_binary_cont_same, wh_binary_cont_same, acc_binary_cont_same,
         all_binary_same, all_cont_same,
         overall_binary_risk, overall_cont_risk) %>%
  arrange(desc(overall_binary_risk))

# Print the main typology table
print(region_typology)

# Print the summary table
print(region_summary, n = Inf, width = Inf)


##===============================================##
##     6. PREVALENCE ESTIMATION OF OUTCOMES      ##
##===============================================##


# For the prevalence side of the analysis we are using standardized mortality ratios. As this is administrative
# data, there is a ton of messiness and likley a lot of hidden bias. We are construcing basic rates of the three outcomes
# by the population, working population, and establishments. We then calculate an expected rate based on the national rate (so the number of establishments by the national rate) 
# which takes into account regional varaition and regulatory differences. Though not perfect and a limitation of the analysis as a true
# "expected" value is hard to quantify without very intense modeling. SMR is then just the observed rate divided by the expected rate. 

# However, gine the rarety of some of the typs of violations and potential division by numbers between 0 and 1 due to averaging, we are both limiting the accounted for zip codes as those with expected counts of at least 1,
# and also smoothing values via Empirical Bayes smoothing. I guess more precicely, we are using the Poisson-Gamma model, in which SMR calculation includes weights for the "prior" count (our alpha), and
# the prior expectation (the beta). I followed Acosta, L. A., Grant Morrison, Angela Li, Karina. (n.d.). Chapter 5 Rate Mapping to work through the calculation.
# Confidence intervalues come via bootstrapping. The process is pretty simple, we resample with replacement between 1 and a theoretical smoothed mean and recalculated EB SMRs and then take the 2.5th and 97.5th percentiles as the bounds.


#-------Getting the basic rates and national rates
model_data <- model_data %>%
  mutate(
    vio_per_100_estab = ifelse(ESTAB_TOTA > 0, (vio_overal / ESTAB_TOTA) * 100, NA),
    wh_per_100_estab = ifelse(ESTAB_TOTA > 0, (wh_overa_5 / ESTAB_TOTA) * 100, NA),
    acc_per_100_estab = ifelse(ESTAB_TOTA > 0, (acc_overal / ESTAB_TOTA) * 100, NA)
  )

national_rates <- model_data %>%
  filter(!is.na(vio_per_100_estab)) %>%
  summarise(
    national_vio_per_100_estab = weighted.mean(vio_per_100_estab, ESTAB_TOTA, na.rm = TRUE),
    national_wh_per_100_estab = weighted.mean(wh_per_100_estab, ESTAB_TOTA, na.rm = TRUE),
    national_acc_per_100_estab = weighted.mean(acc_per_100_estab, ESTAB_TOTA, na.rm = TRUE)
)
print(national_rates)

#-----Expected rates, SMRs, and CIs
model_data <- model_data %>%
  mutate(
    expected_vio_estab = (ESTAB_TOTA * national_rates$national_vio_per_100_estab) / 100,
    expected_wh_estab = (ESTAB_TOTA * national_rates$national_wh_per_100_estab) / 100,
    expected_acc_estab = (ESTAB_TOTA * national_rates$national_acc_per_100_estab) / 100,
    SMR_vio_estab = ifelse(expected_vio_estab > 0, vio_overal / expected_vio_estab, NA),
    SMR_wh_estab = ifelse(expected_wh_estab > 0, wh_overa_5 / expected_wh_estab, NA),
    SMR_acc_estab = ifelse(expected_acc_estab > 0, acc_overal / expected_acc_estab, NA)
  )
summary(as.data.frame(model_data)[c("SMR_vio_estab", "SMR_wh_estab", "SMR_acc_estab" )])

#-------Fixed Empirical Bayes SMRs and Confidence Intervals

#--Global mean
global_vio_SMR_estab <- sum(model_data$vio_overal, na.rm = TRUE) / sum(model_data$expected_vio_estab, na.rm = TRUE)
global_wh_SMR_estab <- sum(model_data$wh_overa_5, na.rm = TRUE) / sum(model_data$expected_wh_estab, na.rm = TRUE)
global_acc_SMR_estab <- sum(model_data$acc_overal, na.rm = TRUE) / sum(model_data$expected_acc_estab, na.rm = TRUE)

#--Variance
n_vio_estab <- sum(!is.na(model_data$SMR_vio_estab))
numerator_vio_estab <- sum(model_data$expected_vio_estab * (model_data$SMR_vio_estab - global_vio_SMR_estab)^2, na.rm = TRUE)
denominator_vio_estab <- sum(model_data$expected_vio_estab, na.rm = TRUE) - global_vio_SMR_estab * (sum(model_data$expected_vio_estab, na.rm = TRUE) / n_vio_estab)
variance_vio_estab <- max(numerator_vio_estab / denominator_vio_estab, 0)
n_wh_estab <- sum(!is.na(model_data$SMR_wh_estab))
numerator_wh_estab <- sum(model_data$expected_wh_estab * (model_data$SMR_wh_estab - global_wh_SMR_estab)^2, na.rm = TRUE)
denominator_wh_estab <- sum(model_data$expected_wh_estab, na.rm = TRUE) - global_wh_SMR_estab * (sum(model_data$expected_wh_estab, na.rm = TRUE) / n_wh_estab)
variance_wh_estab <- max(numerator_wh_estab / denominator_wh_estab, 0)
n_acc_estab <- sum(!is.na(model_data$SMR_acc_estab))
numerator_acc_estab <- sum(model_data$expected_acc_estab * (model_data$SMR_acc_estab - global_acc_SMR_estab)^2, na.rm = TRUE)
denominator_acc_estab <- sum(model_data$expected_acc_estab, na.rm = TRUE) - global_acc_SMR_estab * (sum(model_data$expected_acc_estab, na.rm = TRUE) / n_acc_estab)
variance_acc_estab <- max(numerator_acc_estab / denominator_acc_estab, 0)

#--Alpha parameter
alpha_vio_estab <- global_vio_SMR_estab^2 / variance_vio_estab 
alpha_wh_estab <- global_wh_SMR_estab^2 / variance_wh_estab
alpha_acc_estab <- global_acc_SMR_estab^2 / variance_acc_estab

#--Beta parameter
beta_vio_estab <- global_vio_SMR_estab / variance_vio_estab
beta_wh_estab <- global_wh_SMR_estab / variance_wh_estab
beta_acc_estab <- global_acc_SMR_estab / variance_acc_estab

filtered_smr_data <- model_data %>%
  filter(expected_vio_estab >= 1.0 | expected_wh_estab >= 1.0 | expected_acc_estab >= 1.0) %>% 
  mutate(
    EB_vio_SMR_estab = ifelse(expected_vio_estab >= 1, (vio_overal + alpha_vio_estab) / (expected_vio_estab + beta_vio_estab), NA),
    EB_wh_SMR_estab = ifelse(expected_wh_estab >= 1, (wh_overa_5 + alpha_wh_estab) / (expected_wh_estab + beta_wh_estab), NA),
    EB_acc_SMR_estab = ifelse(expected_acc_estab >= 1, (acc_overal + alpha_acc_estab) / (expected_acc_estab + beta_acc_estab), NA)
  )

#--------Poisson distributed bootstrap standard errors
n_boot <- 500

boot_vio_estab <- matrix(NA, nrow = nrow(filtered_smr_data), ncol = n_boot)
boot_wh_estab <- matrix(NA, nrow = nrow(filtered_smr_data), ncol = n_boot)
boot_acc_estab <- matrix(NA, nrow = nrow(filtered_smr_data), ncol = n_boot)

for (i in 1:nrow(filtered_smr_data)) {
  if (!is.na(filtered_smr_data$EB_vio_SMR_estab[i])) {
    for (b in 1:n_boot) {
      sim_vio <- rpois(1, filtered_smr_data$expected_vio_estab[i] * filtered_smr_data$EB_vio_SMR_estab[i])
      boot_vio_estab[i, b] <- (sim_vio + alpha_vio_estab) / 
                              (filtered_smr_data$expected_vio_estab[i] + beta_vio_estab)
    }
  }
  if (!is.na(filtered_smr_data$EB_wh_SMR_estab[i])) {
    for (b in 1:n_boot) {
      sim_wh <- rpois(1, filtered_smr_data$expected_wh_estab[i] * filtered_smr_data$EB_wh_SMR_estab[i])
      boot_wh_estab[i, b] <- (sim_wh + alpha_wh_estab) / 
                             (filtered_smr_data$expected_wh_estab[i] + beta_wh_estab)
    }
  }
  if (!is.na(filtered_smr_data$EB_acc_SMR_estab[i])) {
    for (b in 1:n_boot) {
      sim_acc <- rpois(1, filtered_smr_data$expected_acc_estab[i] * filtered_smr_data$EB_acc_SMR_estab[i])
      boot_acc_estab[i, b] <- (sim_acc + alpha_acc_estab) / 
                              (filtered_smr_data$expected_acc_estab[i] + beta_acc_estab)
    }
  }
}

filtered_smr_data <- filtered_smr_data %>%
  mutate(
    EB_vio_lower_estab = ifelse(!is.na(EB_vio_SMR_estab), apply(boot_vio_estab, 1, quantile, 0.025, na.rm = TRUE), NA),
    EB_vio_upper_estab = ifelse(!is.na(EB_vio_SMR_estab), apply(boot_vio_estab, 1, quantile, 0.975, na.rm = TRUE), NA),
    EB_wh_lower_estab = ifelse(!is.na(EB_wh_SMR_estab), apply(boot_wh_estab, 1, quantile, 0.025, na.rm = TRUE), NA),
    EB_wh_upper_estab = ifelse(!is.na(EB_wh_SMR_estab), apply(boot_wh_estab, 1, quantile, 0.975, na.rm = TRUE), NA),
    EB_acc_lower_estab = ifelse(!is.na(EB_acc_SMR_estab), apply(boot_acc_estab, 1, quantile, 0.025, na.rm = TRUE), NA),
    EB_acc_upper_estab = ifelse(!is.na(EB_acc_SMR_estab),  apply(boot_acc_estab, 1, quantile, 0.975, na.rm = TRUE), NA),
    EB_vio_sig_estab = (EB_vio_lower_estab > 1.0) | (EB_vio_upper_estab < 1.0),
    EB_wh_sig_estab = (EB_wh_lower_estab > 1.0) | (EB_wh_upper_estab < 1.0),
    EB_acc_sig_estab = (EB_acc_lower_estab > 1.0) | (EB_acc_upper_estab < 1.0)
  )

print(paste("Establishments - Violations:", sum(filtered_smr_data$EB_vio_sig_estab, na.rm = TRUE)))
print(paste("Establishments - Workplace Hazards:", sum(filtered_smr_data$EB_wh_sig_estab, na.rm = TRUE)))
print(paste("Establishments - Accidents:", sum(filtered_smr_data$EB_acc_sig_estab, na.rm = TRUE)))

# #-------Compare Raw vs Empirical Bayes Results - Smoothing should not be that crazy as we remove zip codes with less than 1 of any given denom
comparison_estab <- filtered_smr_data %>%
 st_drop_geometry() %>%
 group_by(region_8) %>%
 summarise(
   regional_SMR_vio = weighted.mean(SMR_vio_estab, w = expected_vio_estab, na.rm = TRUE),
   regional_EB_vio = weighted.mean(EB_vio_SMR_estab, w = expected_vio_estab, na.rm = TRUE),
   regional_SMR_wh = weighted.mean(SMR_wh_estab, w = expected_wh_estab, na.rm = TRUE),
   regional_EB_wh = weighted.mean(EB_wh_SMR_estab, w = expected_wh_estab, na.rm = TRUE),
   regional_SMR_acc = weighted.mean(SMR_acc_estab, w = expected_acc_estab, na.rm = TRUE),
   regional_EB_acc = weighted.mean(EB_acc_SMR_estab, w = expected_acc_estab, na.rm = TRUE),
   .groups = 'drop'
 ) %>%
 arrange(desc(regional_SMR_vio))

region_summary$regional_EB_vio_estab <- comparison_estab$regional_EB_vio[match(region_summary$region, comparison_estab$region_8)]
region_summary$regional_EB_wh_estab <- comparison_estab$regional_EB_wh[match(region_summary$region, comparison_estab$region_8)]
region_summary$regional_EB_acc_estab <- comparison_estab$regional_EB_acc[match(region_summary$region, comparison_estab$region_8)]

regional_cis_estab <- filtered_smr_data %>%
  st_drop_geometry() %>%
  group_by(region_8) %>%
  summarise(
    regional_EB_vio_lower_estab = weighted.mean(EB_vio_lower_estab, w = expected_vio_estab, na.rm = TRUE),
    regional_EB_vio_upper_estab = weighted.mean(EB_vio_upper_estab, w = expected_vio_estab, na.rm = TRUE),
    regional_EB_wh_lower_estab = weighted.mean(EB_wh_lower_estab, w = expected_wh_estab, na.rm = TRUE),
    regional_EB_wh_upper_estab = weighted.mean(EB_wh_upper_estab, w = expected_wh_estab, na.rm = TRUE),
    regional_EB_acc_lower_estab = weighted.mean(EB_acc_lower_estab, w = expected_acc_estab, na.rm = TRUE),
    regional_EB_acc_upper_estab = weighted.mean(EB_acc_upper_estab, w = expected_acc_estab, na.rm = TRUE),
    .groups = 'drop'
  )

region_summary$regional_EB_vio_lower_estab <- regional_cis_estab$regional_EB_vio_lower_estab[match(region_summary$region, regional_cis_estab$region_8)]
region_summary$regional_EB_vio_upper_estab <- regional_cis_estab$regional_EB_vio_upper_estab[match(region_summary$region, regional_cis_estab$region_8)]
region_summary$regional_EB_wh_lower_estab <- regional_cis_estab$regional_EB_wh_lower_estab[match(region_summary$region, regional_cis_estab$region_8)]
region_summary$regional_EB_wh_upper_estab <- regional_cis_estab$regional_EB_wh_upper_estab[match(region_summary$region, regional_cis_estab$region_8)]
region_summary$regional_EB_acc_lower_estab <- regional_cis_estab$regional_EB_acc_lower_estab[match(region_summary$region, regional_cis_estab$region_8)]
region_summary$regional_EB_acc_upper_estab <- regional_cis_estab$regional_EB_acc_upper_estab[match(region_summary$region, regional_cis_estab$region_8)]

region_summary <- region_summary %>%
  mutate(
    regional_EB_vio_width_estab = regional_EB_vio_upper_estab - regional_EB_vio_lower_estab,
    regional_EB_wh_width_estab = regional_EB_wh_upper_estab - regional_EB_wh_lower_estab,
    regional_EB_acc_width_estab = regional_EB_acc_upper_estab - regional_EB_acc_lower_estab
  )

print("Updated region_summary with confidence intervals:")
print(region_summary, n = Inf, width = Inf)

naics_lookup <- c(
  "11" = "Agric", "21" = "Mining", "22" = "Util", "23" = "Const",
  "31" = "Manuf", "42" = "Whole", "44" = "Retail", "48" = "Trans",
  "51" = "Info", "52" = "Finance", "53" = "RealEst", "54" = "Prof",
  "55" = "Mgmt", "56" = "Admin", "61" = "Educ", "62" = "Health",
  "71" = "Arts", "72" = "Accom", "81" = "Other"
)

regional_summary_formatted <- region_summary %>%
  mutate(
    `Top Binary Industries (Vio/WH/Acc)` = paste(
      naics_lookup[as.character(top_labor_binary_industry)],
      naics_lookup[as.character(top_wh_binary_industry)],
      naics_lookup[as.character(top_acc_binary_industry)],
      sep = "/"
    ),
    `Top Continuous Industries (Vio/WH/Acc)` = paste(
      naics_lookup[as.character(top_labor_cont_industry)],
      naics_lookup[as.character(top_wh_cont_industry)],
      naics_lookup[as.character(top_acc_cont_industry)],
      sep = "/"
    ),
    `Overall Risk (Binary/Continuous)` = paste(
      round(overall_binary_risk, 2),
      round(overall_cont_risk, 2),
      sep = "/"
    ),
    `EB SMR Violations` = paste0(
      round(regional_EB_vio_estab, 2),
      " (", round(regional_EB_vio_lower_estab, 2),
      "-", round(regional_EB_vio_upper_estab, 2), ")"
    ),
    `EB SMR Wage & Hour` = paste0(
      round(regional_EB_wh_estab, 2),
      " (", round(regional_EB_wh_lower_estab, 2),
      "-", round(regional_EB_wh_upper_estab, 2), ")"
    ),
    `EB SMR Accidents` = paste0(
      round(regional_EB_acc_estab, 2),
      " (", round(regional_EB_acc_lower_estab, 2),
      "-", round(regional_EB_acc_upper_estab, 2), ")"
    )
  ) %>%
  dplyr::select(
    Region = region,
    `Top Binary Industries (Vio/WH/Acc)`,
    `Top Continuous Industries (Vio/WH/Acc)`,
    `Overall Risk (Binary/Continuous)`,
    `EB SMR Violations`,
    `EB SMR Wage & Hour`,
    `EB SMR Accidents`
  )

print(regional_summary_formatted, n = Inf, width = Inf)

write.csv(regional_summary_formatted, file = file.path(output_folder, "region_summary_with_EB_SMRs_and_CIs.csv"), row.names = FALSE)

#-------Visuals

# Join the SMRs with the spatial data to map them
smr_lookup <- filtered_smr_data %>%
  st_drop_geometry() %>%
  dplyr::select(zip_id, starts_with("EB_"), starts_with("SMR_"))

model_data_with_smrs <- model_data %>%
  mutate(zip_id = as.character(zip_id)) %>%
  left_join(smr_lookup %>% mutate(zip_id = as.character(zip_id)), by = "zip_id")

model_data_with_smrs <- st_as_sf(model_data_with_smrs, sf_column_name = "geometry")

contiguous_states <- setdiff(levels(model_data_with_smrs$state), c("AK", "HI", "PR", "DC"))

model_data_contiguous <- model_data_with_smrs %>%
  subset(state %in% contiguous_states) %>%
  dplyr::filter(!st_is_empty(geometry)) %>%
  st_transform(5070)

# Build region boundaries from contiguous state data
region_boundaries <- model_data_contiguous %>%
  sf::st_drop_geometry() %>%
  dplyr::distinct(state, region_8) %>%
  dplyr::left_join(tigris::states(cb = TRUE), by = c("state" = "STUSPS")) %>%
  sf::st_as_sf() %>%
  dplyr::group_by(region_8) %>%
  dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop")%>%
  st_transform(5070)
plot(region_boundaries)

vio_quantiles <- quantile(model_data_contiguous$EB_vio_SMR_estab, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), na.rm = TRUE)
vio_smr_plot <- ggplot(model_data_contiguous) +
  geom_sf(aes(fill = EB_vio_SMR_estab, geometry = geometry), color = NA) +
  geom_sf(data = region_boundaries, aes(geometry = geometry), fill = NA, color = "black", linewidth = 0.4) +
  scale_fill_distiller(
    palette = "Blues", na.value = "grey90", direction = 1,
    breaks = vio_quantiles, limits = c(0, 2),
    guide = guide_colorbar(label = TRUE, title = NULL, label.theme = element_text(size = 8, angle = 45), barwidth = unit(10, "cm"), barheight = unit(0.6, "cm"))) +
  theme_void() +
  theme(legend.position = 'bottom', text = element_text(family = "Serif")) +
  labs(title = "Violations SMR (Establishments)")

wh_quantiles <- quantile(model_data_contiguous$EB_wh_SMR_estab, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), na.rm = TRUE)
wage_smr_plot <- ggplot(model_data_contiguous) +
  geom_sf(aes(fill = EB_wh_SMR_estab, geometry = geometry), color = "darkgrey", linewidth = .002) +
  geom_sf(data = region_boundaries, aes(geometry = geometry), fill = NA, color = "black", linewidth = 0.4) +
  scale_fill_distiller(
    palette = "Reds", na.value = "grey90", direction = 1,
    breaks = wh_quantiles, limits = c(0, 2),
    guide = guide_colorbar(label = TRUE, title = NULL, label.theme = element_text(size = 8, angle = 45), barwidth = unit(10, "cm"), barheight = unit(0.6, "cm"))) +
  theme_void() +
  theme(legend.position = 'bottom', text = element_text(family = "Serif")) +
  labs(title = "Wage & Hour SMR (Establishments)")

acc_quantiles <- quantile(model_data_contiguous$EB_acc_SMR_estab, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), na.rm = TRUE)
acc_smr_plot <- ggplot(model_data_contiguous) +
  geom_sf(aes(fill = EB_acc_SMR_estab, geometry = geometry), color = "darkgrey", linewidth = .002) +
  geom_sf(data = region_boundaries, aes(geometry = geometry), fill = NA, color = "black", linewidth = 0.4) +
  scale_fill_distiller(
    palette = "Greens", na.value = "grey90", direction = 1,
    breaks = acc_quantiles, limits = c(0, 2),
    guide = guide_colorbar(label = TRUE, title = NULL,label.theme = element_text(size = 8, angle = 45), barwidth = unit(10, "cm"), barheight = unit(0.6, "cm"))) +
  theme_void() +
  theme(legend.position = 'bottom', text = element_text(family = "Serif")) +
  labs(title = "Accident Rate SMR (Establishments)")

top <- plot_grid(vio_smr_plot, wage_smr_plot, ncol = 2)
bottom <- plot_grid(NULL, acc_smr_plot, NULL,ncol = 3, rel_widths = c(1, 2, 1)) 
smr_combined_plot <- plot_grid(top, bottom, ncol = 1, rel_heights = c(1, 1))

ggsave(file.path(output_folder, "smr_maps.png"), plot = smr_combined_plot, width = 10, height = 8, dpi = 300)

#-----Significance maps
model_data_contiguous <- model_data_contiguous %>%
  dplyr::mutate(
    vio_sig_cat = case_when(
      is.na(EB_vio_lower_estab) | is.na(EB_vio_upper_estab) ~ "No Data",
      EB_vio_lower_estab > 1.0 ~ "High",
      EB_vio_upper_estab < 1.0 ~ "Low",
      TRUE ~ "CI Includes 1"
    ),
    wh_sig_cat = case_when(
      is.na(EB_wh_lower_estab) | is.na(EB_wh_upper_estab) ~ "No Data",
      EB_wh_lower_estab > 1.0 ~ "High",
      EB_wh_upper_estab < 1.0 ~ "Low",
      TRUE ~ "CI Includes 1"
    ),
    acc_sig_cat = case_when(
      is.na(EB_acc_lower_estab) | is.na(EB_acc_upper_estab) ~ "No Data",
      EB_acc_lower_estab > 1.0 ~ "High",
      EB_acc_upper_estab < 1.0 ~ "Low",
      TRUE ~ "CI Includes 1"
    )
  )

vio_sig_map <- ggplot(model_data_contiguous) +
  geom_sf(aes(fill = vio_sig_cat, geometry = geometry), color = "darkgrey", linewidth = .002) +
  geom_sf(data = region_boundaries, aes(geometry = geometry), fill = NA, color = "black", linewidth = 0.4) +
  scale_fill_manual(
    name = "CI Excludes 1",
    values = c("High" = "#08519c", "CI Includes 1" = "#f7f7f7", "Low" = "#c6dbef", "No Data" = "grey90"),
    na.value = "grey90",
    guide = guide_legend(title = NULL, label.theme = element_text(size = 8), keywidth = unit(1.5, "cm"), keyheight = unit(0.5, "cm"), nrow = 1)) +
  theme_void() +
  theme(legend.position = "bottom", text = element_text(family = "Serif")) +
  labs(title = "Violations CI Excludes 1")

wh_sig_map <- ggplot(model_data_contiguous) +
  geom_sf(aes(fill = wh_sig_cat, geometry = geometry), color = "darkgrey", linewidth = .002) +
  geom_sf(data = region_boundaries, aes(geometry = geometry), fill = NA, color = "black", linewidth = 0.4) +
  scale_fill_manual(
    name = "CI Excludes 1",
    values = c("High" = "#a50f15", "CI Includes 1" = "#f7f7f7", "Low" = "#fcbba1", "No Data" = "grey90"),
    na.value = "grey90",
    guide = guide_legend(title = NULL, label.theme = element_text(size = 8), keywidth = unit(1.5, "cm"), keyheight = unit(0.5, "cm"), nrow = 1)) +
  theme_void() +
  theme(legend.position = "bottom", text = element_text(family = "Serif")) +
  labs(title = "Wage & Hour CI Excludes 1")

acc_sig_map <- ggplot(model_data_contiguous) +
  geom_sf(aes(fill = acc_sig_cat, geometry = geometry), color = "darkgrey", linewidth = .002) +
  geom_sf(data = region_boundaries, aes(geometry = geometry), fill = NA, color = "black", linewidth = 0.4) +
  scale_fill_manual(
    name = "CI Excludes 1",
    values = c("High" = "#238b45", "CI Includes 1" = "#f7f7f7", "Low" = "#c7e9c0", "No Data" = "grey90"),
    na.value = "grey90",
    guide = guide_legend(title = NULL, label.theme = element_text(size = 8), keywidth = unit(1.5, "cm"), keyheight = unit(0.5, "cm"), nrow = 1)) +
  theme_void() +
  theme(legend.position = "bottom", text = element_text(family = "Serif")) +
  labs(title = "Accidents CI Excludes 1")

top_sig <- plot_grid(vio_sig_map, wh_sig_map, ncol = 2)
bottom_sig <- plot_grid(NULL, acc_sig_map, NULL, ncol = 3, rel_widths = c(1, 2, 1))
sig_combined_plot <- plot_grid(top_sig, bottom_sig, ncol = 1, rel_heights = c(1, 1))

ggsave(filename = file.path(output_folder, "smr_sig_maps.png"), plot = sig_combined_plot, width = 10, height = 8, dpi = 300)

##==================================================##
##     7. REGIONAL VARIATION IN INDUSTRY RISK       ##
##==================================================##


# This is really the start of analyses, looking at the complex bivariate relationships between
# industry presence, region, and OSHA violations, wage violations, and accident rates

# This section examines how industry-specific risks vary geographically across regions
#     Specifically, we take the state effects
#     We then run simple interction models with only a given industry and region
#         -Largely just an exploratory step to better understand how the mix of industy and region affects the presence of OSHA and wage violations
#     We then simply visualize the results of these models
#         1. We create a heatmap of the interaction effects
#         2. We create a barplot of the interaction effects (e.g., essentially just ropes and ladder but bars)
#         3. We calculate the variance of the interaction effects by industry and region
#         4. We create a correlation plot of the interaction effectsn(e.g., within each region, how do the OSHA and wage violations correlate among industries)

# ---Running interaction models for industry and region
#     There are four things that change based on the models
#         1. The outcome (OSHA violations, wage violations, accident rate)
#         2. The family of the model (e.g., binomial for binary, lognormal for continuous)
#         3. The zip_id variable (e.g., zip_id for binary, zip_id_pos for continuous)
#         4. The spatial adjacancy network (e.g., model_adj_binary for binary, model_adj_labor_pos for continuous)

# The three binary first part of the hurdle models (Takes a lot of time to run!)


osha_binary_interactions <- run_interaction_models(industry_vars, "labor_binary", model_data, "binary")
wh_binary_interactions <- run_interaction_models(industry_vars, "wh_binary", model_data, "binary")
acc_binary_interactions <- run_interaction_models(industry_vars, "acc_binary", model_data, "binary")

# The second part of the hurdle models (Also takes a lot of time to run!)
osha_cont_interactions <- run_interaction_models(industry_vars, "vio_rate", labor_positive, "labor_continuous")
wh_cont_interactions <- run_interaction_models(industry_vars, "wh_rate_pe", wh_positive, "wh_continuous")
acc_cont_interactions <- run_interaction_models(industry_vars, "acc_rate", acc_positive, "acc_continuous")

# Binding by outcome
osha_interactions <- bind_rows(
  osha_binary_interactions %>% mutate(model_type = "binary"),
  osha_cont_interactions %>% mutate(model_type = "continuous")
)

wh_interactions <- bind_rows(
  wh_binary_interactions %>% mutate(model_type = "binary"),
  wh_cont_interactions %>% mutate(model_type = "continuous")
) 

acc_interactions <- bind_rows(
  acc_binary_interactions %>% mutate(model_type = "binary"),
  acc_cont_interactions %>% mutate(model_type = "continuous")
)

binary_interactions <- bind_rows(
  osha_binary_interactions %>% mutate(model_type = "binary"),
  wh_binary_interactions %>% mutate(model_type = "binary"),
  acc_binary_interactions %>% mutate(model_type = "binary")
)

continuous_interactions <- bind_rows(
  osha_cont_interactions %>% mutate(model_type = "continuous"),
  wh_cont_interactions %>% mutate(model_type = "continuous"),
  acc_cont_interactions %>% mutate(model_type = "continuous")
)

lim <- .1
full_interactions <- bind_rows(
  osha_interactions ,
  wh_interactions ,
  acc_interactions 
)

industry_names <- c(
  "ESTAB_11_P" = "Agriculture",
  "ESTAB_21_P" = "Mining",
  "ESTAB_22_P" = "Utilities",
  "ESTAB_23_P" = "Construction",
  "ESTAB_31_1" = "Manufacturing",
  "ESTAB_42_P" = "Wholesale Trade",
  "ESTAB_44_1" = "Retail Trade",
  "ESTAB_48_1" = "Transportation",
  "ESTAB_51_P" = "Information",
  "ESTAB_52_P" = "Finance & Insurance",
  "ESTAB_53_P" = "Real Estate",
  "ESTAB_54_P" = "Professional Services",
  "ESTAB_55_P" = "Management",
  "ESTAB_56_P" = "Administrative Services",
  "ESTAB_61_P" = "Education",
  "ESTAB_62_P" = "Healthcare",
  "ESTAB_71_P" = "Arts & Entertainment",
  "ESTAB_72_P" = "Accommodation & Food",
  "ESTAB_81_P" = "Other Services"
)

binary_interactions$industry_name <- industry_names[binary_interactions$industry]
binary_interactions <- binary_interactions %>%
  mutate(
    direction = case_when(
      coefficient > 0 & significant ~ "Positive",
      coefficient < 0 & significant ~ "Negative",
      TRUE ~ "Null"
    ),
    direction = factor(direction, levels = c("Negative", "Null", "Positive")),
    plot_y = case_when(
      direction == "Positive" ~ 1,
      direction == "Negative" ~ -1,
      TRUE ~ 0
    ),
    symbol = case_when(
      direction == "Positive" ~ "+", 
      direction == "Negative" ~ "\u2212", 
      TRUE ~ "\u2022"  
    ),
    outcome = case_when(
      outcome == "labor_binary" ~ "OSHA Violations",
      outcome == "wh_binary" ~ "Wage & Hour Violations",
      outcome == "acc_binary" ~ "Accidents"
    ),
    shade = as.integer(factor(industry_name)) %% 2 == 0
  )
shade_df <-binary_interactions %>% distinct(industry_name, shade) %>% filter(shade)
region_colors <- c(
  "New England" = "#1D70A2",
  "Mid-Atlantic" = "#558564",
  "Southeast" = "#8E3B46",
  "Plains" = "#C191A1",
  "Mountain West" = "#02394A",
  "Southwest" = "#424B54",
  "Great Lakes" = "#CF9893",
  "Pacific" = "#81C14B"
)
binary_interactions$region <- factor(binary_interactions$region, levels = names(region_colors))
test <- ggplot(binary_interactions, aes(x = region, y = plot_y, color = region)) +
  geom_rect(data = shade_df, aes(x = NULL, y = NULL, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "grey90",  inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 0.3) +
  geom_text(aes(label = symbol), position = position_jitter(width = 0.05, height = 0.05, seed = 123), size = 5, fontface = "bold") +
  scale_color_manual(values = region_colors, guide = "none") +
  scale_shape_identity() +
  scale_color_manual(values = region_colors, guide = "none") +
  scale_fill_manual(values = region_colors, na.value = "white", guide = "none") +
  scale_y_continuous(limits = c(-1.5, 1.5), breaks = c(-1, 0, 1), labels = c("Negative", "Null", "Positive"), expand = expansion(mult = 0.2), position = "right") +
  facet_grid(rows = vars(industry_name), cols = vars(outcome), scales = "fixed", switch = "y") +
  labs(x = NULL, y = NULL, title = "Binary WAV Varaible Outcomes: Industry-Region Interactions", subtitle = "Above line = Positive and CrI > 0 | At line = Null | Below line = Negative and CrI < 0") +
  theme_minimal(base_size = 6) +
  theme(
    text = element_text(family = "Serif"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 0, size = 10),
    strip.text.x = element_text(size = 8, face = "bold"),
    panel.spacing.x = unit(4, "pt"),
    panel.spacing.y = unit(4, "pt"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    legend.position = "none"
  )
ggsave(filename = file.path(output_folder, "binary_interactions_plot_test.png"), plot = test, width = 13, height = 8, dpi = 400)


continuous_interactions$industry_name <- industry_names[continuous_interactions$industry]
continuous_interactions <- continuous_interactions %>%
  mutate(
    direction = case_when(
      coefficient > 0 & significant ~ "Positive",
      coefficient < 0 & significant ~ "Negative",
      TRUE ~ "Null"
    ),
    direction = factor(direction, levels = c("Negative", "Null", "Positive")),
    plot_y = case_when(
      direction == "Positive" ~ 1,
      direction == "Negative" ~ -1,
      TRUE ~ 0
    ),
    symbol = case_when(
      direction == "Positive" ~ "+", 
      direction == "Negative" ~ "\u2212", 
      TRUE ~ "\u2022"  
    ),
    outcome = case_when(
      outcome == "vio_rate" ~ "OSHA Violations",
      outcome == "wh_rate_pe" ~ "Wage & Hour Violations",
      outcome == "acc_rate" ~ "Accidents"
    ),
    shade = as.integer(factor(industry_name)) %% 2 == 0
  )
continuous_interactions$region <- factor(continuous_interactions$region, levels = names(region_colors))
test_2 <- ggplot(continuous_interactions, aes(x = region, y = plot_y, color = region)) +
  geom_rect(data = shade_df, aes(x = NULL, y = NULL, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "grey90",  inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 0.3) +
  geom_text(aes(label = symbol), position = position_jitter(width = 0.05, height = 0.05, seed = 123), size = 5, fontface = "bold") +
  scale_color_manual(values = region_colors, guide = "none") +
  scale_shape_identity() +
  scale_color_manual(values = region_colors, guide = "none") +
  scale_fill_manual(values = region_colors, na.value = "white", guide = "none") +
  scale_y_continuous(limits = c(-1.5, 1.5), breaks = c(-1, 0, 1), labels = c("Negative", "Null", "Positive"), expand = expansion(mult = 0.2), position = "right") +
  facet_grid(rows = vars(industry_name), cols = vars(outcome), scales = "fixed", switch = "y") +
  labs(x = NULL, y = NULL, title = "Continuous WAV Varaible Outcomes: Industry-Region Interactions", subtitle = "Above line = Positive and CrI > 0 | At line = Null | Below line = Negative and CrI < 0") +
  theme_minimal(base_size = 6) +
  theme(
    text = element_text(family = "Serif"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 0, size = 10),
    strip.text.x = element_text(size = 8, face = "bold"),
    panel.spacing.x = unit(4, "pt"),
    panel.spacing.y = unit(4, "pt"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    legend.position = "none"
  )
ggsave(filename = file.path(output_folder, "continuous_interactions_plot_test.png"), plot = test_2, width = 12, height = 8, dpi = 300)


legend_df <- data.frame(
  region = factor(names(region_colors), levels = names(region_colors))
)
legend_df$x <- as.numeric(legend_df$region)
n <- length(region_colors)


legend <- ggplot() +
  annotate("rect", xmin = 0.3, xmax = n + 0.7, ymin = -0.6, ymax = 0.6, fill = "grey90", color = NA) +
  annotate("segment", x = 0.3, xend = n + 0.7, y = 0, yend = 0, color = "grey40", linewidth = 0.3, linetype = "dashed") +
  annotate("rect", xmin = 0.3, xmax = n + 0.7, ymin = -0.6, ymax = 0.6, fill = NA, color = "black", linewidth = 0.3) +
  geom_text(data = legend_df, aes(x = x, y = 0, label = "\u2022", color = region), size = 10, fontface = "bold") +
  geom_text(data = legend_df, aes(x = x, y = -0.7, label = region, color = region), size = 3, fontface = "bold", angle = 45, family = "Serif", hjust = 1) +
  annotate("text", x = (n+1)/2, y = -1.15, label = "Notes: Industries appear in same order as in the main plots. \n Symbols designates direction of effect (+/\u2212) and color designates region", size = 4, fontface = "bold", family = "Serif") +
  scale_color_manual(values = region_colors, guide = "none") +
  scale_y_continuous(limits = c(-1.2, 0.8)) +
  scale_x_continuous(limits = c(0, n + 1), breaks = NULL) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(legend.position = "none", text = element_text(family = "Serif"), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)) 
ggsave(filename = file.path(output_folder, "binary_interactions_legend.png"), plot = legend, width = 5, height = 4, dpi = 400)


osha_interactions_table <- create_model_table(osha_interactions)
osha_interactions_table$industry_name <- industry_names[osha_interactions_table$industry]
write.csv(osha_interactions_table, file = file.path(output_folder, "osha_interactions_table.csv"), row.names = FALSE)


wh_interactions_table <- create_model_table(wh_interactions)
wh_interactions_table$industry_name <- industry_names[wh_interactions_table$industry]
write.csv(wh_interactions_table, file = file.path(output_folder, "wh_interactions_table.csv"), row.names = FALSE)


acc_interactions_table <- create_model_table(acc_interactions)
acc_interactions_table$industry_name <- industry_names[acc_interactions_table$industry]
write.csv(acc_interactions_table, file = file.path(output_folder, "acc_interactions_table.csv"), row.names = FALSE)































