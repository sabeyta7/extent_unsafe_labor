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
    ESTAB_62_P, ESTAB_71_P, ESTAB_72_P, ESTAB_TOTA,

    # Geographic identifiers
    state, region_8,

    # Inspection rate
    insp_rate,

    # Geometry
    geometry,

    # Sensitivity variables

    # Natural breaks
    wh_overa_4, wh_rate_na,
    # Two way cap
    wh_overa_6, wh_rate__1
  ) %>%
  na.omit()

naics_descriptions <- c(
    "11" = "Agriculture",
    "21" = "Mining",
    "22" = "Utilities",
    "23" = "Construction", 
    "31" = "Manufacturing",
    "42" = "Wholesale Trade",
    "44" = "Retail Trade",
    "48" = "Transportation",
    "51" = "Information",
    "52" = "Finance & Insurance",
    "53" = "Real Estate",
    "54" = "Professional Services",
    "55" = "Management",
    "56" = "Administrative Services",
    "61" = "Education",
    "62" = "Healthcare",
    "71" = "Arts & Entertainment",
    "72" = "Accommodation & Food")

industry_lookup <- c(
  "11" = "Agric",
  "21" = "Mining",
  "22" = "Util",
  "23" = "Const",
  "31" = "Manuf",
  "42" = "Whole",
  "44" = "Retail",
  "48" = "Trans",
  "51" = "Info",
  "52" = "Finance",
  "53" = "RealEst",
  "54" = "Prof",
  "55" = "Mgmt",
  "56" = "Admin",
  "61" = "Educ",
  "62" = "Health",
  "71" = "Arts",
  "72" = "Accom"
)

##===============================================##
##     SENSITIVITY ANALYSIS - NATURAL 51        ##
##===============================================##

# Binary and positive subsets
model_data$wh_binary_na <- as.numeric(model_data$wh_rate_na > 0)
model_data$wh_binary_tw <- as.numeric(model_data$wh_rate__1 > 0)

wh_positive_na <- model_data %>% filter(wh_rate_na > 0)
wh_positive_tw <- model_data %>% filter(wh_rate__1 > 0)

##===============================================##
##     REGIONAL DESCRIPTIVES                     ##
##===============================================##

regional_descriptives_na <- model_data %>%
  st_drop_geometry() %>%
  group_by(region_8) %>%
  summarise(
    total_wh_violations = sum(wh_overa_4, na.rm = TRUE),
    total_establishments = sum(ESTAB_TOTA, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(regional_wh_rate = total_wh_violations / total_establishments * 100)

regional_descriptives_tw <- model_data %>%
  st_drop_geometry() %>%
  group_by(region_8) %>%
  summarise(
    total_wh_violations = sum(wh_overa_6, na.rm = TRUE),
    total_establishments = sum(ESTAB_TOTA, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(regional_wh_rate = total_wh_violations / total_establishments * 100)

##===============================================##
##     TOP INDUSTRIES                            ##
##===============================================##

industry_cols <- grep("^ESTAB_.*_P$|^ESTAB_.*_1$", names(model_data), value = TRUE)

national_wh_volume_na <- mean(wh_positive_na$wh_rate_na, na.rm = TRUE)
national_wh_volume_tw <- mean(wh_positive_tw$wh_rate__1, na.rm = TRUE)

results_na <- list()
results_tw <- list()

for (region in unique(model_data$region_8)) {
  region_data <- model_data %>% filter(region_8 == region)
  wh_positive_region_na <- wh_positive_na %>% filter(region_8 == region)
  wh_positive_region_tw <- wh_positive_tw %>% filter(region_8 == region)

  for (industry in industry_cols) {
    ind_name <- industry_abbrev[naics_descriptions[gsub("ESTAB_|_P|_1", "", industry)]]

    wh_binary_cor_na <- cor(region_data[[industry]], region_data$wh_binary_na,
                             use = "pairwise.complete.obs")
    wh_binary_cor_tw <- cor(region_data[[industry]], region_data$wh_binary_tw,
                             use = "pairwise.complete.obs")

    wh_cont_cor_na <- NA
    wh_cont_cor_tw <- NA

    if (nrow(wh_positive_region_na) > 5) {
      wh_cont_cor_na <- cor(wh_positive_region_na[[industry]],
                             wh_positive_region_na$wh_rate_na,
                             use = "pairwise.complete.obs")
    }

    if (nrow(wh_positive_region_tw) > 5) {
      wh_cont_cor_tw <- cor(wh_positive_region_tw[[industry]],
                             wh_positive_region_tw$wh_rate__1,
                             use = "pairwise.complete.obs")
    }

    results_na[[paste(region, ind_name, sep = "_")]] <- data.frame(
      region = region,
      industry = ind_name,
      wh_binary_correlation = wh_binary_cor_na,
      wh_cont_correlation = wh_cont_cor_na,
      relative_wh_volume = mean(wh_positive_region_na$wh_rate_na, na.rm = TRUE) /
                           national_wh_volume_na
    )

    results_tw[[paste(region, ind_name, sep = "_")]] <- data.frame(
      region = region,
      industry = ind_name,
      wh_binary_correlation = wh_binary_cor_tw,
      wh_cont_correlation = wh_cont_cor_tw,
      relative_wh_volume = mean(wh_positive_region_tw$wh_rate__1, na.rm = TRUE) /
                           national_wh_volume_tw
    )
  }
}

relative_risk_na <- bind_rows(results_na)
relative_risk_tw <- bind_rows(results_tw)

top_binary_wh_na <- relative_risk_na %>%
  group_by(region) %>%
  arrange(desc(wh_binary_correlation)) %>%
  slice_head(n = 1) %>%
  ungroup()

top_cont_wh_na <- relative_risk_na %>%
  group_by(region) %>%
  filter(!is.na(wh_cont_correlation)) %>%
  arrange(desc(wh_cont_correlation)) %>%
  slice_head(n = 1) %>%
  ungroup()

top_binary_wh_tw <- relative_risk_tw %>%
  group_by(region) %>%
  arrange(desc(wh_binary_correlation)) %>%
  slice_head(n = 1) %>%
  ungroup()

top_cont_wh_tw <- relative_risk_tw %>%
  group_by(region) %>%
  filter(!is.na(wh_cont_correlation)) %>%
  arrange(desc(wh_cont_correlation)) %>%
  slice_head(n = 1) %>%
  ungroup()

##===============================================##
##     EB SMR                                    ##
##===============================================##

national_wh_per_100_estab_na <- (sum(model_data$wh_overa_4, na.rm = TRUE) /
                                  sum(model_data$ESTAB_TOTA, na.rm = TRUE)) * 100
model_data$expected_wh_na <- (model_data$ESTAB_TOTA * national_wh_per_100_estab_na) / 100

national_wh_per_100_estab_tw <- (sum(model_data$wh_overa_6, na.rm = TRUE) /
                                  sum(model_data$ESTAB_TOTA, na.rm = TRUE)) * 100
model_data$expected_wh_tw <- (model_data$ESTAB_TOTA * national_wh_per_100_estab_tw) / 100

global_wh_SMR_na <- sum(model_data$wh_overa_4, na.rm = TRUE) /
                    sum(model_data$expected_wh_na, na.rm = TRUE)

global_wh_SMR_tw <- sum(model_data$wh_overa_6, na.rm = TRUE) /
                    sum(model_data$expected_wh_tw, na.rm = TRUE)

model_data$SMR_wh_na <- ifelse(model_data$expected_wh_na > 0,
                                model_data$wh_overa_4 / model_data$expected_wh_na, NA)

model_data$SMR_wh_tw <- ifelse(model_data$expected_wh_tw > 0,
                                model_data$wh_overa_6 / model_data$expected_wh_tw, NA)

n_wh_na <- sum(!is.na(model_data$SMR_wh_na))
numerator_wh_na <- sum(model_data$expected_wh_na *
                       (model_data$SMR_wh_na - global_wh_SMR_na)^2, na.rm = TRUE)
denominator_wh_na <- sum(model_data$expected_wh_na, na.rm = TRUE) -
                     global_wh_SMR_na * (sum(model_data$expected_wh_na, na.rm = TRUE) / n_wh_na)
variance_wh_na <- max(numerator_wh_na / denominator_wh_na, 0)

alpha_wh_na <- global_wh_SMR_na^2 / variance_wh_na
beta_wh_na <- global_wh_SMR_na / variance_wh_na

n_wh_tw <- sum(!is.na(model_data$SMR_wh_tw))
numerator_wh_tw <- sum(model_data$expected_wh_tw *
                       (model_data$SMR_wh_tw - global_wh_SMR_tw)^2, na.rm = TRUE)
denominator_wh_tw <- sum(model_data$expected_wh_tw, na.rm = TRUE) -
                     global_wh_SMR_tw * (sum(model_data$expected_wh_tw, na.rm = TRUE) / n_wh_tw)
variance_wh_tw <- max(numerator_wh_tw / denominator_wh_tw, 0)

alpha_wh_tw <- global_wh_SMR_tw^2 / variance_wh_tw
beta_wh_tw <- global_wh_SMR_tw / variance_wh_tw

filtered_na <- model_data %>%
  filter(expected_wh_na >= 1) %>%
  mutate(EB_wh_SMR_na = (wh_overa_4 + alpha_wh_na) / (expected_wh_na + beta_wh_na))

filtered_tw <- model_data %>%
  filter(expected_wh_tw >= 1) %>%
  mutate(EB_wh_SMR_tw = (wh_overa_6 + alpha_wh_tw) / (expected_wh_tw + beta_wh_tw))

n_boot <- 500

boot_wh_na <- matrix(NA, nrow = nrow(filtered_na), ncol = n_boot)
for (i in 1:nrow(filtered_na)) {
  if (!is.na(filtered_na$EB_wh_SMR_na[i])) {
    for (b in 1:n_boot) {
      sim <- rpois(1, filtered_na$expected_wh_na[i] * filtered_na$EB_wh_SMR_na[i])
      boot_wh_na[i, b] <- (sim + alpha_wh_na) / (filtered_na$expected_wh_na[i] + beta_wh_na)
    }
  }
}

boot_wh_tw <- matrix(NA, nrow = nrow(filtered_tw), ncol = n_boot)
for (i in 1:nrow(filtered_tw)) {
  if (!is.na(filtered_tw$EB_wh_SMR_tw[i])) {
    for (b in 1:n_boot) {
      sim <- rpois(1, filtered_tw$expected_wh_tw[i] * filtered_tw$EB_wh_SMR_tw[i])
      boot_wh_tw[i, b] <- (sim + alpha_wh_tw) / (filtered_tw$expected_wh_tw[i] + beta_wh_tw)
    }
  }
}

filtered_na <- filtered_na %>%
  mutate(
    EB_wh_lower_na = apply(boot_wh_na, 1, quantile, 0.025, na.rm = TRUE),
    EB_wh_upper_na = apply(boot_wh_na, 1, quantile, 0.975, na.rm = TRUE)
  )

filtered_tw <- filtered_tw %>%
  mutate(
    EB_wh_lower_tw = apply(boot_wh_tw, 1, quantile, 0.025, na.rm = TRUE),
    EB_wh_upper_tw = apply(boot_wh_tw, 1, quantile, 0.975, na.rm = TRUE)
  )

eb_smr_na <- filtered_na %>%
  st_drop_geometry() %>%
  group_by(region_8) %>%
  summarise(
    EB_wh_SMR = mean(EB_wh_SMR_na, na.rm = TRUE),
    EB_wh_lower = mean(EB_wh_lower_na, na.rm = TRUE),
    EB_wh_upper = mean(EB_wh_upper_na, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(EB_SMR_label = paste0(round(EB_wh_SMR, 2),
                                " (", round(EB_wh_lower, 2),
                                "-", round(EB_wh_upper, 2), ")"))

eb_smr_tw <- filtered_tw %>%
  st_drop_geometry() %>%
  group_by(region_8) %>%
  summarise(
    EB_wh_SMR = mean(EB_wh_SMR_tw, na.rm = TRUE),
    EB_wh_lower = mean(EB_wh_lower_tw, na.rm = TRUE),
    EB_wh_upper = mean(EB_wh_upper_tw, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(EB_SMR_label = paste0(round(EB_wh_SMR, 2),
                                " (", round(EB_wh_lower, 2),
                                "-", round(EB_wh_upper, 2), ")"))

##===============================================##
##     COMPILE TABLE                             ##
##===============================================##


region_order <- c("Pacific", "Southwest", "Southeast", "Mid-Atlantic",
                  "Great Lakes", "Mountain West", "New England", "Plains")

sensitivity_table_wide <- bind_rows(
  
  regional_descriptives_na %>%
    left_join(top_binary_wh_na %>% dplyr::select(region, industry), by = c("region_8" = "region")) %>%
    left_join(top_cont_wh_na %>% dplyr::select(region, industry), by = c("region_8" = "region"), suffix = c("_bin", "_cont")) %>%
    left_join(relative_risk_na %>% group_by(region) %>% summarise(relative_wh_volume = mean(relative_wh_volume, na.rm = TRUE)), by = c("region_8" = "region")) %>%
    left_join(eb_smr_na, by = "region_8") %>%
    mutate(censoring_method = "Natural 51"),
  
  regional_descriptives_tw %>%
    left_join(top_binary_wh_tw %>% dplyr::select(region, industry), by = c("region_8" = "region")) %>%
    left_join(top_cont_wh_tw %>% dplyr::select(region, industry), by = c("region_8" = "region"), suffix = c("_bin", "_cont")) %>%
    left_join(relative_risk_tw %>% group_by(region) %>% summarise(relative_wh_volume = mean(relative_wh_volume, na.rm = TRUE)), by = c("region_8" = "region")) %>%
    left_join(eb_smr_tw, by = "region_8") %>%
    mutate(censoring_method = "Double Rate Cap")
  
) %>%
  mutate(region_8 = factor(region_8, levels = region_order)) %>%
  arrange(metric = "Regional Rate", region_8) %>%
  dplyr::select(censoring_method, region_8, 
                regional_wh_rate, industry_bin, industry_cont, 
                relative_wh_volume, EB_SMR_label) %>%
  mutate(
    regional_wh_rate = round(regional_wh_rate, 2),
    relative_wh_volume = round(relative_wh_volume, 2)
  ) %>%
  pivot_longer(cols = c(regional_wh_rate, industry_bin, industry_cont, 
                        relative_wh_volume, EB_SMR_label),
               names_to = "metric",
               values_to = "value",
               values_transform = list(value = as.character)) %>%
  mutate(metric = factor(metric,
                         levels = c("regional_wh_rate", "industry_bin",
                                    "industry_cont", "relative_wh_volume",
                                    "EB_SMR_label"),
                         labels = c("Regional Rate", "Top WH Binary",
                                    "Top WH Continuous", "Overall Risk (Cont)",
                                    "EB SMR"))) %>%
  pivot_wider(names_from = region_8, values_from = value) %>%
  arrange(metric, censoring_method) %>%
  dplyr::select(metric, censoring_method,
                Pacific, Southwest, Southeast, `Mid-Atlantic`,
                `Great Lakes`, `Mountain West`, `New England`, Plains)


print(sensitivity_table_wide, n = Inf, width = Inf)
write.csv(sensitivity_table_wide, file.path(output_folder, "sensitivity_analysis_table.csv"), row.names = FALSE )
