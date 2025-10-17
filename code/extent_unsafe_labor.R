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

# Load shapefile
spatial_data <- st_read(file.path(data_folder, "labor_data.shp"))

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
#    - `vio_overal`, `wh_overa_4`, `acc_overal`, and `insp_overa` are raw counts of violations, wage/hour, accidents, and inspections.
#       -Sub note: `wh_overa_4` is the wage/hour violations capped at the 99th percentile (max 85 for a given zip year). Other versions include the raw count which is just `wh_overal`, 
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
    vio_overal, wh_overa_4, acc_overal, insp_overa, 

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
    geometry
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
      TRUE ~ industry_var
    )
  )

industry_summary <- industry_data_long %>%
  group_by(industry_name) %>%
  summarize(
    mean = mean(percentage, na.rm = TRUE),
    median = median(percentage, na.rm = TRUE),
    sd = sd(percentage, na.rm = TRUE),
    min = min(percentage, na.rm = TRUE),
    max = max(percentage, na.rm = TRUE),
    q25 = quantile(percentage, 0.25, na.rm = TRUE),
    q75 = quantile(percentage, 0.75, na.rm = TRUE)
  ) %>%
  arrange(desc(mean))

all_industries_hist <- ggplot(industry_data_long, aes(x = percentage)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.7) +
  facet_wrap(~ industry_name, scales = "free_y", ncol = 4) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "lightgray", color = NA),
    strip.text = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 9)
  ) +
  labs(
    title = "Distribution of Industry Percentages Across Zip Codes",
    subtitle = "Showing percentage of establishments in each industry",
    x = "Percentage of Establishments",
    y = "Frequency"
  )

ggsave(
  filename = file.path(output_folder, "all_industries_histogram.png"),
  plot = all_industries_hist,
  width = 10,
  height = 8,
  dpi = 300
)

top_industries <- c("ESTAB_23_P", "ESTAB_31_1", "ESTAB_11_P", "ESTAB_48_1", "ESTAB_72_P")
top_industry_names <- industry_data_long %>%
  filter(industry_var %in% top_industries) %>%
  dplyr::select(industry_var, industry_name) %>%
  distinct() %>%
  pull(industry_name)

top_industries_data <- industry_data_long %>%
  filter(industry_var %in% top_industries)

top_industries_hist <- ggplot(top_industries_data, aes(x = percentage)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.7) +
  geom_vline(aes(xintercept = mean), 
             data = industry_summary %>% filter(industry_name %in% top_industry_names),
             color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = median), 
             data = industry_summary %>% filter(industry_name %in% top_industry_names),
             color = "darkgreen", linetype = "solid", linewidth = 1) +
  facet_wrap(~ industry_name, scales = "free", ncol = 2) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "lightgray", color = NA),
    strip.text = element_text(size = 10, face = "bold")
  ) +
  labs(
    title = "Distribution of Top Risk Industries",
    subtitle = "Red dashed line = mean, Green solid line = median",
    x = "Percentage of Establishments",
    y = "Frequency"
  )

industry_density_plot <- ggplot(top_industries_data, aes(x = percentage, fill = industry_name)) +
  geom_density(alpha = 0.4) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  labs(
    title = "Density Comparison of Top Risk Industries",
    x = "Percentage of Establishments",
    y = "Density",
    fill = "Industry"
  )

regional_descriptives <- model_data %>%
  st_drop_geometry() %>%
  group_by(region_8) %>%
  summarise(
    n_zip_codes = n(),
    
    total_violations = sum(vio_overal, na.rm = TRUE),
    total_wh_violations = sum(wh_overa_4, na.rm = TRUE),
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

# Create industry-only table for easier viewing
industry_breakdown <- regional_descriptives %>%
  dplyr::select(region_8, agriculture_per_1k:accommodation_per_1k) %>%
  arrange(region_8)

us_total <- model_data %>%
  st_drop_geometry() %>%
  summarise(
    region_8 = "US_TOTAL",
    n_zip_codes = n(),
    
    total_violations = sum(vio_overal, na.rm = TRUE),
    total_wh_violations = sum(wh_overa_4, na.rm = TRUE),
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
    accommodation_per_1k = weighted.mean(ESTAB_72_P, census_wor, na.rm = TRUE)
  ) %>%
  mutate(
    regional_violation_rate = total_violations / total_establishments * 100,
    regional_wh_rate = total_wh_violations / total_establishments * 100,
    regional_accident_rate = total_accidents / total_establishments * 100,
    regional_inspection_rate = total_inspections / total_establishments * 100
  )

regional_with_us <- bind_rows(regional_descriptives, us_total)

regional_descriptives_t <- regional_with_us %>%
  column_to_rownames("region_8") %>%
  t() %>%
  as.data.frame()
print(round(regional_descriptives_t, 2))
write.csv(regional_descriptives_t, file.path(output_folder, "regional_descriptives.csv"), row.names = FALSE)

##=============================================##
##   3. OUTCOMES BY REGION FOR EACH INDUSTRY   ##
##=============================================##

'''
This section is more descriptive information but looking at the rates of OSHA violations, wage violations, and accident rates in each 
region, broken down by industry. It is just plotting the rate and volume histograms.
'''

# Generating plots and saving them to the output folder
agriculture_plot <- plot_industry_rates_by_region(model_data, labor_positive, wh_positive, acc_positive, "ESTAB_11_P")
ggsave(filename = file.path(output_folder, "agriculture_by_region_descriptive_plot.png"), plot = agriculture_plot, width = 10, height = 8, dpi = 300)

mining_plot <- plot_industry_rates_by_region(model_data, labor_positive, wh_positive, acc_positive, "ESTAB_21_P")
ggsave(filename = file.path(output_folder, "mining_by_region_descriptive_plot.png"), plot = mining_plot, width = 10, height = 8, dpi = 300)

utilities_plot <- plot_industry_rates_by_region(model_data, labor_positive, wh_positive, acc_positive, "ESTAB_22_P")
ggsave(filename = file.path(output_folder, "utilities_by_region_descriptive_plot.png"), plot = utilities_plot, width = 10, height = 8, dpi = 300)

construction_plot <- plot_industry_rates_by_region(model_data, labor_positive, wh_positive, acc_positive, "ESTAB_23_P")
ggsave(filename = file.path(output_folder, "construction_by_region_descriptive_plot.png"), plot = construction_plot, width = 10, height = 8, dpi = 300)

manufacturing_plot <- plot_industry_rates_by_region(model_data, labor_positive, wh_positive, acc_positive, "ESTAB_31_1")
ggsave(filename = file.path(output_folder, "manufacturing_by_region_descriptive_plot.png"), plot = manufacturing_plot, width = 10, height = 8, dpi = 300)

wholesale_trade_plot <- plot_industry_rates_by_region(model_data, labor_positive, wh_positive, acc_positive, "ESTAB_42_P")
ggsave(filename = file.path(output_folder, "wholesale_trade_by_region_descriptive_plot.png"), plot = wholesale_trade_plot, width = 10, height = 8, dpi = 300)

retail_trade_plot <- plot_industry_rates_by_region(model_data, labor_positive, wh_positive, acc_positive, "ESTAB_44_1")
ggsave(filename = file.path(output_folder, "retail_trade_by_region_descriptive_plot.png"), plot = retail_trade_plot, width = 10, height = 8, dpi = 300)

transportation_plot <- plot_industry_rates_by_region(model_data, labor_positive, wh_positive, acc_positive, "ESTAB_48_1")
ggsave(filename = file.path(output_folder, "transportation_by_region_descriptive_plot.png"), plot = transportation_plot, width = 10, height = 8, dpi = 300)

information_plot <- plot_industry_rates_by_region(model_data, labor_positive, wh_positive, acc_positive, "ESTAB_51_P")
ggsave(filename = file.path(output_folder, "information_by_region_descriptive_plot.png"), plot = information_plot, width = 10, height = 8, dpi = 300)

finance_plot <- plot_industry_rates_by_region(model_data, labor_positive, wh_positive, acc_positive, "ESTAB_52_P")
ggsave(filename = file.path(output_folder, "finance_by_region_descriptive_plot.png"), plot = finance_plot, width = 10, height = 8, dpi = 300)

real_estate_plot <- plot_industry_rates_by_region(model_data, labor_positive, wh_positive, acc_positive, "ESTAB_53_P")
ggsave(filename = file.path(output_folder, "real_estate_by_region_descriptive_plot.png"), plot = real_estate_plot, width = 10, height = 8, dpi = 300)

professional_services_plot <- plot_industry_rates_by_region(model_data, labor_positive, wh_positive, acc_positive, "ESTAB_54_P")
ggsave(filename = file.path(output_folder, "professional_services_by_region_descriptive_plot.png"), plot = professional_services_plot, width = 10, height = 8, dpi = 300)

management_plot <- plot_industry_rates_by_region(model_data, labor_positive, wh_positive, acc_positive, "ESTAB_55_P")
ggsave(filename = file.path(output_folder, "management_by_region_descriptive_plot.png"), plot = management_plot, width = 10, height = 8, dpi = 300)

administrative_services_plot <- plot_industry_rates_by_region(model_data, labor_positive, wh_positive, acc_positive, "ESTAB_56_P")
ggsave(filename = file.path(output_folder, "administrative_services_by_region_descriptive_plot.png"), plot = administrative_services_plot, width = 10, height = 8, dpi = 300)

education_plot <- plot_industry_rates_by_region(model_data, labor_positive, wh_positive, acc_positive, "ESTAB_61_P")
ggsave(filename = file.path(output_folder, "education_by_region_descriptive_plot.png"), plot = education_plot, width = 10, height = 8, dpi = 300)

healthcare_plot <- plot_industry_rates_by_region(model_data, labor_positive, wh_positive, acc_positive, "ESTAB_62_P")
ggsave(filename = file.path(output_folder, "healthcare_by_region_descriptive_plot.png"), plot = healthcare_plot, width = 10, height = 8, dpi = 300)

arts_plot <- plot_industry_rates_by_region(model_data, labor_positive, wh_positive, acc_positive, "ESTAB_71_P")
ggsave(filename = file.path(output_folder, "arts_by_region_descriptive_plot.png"), plot = arts_plot, width = 10, height = 8, dpi = 300)

accommodation_plot <- plot_industry_rates_by_region(model_data, labor_positive, wh_positive, acc_positive, "ESTAB_72_P")
ggsave(filename = file.path(output_folder, "accommodation_by_region_descriptive_plot.png"), plot = accommodation_plot, width = 10, height = 8, dpi = 300)

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


top_binary_labor_by_region <- relative_risk %>%
  group_by(region) %>%
  arrange(desc(labor_binary_correlation)) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  arrange(region, desc(labor_binary_correlation))

top_binary_wh_by_region <- relative_risk %>%
  group_by(region) %>%
  arrange(desc(wh_binary_correlation)) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  arrange(region, desc(wh_binary_correlation))

top_binary_acc_by_region <- relative_risk %>%
  group_by(region) %>%
  arrange(desc(acc_binary_correlation)) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  arrange(region, desc(acc_binary_correlation))

# Identify top risk industries by region for continuous outcomes
top_cont_labor_by_region <- relative_risk %>%
  group_by(region) %>%
  filter(!is.na(labor_cont_correlation)) %>%
  arrange(desc(labor_cont_correlation)) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  arrange(region, desc(labor_cont_correlation))

top_cont_wh_by_region <- relative_risk %>%
  group_by(region) %>%
  filter(!is.na(wh_cont_correlation)) %>%
  arrange(desc(wh_cont_correlation)) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  arrange(region, desc(wh_cont_correlation))

top_cont_acc_by_region <- relative_risk %>%
  group_by(region) %>%
  filter(!is.na(acc_cont_correlation)) %>%
  arrange(desc(acc_cont_correlation)) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  arrange(region, desc(acc_cont_correlation))

print(top_binary_labor_by_region, n = Inf, width = Inf)
print(top_binary_wh_by_region, n = Inf, width = Inf)
print(top_binary_acc_by_region, n = Inf, width = Inf)
print(top_cont_labor_by_region, n = Inf, width = Inf)
print(top_cont_wh_by_region, n = Inf, width = Inf)
print(top_cont_acc_by_region, n = Inf, width = Inf)


# Create correlation matrices for all outcomes
labor_binary_mat <- create_correlation_matrix(relative_risk, "labor_binary_correlation")
wh_binary_mat <- create_correlation_matrix(relative_risk, "wh_binary_correlation")
acc_binary_mat <- create_correlation_matrix(relative_risk, "acc_binary_correlation")

labor_cont_mat <- create_correlation_matrix(relative_risk, "labor_cont_correlation")
wh_cont_mat <- create_correlation_matrix(relative_risk, "wh_cont_correlation")
acc_cont_mat <- create_correlation_matrix(relative_risk, "acc_cont_correlation")


# Create heatmaps for all outcomes
labor_binary_heatmap <- create_heatmap(labor_binary_mat, "OSHA Violations - Regional Industry Risk (Occurrence)")
wh_binary_heatmap <- create_heatmap(wh_binary_mat, "Wage & Hour Violations - Regional Industry Risk (Occurrence)")
acc_binary_heatmap <- create_heatmap(acc_binary_mat, "Accident Rate - Regional Industry Risk (Occurrence)")

labor_cont_heatmap <- create_heatmap(labor_cont_mat, "OSHA Violations - Regional Industry Risk (Severity)")
wh_cont_heatmap <- create_heatmap(wh_cont_mat, "Wage & Hour Violations - Regional Industry Risk (Severity)")
acc_cont_heatmap <- create_heatmap(acc_cont_mat, "Accident Rate - Regional Industry Risk (Severity)")

png(file.path(output_folder, "binary_correaltion_heatmap.png"), width = 800, height = 1200)
par(mfrow = c(3, 1))
create_heatmap(labor_binary_mat, "OSHA Violations (Occurrence)")
create_heatmap(wh_binary_mat, "Wage & Hour Violations (Occurrence)")
create_heatmap(acc_binary_mat, "Accident Rate (Occurrence)")
par(mfrow = c(1, 1))
dev.off()

png(file.path(output_folder, "continuous_correlation_heatmap.png"), width = 800, height = 1200)
par(mfrow = c(3, 1))
create_heatmap(labor_cont_mat, "OSHA Violations (Severity)")
create_heatmap(wh_cont_mat, "Wage & Hour Violations (Severity)")
create_heatmap(acc_cont_mat, "Accident Rate (Severity)")
par(mfrow = c(1, 1))
dev.off()

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
    
    # Binary correlation statistics
    mean_labor_binary_cor <- mean(region_data$labor_binary_correlation, na.rm = TRUE)
    mean_wh_binary_cor <- mean(region_data$wh_binary_correlation, na.rm = TRUE)
    mean_acc_binary_cor <- mean(region_data$acc_binary_correlation, na.rm = TRUE)
    
    max_labor_binary_cor <- max(region_data$labor_binary_correlation, na.rm = TRUE)
    max_wh_binary_cor <- max(region_data$wh_binary_correlation, na.rm = TRUE)
    max_acc_binary_cor <- max(region_data$acc_binary_correlation, na.rm = TRUE)
    
    # Continuous correlation statistics
    mean_labor_cont_cor <- mean(region_data$labor_cont_correlation, na.rm = TRUE)
    mean_wh_cont_cor <- mean(region_data$wh_cont_correlation, na.rm = TRUE)
    mean_acc_cont_cor <- mean(region_data$acc_cont_correlation, na.rm = TRUE)
    
    max_labor_cont_cor <- max(region_data$labor_cont_correlation, na.rm = TRUE)
    max_wh_cont_cor <- max(region_data$wh_cont_correlation, na.rm = TRUE)
    max_acc_cont_cor <- max(region_data$acc_cont_correlation, na.rm = TRUE)
    
    # Binary concentration indices
    labor_binary_concentration <- if (mean_labor_binary_cor == 0 || is.infinite(mean_labor_binary_cor)) NA else max_labor_binary_cor / mean_labor_binary_cor
    wh_binary_concentration <- if (mean_wh_binary_cor == 0 || is.infinite(mean_wh_binary_cor)) NA else max_wh_binary_cor / mean_wh_binary_cor
    acc_binary_concentration <- if (mean_acc_binary_cor == 0 || is.infinite(mean_acc_binary_cor)) NA else max_acc_binary_cor / mean_acc_binary_cor
    
    # Continuous concentration indices
    labor_cont_concentration <- if (mean_labor_cont_cor == 0 || is.infinite(mean_labor_cont_cor)) NA else max_labor_cont_cor / mean_labor_cont_cor
    wh_cont_concentration <- if (mean_wh_cont_cor == 0 || is.infinite(mean_wh_cont_cor)) NA else max_wh_cont_cor / mean_wh_cont_cor
    acc_cont_concentration <- if (mean_acc_cont_cor == 0 || is.infinite(mean_acc_cont_cor)) NA else max_acc_cont_cor / mean_acc_cont_cor
    
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
      
      # Binary concentration indices
      labor_binary_concentration = labor_binary_concentration,
      wh_binary_concentration = wh_binary_concentration,
      acc_binary_concentration = acc_binary_concentration,
      
      # Continuous concentration indices
      labor_cont_concentration = labor_cont_concentration,
      wh_cont_concentration = wh_cont_concentration,
      acc_cont_concentration = acc_cont_concentration,
      
      # Overall risk measures
      overall_binary_risk = overall_binary_risk,
      overall_cont_risk = overall_cont_risk,
      
      # Maximum correlation values
      max_labor_binary_cor = max_labor_binary_cor,
      max_wh_binary_cor = max_wh_binary_cor,
      max_acc_binary_cor = max_acc_binary_cor,
      max_labor_cont_cor = max_labor_cont_cor,
      max_wh_cont_cor = max_wh_cont_cor,
      max_acc_cont_cor = max_acc_cont_cor
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

##==================================================##
##          5. CREATE BIVARIATE RISK MAPS           ##
##==================================================##

'''
This section creates maps showing relationship between industry concentration and violation risk. The overall process
is pretty simple: Take industry concentration and outcome rates, break those into tertiles, by which I mean, a small low
cutoff for each (default is .001) and then divide the rest based on 50th percentile to get low(or none), medium, and high. We then
just descriptively assign categorical values to each tertile combination and plot the results. So zip code A is high in specified
industry concentration, but medium in outcome, it gets the high-medium category. 
'''

# Create maps for key industries
manufacturing_maps <- create_all_outcome_maps_update(model_data, "ESTAB_31_1")
construction_maps <- create_all_outcome_maps_update(model_data, "ESTAB_23_P")
agriculture_maps <- create_all_outcome_maps_update(model_data, "ESTAB_11_P")
retail_maps <- create_all_outcome_maps_update(model_data, "ESTAB_44_1")
transportation_maps <- create_all_outcome_maps_update(model_data, "ESTAB_48_1")
accommodation_maps <- create_all_outcome_maps_update(model_data, "ESTAB_72_P")
mining_maps <- create_all_outcome_maps_update(model_data, "ESTAB_21_P")
utilities_maps <- create_all_outcome_maps_update(model_data, "ESTAB_22_P")
wholesale_maps <- create_all_outcome_maps_update(model_data, "ESTAB_42_P")
information_maps <- create_all_outcome_maps_update(model_data, "ESTAB_51_P")
finance_maps <- create_all_outcome_maps_update(model_data, "ESTAB_52_P")
professional_maps <- create_all_outcome_maps_update(model_data, "ESTAB_54_P")
management_maps <- create_all_outcome_maps_update(model_data, "ESTAB_55_P")
administrative_maps <- create_all_outcome_maps_update(model_data, "ESTAB_56_P")
education_maps <- create_all_outcome_maps_update(model_data, "ESTAB_61_P")
healthcare_maps <- create_all_outcome_maps_update(model_data, "ESTAB_62_P")
arts_maps <- create_all_outcome_maps_update(model_data, "ESTAB_71_P")
real_estate_maps <- create_all_outcome_maps_update(model_data, "ESTAB_53_P")

# Let's look at...retail
retail_maps$vio
retail_maps$wage
retail_maps$acc

# If ya wanna view stacked. It's dumb. There's too much going on, but We can I guess...
(retail_maps$vio)/
(retail_maps$wage)/
(retail_maps$acc)


# Save an example map. There are 54 in total, so....no. Just one to hard save
ggsave(filename = file.path(output_folder, "example_retail_map.png"), plot = retail_maps$viom, width = 10, height = 8, dpi = 300)


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
    wh_per_100_estab = ifelse(ESTAB_TOTA > 0, (wh_overa_4 / ESTAB_TOTA) * 100, NA),
    acc_per_100_estab = ifelse(ESTAB_TOTA > 0, (acc_overal / ESTAB_TOTA) * 100, NA),
    vio_per_1000_workers = ifelse(census_wor > 0, (vio_overal / census_wor) * 1000, NA),
    wh_per_1000_workers = ifelse(census_wor > 0, (wh_overa_4 / census_wor) * 1000, NA),
    acc_per_1000_workers = ifelse(census_wor > 0, (acc_overal / census_wor) * 1000, NA),
    vio_per_10k_pop = ifelse(census_pop > 0, (vio_overal / census_pop) * 10000, NA),
    wh_per_10k_pop = ifelse(census_pop > 0, (wh_overa_4 / census_pop) * 10000, NA),
    acc_per_10k_pop = ifelse(census_pop > 0, (acc_overal / census_pop) * 10000, NA),
    establishments_per_1000_workers = ifelse(census_wor > 0, (ESTAB_TOTA / census_wor) * 1000, NA),
    establishments_per_1000_pop = ifelse(census_pop > 0, (ESTAB_TOTA / census_pop) * 1000, NA),
    insp_per_100_estab = (insp_overa / ESTAB_TOTA) * 100,
    insp_per_1000_workers = ifelse(census_wor > 0, (insp_overa / census_wor) * 1000, NA)
  )

national_rates <- model_data %>%
  filter(!is.na(vio_per_100_estab)) %>%
  summarise(
    national_vio_per_100_estab = weighted.mean(vio_per_100_estab, ESTAB_TOTA, na.rm = TRUE),
    national_wh_per_100_estab = weighted.mean(wh_per_100_estab, ESTAB_TOTA, na.rm = TRUE),
    national_acc_per_100_estab = weighted.mean(acc_per_100_estab, ESTAB_TOTA, na.rm = TRUE),
    national_vio_per_1000_workers = weighted.mean(vio_per_1000_workers, census_wor, na.rm = TRUE),
    national_wh_per_1000_workers = weighted.mean(wh_per_1000_workers, census_wor, na.rm = TRUE),
    national_acc_per_1000_workers = weighted.mean(acc_per_1000_workers, census_wor, na.rm = TRUE),
    national_vio_per_10k_pop = weighted.mean(vio_per_10k_pop, census_pop, na.rm = TRUE),
    national_wh_per_10k_pop = weighted.mean(wh_per_10k_pop, census_pop, na.rm = TRUE),
    national_acc_per_10k_pop = weighted.mean(acc_per_10k_pop, census_pop, na.rm = TRUE)  )
print(national_rates)

#-----Expected rates, SMRs, and CIs
model_data <- model_data %>%
  mutate(
    expected_vio_estab = (ESTAB_TOTA * national_rates$national_vio_per_100_estab) / 100,
    expected_wh_estab = (ESTAB_TOTA * national_rates$national_wh_per_100_estab) / 100,
    expected_acc_estab = (ESTAB_TOTA * national_rates$national_acc_per_100_estab) / 100,
    expected_vio_workers = (census_wor * national_rates$national_vio_per_1000_workers) / 1000,
    expected_wh_workers = (census_wor * national_rates$national_wh_per_1000_workers) / 1000,
    expected_acc_workers = (census_wor * national_rates$national_acc_per_1000_workers) / 1000,
    expected_vio_pop = (census_pop * national_rates$national_vio_per_10k_pop) / 10000,
    expected_wh_pop = (census_pop * national_rates$national_wh_per_10k_pop) / 10000,
    expected_acc_pop = (census_pop * national_rates$national_acc_per_10k_pop) / 10000,
    SMR_vio_estab = ifelse(expected_vio_estab > 0, vio_overal / expected_vio_estab, NA),
    SMR_wh_estab = ifelse(expected_wh_estab > 0, wh_overa_4 / expected_wh_estab, NA),
    SMR_acc_estab = ifelse(expected_acc_estab > 0, acc_overal / expected_acc_estab, NA),
    SMR_vio_workers = ifelse(expected_vio_workers > 0, vio_overal / expected_vio_workers, NA),
    SMR_wh_workers = ifelse(expected_wh_workers > 0, wh_overa_4 / expected_wh_workers, NA),
    SMR_acc_workers = ifelse(expected_acc_workers > 0, acc_overal / expected_acc_workers, NA),
    SMR_vio_pop = ifelse(expected_vio_pop > 0, vio_overal / expected_vio_pop, NA),
    SMR_wh_pop = ifelse(expected_wh_pop > 0, wh_overa_4 / expected_wh_pop, NA),
    SMR_acc_pop = ifelse(expected_acc_pop > 0, acc_overal / expected_acc_pop, NA)
  )
summary(as.data.frame(model_data)[c("SMR_vio_estab", "SMR_wh_estab", "SMR_acc_estab", 
                                    "SMR_vio_workers", "SMR_wh_workers", "SMR_acc_workers",
                                    "SMR_vio_pop", "SMR_wh_pop", "SMR_acc_pop")])

#-------Fixed Empirical Bayes SMRs and Confidence Intervals
global_vio_SMR_estab <- sum(model_data$vio_overal, na.rm = TRUE) / sum(model_data$expected_vio_estab, na.rm = TRUE)
global_wh_SMR_estab <- sum(model_data$wh_overa_4, na.rm = TRUE) / sum(model_data$expected_wh_estab, na.rm = TRUE)
global_acc_SMR_estab <- sum(model_data$acc_overal, na.rm = TRUE) / sum(model_data$expected_acc_estab, na.rm = TRUE)

global_vio_SMR_workers <- sum(model_data$vio_overal, na.rm = TRUE) / sum(model_data$expected_vio_workers, na.rm = TRUE)
global_wh_SMR_workers <- sum(model_data$wh_overa_4, na.rm = TRUE) / sum(model_data$expected_wh_workers, na.rm = TRUE)
global_acc_SMR_workers <- sum(model_data$acc_overal, na.rm = TRUE) / sum(model_data$expected_acc_workers, na.rm = TRUE)

global_vio_SMR_pop <- sum(model_data$vio_overal, na.rm = TRUE) / sum(model_data$expected_vio_pop, na.rm = TRUE)
global_wh_SMR_pop <- sum(model_data$wh_overa_4, na.rm = TRUE) / sum(model_data$expected_wh_pop, na.rm = TRUE)
global_acc_SMR_pop <- sum(model_data$acc_overal, na.rm = TRUE) / sum(model_data$expected_acc_pop, na.rm = TRUE)

variance_vio_estab <- var(model_data$SMR_vio_estab, na.rm = TRUE)
variance_wh_estab <- var(model_data$SMR_wh_estab, na.rm = TRUE)
variance_acc_estab <- var(model_data$SMR_acc_estab, na.rm = TRUE)
variance_vio_workers <- var(model_data$SMR_vio_workers, na.rm = TRUE)
variance_wh_workers <- var(model_data$SMR_wh_workers, na.rm = TRUE)
variance_acc_workers <- var(model_data$SMR_acc_workers, na.rm = TRUE)
variance_vio_pop <- var(model_data$SMR_vio_pop, na.rm = TRUE)
variance_wh_pop <- var(model_data$SMR_wh_pop, na.rm = TRUE)
variance_acc_pop <- var(model_data$SMR_acc_pop, na.rm = TRUE)

# Calculate variance for each denominator separately
variance_vio_estab <- var(model_data$SMR_vio_estab, na.rm = TRUE)
variance_wh_estab <- var(model_data$SMR_wh_estab, na.rm = TRUE)
variance_acc_estab <- var(model_data$SMR_acc_estab, na.rm = TRUE)
variance_vio_workers <- var(model_data$SMR_vio_workers, na.rm = TRUE)
variance_wh_workers <- var(model_data$SMR_wh_workers, na.rm = TRUE)
variance_acc_workers <- var(model_data$SMR_acc_workers, na.rm = TRUE)
variance_vio_pop <- var(model_data$SMR_vio_pop, na.rm = TRUE)
variance_wh_pop <- var(model_data$SMR_wh_pop, na.rm = TRUE)
variance_acc_pop <- var(model_data$SMR_acc_pop, na.rm = TRUE)

alpha_vio_estab <- global_vio_SMR_estab^2 / variance_vio_estab 
alpha_wh_estab <- global_wh_SMR_estab^2 / variance_wh_estab
alpha_acc_estab <- global_acc_SMR_estab^2 / variance_acc_estab
alpha_vio_workers <- global_vio_SMR_workers^2 / variance_vio_workers
alpha_wh_workers <- global_wh_SMR_workers^2 / variance_wh_workers
alpha_acc_workers <- global_acc_SMR_workers^2 / variance_acc_workers
alpha_vio_pop <- global_vio_SMR_pop^2 / variance_vio_pop
alpha_wh_pop <- global_wh_SMR_pop^2 / variance_wh_pop
alpha_acc_pop <- global_acc_SMR_pop^2 / variance_acc_pop

beta_vio_estab <- global_vio_SMR_estab / variance_vio_estab
beta_wh_estab <- global_wh_SMR_estab / variance_wh_estab
beta_acc_estab <- global_acc_SMR_estab / variance_acc_estab
beta_vio_workers <- global_vio_SMR_workers / variance_vio_workers
beta_wh_workers <- global_wh_SMR_workers / variance_wh_workers
beta_acc_workers <- global_acc_SMR_workers / variance_acc_workers
beta_vio_pop <- global_vio_SMR_pop / variance_vio_pop
beta_wh_pop <- global_wh_SMR_pop / variance_wh_pop
beta_acc_pop <- global_acc_SMR_pop / variance_acc_pop

filtered_smr_data <- model_data %>%
  filter(expected_vio_estab >= 1.0 | expected_wh_estab >= 1.0 | expected_acc_estab >= 1.0 |
         expected_vio_workers >= 1.0 | expected_wh_workers >= 1.0 | expected_acc_workers >= 1.0 |
         expected_vio_pop >= 1.0 | expected_wh_pop >= 1.0 | expected_acc_pop >= 1.0) %>%
  mutate(
    EB_vio_SMR_estab = ifelse(expected_vio_estab >= 1, (vio_overal + alpha_vio_estab) / (expected_vio_estab + beta_vio_estab), NA),
    EB_wh_SMR_estab = ifelse(expected_wh_estab >= 1, (wh_overa_4 + alpha_wh_estab) / (expected_wh_estab + beta_wh_estab), NA),
    EB_acc_SMR_estab = ifelse(expected_acc_estab >= 1, (acc_overal + alpha_acc_estab) / (expected_acc_estab + beta_acc_estab), NA),
    
    EB_vio_SMR_workers = ifelse(expected_vio_workers >= 1, (vio_overal + alpha_vio_workers) / (expected_vio_workers + beta_vio_workers), NA),
    EB_wh_SMR_workers = ifelse(expected_wh_workers >= 1, (wh_overa_4 + alpha_wh_workers) / (expected_wh_workers + beta_wh_workers), NA),
    EB_acc_SMR_workers = ifelse(expected_acc_workers >= 1, (acc_overal + alpha_acc_workers) / (expected_acc_workers + beta_acc_workers), NA),
    
    EB_vio_SMR_pop = ifelse(expected_vio_pop >= 1, (vio_overal + alpha_vio_pop) / (expected_vio_pop + beta_vio_pop), NA),
    EB_wh_SMR_pop = ifelse(expected_wh_pop >= 1, (wh_overa_4 + alpha_wh_pop) / (expected_wh_pop + beta_wh_pop), NA),
    EB_acc_SMR_pop = ifelse(expected_acc_pop >= 1, (acc_overal + alpha_acc_pop) / (expected_acc_pop + beta_acc_pop), NA)
  )

#--------Poisson distributed bootstrap standard errors
n_boot <- 500

boot_vio_estab <- matrix(NA, nrow = nrow(filtered_smr_data), ncol = n_boot)
boot_wh_estab <- matrix(NA, nrow = nrow(filtered_smr_data), ncol = n_boot)
boot_acc_estab <- matrix(NA, nrow = nrow(filtered_smr_data), ncol = n_boot)

boot_vio_workers <- matrix(NA, nrow = nrow(filtered_smr_data), ncol = n_boot)
boot_wh_workers <- matrix(NA, nrow = nrow(filtered_smr_data), ncol = n_boot)
boot_acc_workers <- matrix(NA, nrow = nrow(filtered_smr_data), ncol = n_boot)

boot_vio_pop <- matrix(NA, nrow = nrow(filtered_smr_data), ncol = n_boot)
boot_wh_pop <- matrix(NA, nrow = nrow(filtered_smr_data), ncol = n_boot)
boot_acc_pop <- matrix(NA, nrow = nrow(filtered_smr_data), ncol = n_boot)

for (i in 1:nrow(filtered_smr_data)) {
  
  # Establishments denom
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
  
  # Worker population as denom
  if (!is.na(filtered_smr_data$EB_vio_SMR_workers[i])) {
    for (b in 1:n_boot) {
      sim_vio <- rpois(1, filtered_smr_data$expected_vio_workers[i] * filtered_smr_data$EB_vio_SMR_workers[i])
      boot_vio_workers[i, b] <- (sim_vio + alpha_vio_workers) / 
                                (filtered_smr_data$expected_vio_workers[i] + beta_vio_workers)
    }
  }
  if (!is.na(filtered_smr_data$EB_wh_SMR_workers[i])) {
    for (b in 1:n_boot) {
      sim_wh <- rpois(1, filtered_smr_data$expected_wh_workers[i] * filtered_smr_data$EB_wh_SMR_workers[i])
      boot_wh_workers[i, b] <- (sim_wh + alpha_wh_workers) / 
                               (filtered_smr_data$expected_wh_workers[i] + beta_wh_workers)
    }
  }
  if (!is.na(filtered_smr_data$EB_acc_SMR_workers[i])) {
    for (b in 1:n_boot) {
      sim_acc <- rpois(1, filtered_smr_data$expected_acc_workers[i] * filtered_smr_data$EB_acc_SMR_workers[i])
      boot_acc_workers[i, b] <- (sim_acc + alpha_acc_workers) / 
                                (filtered_smr_data$expected_acc_workers[i] + beta_acc_workers)
    }
  }
  
  # Pop as denom
  if (!is.na(filtered_smr_data$EB_vio_SMR_pop[i])) {
    for (b in 1:n_boot) {
      sim_vio <- rpois(1, filtered_smr_data$expected_vio_pop[i] * filtered_smr_data$EB_vio_SMR_pop[i])
      boot_vio_pop[i, b] <- (sim_vio + alpha_vio_pop) / 
                            (filtered_smr_data$expected_vio_pop[i] + beta_vio_pop)
    }
  }
  if (!is.na(filtered_smr_data$EB_wh_SMR_pop[i])) {
    for (b in 1:n_boot) {
      sim_wh <- rpois(1, filtered_smr_data$expected_wh_pop[i] * filtered_smr_data$EB_wh_SMR_pop[i])
      boot_wh_pop[i, b] <- (sim_wh + alpha_wh_pop) / 
                           (filtered_smr_data$expected_wh_pop[i] + beta_wh_pop)
    }
  }
  if (!is.na(filtered_smr_data$EB_acc_SMR_pop[i])) {
    for (b in 1:n_boot) {
      sim_acc <- rpois(1, filtered_smr_data$expected_acc_pop[i] * filtered_smr_data$EB_acc_SMR_pop[i])
      boot_acc_pop[i, b] <- (sim_acc + alpha_acc_pop) / 
                            (filtered_smr_data$expected_acc_pop[i] + beta_acc_pop)
    }
  }
}

filtered_smr_data <- filtered_smr_data %>%
  mutate(
    # Estabs
    EB_vio_lower_estab = ifelse(!is.na(EB_vio_SMR_estab), 
                                apply(boot_vio_estab, 1, quantile, 0.025, na.rm = TRUE), NA),
    EB_vio_upper_estab = ifelse(!is.na(EB_vio_SMR_estab), 
                                apply(boot_vio_estab, 1, quantile, 0.975, na.rm = TRUE), NA),
    EB_wh_lower_estab = ifelse(!is.na(EB_wh_SMR_estab), 
                               apply(boot_wh_estab, 1, quantile, 0.025, na.rm = TRUE), NA),
    EB_wh_upper_estab = ifelse(!is.na(EB_wh_SMR_estab), 
                               apply(boot_wh_estab, 1, quantile, 0.975, na.rm = TRUE), NA),
    EB_acc_lower_estab = ifelse(!is.na(EB_acc_SMR_estab), 
                                apply(boot_acc_estab, 1, quantile, 0.025, na.rm = TRUE), NA),
    EB_acc_upper_estab = ifelse(!is.na(EB_acc_SMR_estab), 
                                apply(boot_acc_estab, 1, quantile, 0.975, na.rm = TRUE), NA),
    
    # Worker
    EB_vio_lower_workers = ifelse(!is.na(EB_vio_SMR_workers), 
                                  apply(boot_vio_workers, 1, quantile, 0.025, na.rm = TRUE), NA),
    EB_vio_upper_workers = ifelse(!is.na(EB_vio_SMR_workers), 
                                  apply(boot_vio_workers, 1, quantile, 0.975, na.rm = TRUE), NA),
    EB_wh_lower_workers = ifelse(!is.na(EB_wh_SMR_workers), 
                                 apply(boot_wh_workers, 1, quantile, 0.025, na.rm = TRUE), NA),
    EB_wh_upper_workers = ifelse(!is.na(EB_wh_SMR_workers), 
                                 apply(boot_wh_workers, 1, quantile, 0.975, na.rm = TRUE), NA),
    EB_acc_lower_workers = ifelse(!is.na(EB_acc_SMR_workers), 
                                  apply(boot_acc_workers, 1, quantile, 0.025, na.rm = TRUE), NA),
    EB_acc_upper_workers = ifelse(!is.na(EB_acc_SMR_workers), 
                                  apply(boot_acc_workers, 1, quantile, 0.975, na.rm = TRUE), NA),
    
    # Pop
    EB_vio_lower_pop = ifelse(!is.na(EB_vio_SMR_pop), 
                              apply(boot_vio_pop, 1, quantile, 0.025, na.rm = TRUE), NA),
    EB_vio_upper_pop = ifelse(!is.na(EB_vio_SMR_pop), 
                              apply(boot_vio_pop, 1, quantile, 0.975, na.rm = TRUE), NA),
    EB_wh_lower_pop = ifelse(!is.na(EB_wh_SMR_pop), 
                             apply(boot_wh_pop, 1, quantile, 0.025, na.rm = TRUE), NA),
    EB_wh_upper_pop = ifelse(!is.na(EB_wh_SMR_pop), 
                             apply(boot_wh_pop, 1, quantile, 0.975, na.rm = TRUE), NA),
    EB_acc_lower_pop = ifelse(!is.na(EB_acc_SMR_pop), 
                              apply(boot_acc_pop, 1, quantile, 0.025, na.rm = TRUE), NA),
    EB_acc_upper_pop = ifelse(!is.na(EB_acc_SMR_pop), 
                              apply(boot_acc_pop, 1, quantile, 0.975, na.rm = TRUE), NA),
    
    # Just determining the zip codes that have an SMR CI fully above or below 1
    EB_vio_sig_estab = (EB_vio_lower_estab > 1.0) | (EB_vio_upper_estab < 1.0),
    EB_wh_sig_estab = (EB_wh_lower_estab > 1.0) | (EB_wh_upper_estab < 1.0),
    EB_acc_sig_estab = (EB_acc_lower_estab > 1.0) | (EB_acc_upper_estab < 1.0),
    EB_vio_sig_workers = (EB_vio_lower_workers > 1.0) | (EB_vio_upper_workers < 1.0),
    EB_wh_sig_workers = (EB_wh_lower_workers > 1.0) | (EB_wh_upper_workers < 1.0),
    EB_acc_sig_workers = (EB_acc_lower_workers > 1.0) | (EB_acc_upper_workers < 1.0),
    EB_vio_sig_pop = (EB_vio_lower_pop > 1.0) | (EB_vio_upper_pop < 1.0),
    EB_wh_sig_pop = (EB_wh_lower_pop > 1.0) | (EB_wh_upper_pop < 1.0),
    EB_acc_sig_pop = (EB_acc_lower_pop > 1.0) | (EB_acc_upper_pop < 1.0)
  )

print(paste("Establishments - Violations:", sum(filtered_smr_data$EB_vio_sig_estab, na.rm = TRUE)))
print(paste("Establishments - Workplace Hazards:", sum(filtered_smr_data$EB_wh_sig_estab, na.rm = TRUE)))
print(paste("Establishments - Accidents:", sum(filtered_smr_data$EB_acc_sig_estab, na.rm = TRUE)))

print(paste("Workers - Violations:", sum(filtered_smr_data$EB_vio_sig_workers, na.rm = TRUE)))
print(paste("Workers - Workplace Hazards:", sum(filtered_smr_data$EB_wh_sig_workers, na.rm = TRUE)))
print(paste("Workers - Accidents:", sum(filtered_smr_data$EB_acc_sig_workers, na.rm = TRUE)))

print(paste("Population - Violations:", sum(filtered_smr_data$EB_vio_sig_pop, na.rm = TRUE)))
print(paste("Population - Workplace Hazards:", sum(filtered_smr_data$EB_wh_sig_pop, na.rm = TRUE)))
print(paste("Population - Accidents:", sum(filtered_smr_data$EB_acc_sig_pop, na.rm = TRUE)))


#-------Compare Raw vs Empirical Bayes Results - Smoothing should not be that crazy as we remove zip codes with less than 1 of any given denom
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

comparison_workers <- filtered_smr_data %>%
 st_drop_geometry() %>%
 group_by(region_8) %>%
 summarise(
   regional_SMR_vio = weighted.mean(SMR_vio_workers, w = expected_vio_workers, na.rm = TRUE),
   regional_EB_vio = weighted.mean(EB_vio_SMR_workers, w = expected_vio_workers, na.rm = TRUE),
   
   regional_SMR_wh = weighted.mean(SMR_wh_workers, w = expected_wh_workers, na.rm = TRUE),
   regional_EB_wh = weighted.mean(EB_wh_SMR_workers, w = expected_wh_workers, na.rm = TRUE),
   
   regional_SMR_acc = weighted.mean(SMR_acc_workers, w = expected_acc_workers, na.rm = TRUE),
   regional_EB_acc = weighted.mean(EB_acc_SMR_workers, w = expected_acc_workers, na.rm = TRUE),
   .groups = 'drop'
 ) %>%
 arrange(desc(regional_SMR_vio))

comparison_pop <- filtered_smr_data %>%
 st_drop_geometry() %>%
 group_by(region_8) %>%
 summarise(
   regional_SMR_vio = weighted.mean(SMR_vio_pop, w = expected_vio_pop, na.rm = TRUE),
   regional_EB_vio = weighted.mean(EB_vio_SMR_pop, w = expected_vio_pop, na.rm = TRUE),
   
   regional_SMR_wh = weighted.mean(SMR_wh_pop, w = expected_wh_pop, na.rm = TRUE),
   regional_EB_wh = weighted.mean(EB_wh_SMR_pop, w = expected_wh_pop, na.rm = TRUE),
   
   regional_SMR_acc = weighted.mean(SMR_acc_pop, w = expected_acc_pop, na.rm = TRUE),
   regional_EB_acc = weighted.mean(EB_acc_SMR_pop, w = expected_acc_pop, na.rm = TRUE),
   .groups = 'drop'
 ) %>%
 arrange(desc(regional_SMR_vio))

print(comparison_estab, n = Inf, width = Inf)
print(comparison_workers, n = Inf, width = Inf)
print(comparison_pop, n = Inf, width = Inf)

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

write.csv(region_summary, file = file.path(output_folder, "region_summary_with_EB_SMRs_and_CIs.csv"), row.names = FALSE)

#-------Visuals

# Create SMR lookup and join to full data
smr_lookup <- filtered_smr_data %>%
  st_drop_geometry() %>%
  dplyr::select(zip_id, starts_with("EB_"), starts_with("SMR_"))

model_data_with_smrs <- model_data %>%
  mutate(zip_id = as.character(zip_id)) %>%
  left_join(smr_lookup %>% mutate(zip_id = as.character(zip_id)), by = "zip_id")

# Filter to contiguous states (use the working approach)
contiguous_states <- setdiff(levels(model_data_with_smrs$state), c("AK", "HI", "PR", "DC"))
model_data_contiguous <- model_data_with_smrs[model_data_with_smrs$state %in% contiguous_states, ]

# Add significance categories to the data
model_data_contiguous <- model_data_contiguous %>%
  dplyr::mutate(
    vio_sig_cat = case_when(
      is.na(EB_vio_lower_estab) | is.na(EB_vio_upper_estab) ~ "No Data",
      EB_vio_lower_estab > 1.0 ~ "High",
      EB_vio_upper_estab < 1.0 ~ "Low",
      TRUE ~ "Non-Sig"
    ),
    wh_sig_cat = case_when(
      is.na(EB_wh_lower_estab) | is.na(EB_wh_upper_estab) ~ "No Data",
      EB_wh_lower_estab > 1.0 ~ "High",
      EB_wh_upper_estab < 1.0 ~ "Low",
      TRUE ~ "Non-Sig"
    ),
    acc_sig_cat = case_when(
      is.na(EB_acc_lower_estab) | is.na(EB_acc_upper_estab) ~ "No Data",
      EB_acc_lower_estab > 1.0 ~ "High",
      EB_acc_upper_estab < 1.0 ~ "Low",
      TRUE ~ "Non-Sig"
    )
  )

model_data_contiguous <- st_as_sf(model_data_contiguous)

vio_quantiles <- quantile(model_data_with_smrs$EB_vio_SMR_estab, 
                         probs = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), na.rm = TRUE)

wh_quantiles <- quantile(model_data_with_smrs$EB_wh_SMR_estab, 
                        probs = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), na.rm = TRUE)
acc_quantiles <- quantile(model_data_with_smrs$EB_acc_SMR_estab, 
                         probs = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), na.rm = TRUE)

vio_smr_plot <- ggplot(model_data_contiguous) +
  geom_sf(aes(fill = EB_vio_SMR_estab), color = "darkgrey", linewidth = .002) +
  scale_fill_distiller(
    palette = "Blues", na.value = "grey90", direction = 1,
    breaks = vio_quantiles, limits = c(0, 2),
    guide = guide_colorbar(
      label = TRUE, title = NULL,
      label.theme = element_text(size = 8, angle = 45),
      barwidth = unit(10, "cm"), barheight = unit(0.6, "cm")  # Make color bar longer
    )
  ) +
  theme_void() +
  theme(legend.position = 'bottom') +
  labs(title = "Violations SMR (Establishments)")

wage_smr_plot <- ggplot(model_data_contiguous) +
  geom_sf(aes(fill = EB_wh_SMR_estab), color = "darkgrey", linewidth = .002) +
  scale_fill_distiller(
    palette = "Reds", na.value = "grey90", direction = 1,
    breaks = wh_quantiles, limits = c(0, 2),
    guide = guide_colorbar(
      label = TRUE, title = NULL,
      label.theme = element_text(size = 8, angle = 45),
      barwidth = unit(10, "cm"), barheight = unit(0.6, "cm")  # Make color bar longer
    )
  ) +
  theme_void() +
  theme(legend.position = 'bottom') +
  labs(title = "Wage & Hour SMR (Establishments)")

acc_smr_plot <- ggplot(model_data_contiguous) +
  geom_sf(aes(fill = EB_acc_SMR_estab), color = "darkgrey", linewidth = .002) +
  scale_fill_distiller(
    palette = "Greens", na.value = "grey90", direction = 1,
    breaks = acc_quantiles, limits = c(0, 2),
    guide = guide_colorbar(
      label = TRUE, title = NULL,
      label.theme = element_text(size = 8, angle = 45),
      barwidth = unit(10, "cm"), barheight = unit(0.6, "cm")  # Make color bar longer
    )
  ) +
  theme_void() +
  theme(legend.position = 'bottom') +
  labs(title = "Accident Rate SMR (Establishments)")

top <- plot_grid(vio_smr_plot, wage_smr_plot, ncol = 2)

bottom <- plot_grid(NULL, acc_smr_plot, NULL,
                    ncol = 3, rel_widths = c(1, 2, 1)) 

smr_combined_plot <- plot_grid(top, bottom, ncol = 1, rel_heights = c(1, 1))

ggsave(
  file.path(output_folder, "smr_maps.png"),
  plot = smr_combined_plot,
  width = 10, height = 8, dpi = 300
)

vio_sig_map <- ggplot(model_data_contiguous) +
  geom_sf(aes(fill = vio_sig_cat), color = "darkgrey", linewidth = .002) +
  scale_fill_manual(
    name = "Significance",
    values = c("High" = "#08519c",    # Dark blue
               "Non-Sig" = "#f7f7f7",       # Neutral grey/white
               "Low" = "#c6dbef",     # Light blue
               "No Data" = "grey90"),
    na.value = "grey90",
    guide = guide_legend(
      title = NULL,
      label.theme = element_text(size = 8),
      keywidth = unit(1.5, "cm"),
      keyheight = unit(0.5, "cm"),
      nrow = 1
    )
  ) +
  theme_void() +
  theme(legend.position = "bottom") +
  labs(title = "Violations Significance")

wh_sig_map <- ggplot(model_data_contiguous) +
  geom_sf(aes(fill = wh_sig_cat), color = "darkgrey", linewidth = .002) +
  scale_fill_manual(
    name = "Significance",
    values = c("High" = "#a50f15",    # Dark red
               "Non-Sig" = "#f7f7f7",       # Neutral grey/white
               "Low" = "#fcbba1",     # Light red
               "No Data" = "grey90"),
    na.value = "grey90",
    guide = guide_legend(
      title = NULL,
      label.theme = element_text(size = 8),
      keywidth = unit(1.5, "cm"),
      keyheight = unit(0.5, "cm"),
      nrow = 1
    )
  ) +
  theme_void() +
  theme(legend.position = "bottom") +
  labs(title = "Wage & Hour Significance")

acc_sig_map <- ggplot(model_data_contiguous) +
  geom_sf(aes(fill = acc_sig_cat), color = "darkgrey", linewidth = .002) +
  scale_fill_manual(
    name = "Significance",
    values = c("High" = "#238b45",    # Dark green
               "Non-Sig" = "#f7f7f7",       # Neutral grey/white
               "Low" = "#c7e9c0",     # Light green
               "No Data" = "grey90"),
    na.value = "grey90",
    guide = guide_legend(
      title = NULL,
      label.theme = element_text(size = 8),
      keywidth = unit(1.5, "cm"),
      keyheight = unit(0.5, "cm"),
      nrow = 1
    )
  ) +
  theme_void() +
  theme(legend.position = "bottom") +
  labs(title = "Accidents Significance")

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

binary_interactions <- bind_rows(
  osha_binary_interactions %>% mutate(violation_type = "OSHA Violations", model_type = "binary"),
  wh_binary_interactions %>% mutate(violation_type = "Wage & Hour Violations", model_type = "binary"),
  acc_binary_interactions %>% mutate(violation_type = "Accident Rate", model_type = "binary")
)

continuous_interactions <- bind_rows(
  osha_cont_interactions %>% mutate(violation_type = "OSHA Violations", model_type = "continuous"),
  wh_cont_interactions %>% mutate(violation_type = "Wage & Hour Violations", model_type = "continuous"),
  acc_cont_interactions %>% mutate(violation_type = "Accident Rate", model_type = "continuous")
)

#----Visualizing the interaction models
industry_region_interactions <- bind_rows(binary_interactions, continuous_interactions)


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
  "ESTAB_72_P" = "Accommodation & Food"
)

industry_region_interactions$industry_name <- industry_names[industry_region_interactions$industry]

tables_list <- list(
  osha_binary = create_model_table(industry_region_interactions, "OSHA Violations", "binary"),
  wh_binary = create_model_table(industry_region_interactions, "Wage & Hour Violations", "binary"),
  acc_binary = create_model_table(industry_region_interactions, "Accident Rate", "binary"),
  osha_cont = create_model_table(industry_region_interactions, "OSHA Violations", "continuous"),
  wh_cont = create_model_table(industry_region_interactions, "Wage & Hour Violations", "continuous"),
  acc_cont = create_model_table(industry_region_interactions, "Accident Rate", "continuous")
)

table_names <- c(
  "OSHA Violations (Binary)",
  "Wage & Hour Violations (Binary)",
  "Accident Rate (Binary)",
  "OSHA Violations (Continuous)",
  "Wage & Hour Violations (Continuous)",
  "Accident Rate (Continuous)"
)

for(i in 1:length(tables_list)) {
  tables_list[[i]]$model_type <- table_names[i]
}

combined_table <- do.call(rbind, tables_list)

big_table <- kable(combined_table, 
                   format = "html", 
                   caption = "Industry-Region Interactions") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = TRUE) %>%
  add_header_above(c(" " = 2, "Regions" = ncol(combined_table) - 1)) %>%
  column_spec(1, bold = TRUE) %>%
  footnote(symbol = "* indicates statistical significance at p < 0.05") %>%
  cat(output_folder, file="interaction_table.html")

# Making the final table
industry_abbrev <- c(
  "Accommodation & Food" = "Accom",
  "Administrative Services" = "Admin", 
  "Agriculture" = "Agric",
  "Arts & Entertainment" = "Arts",
  "Construction" = "Const",
  "Education" = "Educ",
  "Finance & Insurance" = "Finance",
  "Healthcare" = "Health",
  "Information" = "Info",
  "Management" = "Mgmt",
  "Manufacturing" = "Manuf",
  "Mining" = "Mining",
  "Professional Services" = "Prof",
  "Real Estate" = "RealEst",
  "Retail Trade" = "Retail",
  "Transportation" = "Trans",
  "Utilities" = "Util",
  "Wholesale Trade" = "Whole"
)
models <- unique(combined_table$model_type)
model_abbrev <- c(
  "OSHA Violations (Binary)" = "OSHA Binary Model",
  "Wage & Hour Violations (Binary)" = "WH Binary Model",
  "Accident Rate (Binary)" = "Accident Binary Model",
  "OSHA Violations (Continuous)" = "OSHA Continuous Model",
  "Wage & Hour Violations (Continuous)" = "WH Continuous Model",
  "Accident Rate (Continuous)" = "Accident Continuous Model"
)
all_cols <- names(combined_table)
regions <- setdiff(all_cols, "model_type")
all_results <- list()
for(model in models) {
  model_data <- combined_table[combined_table$model_type == model, ]
  for(region in regions) {
    industries <- rownames(model_data)
    coefs_raw <- model_data[[region]]
    industries_clean <- gsub("^[^.]+\\.", "", industries)
    coefs_numeric <- as.numeric(gsub("\\*", "", coefs_raw))
    has_sig <- grepl("\\*", coefs_raw)
    df <- data.frame(
      industry = industries_clean,
      coef = coefs_numeric,
      sig = has_sig,
      coef_original = coefs_raw,  # Keep original with asterisks
      stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$coef), ]
    if(nrow(df) == 0) next
    mpi_idx <- which.max(df$coef)
    mpi_ind <- df$industry[mpi_idx]
    mpi_coef <- df$coef[mpi_idx]
    mpi_sig <- df$sig[mpi_idx]
    mpi_abbr <- ifelse(mpi_ind %in% names(industry_abbrev), 
                       industry_abbrev[mpi_ind], 
                       mpi_ind)
    mpi_str <- sprintf("%s (%.2f%s)", mpi_abbr, mpi_coef, ifelse(mpi_sig, "*", ""))
    mni_idx <- which.min(df$coef)
    mni_ind <- df$industry[mni_idx]
    mni_coef <- df$coef[mni_idx]
    mni_sig <- df$sig[mni_idx]
    mni_abbr <- ifelse(mni_ind %in% names(industry_abbrev), 
                       industry_abbrev[mni_ind], 
                       mni_ind)
    mni_str <- sprintf("%s (%.2f%s)", mni_abbr, mni_coef, ifelse(mni_sig, "*", ""))
    lri_idx <- which.max(abs(df$coef))
    lri_ind <- df$industry[lri_idx]
    lri_coef <- df$coef[lri_idx]
    lri_sig <- df$sig[lri_idx]
    lri_abbr <- ifelse(lri_ind %in% names(industry_abbrev), 
                       industry_abbrev[lri_ind], 
                       lri_ind)
    lri_str <- sprintf("%s (%.2f%s)", lri_abbr, lri_coef, ifelse(lri_sig, "*", ""))
    nri_idx <- which.min(abs(df$coef))
    nri_ind <- df$industry[nri_idx]
    nri_coef <- df$coef[nri_idx]
    nri_sig <- df$sig[nri_idx]
    nri_abbr <- ifelse(nri_ind %in% names(industry_abbrev), 
                       industry_abbrev[nri_ind], 
                       nri_ind)
    nri_str <- sprintf("%s (%.2f%s)", nri_abbr, nri_coef, ifelse(nri_sig, "*", ""))
    sig_df <- df[df$sig == TRUE, ]
    if(nrow(sig_df) > 0) {
      msi_idx <- which.max(abs(sig_df$coef))
      msi_ind <- sig_df$industry[msi_idx]
      msi_coef <- sig_df$coef[msi_idx]
      msi_abbr <- ifelse(msi_ind %in% names(industry_abbrev), 
                         industry_abbrev[msi_ind], 
                         msi_ind)
      msi_dir <- ifelse(msi_coef > 0, "Pos", "Neg")
      msi_str <- sprintf("%s (%s)", msi_abbr, msi_dir)
    } else {
      msi_str <- "None ()"
    }
    key <- paste(model, region, sep = "___")
    model_abbr_name <- ifelse(model %in% names(model_abbrev), 
                               model_abbrev[model], 
                               model)
    all_results[[key]] <- list(
      Model = model_abbr_name,
      Region = region,
      MPI = mpi_str,
      MNI = mni_str,
      LRI = lri_str,
      NRI = nri_str,
      MSI = msi_str
    )
  }
}
results_df <- do.call(rbind, lapply(all_results, function(x) {
  data.frame(
    Model = x$Model,
    Region = x$Region,
    MPI = x$MPI,
    MNI = x$MNI,
    LRI = x$LRI,
    NRI = x$NRI,
    MSI = x$MSI,
    stringsAsFactors = FALSE
  )
}))
rownames(results_df) <- NULL
metrics <- c("MPI", "MNI", "LRI", "NRI", "MSI")
final_table3 <- data.frame()
unique_models <- unique(results_df$Model)
for(model in unique_models) {
  model_data <- results_df[results_df$Model == model, ]
  for(i in 1:length(metrics)) {
    row_data <- data.frame(
      Model = model,
      Metric = metrics[i],
      stringsAsFactors = FALSE
    )
    for(region in regions) {
      region_data <- model_data[model_data$Region == region, ]
      if(nrow(region_data) > 0) {
        row_data[[region]] <- region_data[[metrics[i]]]
      } else {
        row_data[[region]] <- NA
      }
    }
    
    final_table3 <- rbind(final_table3, row_data)
  }
}
rownames(final_table3) <- NULL
print(final_table3)
write.csv(final_table3, file.path(output_folder, "final_interaction_table.csv"), row.names = FALSE)

# Create all six interaction plots
osha_binary_interaction_plot <- create_interaction_plot(industry_region_interactions, "OSHA Violations", "binary")
wh_binary_interaction_plot <- create_interaction_plot(industry_region_interactions, "Wage & Hour Violations", "binary")
acc_binary_interaction_plot <- create_interaction_plot(industry_region_interactions, "Accident Rate", "binary")

osha_cont_interaction_plot <- create_interaction_plot(industry_region_interactions, "OSHA Violations", "continuous")
wh_cont_interaction_plot <- create_interaction_plot(industry_region_interactions, "Wage & Hour Violations", "continuous")
acc_cont_interaction_plot <- create_interaction_plot(industry_region_interactions, "Accident Rate", "continuous")

ggsave(filename = file.path(output_folder, "osha_binary_interaction_plot.png"), plot = osha_binary_interaction_plot, width = 10, height = 8, dpi = 300)
ggsave(filename = file.path(output_folder, "wh_binary_interaction_plot.png"), plot = wh_binary_interaction_plot, width = 10, height = 8, dpi = 300)
ggsave(filename = file.path(output_folder, "acc_binary_interaction_plot.png"), plot = acc_binary_interaction_plot, width = 10, height = 8, dpi = 300)

ggsave(filename = file.path(output_folder, "osha_cont_interaction_plot.png"), plot = osha_cont_interaction_plot, width = 10, height = 8, dpi = 300)
ggsave(filename = file.path(output_folder, "wh_cont_interaction_plot.png"), plot = wh_cont_interaction_plot, width = 10, height = 8, dpi = 300)
ggsave(filename = file.path(output_folder, "acc_cont_interaction_plot.png"), plot = acc_cont_interaction_plot, width = 10, height = 8, dpi = 300)


# Create heatmaps for all outcomes and both model types
osha_binary_heatmap <- create_interaction_heatmap(industry_region_interactions, "OSHA Violations", "binary")
wh_binary_heatmap <- create_interaction_heatmap(industry_region_interactions, "Wage & Hour Violations", "binary")
acc_binary_heatmap <- create_interaction_heatmap(industry_region_interactions, "Accident Rate", "binary")

osha_cont_heatmap <- create_interaction_heatmap(industry_region_interactions, "OSHA Violations", "continuous")
wh_cont_heatmap <- create_interaction_heatmap(industry_region_interactions, "Wage & Hour Violations", "continuous")
acc_cont_heatmap <- create_interaction_heatmap(industry_region_interactions, "Accident Rate", "continuous")









#---------Below is an increidbly stupid idea.

# Step 1: Calculate population density and categorize
filtered_smr_data <- filtered_smr_data %>%
  mutate(
    area_sq_m = as.numeric(st_area(geometry)),
    area_sq_km = area_sq_m / 1000000,
    pop_density_per_sqkm = ifelse(area_sq_km > 0, census_pop / area_sq_km, NA),
    density_category = case_when(
      pop_density_per_sqkm < 100 ~ "rural",      # Less than 100 people/sq km
      pop_density_per_sqkm < 500 ~ "suburban",   # 100-500 people/sq km  
      TRUE ~ "urban"                             # More than 500 people/sq km
    )
  )
print(table(filtered_smr_data$density_category, useNA = "ifany"))

# Step 2: Create spatial blocks within state-density groups
all_blocks <- list()
block_id <- 1
state_density_combos <- filtered_smr_data %>%
  st_drop_geometry() %>%
  group_by(state, density_category) %>%
  summarise(n_zips = n(), .groups = "drop") %>%
  arrange(desc(n_zips))

for(i in 1:nrow(state_density_combos)) {
  current_state <- state_density_combos$state[i]
  current_density <- state_density_combos$density_category[i]
  n_zips_in_group <- state_density_combos$n_zips[i]
  group_data <- filtered_smr_data %>% 
    filter(state == current_state, density_category == current_density)
  
  if(nrow(group_data) == 0) next
  
  base_size <- switch(current_density,
    "rural" = 12,      # Rural areas can have larger blocks
    "suburban" = 8,    # Suburban moderate blocks  
    "urban" = 6        # Urban areas smaller blocks
  )
  large_states <- c("CA", "TX", "NY", "PA", "FL", "OH", "MI", "IL", "NC", "NJ")
  small_states <- c("WV", "NH", "NE", "UT", "NV", "NM", "ID", "MT", "SD", "VT", 
                    "HI", "RI", "AK", "ND", "DE", "WY", "DC")
  if(current_state %in% large_states) {
    adjustment <- 0.8  # Smaller blocks for large states
  } else if(current_state %in% small_states) {
    adjustment <- 1.5  # Larger blocks for small states
  } else {
    adjustment <- 1.0  # No adjustment for medium states
  }
  target_block_size <- max(3, round(base_size * adjustment))
  if(n_zips_in_group <= target_block_size) {
    all_blocks[[block_id]] <- list(
      zip_codes = as.character(group_data$zip_id),  # Ensure character
      state = current_state,
      density = current_density,
      block_size = n_zips_in_group
    )
    block_id <- block_id + 1
  } else {
    n_blocks <- ceiling(n_zips_in_group / target_block_size)
    n_blocks <- min(n_blocks, n_zips_in_group)
    coords <- st_coordinates(st_centroid(group_data$geometry))
    tryCatch({
      set.seed(123)
      clusters <- kmeans(coords, centers = n_blocks, nstart = 20, iter.max = 100)
      for(cluster_num in 1:n_blocks) {
        cluster_indices <- which(clusters$cluster == cluster_num)
        if(length(cluster_indices) > 0) {
          cluster_zips <- as.character(group_data$zip_id[cluster_indices])  # Ensure character
          all_blocks[[block_id]] <- list(
            zip_codes = cluster_zips,
            state = current_state,
            density = current_density,
            block_size = length(cluster_zips)
          )
          block_id <- block_id + 1
        }
      }
      
    }, error = function(e) {
      all_blocks[[block_id]] <- list(
        zip_codes = as.character(group_data$zip_id),
        state = current_state,
        density = current_density,
        block_size = n_zips_in_group
      )
      block_id <- block_id + 1
    })
  }
}

block_sizes <- sapply(all_blocks, function(x) x$block_size)
print(summary(block_sizes))
all_zips_in_blocks <- unique(unlist(lapply(all_blocks, function(x) x$zip_codes)))
missing_zips <- setdiff(as.character(filtered_smr_data$zip_id), all_zips_in_blocks)
print(length(missing_zips))


# Step 3: Perform spatial bootstrap to estimate overall SMRs and confidence intervals
n_bootstrap <- 500
bootstrap_overall_smrs <- list()
combinations <- data.frame(
  outcome = c("vio", "wh", "acc", "vio", "wh", "acc", "vio", "wh", "acc"),
  denominator = c("estab", "estab", "estab", "workers", "workers", "workers", "pop", "pop", "pop"),
  stringsAsFactors = FALSE
)
combinations$combo_name <- paste(combinations$outcome, combinations$denominator, sep = "_")
combinations$description <- c(
  "Violations per 100 establishments",
  "Workplace hazards per 100 establishments", 
  "Accidents per 100 establishments",
  "Violations per 1000 workers",
  "Workplace hazards per 1000 workers",
  "Accidents per 1000 workers", 
  "Violations per 10000 population",
  "Workplace hazards per 10000 population",
  "Accidents per 10000 population"
)
get_columns <- function(outcome, denominator) {
  observed_col <- paste0(outcome, "_overal")
  if(outcome == "wh") observed_col <- "wh_overa_4"
  
  expected_col <- paste0("expected_", outcome, "_", denominator)
  
  return(list(observed = observed_col, expected = expected_col))
}

original_overall_smrs <- list()

for(i in 1:nrow(combinations)) {
  outcome <- combinations$outcome[i]
  denominator <- combinations$denominator[i]
  combo_name <- combinations$combo_name[i]
  description <- combinations$description[i]
  
  cols <- get_columns(outcome, denominator)
  
  if(all(c(cols$observed, cols$expected) %in% colnames(filtered_smr_data))) {
    total_observed <- sum(filtered_smr_data[[cols$observed]], na.rm = TRUE)
    total_expected <- sum(filtered_smr_data[[cols$expected]], na.rm = TRUE)
    
    if(total_expected > 0) {
      original_smr <- total_observed / total_expected
      original_overall_smrs[[combo_name]] <- original_smr
    }
  }
}

for(b in 1:n_bootstrap) {
  
  sampled_block_ids <- sample(1:length(all_blocks), size = length(all_blocks), replace = TRUE)
  
  bootstrap_zip_sample <- c()
  for(block_id in sampled_block_ids) {
    block_zips <- all_blocks[[block_id]]$zip_codes
    bootstrap_zip_sample <- c(bootstrap_zip_sample, block_zips)
  }

  bootstrap_data <- data.frame()
  for(zip_code in bootstrap_zip_sample) {
    zip_row <- filtered_smr_data[filtered_smr_data$zip_id == zip_code, ]
    if(nrow(zip_row) > 0) {
      bootstrap_data <- rbind(bootstrap_data, zip_row)
    }
  }

  bootstrap_smrs_this_iter <- list()
  
  for(i in 1:nrow(combinations)) {
    outcome <- combinations$outcome[i]
    denominator <- combinations$denominator[i]
    combo_name <- combinations$combo_name[i]
    
    cols <- get_columns(outcome, denominator)
    
    if(all(c(cols$observed, cols$expected) %in% colnames(bootstrap_data))) {
      total_observed <- sum(bootstrap_data[[cols$observed]], na.rm = TRUE)
      total_expected <- sum(bootstrap_data[[cols$expected]], na.rm = TRUE)
      
      if(total_expected > 0) {
        bootstrap_smr <- total_observed / total_expected
        bootstrap_smrs_this_iter[[combo_name]] <- bootstrap_smr
      }
    }
  }
  bootstrap_overall_smrs[[b]] <- bootstrap_smrs_this_iter
}

final_results <- data.frame(
  combination = character(),
  description = character(),
  original_smr = numeric(),
  bootstrap_mean = numeric(),
  bootstrap_sd = numeric(),
  lower_95_ci = numeric(),
  upper_95_ci = numeric(),
  ci_width = numeric(),
  cv = numeric(),
  n_bootstrap = integer(),
  stringsAsFactors = FALSE
)

for(combo_name in names(original_overall_smrs)) {
  bootstrap_values <- sapply(bootstrap_overall_smrs, function(x) x[[combo_name]])
  bootstrap_values <- bootstrap_values[!is.na(bootstrap_values)]
  
  if(length(bootstrap_values) >= 100) {  # Need reasonable number for CI
    original_smr <- original_overall_smrs[[combo_name]]
    description <- combinations$description[combinations$combo_name == combo_name]
    
    bootstrap_mean <- mean(bootstrap_values)
    bootstrap_sd <- sd(bootstrap_values)
    lower_ci <- quantile(bootstrap_values, 0.025)
    upper_ci <- quantile(bootstrap_values, 0.975)
    ci_width <- upper_ci - lower_ci
    cv <- bootstrap_sd / bootstrap_mean
    
    result_row <- data.frame(
      combination = combo_name,
      description = description,
      original_smr = original_smr,
      bootstrap_mean = bootstrap_mean,
      bootstrap_sd = bootstrap_sd,
      lower_95_ci = lower_ci,
      upper_95_ci = upper_ci,
      ci_width = ci_width,
      cv = cv,
      n_bootstrap = length(bootstrap_values),
      stringsAsFactors = FALSE
    )
    
    final_results <- rbind(final_results, result_row)
  }
}

final_results_with_context <- merge(combinations, final_results, by.x = "combo_name", by.y = "combination")
write.csv(final_results_with_context, file.path(output_folder, "spatial_bootstrap_overall_smrs.csv"), row.names = FALSE)



