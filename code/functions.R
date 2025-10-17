##===================================##
##             Functions             ##
##===================================##

#----Interactions
run_interaction_models <- function(industry_vars, outcome_var, model_data, outcome_type = "binary") {

  # This function runs the interaction models at the end of the analysis. It takes in:
  # ---------------
  # industry_vars: A vector of industry variable names to include in the model iterations.
  # outcome_var: The name of the outcome variable to model.
  # model_data: The data, which will either be the main model data or one of the positive subsets. Remember that the data subset should probs match the outcome!
  # outcome_type: The type of outcome variable, which can be "binary", "labor_continuous", "wh_continuous", or "acc_continuous".
  # ---------------

  # This function filters the data based on the outcome and data specified, and then runs a series of regression models in which each
  # industry is independently replaced and interacted with region, while also controlling for inspection rate. This model also includes 
  # spatial random effects via BYM specification. If the outcome variable is binary, it uses a binomial family; if it is continuous, it uses a lognormal family.
  # I want to use latex notation here to outline the model, but this is not Julia, and while the hipster in me wants to say something like 
  # "oh well, Julia is actually my preferred language of choice", it is limited and does kinda suck, but being able to use in-line latex 
  # notation is quite sweet. Anyway, you get it, the model is not specifying that much.

  # The function returns a data frame with the interaction coefficients, confidence intervals, "significance" (significance does not really even make sense here, because, you know,
  # Bayesian. But if no "*" then no publish. Go academia.) for each industry-region interaction, and then the total effect of the main effect
  # plus the interaction coef.

  results <- list()
  
  if (outcome_type == "binary") {
    family_type <- "binomial"
    zip_id <- "zip_id"
    model_adj <- "model_adj_binary"
  } else if (outcome_type == "labor_continuous") {
    family_type <- "lognormal"
    zip_id <- "zip_id_pos"
    model_adj <- "model_adj_labor_pos"
  } else if (outcome_type == "wh_continuous") {
    family_type <- "lognormal"
    zip_id <- "zip_id_pos"
    model_adj <- "model_adj_wh_pos"
  } else if (outcome_type == "acc_continuous") {
    family_type <- "lognormal"
    zip_id <- "zip_id_pos"
    model_adj <- "model_adj_acc_pos"
  } else {
    stop("Invalid outcome_type")
  }
  
  total_n <- nrow(model_data)
  
  region_counts <- model_data %>%
    group_by(region_8) %>%
    summarize(
      region_n = n(),
      outcome_n = sum(!!sym(outcome_var) > 0, na.rm = TRUE)
    )
  
  for(industry_var in industry_vars) {
    model_formula <- as.formula(paste0(
      outcome_var, " ~ ", industry_var, " * region_8 + insp_rate + ", 
      "f(", zip_id, ", model = 'bym', graph = ", model_adj, ")"
    ))
    
    interaction_model <- inla(
      model_formula,
      data = model_data,
      family = family_type,
      control.predictor = list(compute = TRUE),
      control.compute = list(dic = TRUE, waic = TRUE)
    )
    
    fixed_effects <- interaction_model$summary.fixed
    interaction_pattern <- paste0(industry_var, ":region_8")
    interaction_vars <- rownames(fixed_effects)[grep(interaction_pattern, rownames(fixed_effects))]
    
    if(length(interaction_vars) > 0) {
      interactions <- data.frame(
        industry = industry_var,
        interaction = interaction_vars,
        coefficient = fixed_effects[interaction_vars, "mean"],
        lower_ci = fixed_effects[interaction_vars, "0.025quant"],
        upper_ci = fixed_effects[interaction_vars, "0.975quant"],
        significant = (fixed_effects[interaction_vars, "0.025quant"] > 0) | 
                   (fixed_effects[interaction_vars, "0.975quant"] < 0),
        outcome = outcome_var
      )
      
      interactions$region <- gsub(paste0(industry_var, ":region_8"), "", interactions$interaction)
      main_effect <- fixed_effects[industry_var, "mean"]
      interactions$total_effect <- main_effect + interactions$coefficient
      
      results[[industry_var]] <- interactions
    }
  }
  
  combined_results <- bind_rows(results)
  
  attr(combined_results, "sample_size") <- list(
    total = total_n,
    outcome_var = outcome_var,
    outcome_type = outcome_type
  )
  
  return(combined_results)
}

#------Making interaction table
create_model_table <- function(interaction_data, violation_type, model_type){

  # This function creates a table of the interaction results for a given violation type and model type. It takes in:
  # ---------------
  # interaction_data: The data frame containing the interaction results from the `run_interaction_models` function.
  # violation_type: The type of violation to filter the data by (e.g., "labor", "wh", "acc").
  # model_type: The type of model to filter the data by (e.g., "binary", "labor_continuous", "wh_continuous", "acc_continuous").
  # ---------------
  # The function filters the interaction data based on the specified violation type and model type, then reshapes the data to create a wide-format table.
  # It returns a data frame where rows are industries, columns are regions, and cell values are the total effects with significance indicated by an asterisk.

  subset_data <- interaction_data %>%
    filter(violation_type == !!violation_type, model_type == !!model_type)

  subset_data <- subset_data%>%mutate(coef_sig = ifelse(significant == TRUE, paste0(as.character(round(total_effect, 2)), "*"), as.character(round(total_effect, 2))))

  wide_data <- subset_data %>%
    dplyr::select(industry_name, region, coef_sig) %>%
    pivot_wider(names_from = region, values_from = coef_sig)

  table_data <- as.data.frame(wide_data[,-1])
  rownames(table_data)<-wide_data$industry_name

  return(table_data)
}

#----Heatmap for correlations
create_interaction_heatmap <- function(interaction_data, violation_type, model_type) {
  
  # This function creates a heatmap of the interaction results for a given violation type and model type. It takes in:
  # ---------------
  # interaction_data: The data frame containing the interaction results from the `run_interaction_models` function.
  # violation_type: The type of violation to filter the data by (e.g., "labor", "wh", "acc").
  # model_type: The type of model to filter the data by (e.g., "binary", "labor_continuous", "wh_continuous", "acc_continuous").
  # ---------------
  # The function filters the interaction data based on the specified violation type and model type, then reshapes the data to create a matrix suitable for heatmap visualization.
  # It returns a heatmap where rows are industries, columns are regions, and cell values are the total effects. It is, you know, pretty much what happens
  # above, but in a heatmap format.
  
  subset_data <- interaction_data %>%
    filter(violation_type == !!violation_type, model_type == !!model_type)
  
  wide_data <- subset_data %>%
    dplyr::select(industry_name, region, total_effect) %>%
    pivot_wider(names_from = region, values_from = total_effect)
  
  mat_data <- as.matrix(wide_data[, -1])
  rownames(mat_data) <- wide_data$industry_name
  
  heatmap(mat_data, 
         col = colorRampPalette(c("blue", "white", "red"))(100),
         main = paste0(violation_type, " - ", 
                      ifelse(model_type == "binary", "Presence", "Volume"), 
                      " by Region"),
         xlab = "Region", 
         ylab = "Industry",
         scale = "none")
  
  return(mat_data)
}

#----Interaction plot
create_interaction_plot <- function(interaction_data, violation_type, model_type) {

  # This function creates a bar plot of the interaction results for a given violation type and model type. It takes in:
  # ---------------
  # interaction_data: The data frame containing the interaction results from the `run_interaction_models` function.
  # violation_type: The type of violation to filter the data by (e.g., "labor", "wh", "acc").
  # model_type: The type of model to filter the data by (e.g., "binary", "labor_continuous", "wh_continuous", "acc_continuous").
  # ---------------
  # The function filters the interaction data based on the specified violation type and model type, then creates a bar plot with error bars.
  # It returns a ggplot object where the x-axis is regions, y-axis is interaction coefficients, and bars are colored by significance.

  subset_data <- interaction_data %>%
    filter(violation_type == !!violation_type, model_type == !!model_type)
  
  subset_data$significance <- ifelse(subset_data$significant, "Significant", "Not Significant")
  
  plot_title <- paste0(violation_type, " - ", 
                      ifelse(model_type == "binary", "Presence", "Volume"), 
                      " Industry-Region Effects")
  
  ggplot(subset_data, aes(x = reorder(region, coefficient), y = coefficient, fill = significance)) +
    geom_col() +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
    facet_wrap(~ industry_name, scales = "free_y") +
    coord_flip() +
    scale_fill_manual(values = c("Significant" = "darkred", "Not Significant" = "gray70")) +
    theme_minimal() +
    labs(title = plot_title,
         subtitle = "Effect of industry concentration varies by region",
         x = "Region",
         y = "Interaction Coefficient",
         fill = "Statistical\nSignificance")
}

#---Correlation Matrix Heatmap
create_heatmap <- function(mat, title) {
  if(nrow(mat) <= 20 && ncol(mat) <= 20) {
    corrplot(mat, is.corr = FALSE, method = "color",
             tl.col = "black", tl.srt = 45, addCoef.col = "black",
             title = title,
             mar = c(0,0,1,0))
  } else {
    heatmap(mat,
            col = colorRampPalette(c("blue", "white", "red"))(100),
            main = title,
            xlab = "Industry", 
            ylab = "Region",
            scale = "none")
  }
  
  return(mat)
}

#----Function to take in outcomes and produce bivariate risk maps
industry_risk_overlap_map <- function(gdf, industry_var, 
                                     outcome_x = "vio_rate", 
                                     x_label = "OSHA Violations",
                                     industry_label = "Industry",
                                     zero_threshold = 0.001) {

  # This function creates bivariate risk maps comparing an industry concentration variable with an outcome variable. It takes in:
  # ---------------
  # gdf: A spatial data frame containing the geographic data and variables.
  # industry_var: The name of the industry variable to compare. Here it is the rate of establishment per 1000 workers in a given industry.
  # outcome_x: The name of the outcome variable to compare. Default is "vio_rate" for OSHA violations.
  # x_label: The label for the x-axis in the legend. Default is "OSHA Violations".
  # industry_label: The label for the y-axis in the legend. Default is "Industry".
  # zero_threshold: A threshold value to define "Low" category. Default is 0.001.
  # ---------------
  # The function categorizes both the industry and outcome variables into tertiles (Low, Medium, High) based on the specified threshold and quantiles.
  # It then creates a bivariate category by combining these tertiles and generates a map with a corresponding legend. Here, tertiles are low when 
  # the value is less than or equal to the zero_threshold, medium when it is between the zero_threshold and the median of non-zero values, and high when it is above the median of non-zero values.
  # It then simply applies a categorical color scheme to the bivariate categories and plots the map. The correlation here is merely descirptive and not actual correlation coefficients.
                                 
  ind_vals <- gdf[[industry_var]]
  ind_vals <- as.numeric(as.character(ind_vals))
  
  zero_threshold <- zero_threshold

  gdf$ind_tertile <- case_when(
    ind_vals <= zero_threshold ~ "Low",
    ind_vals <= quantile(ind_vals[ind_vals > zero_threshold], .5, na.rm = TRUE) ~ "Medium",
    TRUE ~ "High"
  )

  x_vals <- gdf[[outcome_x]]
  
  gdf$x_tertile <- case_when(
    x_vals <= zero_threshold ~ "Low",
    x_vals <= quantile(x_vals[x_vals > zero_threshold], .5, na.rm = TRUE) ~ "Medium",
    TRUE ~ "High"
  )

  gdf$ind_tertile[is.na(gdf$ind_tertile)] <- "Low"
  gdf$x_tertile[is.na(gdf$x_tertile)] <- "Low"
  
  gdf$bivariate_cat <- paste(gdf$x_tertile, gdf$ind_tertile, sep = "-")
  
  bivariate_colors <- c(
    "Low-Low" = "#f7f7f7",
    "Low-Medium" = "#dfc8e7",
    "Low-High" = "#9970ab",
    "Medium-Low" = "#c4e6e1",
    "Medium-Medium" = "#9ebcda",
    "Medium-High" = "#6059a9",
    "High-Low" = "#67c6c1",
    "High-Medium" = "#5395b4",
    "High-High" = "#3c4e8d"
  )
  
  ind_name <- gsub("ESTAB_|_P|_1", "", industry_var)
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
    "72" = "Accommodation & Food"
  )
  industry_name <- naics_descriptions[ind_name]
  if(is.na(industry_name)) industry_name <- ind_name
  
  legend_data <- expand.grid(
    x = c("Low", "Medium", "High"),
    y = c("Low", "Medium", "High")
  )
  legend_data$bivariate_cat <- paste(legend_data$x, legend_data$y, sep = "-")
  
  legend_plot <- ggplot() +
    geom_tile(data = legend_data, aes(x = x, y = y, fill = bivariate_cat)) +
    scale_fill_manual(values = bivariate_colors, guide = "none") +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ) +
    labs(
      x = x_label,
      y = industry_label
    )
  
  map <- ggplot(gdf) +
    geom_sf(aes(fill = bivariate_cat), color = NA) +
    scale_fill_manual(values = bivariate_colors, guide = "none") +
    theme_minimal() +
    labs(title = paste("Bivariate Risk Map:", industry_name),
         subtitle = paste(x_label, "vs.", industry_label)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  map_with_legend <- cowplot::plot_grid(
    map, 
    legend_plot, 
    ncol = 2,
    rel_widths = c(3, 1)
  )
  
  return(map_with_legend)
}


#----Function to combine all the bivariate maps
create_all_outcome_maps_update <- function(gdf, industry_var, 
                                           zero_threshold = 0.001) {

  # This function creates all three bivariate risk maps for a given industry variable. It takes in:
  # ---------------
  # gdf: A spatial data frame containing the geographic data and variables.
  # industry_var: The name of the industry variable to compare. Here it is the rate of establishment per 1000 workers in a given industry.
  # zero_threshold: A threshold value to define "Low" category. Default is 0.001.
  # ---------------
  # The function calls the `industry_risk_overlap_map` function three times to create maps comparing:
  # 1. OSHA Violations vs. Industry Concentration
  # 2. Wage & Hour Violations vs. Industry Concentration
  # 3. Accident Rate vs. Industry Concentration

  
  # OSHA vs. Wage & Hour
  vio_ind_map <- industry_risk_overlap_map(
    gdf, industry_var, 
    outcome_x = "vio_rate", 
    x_label = "OSHA Violations",
    industry_label = "Industry"
  )
  
  # OSHA vs. Accident
  wage_ind_map <- industry_risk_overlap_map(
    gdf, industry_var, 
    outcome_x = "wh_rate_pe", 
    x_label = "OSHA Violations",
    industry_label = "Industry"
  )
  
  # Wage & Hour vs. Accident
  acc_ind_map <- industry_risk_overlap_map(
    gdf, industry_var, 
    outcome_x = "acc_rate", 
    x_label = "Wage & Hour Violations",
    industry_label = "Industry"
  )
  
  # Get industry name for the title
  ind_name <- gsub("ESTAB_|_P|_1", "", industry_var)
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
    "72" = "Accommodation & Food"
  )
  industry_name <- naics_descriptions[ind_name]
  if(is.na(industry_name)) industry_name <- ind_name
  
  return(list(
    vio = vio_ind_map,
    wage = wage_ind_map,
    acc = acc_ind_map
  ))
}

#----Function to plot industry rates by region
plot_industry_rates_by_region <- function(model_data, labor_positive, wh_positive, acc_positive, industry_var) {

  # This function creates bar plots of industry rates and outcome rates by region. It takes in:
  # ---------------
  # model_data: The main data frame containing all data.
  # labor_positive: A subset of the data containing only records with positive OSHA violation rates.
  # wh_positive: A subset of the data containing only records with positive Wage & Hour violation rates.
  # acc_positive: A subset of the data containing only records with positive Accident rates.
  # industry_var: The name of the industry variable to analyze.
  # ---------------
  # The function calculates the percentage of establishments in the specified industry and the rates of OSHA violations, Wage & Hour violations, and Accident rates by region.
  # It then creates two bar plots: one for occurrence rates (binary) and another for severity rates (continuous) among places with positive outcomes.
  # The function returns a combined plot with both bar plots side by side, along with sample sizes as an attribute. Not much is going to be used in the actual paper and
  # accomplishes little more than a descriptives table, but bar plots are fun and it is here to look through if industry/region descriptives want to be visualized
  # as opposed to simply tabulated.

    ind_name <- gsub("ESTAB_|_P|_1", "", industry_var)
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
        "72" = "Accommodation & Food"
    )
    industry_name <- naics_descriptions[ind_name]
    if(is.na(industry_name)) industry_name <- ind_name
    
    binary_data <- model_data %>%
        dplyr::select(region_8, !!sym(industry_var), labor_binary, wh_binary, acc_binary) %>%
        group_by(region_8) %>%
        summarize(
            ind_percent = mean(!!sym(industry_var), na.rm = TRUE),
            labor_rate = mean(labor_binary, na.rm = TRUE),
            wh_rate = mean(wh_binary, na.rm = TRUE),
            acc_rate = mean(acc_binary, na.rm = TRUE),
            n = n(),
            labor_n = sum(labor_binary, na.rm = TRUE),
            wh_n = sum(wh_binary, na.rm = TRUE),
            acc_n = sum(acc_binary, na.rm = TRUE)
        )
    
    binary_data_long <- binary_data %>%
        pivot_longer(cols = c(labor_rate, wh_rate, acc_rate), 
                    names_to = "outcome_type", 
                    values_to = "rate") %>%
        mutate(
            outcome_type = case_when(
                outcome_type == "labor_rate" ~ "OSHA Violations",
                outcome_type == "wh_rate" ~ "Wage & Hour Violations",
                outcome_type == "acc_rate" ~ "Accident Rate"
            ),
            count = case_when(
                outcome_type == "OSHA Violations" ~ labor_n,
                outcome_type == "Wage & Hour Violations" ~ wh_n,
                outcome_type == "Accident Rate" ~ acc_n
            ),
            label = paste0(round(rate*100, 1), "% (", count, ")")
        )
    
    binary_plot <- ggplot(binary_data_long, aes(x = reorder(region_8, ind_percent), 
                                              y = rate, fill = outcome_type)) +
        geom_col(position = "dodge") +
        geom_text(aes(label = label), position = position_dodge(width = 0.9), 
                  hjust = -0.1, size = 2.5, angle = 0) +
        coord_flip() +
        scale_fill_manual(values = c("OSHA Violations" = "cyan3", 
                                    "Wage & Hour Violations" = "magenta3",
                                    "Accident Rate" = "yellow3")) +
        theme_minimal() +
        labs(title = "Occurrence Rate by Region",
            subtitle = paste0("Total N=", sum(binary_data$n), 
                             " | Probability of violations/accidents"),
            x = NULL,
            y = "Occurrence Rate",
            fill = "Outcome Type") +
        theme(legend.position = "top")
    
    # Labor violation rates among places with violations
    labor_positive_region <- labor_positive %>%
        dplyr::select(region_8, !!sym(industry_var), vio_rate) %>%
        group_by(region_8) %>%
        summarize(
            labor_rate = mean(vio_rate, na.rm = TRUE),
            n_labor = n()
        )
    
    # Wage & Hour violation rates among places with violations
    wh_positive_region <- wh_positive %>%
        dplyr::select(region_8, !!sym(industry_var), wh_rate_pe) %>%
        group_by(region_8) %>%
        summarize(
            wh_rate = mean(wh_rate_pe, na.rm = TRUE),
            n_wh = n()
        )
    
    # Accident rates among places with accidents
    acc_positive_region <- acc_positive %>%
        dplyr::select(region_8, !!sym(industry_var), acc_rate) %>%
        group_by(region_8) %>%
        summarize(
            acc_rate = mean(acc_rate, na.rm = TRUE),
            n_acc = n()
        )
    
    # Industry concentrations
    industry_by_region <- model_data %>%
        dplyr::select(region_8, !!sym(industry_var)) %>%
        group_by(region_8) %>%
        summarize(
            ind_percent = mean(!!sym(industry_var), na.rm = TRUE),
            n = n()
        )
    
    continuous_data <- industry_by_region %>%
        left_join(labor_positive_region, by = "region_8") %>%
        left_join(wh_positive_region, by = "region_8") %>%
        left_join(acc_positive_region, by = "region_8")
    
    total_labor <- sum(labor_positive_region$n_labor, na.rm = TRUE)
    total_wh <- sum(wh_positive_region$n_wh, na.rm = TRUE)
    total_acc <- sum(acc_positive_region$n_acc, na.rm = TRUE)
    
    # Convert to long format for plotting
    continuous_data_long <- continuous_data %>%
        pivot_longer(cols = c(labor_rate, wh_rate, acc_rate), 
                    names_to = "outcome_type", 
                    values_to = "rate") %>%
        mutate(
            outcome_type = case_when(
                outcome_type == "labor_rate" ~ "OSHA Violations",
                outcome_type == "wh_rate" ~ "Wage & Hour Violations",
                outcome_type == "acc_rate" ~ "Accident Rate"
            ),
            count = case_when(
                outcome_type == "OSHA Violations" ~ n_labor,
                outcome_type == "Wage & Hour Violations" ~ n_wh,
                outcome_type == "Accident Rate" ~ n_acc
            ),
            label = paste0(round(rate, 2), " (", count, ")")
        )
    
    region_order <- levels(reorder(binary_data$region_8, binary_data$ind_percent))
    continuous_data_long$region_8 <- factor(continuous_data_long$region_8, levels = region_order)
    
    continuous_plot <- ggplot(continuous_data_long, aes(x = region_8, 
                                                        y = rate, fill = outcome_type)) +
        geom_col(position = "dodge") +
        geom_text(aes(label = label), position = position_dodge(width = 0.9), 
                  hjust = -0.1, size = 2.5, angle = 0) +
        coord_flip() +
        scale_fill_manual(values = c("OSHA Violations" = "cyan3", 
                                    "Wage & Hour Violations" = "magenta3",
                                    "Accident Rate" = "yellow3")) +
        theme_minimal() +
        labs(title = "Severity Rate by Region",
            subtitle = paste0("Violations: OSHA=", total_labor, 
                             ", W&H=", total_wh, 
                             ", Acc=", total_acc),
            x = "Region",
            y = "Severity Rate",
            fill = "Outcome Type") +
        theme(legend.position = "top")
    
    combined_plot <- cowplot::plot_grid(
        binary_plot, continuous_plot, 
        ncol = 2,
        rel_widths = c(1, 1),
        align = "h"
    )
    
    title <- cowplot::ggdraw() + 
        cowplot::draw_label(
            paste("Regional Risk Profile for", industry_name),
            fontface = "bold",
            x = 0.5,
            hjust = 0.5,
            size = 14
        )
    
    final_plot <- cowplot::plot_grid(
        title, combined_plot,
        ncol = 1,
        rel_heights = c(0.1, 1)
    )
    
    attr(final_plot, "sample_sizes") <- list(
        total = nrow(model_data),
        labor_positive = nrow(labor_positive),
        wh_positive = nrow(wh_positive),
        acc_positive = nrow(acc_positive),
        by_region = binary_data %>% dplyr::select(region_8, n, labor_n, wh_n, acc_n)
    )
    
    return(final_plot)
}

#----Creating a correlation matrix
create_correlation_matrix <- function(data, correlation_column) {

  # This function creates a correlation matrix based on industry concentration and outcome. It takes in:
  # ---------------
  # data: The data frame containing the interaction results from the `run_interaction_models` function.
  # correlation_column: The name of the column to use for the correlation values (e.g., "total_effect").
  # ---------------
  # The function reshapes the data to create a matrix where rows are regions, columns are industries, and cell values are the specified correlation values.
  # It returns a matrix suitable for heatmap visualization.

  required_cols <- c("region", "industry", correlation_column)
  missing_cols <- setdiff(required_cols, names(data))
  
  if(length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  data[[correlation_column]][is.na(data[[correlation_column]])] <- 0
  
  if(length(unique(data$region)) < 2 || length(unique(data$industry)) < 2) {
    stop("Need at least 2 regions and 2 industries for heatmap")
  }
  
  regions <- unique(data$region)
  industries <- unique(data$industry)
  
  mat <- matrix(NA, nrow = length(regions), ncol = length(industries))
  rownames(mat) <- regions
  colnames(mat) <- industries
  
  for(i in 1:nrow(data)) {
    r <- data$region[i]
    ind <- data$industry[i]
    val <- data[[correlation_column]][i]
    
    r_idx <- match(r, regions)
    i_idx <- match(ind, industries)
    
    if(!is.na(r_idx) && !is.na(i_idx)) {
      mat[r_idx, i_idx] <- val
    }
  }
  
  return(mat)
}