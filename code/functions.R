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
  model_data$region_8 <- as.factor(model_data$region_8)
  
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
    group_by(region_8) %>% summarize(region_n = n(),outcome_n = sum(!!sym(outcome_var) > 0, na.rm = TRUE))
  
  for(industry_var in industry_vars) {
    model_formula <- as.formula(paste0(
      outcome_var, " ~ ", industry_var, " * region_8 + insp_rate + ", 
      "f(", zip_id, ", model = 'bym', graph = ", model_adj, ")"
    ))

    region_levels <- levels(model_data$region_8)[-1]
    lincomb_list <- list()
    lincomb_list[["reference_region_effect"]] <- inla.make.lincomb(setNames(list(1), industry_var))
    for (region in region_levels) {
      lincomb_name <- paste0("interaction_", region)
      interaction_term <- paste0(industry_var, ":region_8", region)
      lincomb_list[[lincomb_name]] <- inla.make.lincomb(setNames(list(1,1), c(industry_var, interaction_term)))
    }
    all_lincomb <- do.call(c, lincomb_list)

    interaction_model <- inla(
      model_formula,
      data = model_data,
      family = family_type,
      control.predictor = list(compute = TRUE),
      control.compute = list(dic = TRUE, waic = TRUE),
      lincomb = all_lincomb
    )
    
    lincomb_results <- interaction_model$summary.lincomb.derived
    

    interactions_df <- data.frame(
      industry = industry_var,
      lc_name = rownames(lincomb_results),
      interaction = c(NA, paste0(industry_var, ":region_8", region_levels)),
      region = c(levels(model_data$region_8)[1], region_levels),
      coefficient = lincomb_results[, "mean"],
      lower_ci = lincomb_results[, "0.025quant"],
      upper_ci = lincomb_results[, "0.975quant"],
      outcome = outcome_var
    )

    interactions_df$significant <- (interactions_df$lower_ci > 0) | (interactions_df$upper_ci < 0)
    interactions_df$ref_region <- interactions_df$lc_name == "reference_region_effect"
    results[[industry_var]] <- interactions_df
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
create_model_table <- function(interaction_data){

  # This function creates a table of the interaction results for a given violation type and model type. It takes in:
  # ---------------
  # interaction_data: The data frame containing the interaction results from the `run_interaction_models` function.
  # violation_type: The type of violation to filter the data by (e.g., "labor", "wh", "acc").
  # model_type: The type of model to filter the data by (e.g., "binary", "labor_continuous", "wh_continuous", "acc_continuous").
  # ---------------
  # The function filters the interaction data based on the specified violation type and model type, then reshapes the data to create a wide-format table.
  # It returns a data frame where rows are industries, columns are regions, and cell values are the total effects with significance indicated by an asterisk.

  data_split <- interaction_data%>%mutate(
    cell_str = paste0(round(coefficient, 2), ifelse(significant, " *", ""),
    "\n", "(", round(lower_ci, 2), ", ", round(upper_ci, 2), ")"))

  wide_data <- data_split %>%
    dplyr::select(industry, region, model_type, cell_str) %>%
    pivot_wider(names_from = c(region, model_type), values_from = cell_str)
  

  # full_wide <- bind_rows(coef_wide, ci_wide) %>%
  #   mutate(row = factor(row, levels = c("coef", "ci"))) %>%
  #   arrange(industry, row) %>%
  #   group_by(industry) %>%
  #   mutate(industry = ifelse(row == "coef", industry, "")) %>%
  #   ungroup() %>%
  #   dplyr::select(-row) %>%
  #   as.data.frame()
  return(as.data.frame(wide_data))
}
