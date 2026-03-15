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

