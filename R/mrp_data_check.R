#' Check Data Consistency for MRP Models
#'
#' Validates the consistency between a model formula, survey data, and tabulation
#' tables used in multilevel regression with poststratification (MRP) analysis.
#' This function examines predictor variables to ensure data types match and 
#' identifies missing or additional categories across datasets.
#'
#' @param model A fitted model object or formula containing the variables to check.
#'   The function extracts variable names using \code{all.vars(formula(model))}.
#' @param survey_data A data frame containing the survey data with predictor variables.
#' @param table_list A list of data frames representing tabulation tables, each
#'   containing predictor variables that should match the survey data structure.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{tab}{Integer indicating which table from \code{table_list} (NA for variables not found in any table)}
#'   \item{variable}{Character name of the predictor variable}
#'   \item{k_svy}{Number of unique categories in the survey data}
#'   \item{k_tab}{Number of unique categories in the tabulation table}
#'   \item{k_tab_miss}{Number of categories present in survey but missing from table}
#'   \item{k_tab_add}{Number of categories present in table but missing from survey}
#'   \item{k_tab_NA}{Logical indicating if the table variable contains NA values}
#'   \item{k_svy_NA}{Logical indicating if the survey variable contains NA values}
#'   \item{tab_miss}{List column containing categories missing from the table}
#'   \item{tab_add}{List column containing additional categories in the table}
#' }
#'
#' @details
#' The function performs the following checks:
#' \itemize{
#'   \item Extracts predictor variables from the model formula (excluding random effects notation)
#'   \item Compares data types between survey and tabulation data
#'   \item Counts unique categories for each variable in both datasets
#'   \item Identifies missing and additional categories
#'   \item Checks for presence of NA values
#'   \item Reports variables that appear in the model but not in any tabulation table
#' }
#'
#' For factor variables, the function uses \code{levels()} to determine categories.
#' For other variable types, it uses unique non-NA values.
#'
#' @section Errors:
#' The function will stop with an error if data types don't match between
#' survey data and tabulation tables for the same variable.
#'
#' @examples
#' \dontrun{
#' # Assuming you have a fitted model, survey data, and tabulation tables
#' model <- glmer(outcome ~ age + gender + (1|state), data = survey_data, family = binomial)
#' tables <- list(table1, table2, table3)
#' 
#' # Check data consistency
#' check_results <- mrp_data_check(model, survey_data, tables)
#' print(check_results)
#' 
#' # Identify problematic variables
#' problems <- subset(check_results, k_tab_miss > 0 | k_tab_add > 0)
#' }
#'
#' @seealso \code{\link{formula}}, \code{\link{all.vars}}
#' @export
#' @importFrom stats formula
mrp_data_check <- function(model, survey_data, table_list) {
    
    # Extract variable names from model formula
    formula_vars <- all.vars(formula(model))
    response_var <- formula_vars[1]
    predictor_vars <- formula_vars[-1]
    
    # Remove random effect grouping variables notation
    predictor_vars <- gsub("\\|.*", "", predictor_vars)
    predictor_vars <- unique(trimws(predictor_vars))
    
    # Initialize results data frame
    result_df <- data.frame(
        tab = integer(0),
        variable = character(0),
        k_svy = integer(0),
        k_tab = integer(0),
        k_tab_miss = integer(0),
        k_tab_add = integer(0),
        k_tab_NA = logical(0),
        k_svy_NA = logical(0),
        tab_miss = I(list()),
        tab_add = I(list()),
        stringsAsFactors = FALSE
    )
    
    # Process each table separately
    for (i in seq_along(table_list)) {
        tab_table <- table_list[[i]]
        
        # Check each variable in this table
        for (var in names(tab_table)) {
            
            # Only process predictor variables
            if (var %in% predictor_vars) {
                
                # Get survey variable info
                svy_var <- survey_data[[var]]
                svy_type <- typeof(svy_var)
                
                # Get tab variable info
                tab_var <- tab_table[[var]]
                tab_type <- typeof(tab_var)
                
                # Check type consistency
                if (svy_type != tab_type) {
                    stop(paste("Type mismatch for variable", var, "- Survey:", svy_type, 
                              "vs Tab table", i, ":", tab_type))
                }
                
                # Get survey categories (levels for factors, unique values for others)
                if (is.factor(svy_var)) {
                    svy_categories <- levels(svy_var)
                } else {
                    svy_categories <- sort(unique(svy_var[!is.na(svy_var)]))
                }
                
                # Get tab categories
                if (is.factor(tab_var)) {
                    tab_categories <- levels(tab_var)
                } else {
                    tab_categories <- sort(unique(tab_var[!is.na(tab_var)]))
                }
                
                # Calculate metrics
                k_svy <- length(svy_categories)
                k_tab <- length(tab_categories)
                
                # Check for NA values
                k_svy_NA <- any(is.na(svy_var))
                k_tab_NA <- any(is.na(tab_var))
                
                # Find missing and additional categories
                missing_in_tabs <- setdiff(svy_categories, tab_categories)
                additional_in_tabs <- setdiff(tab_categories, svy_categories)
                
                k_tab_miss <- length(missing_in_tabs)
                k_tab_add <- length(additional_in_tabs)
                
                # Convert to lists for storage in data frame
                tab_miss <- if (k_tab_miss > 0) list(missing_in_tabs) else list(NA)
                tab_add <- if (k_tab_add > 0) list(additional_in_tabs) else list(NA)
                
                # Add to result data frame
                result_df <- rbind(result_df, data.frame(
                    tab = i,
                    variable = var,
                    k_svy = k_svy,
                    k_tab = k_tab,
                    k_tab_miss = k_tab_miss,
                    k_tab_add = k_tab_add,
                    k_tab_NA = k_tab_NA,
                    k_svy_NA = k_svy_NA,
                    tab_miss = I(tab_miss),
                    tab_add = I(tab_add),
                    stringsAsFactors = FALSE
                ))
            }
        }
    }
    
    # Add rows for predictor variables not found in any tabs
    variables_in_tabs <- unique(result_df$variable)
    missing_variables <- setdiff(predictor_vars, variables_in_tabs)
    
    for (var in missing_variables) {
        result_df <- rbind(result_df, data.frame(
            tab = NA,
            variable = var,
            k_svy = NA,
            k_tab = NA,
            k_tab_miss = NA,
            k_tab_add = NA,
            k_tab_NA = NA,
            k_svy_NA = NA,
            tab_miss = I(list(NA)),
            tab_add = I(list(NA)),
            stringsAsFactors = FALSE
        ))
    }
    
    return(result_df)
}

