# Functions for H3 outputs (06b).
# Attach dplyr (or tidyverse); glmmTMB is used through its namespace.

#' Summarise fixed effects, random effects and fit of an H3 model
#'
#' Extracts conditional-model estimates and diagnostics into a long table.
#'
#' @param model A fitted glmmTMB binomial logit model with habitat and cell_id
#'   grouping variables in its model frame.
#' @param model_name Label identifying the model in the output.
#'
#' @return A data frame with model, section, term, estimate, std_error,
#'   conf_low, conf_high, z_value and p_value. Inapplicable fields are NA.
#'
#' @details
#' Fixed effects receive normal-approximation 95% Wald intervals. Separate
#' rows report exponentiated estimates and bounds. Random-effect rows contain
#' variances, SDs and correlations, without confidence intervals. Fit rows
#' include sample size, grouping counts, likelihood, information criteria,
#' residual degrees of freedom, convergence code and Hessian status.
#' The function reports diagnostics without enforcing successful convergence.
#' Likelihood-ratio comparisons are calculated separately in script 06b.
#' Uses dplyr and glmmTMB; no global analysis objects are required.
extract_h3_results <- function(model, model_name) {
        
        model_summary <- summary(model)
        coefficients <- model_summary$coefficients$cond
        
        fixed_results <- data.frame(
                model = model_name,
                section = "Fixed effects",
                term = rownames(coefficients),
                estimate = coefficients[, "Estimate"],
                std_error = coefficients[, "Std. Error"],
                conf_low = coefficients[, "Estimate"] -
                        qnorm(0.975) * coefficients[, "Std. Error"],
                conf_high = coefficients[, "Estimate"] +
                        qnorm(0.975) * coefficients[, "Std. Error"],
                z_value = coefficients[, "z value"],
                p_value = coefficients[, "Pr(>|z|)"],
                row.names = NULL
        )
        
        # Report exponentiated coefficients in separate rows to keep tables narrow.
        odds_results <- fixed_results |>
                transmute(
                        model,
                        section = "Exponentiated fixed effects",
                        term,
                        estimate = exp(estimate),
                        conf_low = exp(conf_low),
                        conf_high = exp(conf_high)
                )
        
        # Extract variances, standard deviations, and correlations.
        random_components <- glmmTMB::VarCorr(model)$cond
        
        random_results <- lapply(names(random_components), function(group) {
                
                component <- random_components[[group]]
                component_sd <- attr(component, "stddev")
                component_cor <- attr(component, "correlation")
                
                results <- bind_rows(
                        data.frame(
                                model = model_name,
                                section = "Random effects",
                                term = paste(
                                        group, names(component_sd), "variance:"
                                ),
                                estimate = unname(component_sd^2)
                        ),
                        data.frame(
                                model = model_name,
                                section = "Random effects",
                                term = paste(
                                        group, names(component_sd), "SD:"
                                ),
                                estimate = unname(component_sd)
                        )
                )
                
                if (length(component_sd) > 1L) {
                        
                        pairs <- which(
                                upper.tri(component_cor),
                                arr.ind = TRUE
                        )
                        
                        results <- bind_rows(
                                results,
                                data.frame(
                                        model = model_name,
                                        section = "Random effects",
                                        term = paste(
                                                group,
                                                rownames(component_cor)[pairs[, 1]],
                                                "vs",
                                                colnames(component_cor)[pairs[, 2]],
                                                "correlation"
                                        ),
                                        estimate = component_cor[pairs]
                                )
                        )
                }
                
                results
        }) |>
                bind_rows()
        
        # Count observations and grouping levels in the fitted model frame.
        model_frame <- model.frame(model)
        
        fit_values <- c(
                "Observations" = nobs(model),
                "Habitat groups" = n_distinct(model_frame$habitat),
                "Cell groups" = n_distinct(model_frame$cell_id),
                "Estimated parameters" = attr(logLik(model), "df"),
                "Log-likelihood" = as.numeric(logLik(model)),
                "Minus twice log-likelihood" = -2 * as.numeric(logLik(model)),
                "AIC" = AIC(model),
                "BIC" = BIC(model),
                "Residual degrees of freedom" = df.residual(model),
                "Convergence code (0 = success)" = model$fit$convergence,
                "Positive-definite Hessian (1 = TRUE)" =
                        as.integer(model$sdr$pdHess)
        )
        
        fit_results <- data.frame(
                model = model_name,
                section = "Model fit and convergence",
                term = names(fit_values),
                estimate = unname(fit_values)
        )
        
        bind_rows(
                fixed_results,
                odds_results,
                random_results,
                fit_results
        )
}

#' Average a ten-percentage-point protection contrast at one occupancy
#'
#' Calculates the mean change in predicted occurrence-loss probability for
#' a protection increase of ten percentage points across a supplied grid.
#'
#' @param occupancy Positive initial habitat occupancy, expressed as cells.
#' @param protection_grid Baseline coverage in units of ten percentage points.
#'   Script 06b uses seq(0, 9, by = 0.1), equivalent to 0 through 90 percent.
#' @param occupancy_reference Mean log occupancy used to centre the predictor.
#' @param beta Named fixed-effect coefficient vector for the interaction model.
#' @param beta_draws Matrix of coefficient draws, with one draw per row and
#'   columns ordered as beta. Draws are generated once in script 06b.
#'
#' @return A one-row data frame with n_cells_t0, estimate, conf_low, conf_high.
#'   Contrasts and bounds are in percentage points.
#'
#' @details
#' Uses the fixed-effects design pa_10 * log_occupancy_c with random effects
#' set to zero. Each baseline grid value receives equal weight. The same
#' coefficient draw is used for baseline and increased protection, preserving
#' their covariance. Bounds are the 2.5th and 97.5th percentiles of simulated
#' mean contrasts. These are pointwise intervals; averaging is over the
#' supplied theoretical grid, not the observed protection distribution.
#' Uses base R and stats; the original contrast calculation is unchanged.
average_h3_contrast <- function(occupancy, protection_grid, occupancy_reference,
                                beta, beta_draws) {
        
        baseline <- data.frame(
                pa_10 = protection_grid,
                log_occupancy_c = log(occupancy) - occupancy_reference
        )
        
        increased <- baseline
        increased$pa_10 <- increased$pa_10 + 1
        
        X0 <- model.matrix(
                ~ pa_10 * log_occupancy_c,
                data = baseline
        )[, names(beta), drop = FALSE]
        
        X1 <- model.matrix(
                ~ pa_10 * log_occupancy_c,
                data = increased
        )[, names(beta), drop = FALSE]
        
        # Express probability differences in percentage points.
        estimate <- 100 * mean(
                plogis(X1 %*% beta) - plogis(X0 %*% beta)
        )
        
        simulated_contrasts <- 100 * colMeans(
                plogis(X1 %*% t(beta_draws)) -
                        plogis(X0 %*% t(beta_draws))
        )
        
        interval <- quantile(
                simulated_contrasts,
                probs = c(0.025, 0.975),
                names = FALSE
        )
        
        data.frame(
                n_cells_t0 = occupancy,
                estimate = estimate,
                conf_low = interval[1],
                conf_high = interval[2]
        )
}
