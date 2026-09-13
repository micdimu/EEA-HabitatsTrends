# Shared functions for H2 main results (05a) and supplementary analyses (05b).
# Attach tidyverse in the calling scripts: dplyr, tidyr, tibble and ggplot2
# functions are used without namespaces. scales is called with ::.
# Function definitions are unchanged from the supplied scripts.

#' Predict H2 responses with a supplied covariance matrix
#'
#' Calculates fitted probabilities and pointwise 95% confidence intervals
#' across the observed protection range of each habitat group.
#'
#' @param model A fitted logit-link GLM using total_pa_land_cov and taxGroup.
#' @param V Coefficient covariance matrix in the same order as coef(model).
#'   Supply the cell-cluster-robust matrix for the main H2 analysis.
#' @param data Model data containing taxGroup and total_pa_land_cov, with
#'   finite coverage values and factor coding consistent with the model.
#'
#' @return A tibble containing taxGroup, total_pa_land_cov, estimate, lower
#'   and upper, with 101 coverage values per group.
#'
#' @details
#' Uses the model matrix and supplied covariance to calculate uncertainty
#' on the logit scale. Bounds use the standard normal 97.5th percentile and
#' are transformed with the inverse logit. The function does not estimate
#' the covariance matrix or apply a gain-loss balance transformation.
predict_h2 <- function(model, V, data) {
        
        # Predict within the observed coverage range of each group.
        newdata <- data |>
                group_by(taxGroup) |>
                summarise(
                        coverage = list(seq(
                                min(total_pa_land_cov),
                                max(total_pa_land_cov),
                                length.out = 101
                        )),
                        .groups = "drop"
                ) |>
                tidyr::unnest_longer(
                        coverage, values_to = "total_pa_land_cov"
                )
        
        X <- model.matrix(
                delete.response(terms(model)),
                data = newdata,
                contrasts.arg = model$contrasts
        )
        
        eta <- drop(X %*% coef(model))
        se <- sqrt(pmax(0, rowSums((X %*% V) * X)))
        
        newdata |>
                mutate(
                        estimate = plogis(eta),
                        lower = plogis(eta - qnorm(0.975) * se),
                        upper = plogis(eta + qnorm(0.975) * se)
                )
}

#' Plot H2 predictions by habitat group
#'
#' Draws fitted curves and confidence ribbons using the three-group palette.
#'
#' @param predictions Data frame containing total_pa_land_cov, taxGroup,
#'   estimate, lower and upper, normally returned by predict_h2.
#' @param y_label Vertical axis label.
#'
#' @return A ggplot object with a bottom legend and percentage response axis.
#'
#' @details
#' Colours follow the ordering of the three taxGroup levels. No model is
#' fitted and no balance transformation is applied. In the main script,
#' balance estimates and bounds are transformed beforehand; the response
#' scale is then replaced by a numeric scale ranging from -1 to 1.
plot_h2 <- function(predictions, y_label) {
        ggplot(
                predictions,
                aes(total_pa_land_cov, estimate, colour = taxGroup, fill = taxGroup)
        ) +
                geom_ribbon(
                        aes(ymin = lower, ymax = upper),
                        alpha = 0.15, colour = NA
                ) +
                geom_line(linewidth = 0.9) +
                scale_colour_manual(values = c("#C79536", "#28734B", "#80639D")) +
                scale_fill_manual(values = c("#C79536", "#28734B", "#80639D")) +
                scale_y_continuous(labels = scales::label_percent()) +
                labs(
                        x = "Protected terrestrial area (%)",
                        y = y_label,
                        colour = NULL,
                        fill = NULL
                ) +
                theme_classic(base_size = 12) +
                theme(legend.position = "bottom")
}

#' Model habitat gain in relation to loss in another group
#'
#' Fits a land-weighted quasibinomial logit model for the occurrence of any
#' gain, tests its loss-by-protection interaction and prepares a plot.
#'
#' @param gain_group Prefix identifying the response group in wide data.
#' @param loss_group Prefix identifying the predictor group in wide data.
#' @param gain_label,loss_label Display labels for the two groups.
#' @param data Wide data containing CellCode, total_pa_land_cov, land_weight,
#'   <gain_group>_gained, <loss_group>_lost and both groups' _initial columns.
#'
#' @return A list containing model, filtered data, predictions, plot and
#'   a one-row test tibble with response, predictor, n_cells,
#'   interaction_F and interaction_p.
#'
#' @details
#' Removes incomplete records and requires finite coverage, positive finite
#' weights and initial presence of the loss group. Initial absence of the
#' gain group remains eligible. Gain and loss are converted to binary states.
#' Both response and predictor must vary in the retained sample.
#'
#' The model includes loss * coverage plus initial richness of both groups.
#' An F test compares it with its additive counterpart on the same sample.
#' Predictions use 150 coverage values within the range shared by both loss
#' categories, holding initial richness at model-sample medians. Pointwise
#' mean-confidence bounds use a residual-df t critical value on the logit
#' scale. Inference is model-based, without a spatial or cluster correction.
#' The interaction P value is unadjusted; Holm adjustment occurs in 05b.
fit_mosaic <- function(gain_group, loss_group, gain_label, loss_label, data) {
        
        d <- data |>
                transmute(
                        CellCode,
                        coverage = total_pa_land_cov,
                        land_weight,
                        gain = as.integer(
                                .data[[paste0(gain_group, "_gained")]] > 0
                        ),
                        loss = factor(
                                as.integer(
                                        .data[[paste0(loss_group, "_lost")]] > 0
                                ),
                                levels = c(0, 1),
                                labels = c("No loss", "Loss")
                        ),
                        initial_gain = .data[[paste0(gain_group, "_initial")]],
                        initial_loss = .data[[paste0(loss_group, "_initial")]]
                ) |>
                # Only required model variables are checked here; beta_jaccard
                # is not included, so its NA values cannot remove double absences.
                drop_na() |>
                filter(
                        is.finite(coverage),
                        is.finite(land_weight),
                        land_weight > 0,
                        initial_loss > 0
                )
        
        # Loss must be possible in the predictor group.
        # Both gain = 0 and gain = 1 are retained when initial_gain = 0.
        
        stopifnot(
                n_distinct(d$loss) == 2L,
                n_distinct(d$gain) == 2L
        )
        
        model <- glm(
                gain ~ loss * coverage + initial_gain + initial_loss,
                family = quasibinomial(),
                weights = land_weight,
                data = d
        )
        
        # Compare nested models fitted to exactly the same observations.
        # These F tests and prediction intervals use model-based dispersion;
        # they do not correct for spatial dependence between cells.
        additive <- update(model, . ~ . - loss:coverage)
        
        interaction_test <- anova(additive, model, test = "F")
        
        # Use the coverage range shared by both loss categories.
        ranges <- d |>
                group_by(loss) |>
                summarise(
                        low = min(coverage),
                        high = max(coverage),
                        .groups = "drop"
                )
        
        low <- max(ranges$low)
        high <- min(ranges$high)
        stopifnot(low < high)
        
        # Hold initial richness at its median within this model's sample.
        # Curves are conditional predictions, not averages over all cells.
        nd <- expand_grid(
                coverage = seq(low, high, length.out = 150),
                loss = factor(
                        c("No loss", "Loss"),
                        levels = c("No loss", "Loss")
                )
        ) |>
                mutate(
                        initial_gain = median(d$initial_gain),
                        initial_loss = median(d$initial_loss)
                )
        
        pr <- predict(model, newdata = nd, type = "link", se.fit = TRUE)
        crit <- qt(0.975, df.residual(model))
        
        nd <- nd |>
                mutate(
                        estimate = plogis(pr$fit),
                        lower = plogis(pr$fit - crit * pr$se.fit),
                        upper = plogis(pr$fit + crit * pr$se.fit)
                )
        
        panel <- ggplot(
                nd,
                aes(coverage, estimate, colour = loss, fill = loss)
        ) +
                geom_ribbon(
                        aes(ymin = lower, ymax = upper),
                        alpha = 0.15, colour = NA
                ) +
                geom_line(linewidth = 0.9) +
                scale_colour_manual(
                        values = c("No loss" = "#547A94", "Loss" = "#C36B35"),
                        drop = FALSE
                ) +
                scale_fill_manual(
                        values = c("No loss" = "#547A94", "Loss" = "#C36B35"),
                        drop = FALSE
                ) +
                scale_x_continuous(
                        limits = c(0, 100),
                        breaks = c(0, 25, 50, 75, 100)
                ) +
                scale_y_continuous(
                        limits = c(0, 1),
                        breaks = c(0, 0.25, 0.5, 0.75, 1),
                        labels = scales::label_percent()
                ) +
                labs(
                        title = paste0(gain_label, " gain"),
                        subtitle = paste0("Predictor: ", tolower(loss_label), " loss"),
                        x = "Protected terrestrial area (%)",
                        y = "Probability of habitat gain",
                        colour = "Predictor-group loss",
                        fill = "Predictor-group loss"
                ) +
                theme_classic(base_size = 11) +
                theme(
                        plot.title = element_text(face = "bold"),
                        legend.position = "bottom"
                )
        
        list(
                model = model,
                data = d,
                predictions = nd,
                plot = panel,
                test = tibble(
                        response = paste(gain_label, "gain"),
                        predictor = paste(loss_label, "loss"),
                        n_cells = nrow(d),
                        interaction_F = interaction_test$F[2],
                        interaction_p = interaction_test$`Pr(>F)`[2]
                )
        )
}

#' Contrast habitat-gain probabilities between loss categories
#'
#' Estimates Loss minus No loss probability differences at selected
#' protected-area coverage levels.
#'
#' @param result A result list returned by fit_mosaic.
#' @param coverage_values Numeric coverage percentages; defaults to 0, 30, 100.
#'
#' @return A tibble containing response, predictor, protection_pct,
#'   initial_gain, initial_loss, gain_no_loss_pct, gain_loss_pct,
#'   difference_pp, lower_pp and upper_pp.
#'
#' @details
#' Coverage values outside the range shared by both loss categories are
#' omitted; the function stops if no requested values remain. Richness is
#' fixed at model-sample medians, matching the plotted predictions.
#' Delta-method SEs use the model covariance, including covariance between
#' the two predictions. Pointwise 95% intervals use a residual-df t critical
#' value and are not adjusted for multiple comparisons or spatial dependence.
#' Probabilities are percentages; contrasts and bounds are percentage points.
contrast_mosaic <- function(result, coverage_values = c(0, 30, 100)) {
        
        model <- result$model
        d <- result$data
        
        # Keep predictions within the coverage range shared by both categories.
        ranges <- d |>
                group_by(loss) |>
                summarise(
                        low = min(coverage),
                        high = max(coverage),
                        .groups = "drop"
                )
        
        coverage_values <- coverage_values[
                coverage_values >= max(ranges$low) &
                        coverage_values <= min(ranges$high)
        ]
        
        stopifnot(length(coverage_values) > 0)
        
        # Use the same initial-richness values as in the figure.
        nd_no <- data.frame(
                coverage = coverage_values,
                loss = factor("No loss", levels = levels(d$loss)),
                initial_gain = median(d$initial_gain),
                initial_loss = median(d$initial_loss)
        )
        
        nd_loss <- nd_no
        nd_loss$loss <- factor("Loss", levels = levels(d$loss))
        
        model_terms <- delete.response(terms(model))
        
        X_no <- model.matrix(
                model_terms, nd_no, contrasts.arg = model$contrasts
        )
        
        X_loss <- model.matrix(
                model_terms, nd_loss, contrasts.arg = model$contrasts
        )
        
        p_no <- plogis(drop(X_no %*% coef(model)))
        p_loss <- plogis(drop(X_loss %*% coef(model)))
        
        # Delta-method gradient for the difference between probabilities.
        # The covariance matrix includes quasibinomial dispersion.
        gradient <- sweep(X_loss, 1, p_loss * (1 - p_loss), "*") -
                sweep(X_no, 1, p_no * (1 - p_no), "*")
        
        V <- vcov(model)
        
        se_difference <- sqrt(pmax(
                0,
                rowSums((gradient %*% V) * gradient)
        ))
        
        difference <- p_loss - p_no
        crit <- qt(0.975, df.residual(model))
        
        # Probabilities are percentages; contrasts are percentage points.
        tibble(
                response = result$test$response,
                predictor = result$test$predictor,
                protection_pct = coverage_values,
                initial_gain = nd_no$initial_gain,
                initial_loss = nd_no$initial_loss,
                gain_no_loss_pct = 100 * p_no,
                gain_loss_pct = 100 * p_loss,
                difference_pp = 100 * difference,
                lower_pp = 100 * (difference - crit * se_difference),
                upper_pp = 100 * (difference + crit * se_difference)
        )
}

#' Model proportional habitat loss in relation to other-group gain
#'
#' Fits a land-weighted quasibinomial logit model for lost versus shared
#' habitats, tests its gain-by-protection interaction and prepares a plot.
#'
#' @param loss_group Prefix identifying the response group in wide data.
#' @param gain_group Prefix identifying the predictor group in wide data.
#' @param loss_label,gain_label Display labels for the two groups.
#' @param data Wide data containing CellCode, coverage, land_weight,
#'   <loss_group>_lost, <loss_group>_shared, <gain_group>_gained and
#'   <gain_group>_initial.
#'
#' @return A list containing model, filtered data, predictions, plot and a
#'   one-row test tibble with response, predictor, n_cells, interaction_F,
#'   interaction_p and dispersion.
#'
#' @details
#' Removes incomplete records and requires finite coverage, positive finite
#' weights and a positive lost-plus-shared denominator. Initial absence of
#' the gain group is allowed. Both gain categories must occur; loss and
#' shared counts must be nonnegative and have positive totals.
#'
#' The response is cbind(lost, shared); predictors are gain * coverage and
#' initial predictor-group richness. An F test compares the interaction
#' model with its additive counterpart. The response denominator supplies
#' trial counts rather than a separate richness effect.
#'
#' Predictions span 150 values within the coverage range shared by both gain
#' categories, at median initial predictor richness. Pointwise 95% bounds
#' use model-based SEs and a residual-df t critical value on the logit scale.
#' Interaction P values are unadjusted; Holm adjustment occurs in 05b.
#' Inference does not correct for spatial dependence between cells.
fit_mosaic_loss <- function(
                loss_group, gain_group, loss_label, gain_label, data
) {
        
        d <- data |>
                transmute(
                        CellCode,
                        coverage,
                        land_weight,
                        lost = .data[[paste0(loss_group, "_lost")]],
                        shared = .data[[paste0(loss_group, "_shared")]],
                        gain = factor(
                                as.integer(
                                        .data[[paste0(gain_group, "_gained")]] > 0
                                ),
                                levels = c(0, 1),
                                labels = c("No gain", "Gain")
                        ),
                        initial_predictor =
                                .data[[paste0(gain_group, "_initial")]]
                ) |>
                drop_na() |>
                filter(
                        is.finite(coverage),
                        is.finite(land_weight),
                        land_weight > 0,
                        lost + shared > 0
                )
        
        # The response group must be initially present for loss to be defined.
        # The predictor group may be initially absent, with or without later gain.
        stopifnot(
                n_distinct(d$gain) == 2L,
                all(d$lost >= 0),
                all(d$shared >= 0),
                sum(d$lost) > 0,
                sum(d$shared) > 0
        )
        
        # Model lost / (lost + shared), with an additional terrestrial-area weight.
        # The denominator supplies initial response-group richness as trial count;
        # it does not estimate a separate effect of richness on loss probability.
        model <- glm(
                cbind(lost, shared) ~
                        gain * coverage + initial_predictor,
                family = quasibinomial(),
                weights = land_weight,
                data = d
        )
        
        additive <- update(model, . ~ . - gain:coverage)
        interaction_test <- anova(additive, model, test = "F")
        
        # Predict only within the coverage range shared by both gain categories.
        ranges <- d |>
                group_by(gain) |>
                summarise(
                        low = min(coverage),
                        high = max(coverage),
                        .groups = "drop"
                )
        
        low <- max(ranges$low)
        high <- min(ranges$high)
        stopifnot(low < high)
        
        nd <- expand_grid(
                coverage = seq(low, high, length.out = 150),
                gain = factor(
                        c("No gain", "Gain"),
                        levels = c("No gain", "Gain")
                )
        ) |>
                mutate(
                        initial_predictor = median(d$initial_predictor)
                )
        
        # Pointwise 95% confidence intervals include quasibinomial dispersion.
        # Inference does not account for spatial dependence between cells.
        pr <- predict(model, newdata = nd, type = "link", se.fit = TRUE)
        crit <- qt(0.975, df.residual(model))
        
        nd <- nd |>
                mutate(
                        estimate = plogis(pr$fit),
                        lower = plogis(pr$fit - crit * pr$se.fit),
                        upper = plogis(pr$fit + crit * pr$se.fit)
                )
        
        panel <- ggplot(
                nd,
                aes(coverage, estimate, colour = gain, fill = gain)
        ) +
                geom_ribbon(
                        aes(ymin = lower, ymax = upper),
                        alpha = 0.15, colour = NA
                ) +
                geom_line(linewidth = 0.9) +
                scale_colour_manual(
                        values = c("No gain" = "#547A94", "Gain" = "#C36B35")
                ) +
                scale_fill_manual(
                        values = c("No gain" = "#547A94", "Gain" = "#C36B35")
                ) +
                scale_x_continuous(
                        limits = c(0, 100),
                        breaks = c(0, 25, 50, 75, 100)
                ) +
                scale_y_continuous(
                        limits = c(0, 1),
                        breaks = c(0, 0.25, 0.5, 0.75, 1),
                        labels = scales::label_percent()
                ) +
                labs(
                        title = paste(loss_label, "loss"),
                        subtitle = paste0(
                                "Predictor: ", tolower(gain_label), " gain"
                        ),
                        x = "Protected terrestrial area (%)",
                        y = "Proportional habitat loss",
                        colour = "Predictor-group gain",
                        fill = "Predictor-group gain"
                ) +
                theme_classic(base_size = 11) +
                theme(
                        plot.title = element_text(face = "bold"),
                        legend.position = "bottom"
                )
        
        list(
                model = model,
                data = d,
                predictions = nd,
                plot = panel,
                test = tibble(
                        response = paste(loss_label, "loss"),
                        predictor = paste(gain_label, "gain"),
                        n_cells = nrow(d),
                        interaction_F = interaction_test$F[2],
                        interaction_p = interaction_test$`Pr(>F)`[2],
                        dispersion = summary(model)$dispersion
                )
        )
}

#' Contrast proportional habitat loss between gain categories
#'
#' Estimates Gain minus No gain differences in predicted loss probability
#' at selected protected-area coverage levels.
#'
#' @param result A result list returned by fit_mosaic_loss.
#' @param coverage_values Numeric coverage percentages; defaults to 0, 30, 100.
#'
#' @return A tibble containing response, predictor, protection_pct,
#'   initial_predictor, loss_no_gain_pct, loss_gain_pct, difference_pp,
#'   lower_pp and upper_pp.
#'
#' @details
#' Retains only coverage values within the range shared by both gain
#' categories and stops if none remain. Initial predictor richness is held
#' at its model-sample median. Delta-method SEs include covariance between
#' predictions and quasibinomial dispersion through the model covariance.
#' Pointwise 95% intervals use a residual-df t critical value without
#' multiplicity or spatial corrections. Probabilities are percentages;
#' differences and confidence bounds are percentage points.
contrast_mosaic_loss <- function(result, coverage_values = c(0, 30, 100)) {
        
        model <- result$model
        d <- result$data
        
        # Restrict contrasts to the coverage range shared by both categories.
        ranges <- d |>
                group_by(gain) |>
                summarise(
                        low = min(coverage),
                        high = max(coverage),
                        .groups = "drop"
                )
        
        coverage_values <- coverage_values[
                coverage_values >= max(ranges$low) &
                        coverage_values <= min(ranges$high)
        ]
        
        stopifnot(length(coverage_values) > 0)
        
        # Hold predictor-group initial richness at the value used in the figure.
        nd_no <- data.frame(
                coverage = coverage_values,
                gain = factor("No gain", levels = levels(d$gain)),
                initial_predictor = median(d$initial_predictor)
        )
        
        nd_gain <- nd_no
        nd_gain$gain <- factor("Gain", levels = levels(d$gain))
        
        model_terms <- delete.response(terms(model))
        
        X_no <- model.matrix(
                model_terms, nd_no, contrasts.arg = model$contrasts
        )
        
        X_gain <- model.matrix(
                model_terms, nd_gain, contrasts.arg = model$contrasts
        )
        
        p_no <- plogis(drop(X_no %*% coef(model)))
        p_gain <- plogis(drop(X_gain %*% coef(model)))
        
        # Delta-method uncertainty accounts for covariance between predictions.
        gradient <- sweep(X_gain, 1, p_gain * (1 - p_gain), "*") -
                sweep(X_no, 1, p_no * (1 - p_no), "*")
        
        V <- vcov(model)
        
        se_difference <- sqrt(pmax(
                0,
                rowSums((gradient %*% V) * gradient)
        ))
        
        difference <- p_gain - p_no
        crit <- qt(0.975, df.residual(model))
        
        tibble(
                response = result$test$response,
                predictor = result$test$predictor,
                protection_pct = coverage_values,
                initial_predictor = nd_no$initial_predictor,
                loss_no_gain_pct = 100 * p_no,
                loss_gain_pct = 100 * p_gain,
                difference_pp = 100 * difference,
                lower_pp = 100 * (difference - crit * se_difference),
                upper_pp = 100 * (difference + crit * se_difference)
        )
}
