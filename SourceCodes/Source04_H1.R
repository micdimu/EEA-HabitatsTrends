# Functions for H1: Figure 2 and model-summary tables.
# Attach dplyr, ggplot2 and scales in the main script.

#' Predict a logit-link model on the response scale
#'
#' Returns fitted probabilities and approximate 95% pointwise confidence
#' intervals for a supplied prediction grid.
#'
#' @param model A fitted logit-link GLM supporting link predictions with SEs.
#' @param model_name Label identifying the model in the returned data.
#' @param new_data A data frame containing all predictors required by model.
#'
#' @return The prediction grid with columns `model`, `fit`, `lower`, `upper`.
#'
#' @details
#' Bounds are calculated as link predictions plus or minus 1.96 standard
#' errors, then transformed with plogis. These are confidence intervals for
#' the fitted mean, not prediction intervals. The logit link is assumed.
#' Requires dplyr to be attached.
get_model_predictions <- function(
                model,
                model_name,
                new_data
) {
        
        prediction <- predict(
                model,
                newdata = new_data,
                type = "link",
                se.fit = TRUE
        )
        
        new_data |>
                mutate(
                        model = model_name,
                        fit = plogis(prediction$fit),
                        lower = plogis(
                                prediction$fit -
                                        1.96 * prediction$se.fit
                        ),
                        upper = plogis(
                                prediction$fit +
                                        1.96 * prediction$se.fit
                        )
                )
}

#' Create explained-deviance labels for the two H1 models
#'
#' Calculates in-sample explained deviance and positions line keys and text
#' labels for the linear and logarithmic models.
#'
#' @param null_model The corresponding intercept-only fitted model.
#' @param linear_model Fitted model with linear protected-area coverage.
#' @param log_model Fitted model with log-transformed protected-area coverage.
#' @param y_range Numeric vector containing the lower and upper panel limits.
#' @param x_start Horizontal starting position for line keys; defaults to 2.
#'
#' @return A two-row data frame containing `model`, `explained_deviance`,
#'   `asterisk`, `label`, `x_start`, `x_end`, `x_text` and `y`.
#'
#' @details
#' Explained deviance is 100 times one minus model deviance divided by null
#' deviance. Models must use comparable data and weights, with nonzero null
#' deviance. An asterisk marks the largest value; ties select the first model.
#' The asterisk does not indicate statistical significance.
#' Requires dplyr to be attached.
make_fit_labels <- function(
                null_model,
                linear_model,
                log_model,
                y_range,
                x_start = 2
) {
        
        null_deviance <- deviance(null_model)
        
        fit_labels <- data.frame(
                model = c("Linear", "Logarithmic"),
                explained_deviance = 100 * c(
                        1 - deviance(linear_model) / null_deviance,
                        1 - deviance(log_model) / null_deviance
                )
        )
        
        # Add an asterisk only to the model with the highest explained deviance
        fit_labels$asterisk <- ""
        fit_labels$asterisk[which.max(
                fit_labels$explained_deviance
        )] <- "*"
        
        # Generate labels and positions within the response scale
        fit_labels |>
                mutate(
                        label = sprintf(
                                "%s (D² = %.3f%%)%s",
                                model,
                                explained_deviance,
                                asterisk
                        ),
                        x_start = x_start,
                        x_end = x_start + 5,
                        x_text = x_start + 6.5,
                        y = min(y_range) +
                                c(0.12, 0.045) * diff(y_range)
                )
}

#' Draw an H1 response panel with observations and fitted models
#'
#' Combines binned observations, fitted curves, confidence ribbons and
#' explained-deviance labels in a common graphical layout.
#'
#' @param binned_data Data frame containing `protection` and `observed`.
#' @param predictions Data frame containing `total_pa_land_cov`, `model`,
#'   `fit`, `lower` and `upper`, as returned by get_model_predictions.
#' @param fit_labels Data frame returned by make_fit_labels.
#' @param y_label Vertical axis label, as text or an expression.
#' @param panel_tag Panel identifier, such as "a" or "b".
#' @param model_colours Named colour vector matching the model labels.
#' @param model_lines Named line-type vector matching the model labels.
#'
#' @return A ggplot object for one biological response.
#'
#' @details
#' The protection axis spans 0 to 100 and the response axis is formatted as
#' percentages. Binned observations are descriptive; the function fits no
#' models. Requires ggplot2 and scales to be attached. Colours and line types
#' are supplied explicitly rather than read from global objects.
make_h1_panel <- function(
                binned_data,
                predictions,
                fit_labels,
                y_label,
                panel_tag,
                model_colours,
                model_lines
) {
        
        ggplot() +
                geom_point(
                        data = binned_data,
                        aes(
                                x = protection,
                                y = observed
                        ),
                        shape = 21,
                        size = 1.5,
                        stroke = 0.6,
                        colour = "black",
                        fill = "white"
                ) +
                geom_ribbon(
                        data = predictions,
                        aes(
                                x = total_pa_land_cov,
                                ymin = lower,
                                ymax = upper,
                                fill = model
                        ),
                        alpha = 0.12,
                        colour = NA
                ) +
                geom_line(
                        data = predictions,
                        aes(
                                x = total_pa_land_cov,
                                y = fit,
                                colour = model,
                                linetype = model
                        ),
                        linewidth = 1.15
                ) +
                geom_segment(
                        data = fit_labels,
                        aes(
                                x = x_start,
                                xend = x_end,
                                y = y,
                                yend = y,
                                colour = model,
                                linetype = model
                        ),
                        linewidth = 1.1,
                        inherit.aes = FALSE,
                        show.legend = FALSE
                ) +
                geom_text(
                        data = fit_labels,
                        aes(
                                x = x_text,
                                y = y,
                                label = label
                        ),
                        colour = "black",
                        hjust = 0,
                        vjust = 0.5,
                        size = 3.1,
                        inherit.aes = FALSE
                ) +
                scale_colour_manual(
                        values = model_colours
                ) +
                scale_fill_manual(
                        values = model_colours
                ) +
                scale_linetype_manual(
                        values = model_lines
                ) +
                scale_x_continuous(
                        limits = c(0, 100),
                        breaks = seq(0, 100, by = 20),
                        expand = expansion(
                                mult = c(0.005, 0.01)
                        )
                ) +
                scale_y_continuous(
                        labels = label_percent(
                                accuracy = 0.1
                        ),
                        expand = expansion(
                                mult = c(0.05, 0.08)
                        )
                ) +
                labs(
                        tag = panel_tag,
                        x = NULL,
                        y = y_label
                ) +
                guides(
                        colour = "none",
                        fill = "none",
                        linetype = "none"
                ) +
                theme_classic(
                        base_size = 8
                ) +
                theme(
                        axis.title.y = element_text(
                                size = 11
                        ),
                        axis.text = element_text(
                                colour = "black"
                        ),
                        plot.tag = element_text(
                                face = "bold",
                                size = 10
                        ),
                        plot.tag.position = c(
                                0.01,
                                0.99
                        ),
                        plot.margin = margin(
                                3,
                                8,
                                3,
                                5.5
                        )
                )
}

#' Summarise the linear and logarithmic models for an H1 response
#'
#' Extracts coefficients, comparisons with the null model, explained
#' deviance and predicted changes across the protection gradient.
#'
#' @param hypothesis Hypothesis identifier, such as "H1a" or "H1b".
#' @param response Descriptive label for the response variable.
#' @param null_model Corresponding intercept-only quasibinomial logit model.
#' @param linear_model Model with a linear total_pa_land_cov effect.
#' @param log_model Model with a log1p(total_pa_land_cov) effect.
#'
#' @return A two-row data frame with hypothesis, response and model labels;
#'   `estimate`, `standard_error`, `F_value`, `p_value`, `residual_deviance`,
#'   `explained_deviance`, `predicted_at_0`, `predicted_at_100`,
#'   `change_percentage_points`, `relative_change` and `best_in_sample`.
#'
#' @details
#' Each alternative is compared separately with the null model using an F
#' test. All models must use the same observations and weights. The second
#' coefficient is taken as the protection effect on the logit scale.
#' Predictions use protection values of 0 and 100 and the inverse logit.
#' Explained deviance and relative change are percentages; endpoint change
#' is in percentage points. All ties for greatest explained deviance are
#' flagged by best_in_sample. Zero baseline predictions or null deviance
#' are not handled specially. Requires dplyr to be attached.
summarise_h1_models <- function(
                hypothesis,
                response,
                null_model,
                linear_model,
                log_model
) {
        
        models <- list(
                linear_model,
                log_model
        )
        
        model_names <- c(
                "Linear",
                "Logarithmic"
        )
        
        tests <- lapply(
                models,
                function(model) {
                        anova(
                                null_model,
                                model,
                                test = "F"
                        )
                }
        )
        
        endpoint_predictions <- lapply(
                models,
                function(model) {
                        plogis(
                                predict(
                                        model,
                                        newdata = data.frame(
                                                total_pa_land_cov = c(
                                                        0,
                                                        100
                                                )
                                        ),
                                        type = "link"
                                )
                        )
                }
        )
        
        data.frame(
                hypothesis = hypothesis,
                response = response,
                model = model_names,
                estimate = vapply(
                        models,
                        function(model) {
                                coef(model)[2]
                        },
                        numeric(1)
                ),
                standard_error = vapply(
                        models,
                        function(model) {
                                coef(summary(model))[2, "Std. Error"]
                        },
                        numeric(1)
                ),
                F_value = vapply(
                        tests,
                        function(test) {
                                test$F[2]
                        },
                        numeric(1)
                ),
                p_value = vapply(
                        tests,
                        function(test) {
                                test$`Pr(>F)`[2]
                        },
                        numeric(1)
                ),
                residual_deviance = vapply(
                        models,
                        deviance,
                        numeric(1)
                ),
                explained_deviance = 100 * vapply(
                        models,
                        function(model) {
                                1 -
                                        deviance(model) /
                                        deviance(null_model)
                        },
                        numeric(1)
                ),
                predicted_at_0 = vapply(
                        endpoint_predictions,
                        function(prediction) {
                                prediction[1]
                        },
                        numeric(1)
                ),
                predicted_at_100 = vapply(
                        endpoint_predictions,
                        function(prediction) {
                                prediction[2]
                        },
                        numeric(1)
                )
        ) |>
                mutate(
                        change_percentage_points =
                                100 * (
                                        predicted_at_100 -
                                                predicted_at_0
                                ),
                        relative_change = 100 *
                                change_percentage_points /
                                (predicted_at_0 * 100),
                        best_in_sample =
                                explained_deviance ==
                                max(explained_deviance)
                )
}
