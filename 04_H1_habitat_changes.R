library(tidyverse)
library(sf)
library(terra)
library(patchwork)
library(ggeffects)
library(mgcv)
library(scales)
library(patchwork)
library(gt)

source("Source.R")

#### Load data ####

temporal_beta <- read.csv("processed/temporal_beta_cell_by_cell.csv") |> 
        mutate(net_change = richness_2013_2018 - richness_2007_2012)

grid_pa <- st_read("data/PAs/derived/europe_10km_protected_area_coverage.gpkg")


grid_T <- grid_pa |>
        right_join(temporal_beta, by = c("CellCode" = "cell_id")) 

#### H1 | Protected-area coverage and habitat change ####

# General hypothesis ------------------------------------------------------
# Greater protected-area coverage should reduce habitat change through time.
# We test two complementary components:
#
# H1a | Habitat retention
# The probability that a habitat present in the first monitoring cycle is lost
# should decrease as protected-area coverage increases.
#
# H1b | Compositional stability
# Temporal Jaccard dissimilarity should decrease as protected-area coverage
# increases, indicating greater stability in habitat composition.
#
# For both responses, we compare:
#   1. a null model with no protection effect;
#   2. a linear effect, representing a constant change along the gradient;
#   3. a logarithmic effect, representing stronger benefits at low protection
#      levels followed by diminishing returns.
#
# Quasibinomial models account for overdispersion. Grid cells are weighted by
# their terrestrial fraction. H1a additionally accounts for initial habitat
# richness through its two-column binomial response; H1b weights each Jaccard
# value by the number of habitats contributing to the union.


#### H1a | Habitat retention ####

# Response: losses among habitats present in the first monitoring cycle.
# In cbind(lost_species, shared_species), lost habitats are the events and
# shared habitats are the retained cases. The model therefore estimates the
# probability that an initially present habitat is lost.

h1a_data <- grid_T |>
        st_drop_geometry() |>
        mutate(
                initial_richness = lost_species + shared_species,
                loss_proportion = lost_species / initial_richness,
                land_weight = land_cov / 100
        ) |>
        filter(
                initial_richness > 0,
                land_weight > 0,
                !is.na(total_pa_land_cov),
                !is.na(lost_species),
                !is.na(shared_species)
        )

# Null model: habitat-loss probability is unrelated to protection
h1a_null_model <- glm(
        cbind(lost_species, shared_species) ~ 1,
        data = h1a_data,
        weights = land_weight,
        family = quasibinomial(link = "logit")
)

# Linear alternative: protection has a constant effect across its gradient
h1a_linear_model <- glm(
        cbind(lost_species, shared_species) ~ total_pa_land_cov,
        data = h1a_data,
        weights = land_weight,
        family = quasibinomial(link = "logit")
)

# Logarithmic alternative: the protection effect is strongest at low coverage
# and progressively weakens as coverage increases
h1a_log_model <- glm(
        cbind(lost_species, shared_species) ~ log1p(total_pa_land_cov),
        data = h1a_data,
        weights = land_weight,
        family = quasibinomial(link = "logit")
)

# Test each protection model separately against the null model.
# The linear and logarithmic alternatives are not nested in one another.

h1a_linear_test <- anova(
        h1a_null_model,
        h1a_linear_model,
        test = "F"
)

h1a_log_test <- anova(
        h1a_null_model,
        h1a_log_model,
        test = "F"
)

h1a_linear_test
h1a_log_test

summary(h1a_linear_model)
summary(h1a_log_model)


#### H1b | Compositional stability ####

# Response: temporal Jaccard dissimilarity, available as beta_jaccard.
# It represents the proportion of the habitat union involved in losses or
# gains. Lower values therefore indicate more stable habitat composition.
#
# beta_weight combines the terrestrial fraction of each grid cell with the
# habitat union on which beta_jaccard is based. A beta value calculated from
# few habitats therefore receives less weight than one based on a richer union.

h1b_data <- grid_T |>
        st_drop_geometry() |>
        mutate(
                union_richness =
                        lost_species +
                        gained_species +
                        shared_species,
                land_weight = land_cov / 100,
                beta_weight = land_weight * union_richness
        ) |>
        filter(
                union_richness > 0,
                beta_weight > 0,
                !is.na(beta_jaccard),
                !is.na(total_pa_land_cov)
        )

# Null model: temporal beta diversity is unrelated to protection
h1b_null_model <- glm(
        beta_jaccard ~ 1,
        data = h1b_data,
        weights = beta_weight,
        family = quasibinomial(link = "logit")
)

# Linear alternative: protection has a constant effect across its gradient
h1b_linear_model <- glm(
        beta_jaccard ~ total_pa_land_cov,
        data = h1b_data,
        weights = beta_weight,
        family = quasibinomial(link = "logit")
)

# Logarithmic alternative: the protection effect is strongest at low coverage
# and progressively weakens as coverage increases
h1b_log_model <- glm(
        beta_jaccard ~ log1p(total_pa_land_cov),
        data = h1b_data,
        weights = beta_weight,
        family = quasibinomial(link = "logit")
)

# Test each protection model separately against the null model

h1b_linear_test <- anova(
        h1b_null_model,
        h1b_linear_model,
        test = "F"
)

h1b_log_test <- anova(
        h1b_null_model,
        h1b_log_model,
        test = "F"
)

h1b_linear_test
h1b_log_test

summary(h1b_linear_model)
summary(h1b_log_model)


#### Figures and outputs ####

# Create output directories when they do not already exist

dir.create(
        "figures",
        recursive = TRUE,
        showWarnings = FALSE
)

dir.create(
        "tables",
        recursive = TRUE,
        showWarnings = FALSE
)


##### Figure 2 | Protected-area coverage and habitat change ####

# Calculate land-weighted observed habitat-loss proportions in
# 5-percentage-point protection bins.
#
# The points are descriptive summaries. Statistical inference comes from
# the quasibinomial models fitted above.

h1a_binned <- h1a_data |>
        mutate(
                protection_bin = cut(
                        total_pa_land_cov,
                        breaks = seq(0, 100, by = 5),
                        include.lowest = TRUE
                )
        ) |>
        group_by(protection_bin) |>
        summarise(
                protection = weighted.mean(
                        total_pa_land_cov,
                        w = land_weight * initial_richness
                ),
                observed = sum(
                        land_weight * lost_species
                ) / sum(
                        land_weight * initial_richness
                ),
                n_cells = n(),
                .groups = "drop"
        )

# Calculate weighted observed temporal beta diversity in
# 5-percentage-point protection bins

h1b_binned <- h1b_data |>
        mutate(
                protection_bin = cut(
                        total_pa_land_cov,
                        breaks = seq(0, 100, by = 5),
                        include.lowest = TRUE
                )
        ) |>
        group_by(protection_bin) |>
        summarise(
                protection = weighted.mean(
                        total_pa_land_cov,
                        w = beta_weight
                ),
                observed = weighted.mean(
                        beta_jaccard,
                        w = beta_weight
                ),
                n_cells = n(),
                .groups = "drop"
        )

# Create a common prediction grid for both response variables

prediction_grid <- data.frame(
        total_pa_land_cov = seq(
                0,
                100,
                length.out = 500
        )
)

# Extract predictions and 95% confidence intervals from a fitted model

get_model_predictions <- function(
                model,
                model_name,
                new_data = prediction_grid
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

# Habitat-loss predictions

h1a_predictions <- bind_rows(
        get_model_predictions(
                h1a_linear_model,
                "Linear"
        ),
        get_model_predictions(
                h1a_log_model,
                "Logarithmic"
        )
)

# Temporal beta-diversity predictions

h1b_predictions <- bind_rows(
        get_model_predictions(
                h1b_linear_model,
                "Linear"
        ),
        get_model_predictions(
                h1b_log_model,
                "Logarithmic"
        )
)

# Common graphical settings

model_colours <- c(
        "Linear" = "#0072B2",
        "Logarithmic" = "#D55E00"
)

model_lines <- c(
        "Linear" = "solid",
        "Logarithmic" = "22"
)

# make_fit_labels() is already defined in Source.R.
#
# It calculates explained deviance directly from the fitted models,
# adds an asterisk to the model with the largest explained deviance,
# and returns the coordinates needed for the line keys and text labels.

h1a_y_range <- range(
        c(
                h1a_binned$observed,
                h1a_predictions$lower,
                h1a_predictions$upper
        ),
        na.rm = TRUE
)

h1b_y_range <- range(
        c(
                h1b_binned$observed,
                h1b_predictions$lower,
                h1b_predictions$upper
        ),
        na.rm = TRUE
)

h1a_fit_labels <- make_fit_labels(
        null_model = h1a_null_model,
        linear_model = h1a_linear_model,
        log_model = h1a_log_model,
        y_range = h1a_y_range
)

h1b_fit_labels <- make_fit_labels(
        null_model = h1b_null_model,
        linear_model = h1b_linear_model,
        log_model = h1b_log_model,
        y_range = h1b_y_range
)

# make_h1_panel is already defined in Source.R.
# Panel a: probability of habitat loss

panel_h1a <- make_h1_panel(
        binned_data = h1a_binned,
        predictions = h1a_predictions,
        fit_labels = h1a_fit_labels,
        y_label = "Probability of habitat loss",
        panel_tag = "a"
)

# Panel b: temporal beta diversity

panel_h1b <- make_h1_panel(
        binned_data = h1b_binned,
        predictions = h1b_predictions,
        fit_labels = h1b_fit_labels,
        y_label = expression(
                "Temporal " * beta * " diversity"
        ),
        panel_tag = "b"
)

# Bottom panel: terrestrial-area-weighted distribution of grid cells.
# The complete x axis is retained in all three panels, while the axis title
# appears only below the bottom panel.

# Bottom panel: distribution of the number of grid cells
# along the protected-area coverage gradient

panel_cells <- ggplot(
        h1a_data,
        aes(x = total_pa_land_cov)
) +
        geom_histogram(
                binwidth = 2.5,
                boundary = 0,
                fill = "grey65",
                colour = "white",
                linewidth = 0.15
        ) +
        scale_x_continuous(
                limits = c(0, 100),
                breaks = seq(0, 100, by = 20),
                expand = expansion(
                        mult = c(0.005, 0.01)
                )
        ) +
        scale_y_continuous(
                labels = label_number(
                        big.mark = ","
                ),
                expand = expansion(
                        mult = c(0, 0.05)
                )
        ) +
        labs(
                x = "Protected terrestrial area (%)",
                y = "Number of grid cells"
        ) +
        theme_classic(
                base_size = 8
        ) +
        theme(
                axis.title.x = element_text(
                        size = 8
                ),
                axis.title.y = element_text(
                        size = 8
                ),
                axis.text = element_text(
                        colour = "black"
                ),
                plot.margin = margin(
                        0,
                        8,
                        5.5,
                        5.5
                )
        )
# Assemble Figure 2

figure_2 <- panel_h1a /
        panel_h1b /
        panel_cells +
        plot_layout(
                heights = c(4, 4, 2)
        )

figure_2

# Export Figure 2

ggsave(
        filename = file.path(
                "figures",
                "Figure_2.tiff"
        ),
        plot = figure_2,
        width = 180,
        height = 190,
        units = "mm",
        dpi = 600,
        device = "tiff",
        compression = "lzw",
        bg = "white"
)


##### Table | Summary of H1 model results ####

# Create one consistent summary for the four candidate models.
#
# Estimates are expressed on the logit scale.
# F tests and P values compare each candidate model with its null model.
# Explained deviance is the percentage of null deviance explained in-sample.
# Predicted changes compare 0% and 100% protected terrestrial area.

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

# Combine H1a and H1b results

h1_model_results <- bind_rows(
        summarise_h1_models(
                hypothesis = "H1a",
                response = "Habitat-loss probability",
                null_model = h1a_null_model,
                linear_model = h1a_linear_model,
                log_model = h1a_log_model
        ),
        summarise_h1_models(
                hypothesis = "H1b",
                response = "Temporal Jaccard dissimilarity",
                null_model = h1b_null_model,
                linear_model = h1b_linear_model,
                log_model = h1b_log_model
        )
)

h1_model_results

# Create a formatted, paper-ready table.
# The asterisk identifies the model with the greatest in-sample explained
# deviance within each response; it is not a statistical-significance symbol.

h1_model_table <- h1_model_results |>
        transmute(
                Hypothesis = hypothesis,
                Response = response,
                Model = paste0(
                        model,
                        if_else(
                                best_in_sample,
                                "*",
                                ""
                        )
                ),
                `Estimate (SE)` = sprintf(
                        "%.4f (%.4f)",
                        estimate,
                        standard_error
                ),
                F = sprintf(
                        "%.2f",
                        F_value
                ),
                P = format.pval(
                        p_value,
                        digits = 3,
                        eps = 0.001
                ),
                `Residual deviance` = sprintf(
                        "%.1f",
                        residual_deviance
                ),
                `Explained deviance (%)` = sprintf(
                        "%.3f",
                        explained_deviance
                ),
                `Predicted at 0%` = sprintf(
                        "%.3f",
                        predicted_at_0
                ),
                `Predicted at 100%` = sprintf(
                        "%.3f",
                        predicted_at_100
                ),
                `Change (percentage points)` = sprintf(
                        "%.2f",
                        change_percentage_points
                ),
                `Relative change (%)` = sprintf(
                        "%.2f",
                        relative_change
                )
        )

h1_model_table

# Print a formatted table when knitr is available

if (requireNamespace("knitr", quietly = TRUE)) {
        
        knitr::kable(
                h1_model_table,
                align = c(
                        "l",
                        "l",
                        "l",
                        rep("r", 9)
                ),
                caption = paste(
                        "Protected-area coverage and habitat change.",
                        "Asterisks identify the model with the greatest",
                        "in-sample explained deviance within each response."
                )
        )
}

##### Publication-ready gt table ####

# Prepare the data used in the formatted table

h1_model_gt_data <- h1_model_results |>
        mutate(
                Hypothesis = case_when(
                        hypothesis == "H1a" ~
                                "H1a | Habitat retention",
                        hypothesis == "H1b" ~
                                "H1b | Compositional stability"
                ),
                Model = paste0(
                        model,
                        if_else(
                                best_in_sample,
                                "*",
                                ""
                        )
                )
        ) |>
        transmute(
                Hypothesis,
                Model,
                Estimate = estimate,
                SE = standard_error,
                F = F_value,
                P = format.pval(
                        p_value,
                        digits = 3,
                        eps = 0.001
                ),
                `Residual deviance` = residual_deviance,
                `Explained deviance` = explained_deviance,
                `Protection 0%` = predicted_at_0,
                `Protection 100%` = predicted_at_100,
                `Change` = change_percentage_points
        )

# Create the formatted gt table

h1_model_gt <- h1_model_gt_data |>
        gt::gt(
                groupname_col = "Hypothesis"
        ) |>
        gt::tab_header(
                title = gt::md(
                        "**Protected-area coverage and habitat change**"
                ),
                subtitle = paste(
                        "Quasibinomial models of habitat-loss probability",
                        "and temporal Jaccard dissimilarity"
                )
        ) |>
        gt::tab_spanner(
                label = "Model coefficient",
                columns = c(
                        Estimate,
                        SE
                )
        ) |>
        gt::tab_spanner(
                label = "Comparison with null model",
                columns = c(
                        F,
                        P
                )
        ) |>
        gt::tab_spanner(
                label = "Model fit",
                columns = c(
                        `Residual deviance`,
                        `Explained deviance`
                )
        ) |>
        gt::tab_spanner(
                label = "Predicted response",
                columns = c(
                        `Protection 0%`,
                        `Protection 100%`,
                        Change
                )
        ) |>
        gt::cols_label(
                Model = "Model",
                Estimate = "Estimate",
                SE = "SE",
                F = gt::html("<i>F</i>"),
                P = gt::html("<i>P</i>"),
                `Residual deviance` = "Residual deviance",
                `Explained deviance` =
                        gt::html("Explained deviance, <i>D</i><sup>2</sup>"),
                `Protection 0%` = "0%",
                `Protection 100%` = "100%",
                Change = "Change (pp)"
        ) |>
        gt::fmt_number(
                columns = c(
                        Estimate,
                        SE
                ),
                decimals = 4
        ) |>
        gt::fmt_number(
                columns = F,
                decimals = 2
        ) |>
        gt::fmt_number(
                columns = `Residual deviance`,
                decimals = 1,
                use_seps = TRUE
        ) |>
        gt::fmt_number(
                columns = `Explained deviance`,
                decimals = 3,
                pattern = "{x}%"
        ) |>
        gt::fmt_percent(
                columns = c(
                        `Protection 0%`,
                        `Protection 100%`
                ),
                decimals = 1
        ) |>
        gt::fmt_number(
                columns = Change,
                decimals = 2
        ) |>
        gt::tab_footnote(
                footnote = paste(
                        "Asterisks identify the model with the greatest",
                        "in-sample explained deviance within each response;",
                        "they do not indicate statistical significance."
                ),
                locations = gt::cells_body(
                        columns = Model,
                        rows = grepl(
                                "\\*$",
                                Model
                        )
                )
        ) |>
        gt::tab_footnote(
                footnote = paste(
                        "Estimates and standard errors are expressed",
                        "on the logit scale."
                ),
                locations = gt::cells_column_labels(
                        columns = c(
                                Estimate,
                                SE
                        )
                )
        ) |>
        gt::tab_footnote(
                footnote = paste(
                        "F tests and P values compare each protection",
                        "model separately with its corresponding null model."
                ),
                locations = gt::cells_column_labels(
                        columns = c(
                                F,
                                P
                        )
                )
        ) |>
        gt::tab_footnote(
                footnote = paste(
                        "Explained deviance is calculated relative to",
                        "the corresponding null model."
                ),
                locations = gt::cells_column_labels(
                        columns = `Explained deviance`
                )
        ) |>
        gt::tab_options(
                table.width = gt::pct(100),
                table.font.names = c(
                        "Arial",
                        "Helvetica",
                        "sans-serif"
                ),
                table.font.size = gt::px(11),
                heading.title.font.size = gt::px(14),
                heading.subtitle.font.size = gt::px(11),
                column_labels.font.weight = "bold",
                row_group.font.weight = "bold",
                row_group.background.color = "#F2F2F2",
                table.border.top.color = "black",
                table.border.top.width = gt::px(1.5),
                table.border.bottom.color = "black",
                table.border.bottom.width = gt::px(1.5),
                column_labels.border.top.color = "black",
                column_labels.border.bottom.color = "black",
                column_labels.border.bottom.width = gt::px(1),
                data_row.padding = gt::px(5),
                source_notes.font.size = gt::px(9)
        )

# Display the table in the Viewer

h1_model_gt

###### table export ######

write.csv(
        h1_model_results,
        file = file.path(
                "tables",
                "Table_H1_model_results.csv"
        ),
        row.names = FALSE
)

# Formatted HTML table:
# open it in a browser and copy-paste it into Word

gt::gtsave(
        data = h1_model_gt,
        filename = "Table_H1_model_results.html",
        path = "tables",
        inline_css = TRUE
)

# Editable Word document

gt::gtsave(
        data = h1_model_gt,
        filename = "Table_H1_model_results.docx",
        path = "tables"
)
