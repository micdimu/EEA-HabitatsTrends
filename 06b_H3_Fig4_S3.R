library(tidyverse)
library(patchwork)
library(glmmTMB)
library(gt)

source("SourceCodes/Source06_H3.R")


#### Load Models ####
# If they are not available run 06a first. This script can then run in a fresh R session.

h3_additive <- readRDS("processed/models/h3_additive.rds")
h3_interaction <- readRDS("processed/models/h3_interaction.rds")
occupancy_reference <- readRDS("processed/models/h3_occupancy_reference.rds")
h3_data <- readRDS("processed/models/h3_analysis_data.rds")

### Graphical output ####

beta <- glmmTMB::fixef(h3_interaction)$cond
beta_vcov <- as.matrix(vcov(h3_interaction)$cond)
beta_vcov <- beta_vcov[names(beta), names(beta), drop = FALSE]

# Recover the interaction P-value from the likelihood-ratio test.
model_comparison <- anova(h3_additive, h3_interaction)
interaction_p <- model_comparison[["Pr(>Chisq)"]][2]

interaction_label <- sprintf(
        "Protection × occupancy: likelihood-ratio P = %.3f",
        interaction_p
)

# Define evenly spaced coordinates on the log10 occupancy scale.
log_occupancy_grid <- seq(
        log10(min(h3_data$n_cells_t0)),
        log10(max(h3_data$n_cells_t0)),
        length.out = 200
)

occupancy_grid <- 10^log_occupancy_grid

# Use the same occupancy axis labels across figures.
occupancy_breaks <- c(1, 10, 100, 1000, 10000)
occupancy_breaks <- occupancy_breaks[
        occupancy_breaks >= min(h3_data$n_cells_t0) &
                occupancy_breaks <= max(h3_data$n_cells_t0)
]


##### Simulate fixed-effect uncertainty ####

# Use the same coefficient draw for both scenarios to preserve covariance.
set.seed(123)

beta_draws <- MASS::mvrnorm(
        n = 2000,
        mu = beta,
        Sigma = beta_vcov
)


##### Predict loss at selected protection levels ####

# Keep habitat occupancy continuous and compare three protection scenarios.
response_grid <- expand.grid(
        n_cells_t0 = occupancy_grid,
        protection = c(10, 50, 90)
)

response_grid <- response_grid |>
        mutate(
                pa_10 = protection / 10,
                log_occupancy_c = log(n_cells_t0) - occupancy_reference,
                protection_label = factor(
                        protection,
                        levels = c(10, 50, 90),
                        labels = c("10%", "50%", "90%")
                )
        )

# Construct the fixed-effect design matrix.
X <- model.matrix(
        ~ pa_10 * log_occupancy_c,
        data = response_grid
)[, names(beta), drop = FALSE]

# Calculate predictions and standard errors on the logit scale.
# Random effects are set to zero.
eta <- drop(X %*% beta)
eta_se <- sqrt(rowSums((X %*% beta_vcov) * X))

# Transform predictions and pointwise confidence limits to percentages.
response_grid <- response_grid |>
        mutate(
                loss_pct = 100 * plogis(eta),
                conf_low = 100 * plogis(eta - qnorm(0.975) * eta_se),
                conf_high = 100 * plogis(eta + qnorm(0.975) * eta_se)
        )


###### Plot loss probability across habitat occupancy ####

protection_colours <- c(
        "10%" = "#B66A40",
        "50%" = "#737373",
        "90%" = "#247A78"
)

# Display one occupancy value per habitat.
habitat_occupancy <- h3_data |>
        distinct(habitat, n_cells_t0)

h3_response_plot <- ggplot(
        response_grid,
        aes(
                x = n_cells_t0,
                y = loss_pct,
                colour = protection_label,
                fill = protection_label
        )
) +
        geom_ribbon(
                aes(ymin = conf_low, ymax = conf_high),
                alpha = 0.10,
                colour = NA,
                show.legend = FALSE
        ) +
        geom_line(linewidth = 1) +
        geom_rug(
                data = habitat_occupancy,
                aes(x = n_cells_t0),
                inherit.aes = FALSE,
                sides = "b",
                colour = "grey40",
                alpha = 0.25
        ) +
        scale_colour_manual(values = protection_colours) +
        scale_fill_manual(values = protection_colours) +
        scale_x_log10(
                breaks = occupancy_breaks,
                labels = scales::label_comma()
        ) +
        scale_y_continuous(
                limits = c(0, NA),
                expand = expansion(mult = c(0, 0.05))
        ) +
        labs(
                x = "Initial habitat occupancy (number of cells)",
                y = "Predicted occurrence loss (%)",
                colour = "Protected-area coverage"
        ) +
        theme_classic(base_size = 12) +
        theme(
                legend.position = "top",
                legend.title = element_text(size = 10),
                legend.text = element_text(size = 10)
        )

##### Calculate the contrast between protection scenarios ######

# Compare 90% versus 10% protection at each occupancy value.
low_protection <- data.frame(
        pa_10 = 1,
        log_occupancy_c = log(occupancy_grid) - occupancy_reference
)

high_protection <- low_protection
high_protection$pa_10 <- 9

X_low <- model.matrix(
        ~ pa_10 * log_occupancy_c,
        data = low_protection
)[, names(beta), drop = FALSE]

X_high <- model.matrix(
        ~ pa_10 * log_occupancy_c,
        data = high_protection
)[, names(beta), drop = FALSE]

# Use the same coefficient draws for both protection scenarios.
contrast_draws <- 100 * (
        plogis(X_high %*% t(beta_draws)) -
                plogis(X_low %*% t(beta_draws))
)

contrast_intervals <- t(apply(
        contrast_draws,
        1,
        quantile,
        probs = c(0.025, 0.975),
        names = FALSE
))

scenario_contrast <- data.frame(
        n_cells_t0 = occupancy_grid,
        estimate = 100 * (
                plogis(drop(X_high %*% beta)) -
                        plogis(drop(X_low %*% beta))
        ),
        conf_low = contrast_intervals[, 1],
        conf_high = contrast_intervals[, 2]
)


###### Plot the protection scenario contrast ####

h3_scenario_plot <- ggplot(
        scenario_contrast,
        aes(x = n_cells_t0, y = estimate)
) +
        geom_hline(
                yintercept = 0,
                linetype = "dashed",
                colour = "grey45",
                linewidth = 0.5
        ) +
        geom_ribbon(
                aes(ymin = conf_low, ymax = conf_high),
                fill = "#247A78",
                alpha = 0.15
        ) +
        geom_line(
                colour = "#247A78",
                linewidth = 1
        ) +
        geom_rug(
                data = habitat_occupancy,
                aes(x = n_cells_t0),
                inherit.aes = FALSE,
                sides = "b",
                colour = "grey40",
                alpha = 0.25
        ) +
        scale_x_log10(
                breaks = occupancy_breaks,
                labels = scales::label_comma()
        ) +
        labs(
                x = "Initial habitat occupancy (number of cells)",
                y = "Difference in loss probability\n(90% − 10% protection; percentage points)"
        ) +
        theme_classic(base_size = 12)


##### Combine response curves and scenario contrast ####

figure_4 <- (
        h3_response_plot / h3_scenario_plot
) +
        plot_layout(ncol = 1, heights = c(1, 1), guides = "collect") +
        plot_annotation(tag_levels = "A") &
        theme(legend.position = "top")

figure_4


##### Predict loss across protection and occupancy ####

prediction_grid <- expand.grid(
        protection = seq(0, 100, by = 1),
        log10_occupancy = log_occupancy_grid
)

prediction_grid <- prediction_grid |>
        mutate(
                n_cells_t0 = 10^log10_occupancy,
                pa_10 = protection / 10,
                log_occupancy_c = log(n_cells_t0) - occupancy_reference
        )

# Calculate predictions with habitat and cell random effects set to zero.
X <- model.matrix(
        ~ pa_10 * log_occupancy_c,
        data = prediction_grid
)[, names(beta), drop = FALSE]

prediction_grid$loss_pct <- 100 * plogis(
        drop(X %*% beta)
)


###### Plot predicted loss probability ####

# Plot on log10 coordinates to obtain a regularly spaced raster.
# Contour lines connect equal predicted loss probabilities.
h3_surface_plot <- ggplot(
        prediction_grid,
        aes(x = protection, y = log10_occupancy)
) +
        geom_raster(aes(fill = loss_pct)) +
        geom_contour(
                aes(z = loss_pct),
                bins = 6,
                colour = "white",
                linewidth = 0.35,
                alpha = 0.7
        ) +
        scale_fill_viridis_c(
                option = "C",
                name = "Predicted loss\nprobability (%)"
        ) +
        scale_x_continuous(
                breaks = seq(0, 100, 20),
                expand = expansion(mult = 0)
        ) +
        scale_y_continuous(
                breaks = log10(occupancy_breaks),
                labels = scales::label_comma()(occupancy_breaks),
                expand = expansion(mult = 0)
        ) +
        labs(
                x = "Protected-area coverage (%)",
                y = "Initial habitat occupancy (number of cells)"
        ) +
        theme_classic(base_size = 12) +
        theme(
                plot.subtitle = element_text(size = 10),
                legend.title = element_text(size = 10)
        )




#### Prepare supplementary contrasts ####

# Average equally over baseline coverage from 0% to 90%.
# These contrasts describe the full theoretical protection gradient.
protection_grid <- seq(0, 9, by = 0.1)


##### Calculate average protection contrasts ####

h3_contrasts <- lapply(
        occupancy_grid,
        average_h3_contrast,
        protection_grid = protection_grid,
        occupancy_reference = occupancy_reference,
        beta = beta,
        beta_draws = beta_draws
) |>
        bind_rows()


###### Plot supplementary protection contrasts ####

h3_contrast_plot <- ggplot(
        h3_contrasts,
        aes(x = n_cells_t0, y = estimate)
) +
        geom_hline(
                yintercept = 0,
                linetype = "dashed",
                colour = "grey45",
                linewidth = 0.5
        ) +
        geom_ribbon(
                aes(ymin = conf_low, ymax = conf_high),
                fill = "#247A78",
                alpha = 0.20
        ) +
        geom_line(
                colour = "#247A78",
                linewidth = 1
        ) +
        geom_rug(
                data = habitat_occupancy,
                aes(x = n_cells_t0),
                inherit.aes = FALSE,
                sides = "b",
                alpha = 0.25
        ) +
        scale_x_log10(
                breaks = occupancy_breaks,
                labels = scales::label_comma()
        ) +
        labs(
                x = "Initial habitat occupancy (number of cells)",
                y = "Change in loss probability\n(+10 points of protection; percentage points)",
                caption = paste(
                        "Contrasts for +10 percentage points of protected-area coverage;",
                        "shading shows pointwise 95% confidence intervals."
                )
        ) +
        theme_classic(base_size = 12) +
        theme(
                plot.subtitle = element_text(size = 10),
                plot.caption = element_text(hjust = 0, size = 9)
        )


###### Combine supplementary panels #####

# Panel A: predicted probability across protection and occupancy.
# Panel B: average +10-percentage-point contrast with pointwise 95% intervals.
# These intervals are not an uncertainty map for panel A.
figure_s_h3 <- (
        h3_surface_plot / h3_contrast_plot
) +
        plot_layout(ncol = 1, heights = c(1, 1), guides = "keep") +
        plot_annotation(tag_levels = "A")

figure_s_h3


##### Save graphical outputs ######

# Create output directories if they do not already exist.
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/supplementary", recursive = TRUE, showWarnings = FALSE)

# Save the main figure as a high-resolution TIFF.
ggsave(
        filename = "output/figures/Figure_4.tiff",
        plot = figure_4,
        device = "tiff",
        width = 180,
        height = 220,
        units = "mm",
        dpi = 600,
        compression = "lzw",
        bg = "white"
)

# Save the supplementary figure using the same export settings.
ggsave(
        filename = "output/supplementary/Figure_S_H3.tiff",
        plot = figure_s_h3,
        device = "tiff",
        width = 260,
        height = 240,
        units = "mm",
        dpi = 600,
        compression = "lzw",
        bg = "white"
)

# Preserve prediction tables, simulation settings, and plotting objects.
saveRDS(
        list(
                response_grid = response_grid,
                scenario_contrast = scenario_contrast,
                prediction_grid = prediction_grid,
                average_contrasts = h3_contrasts,
                occupancy_reference = occupancy_reference,
                interaction_p = interaction_p,
                seed = 123,
                n_draws = nrow(beta_draws),
                baseline_protection_pct = protection_grid * 10,
                random_effects = "Set to zero",
                figure_4 = figure_4,
                figure_s_h3 = figure_s_h3
        ),
        file = "output/figures/H3_graphical_outputs.rds"
)


#### Prepare supplementary tables ####

table_path <- file.path("output", "supplementary", "tables")

dir.create(
        table_path,
        recursive = TRUE,
        showWarnings = FALSE
)


##### Summarise habitat occurrence loss ####

# Calculate descriptive statistics within the H3 analysis sample.
# For a binary response, mean(loss) is the proportion of baseline occurrences lost.
h3_habitat_results <- h3_data |>
        group_by(habitat) |>
        summarise(
                baseline_cells = n(),
                lost_cells = sum(loss),
                maintained_cells = sum(loss == 0L),
                loss_pct = 100 * mean(loss),
                persistence_pct = 100 * mean(loss == 0L),
                .groups = "drop"
        ) |>
        mutate(
                habitat = as.character(habitat),
                relative_loss_rank = min_rank(desc(loss_pct)),
                absolute_loss_rank = min_rank(desc(lost_cells))
        ) |>
        arrange(desc(loss_pct), desc(lost_cells), habitat) |>
        select(
                habitat,
                baseline_cells,
                lost_cells,
                maintained_cells,
                loss_pct,
                persistence_pct,
                relative_loss_rank,
                absolute_loss_rank
        )

# Calculate overall summaries for the table notes.
pooled_loss_pct <- with(
        h3_habitat_results,
        100 * sum(lost_cells) / sum(baseline_cells)
)

mean_habitat_loss <- mean(h3_habitat_results$loss_pct)
median_habitat_loss <- median(h3_habitat_results$loss_pct)
habitat_loss_iqr <- quantile(
        h3_habitat_results$loss_pct,
        probs = c(0.25, 0.75),
        names = FALSE
)


##### Format the habitat summary table ####

h3_habitat_gt <- h3_habitat_results |>
        gt() |>
        tab_header(
                title = "Habitat occurrence loss between reporting periods",
                subtitle = "Descriptive statistics for the H3 analysis sample"
        ) |>
        cols_label(
                habitat = "Habitat code",
                baseline_cells = "Baseline cells",
                lost_cells = "Lost cells",
                maintained_cells = "Maintained cells",
                loss_pct = "Loss (%)",
                persistence_pct = "Persistence (%)",
                relative_loss_rank = "Relative loss rank",
                absolute_loss_rank = "Absolute loss rank"
        ) |>
        fmt_integer(
                columns = c(
                        baseline_cells,
                        lost_cells,
                        maintained_cells,
                        relative_loss_rank,
                        absolute_loss_rank
                )
        ) |>
        fmt_number(
                columns = c(loss_pct, persistence_pct),
                decimals = 2
        ) |>
        tab_source_note(
                source_note = paste(
                        "Baseline: 2007–2012; follow-up: 2013–2018.",
                        "Both periods are aligned after upstream completion of missing records with zero reported habitats.",
                        "Loss denotes a baseline presence no longer reported at follow-up;",
                        "persistence denotes a presence reported in both periods.",
                        "Gains are excluded because the analysis conditions on baseline presence."
                )
        ) |>
        tab_source_note(
                source_note = paste(
                        "Rows are ordered by decreasing percentage loss.",
                        "Rank 1 denotes the highest loss; ties receive the same rank.",
                        "Relative rankings should be interpreted alongside baseline cell counts:",
                        "large percentages may reflect very few initial occurrences."
                )
        ) |>
        tab_source_note(
                source_note = sprintf(
                        paste(
                                "Across habitats, mean loss = %.2f%%;",
                                "median = %.2f%% (IQR %.2f–%.2f%%).",
                                "Pooled loss across all baseline occurrences = %.2f%%.",
                                "The mean across habitats gives each habitat equal weight;",
                                "the pooled proportion weights habitats by baseline occupancy."
                        ),
                        mean_habitat_loss,
                        median_habitat_loss,
                        habitat_loss_iqr[1],
                        habitat_loss_iqr[2],
                        pooled_loss_pct
                )
        ) |>
        tab_options(
                table.font.size = px(10),
                data_row.padding = px(3)
        )


##### Extract detailed model results ####

# Store results in a long table to accommodate coefficients and model summaries.
# Confidence intervals for fixed effects are approximate Wald intervals.


h3_model_results <- bind_rows(
        extract_h3_results(h3_additive, "Additive"),
        extract_h3_results(h3_interaction, "Interaction")
)


##### Add the interaction test ####

model_comparison <- anova(h3_additive, h3_interaction)

comparison_results <- data.frame(
        model = "Model comparison",
        section = "Likelihood-ratio test",
        term = c(
                "Interaction versus additive: chi-square",
                "Difference in degrees of freedom",
                "AIC difference: interaction minus additive"
        ),
        estimate = c(
                model_comparison[["Chisq"]][2],
                model_comparison[["Chi Df"]][2],
                AIC(h3_interaction) - AIC(h3_additive)
        ),
        p_value = c(
                model_comparison[["Pr(>Chisq)"]][2],
                NA_real_,
                NA_real_
        )
)

h3_model_results <- bind_rows(
        h3_model_results,
        comparison_results
) |>
        select(
                model, section, term,
                estimate, std_error,
                conf_low, conf_high,
                z_value, p_value
        )


##### Format the model results table ####

# Combine confidence limits for display while retaining numeric columns in CSV.
h3_model_display <- h3_model_results |>
        mutate(
                table_group = paste(model, section, sep = " — "),
                confidence_interval = ifelse(
                        is.na(conf_low),
                        NA_character_,
                        sprintf("%.4f to %.4f", conf_low, conf_high)
                )
        ) |>
        select(
                table_group, term, estimate,
                std_error, confidence_interval,
                z_value, p_value
        )

h3_model_gt <- h3_model_display |>
        gt(groupname_col = "table_group") |>
        tab_header(
                title = "Models of habitat occurrence loss",
                subtitle = "Binomial mixed models with a logit link"
        ) |>
        cols_label(
                term = "Term or statistic",
                estimate = "Estimate / value",
                std_error = "SE",
                confidence_interval = "95% CI",
                z_value = "z",
                p_value = "P"
        ) |>
        fmt_number(
                columns = c(estimate, std_error),
                decimals = 4
        ) |>
        fmt_number(
                columns = z_value,
                decimals = 3
        ) |>
        fmt(
                columns = p_value,
                fns = function(x) {
                        ifelse(
                                is.na(x),
                                NA_character_,
                                ifelse(
                                        x < 0.001,
                                        "<0.001",
                                        sprintf("%.3f", x)
                                )
                        )
                }
        ) |>
        sub_missing(missing_text = "—") |>
        tab_source_note(
                source_note = paste(
                        "Response: loss = 1 for a baseline occurrence no longer reported",
                        "at follow-up; loss = 0 for a maintained occurrence.",
                        "Both models were fitted to the same baseline habitat–cell records,",
                        "without land-area weights."
                )
        ) |>
        tab_source_note(
                source_note = paste(
                        "Additive fixed effects: pa_10 + log_occupancy_c.",
                        "Interaction fixed effects: pa_10 * log_occupancy_c.",
                        "Both models include (1 + pa_10 | habitat) + (1 | cell_id).",
                        "The habitat intercept and protection slope are correlated random effects."
                )
        ) |>
        tab_source_note(
                source_note = sprintf(
                        paste(
                                "pa_10 = total protected-area coverage / 10.",
                                "log_occupancy_c = ln(initial occupied cells) minus %.6f.",
                                "This centring constant is the mean log occupancy across habitats;",
                                "the corresponding reference occupancy is %.2f cells.",
                                "In the interaction model, the protection coefficient applies",
                                "at this reference occupancy, and the occupancy coefficient",
                                "applies at zero protection."
                        ),
                        occupancy_reference,
                        exp(occupancy_reference)
                )
        ) |>
        tab_source_note(
                source_note = paste(
                        "Fixed-effect estimates are on the log-odds scale.",
                        "Their P-values are two-sided Wald z tests.",
                        "Exponentiated slope coefficients are odds ratios;",
                        "the exponentiated intercept represents baseline odds.",
                        "The exponentiated interaction is a ratio of protection odds ratios",
                        "per one-unit increase in log occupancy, not a probability ratio."
                )
        ) |>
        tab_source_note(
                source_note = paste(
                        "Fixed-effect confidence intervals are approximate 95% Wald intervals;",
                        "exponentiated intervals are obtained by exponentiating their limits.",
                        "Random-effect variances and SDs refer to the logit model.",
                        "No confidence intervals are reported for variance components."
                )
        ) |>
        tab_source_note(
                source_note = paste(
                        "The likelihood-ratio test compares the additive and interaction models.",
                        "A negative AIC difference favours the interaction model.",
                        "Cell random intercepts account for records sharing a cell",
                        "but do not model dependence between neighbouring cells.",
                        "The models describe associations with cell-level protection,",
                        "not causal effects of habitat-level protection."
                )
        ) |>
        tab_options(
                table.font.size = px(10),
                data_row.padding = px(3)
        )


#### Export supplementary tables ####

# Save numeric data, formatted HTML, and editable Word tables.
tables_to_export <- list(
        Table_S_H3_habitat_loss = list(
                data = h3_habitat_results,
                formatted = h3_habitat_gt
        ),
        Table_S_H3_model_results = list(
                data = h3_model_results,
                formatted = h3_model_gt
        )
)

for (table_name in names(tables_to_export)) {
        
        current_table <- tables_to_export[[table_name]]
        
        write.csv(
                current_table$data,
                file = file.path(table_path, paste0(table_name, ".csv")),
                row.names = FALSE
        )
        
        gt::gtsave(
                data = current_table$formatted,
                filename = paste0(table_name, ".html"),
                path = table_path,
                inline_css = TRUE
        )
        
        gt::gtsave(
                data = current_table$formatted,
                filename = paste0(table_name, ".docx"),
                path = table_path
        )
}
