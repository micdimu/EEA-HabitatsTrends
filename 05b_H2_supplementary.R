library(tidyverse)
library(sf)
library(patchwork)
library(gt)

source("SourceCodes/Source05_H2.R")

#### Load data ####

# Run from the project root. Keep double absences (0 -> 0): their
# beta_jaccard is undefined (NA), but their gain/loss counts are valid zeros.
# Filter missing values only for variables required by each analysis.
temporal_beta_group <- read.csv("processed/temporal_group_cell_by_cell.csv", check.names = FALSE) |>
        mutate(net_change = richness_2013_2018 - richness_2007_2012)

grid_pa <- st_read("data/PAs/derived/europe_10km_protected_area_coverage.gpkg")


grid_group <- grid_pa |>
        right_join(temporal_beta_group, by = c("CellCode" = "cell_id")) |> 
        mutate(taxGroup = as.factor(taxGroup))

#### H2b: mosaic analysis ####
# Six pairwise models: gain of one group versus loss of another group.
# Shrubland labels refer exclusively to Sclerophyllous scrub.
# Start from grid_group, not d_loss or d_balance, to retain response-group
# absences in both periods and colonisations from initial absence.
# 1. Prepare long data and then create one row per cell -----------------

mosaic_long <- grid_group |>
        sf::st_drop_geometry() |>
        filter(taxGroup %in% c(
                "Grasslands", "Sclerophyllous scrub", "Forests"
        )) |>
        transmute(
                CellCode,
                total_pa_land_cov,
                land_weight = land_cov / 100,
                group = recode(
                        as.character(taxGroup),
                        "Grasslands" = "grass",
                        "Sclerophyllous scrub" = "scrub",
                        "Forests" = "forest"
                ),
                lost = lost_species,
                gained = gained_species,
                initial = richness_2007_2012
        )

stopifnot(
        anyDuplicated(mosaic_long[c("CellCode", "group")]) == 0L
)


# Rebuild this object on every run: never reuse an older in-memory copy.
# Valid double absences remain zero. Missing records remain NA, because
# missing information must not be interpreted automatically as absence.
mosaic_data <- mosaic_long |>
        pivot_wider(
                id_cols = c(CellCode, total_pa_land_cov, land_weight),
                names_from = group,
                values_from = c(lost, gained, initial),
                names_glue = "{group}_{.value}"
        )

stopifnot(anyDuplicated(mosaic_data$CellCode) == 0L)

# 2. Define the six comparisons in panel order -------------------------

comparisons <- tibble(
        gain_group = c("grass", "grass", "scrub", "scrub", "forest", "forest"),
        loss_group = c("scrub", "forest", "grass", "forest", "grass", "scrub"),
        gain_label = c(
                "Grassland", "Grassland", "Shrubland",
                "Shrubland", "Forest", "Forest"
        ),
        loss_label = c(
                "Shrubland", "Forest", "Grassland",
                "Forest", "Grassland", "Shrubland"
        )
)

# 3. Fit each model and generate predictions ---------------------------



# 4. Run all six comparisons -------------------------------------------

mosaic_results <- lapply(seq_len(nrow(comparisons)), function(i) {
        fit_mosaic(
                gain_group = comparisons$gain_group[i],
                loss_group = comparisons$loss_group[i],
                gain_label = comparisons$gain_label[i],
                loss_label = comparisons$loss_label[i],
                data = mosaic_data
        )
})

names(mosaic_results) <- paste(
        comparisons$gain_group,
        comparisons$loss_group,
        sep = "_gain_vs_"
)

#### Probability contrasts at selected protection levels ####



mosaic_contrasts <- bind_rows(
        lapply(mosaic_results, contrast_mosaic)
)

mosaic_contrasts |>
        mutate(across(where(is.numeric), \(x) round(x, 2))) |>
        print(n = Inf, width = Inf)

# Holm adjustment across the six interaction tests.
interaction_tests <- bind_rows(
        lapply(mosaic_results, function(x) x$test)
) |>
        mutate(interaction_p_holm = p.adjust(interaction_p, method = "holm"))

interaction_tests

# Example: inspect an individual model
summary(mosaic_results$forest_gain_vs_grass$model)

# 5. Assemble three rows and two columns -------------------------------

p_mosaic <- (
        wrap_plots(
                lapply(mosaic_results, function(x) x$plot),
                ncol = 2,
                guides = "collect"
        ) +
                plot_annotation(tag_levels = "a")
) & theme(legend.position = "bottom")

p_mosaic


#### Diagnostics: shrubland gains versus grassland losses ####
# Frequencies and mean predictions below are unweighted; model fitting uses
# terrestrial-fraction weights. Pearson residuals include those weights.
d <- mosaic_results$scrub_gain_vs_grass$data

d |>
        group_by(loss) |>
        summarise(
                n_cells = n(),
                n_gains = sum(gain),
                gain_frequency = mean(gain),
                .groups = "drop"
        )

summary(mosaic_results$scrub_gain_vs_grass$model)$dispersion

m <- mosaic_results$scrub_gain_vs_grass$model

check_scrub <- mosaic_results$scrub_gain_vs_grass$data |>
        mutate(
                predicted = fitted(m),
                pearson_sq = residuals(m, type = "pearson")^2
        )

# Observed gains and predictions by initial shrubland richness
check_scrub |>
        group_by(initial_gain) |>
        summarise(
                n_cells = n(),
                n_gains = sum(gain),
                observed_frequency = mean(gain),
                mean_prediction = mean(predicted),
                pearson_contribution = sum(pearson_sq),
                .groups = "drop"
        )

# Cells contributing most to the estimated dispersion
check_scrub |>
        arrange(desc(pearson_sq)) |>
        select(
                CellCode, initial_gain, initial_loss,
                loss, coverage, gain, predicted, pearson_sq
        ) |>
        head(15)


# Confirm that double absences remain available in the source data.
grid_group |>
        sf::st_drop_geometry() |>
        filter(
                taxGroup == "Sclerophyllous scrub",
                richness_2007_2012 == 0
        ) |>
        count(
                richness_2013_2018,
                name = "n_cells"
        )

# Check that the fitted sample includes non-colonised, initially empty cells.
mosaic_results$scrub_gain_vs_grass$data |>
        filter(initial_gain == 0) |>
        count(gain, name = "n_cells")

#### Summary table: mosaic probability contrasts ####

# Format probability differences and their pointwise confidence intervals.
mosaic_table_data <- mosaic_contrasts |>
        mutate(
                contrast = sprintf(
                        "%+.2f [%+.2f, %+.2f]",
                        difference_pp, lower_pp, upper_pp
                ),
                protection = paste0("pa_", protection_pct)
        ) |>
        select(response, predictor, protection, contrast) |>
        pivot_wider(
                names_from = protection,
                values_from = contrast
        ) |>
        left_join(
                interaction_tests |>
                        select(
                                response, predictor,
                                n_cells, interaction_p_holm
                        ),
                by = c("response", "predictor")
        ) |>
        select(
                response, predictor, n_cells,
                pa_0, pa_30, pa_100,
                interaction_p_holm
        )

mosaic_gt <- mosaic_table_data |>
        gt(groupname_col = "response") |>
        tab_header(
                title = "Habitat gains associated with losses in other groups",
                subtitle = "Predicted probability differences along protected-area coverage"
        ) |>
        cols_label(
                predictor = "Loss predictor",
                n_cells = "Cells",
                pa_0 = "0%",
                pa_30 = "30%",
                pa_100 = "100%",
                interaction_p_holm = "Interaction P"
        ) |>
        tab_spanner(
                label = "Protected terrestrial area",
                columns = c(pa_0, pa_30, pa_100)
        ) |>
        fmt_integer(
                columns = n_cells,
                use_seps = TRUE
        ) |>
        fmt(
                columns = interaction_p_holm,
                fns = function(x) {
                        ifelse(
                                is.na(x), NA_character_,
                                ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
                        )
                }
        ) |>
        cols_align(
                align = "left",
                columns = predictor
        ) |>
        cols_align(
                align = "center",
                columns = c(
                        n_cells, pa_0, pa_30, pa_100,
                        interaction_p_holm
                )
        ) |>
        tab_source_note(
                source_note = paste(
                        "Values are probability differences in percentage points",
                        "(Loss − No loss), with pointwise 95% confidence intervals",
                        "in brackets. Positive values indicate a higher probability",
                        "of habitat gain in cells with loss of the predictor group."
                )
        ) |>
        tab_source_note(
                source_note = paste(
                        "Predictions hold initial richness of both groups at their",
                        "model-specific medians. Only cells with the predictor group",
                        "initially present are included; initial absence of the",
                        "response group is allowed. Shrubland denotes",
                        "Sclerophyllous scrub."
                )
        ) |>
        tab_source_note(
                source_note = paste(
                        "Interaction P values are Holm-adjusted across six tests.",
                        "Confidence intervals are not adjusted for multiple comparisons.",
                        "Inference does not account for spatial dependence between cells."
                )
        ) |>
        tab_options(
                table.font.names = "Arial",
                table.font.size = px(12),
                heading.align = "left",
                row_group.font.weight = "bold",
                data_row.padding = px(8),
                source_notes.font.size = px(10)
        )

mosaic_gt


#### H2c: proportional habitat losses associated with other-group gains ####

# Rebuild wide data including shared habitats, required for loss denominators.
mosaic_loss_data <- grid_group |>
        sf::st_drop_geometry() |>
        filter(taxGroup %in% c(
                "Grasslands", "Sclerophyllous scrub", "Forests"
        )) |>
        transmute(
                CellCode,
                coverage = total_pa_land_cov,
                land_weight = land_cov / 100,
                group = recode(
                        as.character(taxGroup),
                        "Grasslands" = "grass",
                        "Sclerophyllous scrub" = "scrub",
                        "Forests" = "forest"
                ),
                lost = lost_species,
                shared = shared_species,
                gained = gained_species,
                initial = richness_2007_2012
        ) |>
        pivot_wider(
                id_cols = c(CellCode, coverage, land_weight),
                names_from = group,
                values_from = c(lost, shared, gained, initial),
                names_glue = "{group}_{.value}"
        )

stopifnot(anyDuplicated(mosaic_loss_data$CellCode) == 0L)

# Panel order: three response groups in rows, two gain predictors in columns.
loss_comparisons <- tibble(
        loss_group = c("grass", "grass", "scrub", "scrub", "forest", "forest"),
        gain_group = c("scrub", "forest", "grass", "forest", "grass", "scrub"),
        loss_label = c(
                "Grassland", "Grassland", "Shrubland",
                "Shrubland", "Forest", "Forest"
        ),
        gain_label = c(
                "Shrubland", "Forest", "Grassland",
                "Forest", "Grassland", "Shrubland"
        )
)



mosaic_loss_results <- lapply(seq_len(nrow(loss_comparisons)), function(i) {
        fit_mosaic_loss(
                loss_group = loss_comparisons$loss_group[i],
                gain_group = loss_comparisons$gain_group[i],
                loss_label = loss_comparisons$loss_label[i],
                gain_label = loss_comparisons$gain_label[i],
                data = mosaic_loss_data
        )
})

names(mosaic_loss_results) <- paste(
        loss_comparisons$loss_group,
        loss_comparisons$gain_group,
        sep = "_loss_vs_"
)

# Adjust across the six interaction tests in this loss analysis.
loss_interaction_tests <- bind_rows(
        lapply(mosaic_loss_results, function(x) x$test)
) |>
        mutate(
                interaction_p_holm = p.adjust(interaction_p, method = "holm")
        )

loss_interaction_tests

p_mosaic_loss <- (
        wrap_plots(
                lapply(mosaic_loss_results, function(x) x$plot),
                ncol = 2,
                guides = "collect"
        ) +
                plot_annotation(tag_levels = "a")
) & theme(legend.position = "bottom")

p_mosaic_loss

#### Summary table: proportional habitat loss contrasts ####



mosaic_loss_contrasts <- bind_rows(
        lapply(mosaic_loss_results, contrast_mosaic_loss)
)

# Format contrasts as percentage-point differences with 95% confidence intervals.
mosaic_loss_table_data <- mosaic_loss_contrasts |>
        mutate(
                contrast = sprintf(
                        "%+.2f [%+.2f, %+.2f]",
                        difference_pp, lower_pp, upper_pp
                ),
                protection = paste0("pa_", protection_pct)
        ) |>
        select(response, predictor, protection, contrast) |>
        pivot_wider(
                names_from = protection,
                values_from = contrast
        ) |>
        left_join(
                loss_interaction_tests |>
                        select(
                                response, predictor,
                                n_cells, interaction_p_holm
                        ),
                by = c("response", "predictor")
        )

# Retain all requested columns even if a coverage level is outside model support.
for (column in c("pa_0", "pa_30", "pa_100")) {
        if (!column %in% names(mosaic_loss_table_data)) {
                mosaic_loss_table_data[[column]] <- NA_character_
        }
}

mosaic_loss_table_data <- mosaic_loss_table_data |>
        select(
                response, predictor, n_cells,
                pa_0, pa_30, pa_100,
                interaction_p_holm
        )

mosaic_loss_gt <- mosaic_loss_table_data |>
        gt(groupname_col = "response") |>
        tab_header(
                title = "Habitat losses associated with gains in other groups",
                subtitle = "Predicted proportional-loss differences along protected-area coverage"
        ) |>
        cols_label(
                predictor = "Gain predictor",
                n_cells = "Cells",
                pa_0 = "0%",
                pa_30 = "30%",
                pa_100 = "100%",
                interaction_p_holm = "Interaction P"
        ) |>
        tab_spanner(
                label = "Protected terrestrial area",
                columns = c(pa_0, pa_30, pa_100)
        ) |>
        fmt_integer(
                columns = n_cells,
                use_seps = TRUE
        ) |>
        fmt(
                columns = interaction_p_holm,
                fns = function(x) {
                        ifelse(
                                is.na(x), NA_character_,
                                ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
                        )
                }
        ) |>
        sub_missing(
                columns = c(pa_0, pa_30, pa_100),
                missing_text = "Not estimated"
        ) |>
        cols_align(
                align = "left",
                columns = predictor
        ) |>
        cols_align(
                align = "center",
                columns = c(
                        n_cells, pa_0, pa_30, pa_100,
                        interaction_p_holm
                )
        ) |>
        tab_source_note(
                source_note = paste(
                        "Values are differences in predicted proportional habitat",
                        "loss in percentage points (Gain − No gain), with pointwise",
                        "95% confidence intervals in brackets. Positive values indicate",
                        "greater proportional loss in cells with gain of the predictor group."
                )
        ) |>
        tab_source_note(
                source_note = paste(
                        "Loss is defined as lost / (lost + shared). Only cells with",
                        "the response group initially present are included.",
                        "Predictor-group initial richness is held at its model-specific",
                        "median. Shrubland denotes Sclerophyllous scrub."
                )
        ) |>
        tab_source_note(
                source_note = paste(
                        "Interaction P values are Holm-adjusted across the six loss models.",
                        "Confidence intervals are not adjusted for multiple comparisons.",
                        "Inference does not account for spatial dependence between cells.",
                        "Contrasts outside the shared coverage range are not estimated."
                )
        ) |>
        tab_options(
                table.font.names = "Arial",
                table.font.size = px(12),
                heading.align = "left",
                row_group.font.weight = "bold",
                data_row.padding = px(8),
                source_notes.font.size = px(10)
        )

mosaic_loss_gt

# Export the formatted table and the underlying numeric contrasts.
# dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
# 
# gtsave(
#         mosaic_loss_gt,
#         filename = "Table_H2_mosaic_loss_contrasts.html",
#         path = "output/tables"
# )
# 
# write.csv(
#         mosaic_loss_contrasts,
#         file = "output/tables/Table_H2_mosaic_loss_contrasts.csv",
#         row.names = FALSE
# )