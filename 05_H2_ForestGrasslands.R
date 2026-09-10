library(tidyverse)
library(sf)
library(terra)
library(patchwork)
library(ggeffects)
source("Source.R")

#### Load data ####


temporal_beta_group <- read.csv("processed/temporal_group_cell_by_cell.csv") |> 
        mutate(net_change = richness_2013_2018 - richness_2007_2012) |> 
        drop_na()

grid_pa <- st_read("data/PAs/derived/europe_10km_protected_area_coverage.gpkg")


grid_group <- grid_pa |>
        right_join(temporal_beta_group, by = c("CellCode" = "cell_id")) |> 
        mutate(taxGroup = as.factor(taxGroup))

#        filter(taxGroup %in% c("Forests", "Sclerophyllous scrub", "Grasslands")) |> 


glimpse(grid_group)
grid_group$taxGroup |> levels()


lostG_model <- glm(formula = lost_species ~ log(total_pa_cov + 1) * taxGroup , data = grid_group, family = poisson, weights = land_cov)
lostG_model <- glm(formula = gained_species ~ log(total_pa_cov + 1) * taxGroup , data = grid_group, family = gaussian, weights = land_cov)

summary(lostG_model)

ggpredict(lostG_model) |> 
        plot()



# Griglia di valori per la predizione
pred_data <- expand_grid(
        total_pa_cov = seq(
                min(grid_group$total_pa_cov, na.rm = TRUE),
                max(grid_group$total_pa_cov, na.rm = TRUE),
                length.out = 200
        ),
        taxGroup = unique(grid_group$taxGroup)
)

# Predizioni sulla scala del link
pred <- predict(
        lostG_model,
        newdata = pred_data,
        type = "link",
        se.fit = TRUE
)

# Ritorno alla scala originale della risposta
pred_data <- pred_data |>
        mutate(
                fit = exp(pred$fit),
                lower = exp(pred$fit - 1.96 * pred$se.fit),
                upper = exp(pred$fit + 1.96 * pred$se.fit)
        )

# Grafico
ggplot( data = pred_data)+
        geom_ribbon(aes(x = total_pa_cov,
                        y = fit,
                        ymin = lower,
                        ymax = upper,
                        fill = taxGroup,
                        col = NULL
                ),
                alpha = 0.15
        ) +
        geom_line(aes(x = total_pa_cov,
                    y = fit,
                    col = taxGroup),
                linewidth = 1
        ) +
        labs(
                x = "Protected area coverage (%)",
                y = "Number of gained Habitat",
                colour = "Taxonomic group",
                fill = "Taxonomic group"
        ) +
        theme_classic() +
        theme(legend.position = "top"
        )+
        geom_hline(yintercept = 0)

### testare se l'aumento delle foreste o delle shrub è dovuto a una diminuzione delle praterie





library(dplyr)
library(ggplot2)
library(sandwich)
library(lmtest)

# 1. Dati e riepilogo per TUTTI i gruppi -----------------------------

dat <- sf::st_drop_geometry(grid_group) |>
        filter(
                !is.na(CellCode),
                !is.na(taxGroup),
                is.finite(total_pa_land_cov),
                is.finite(land_cov),
                land_cov > 0,
                !is.na(lost_species),
                !is.na(shared_species),
                !is.na(gained_species)
        ) |>
        mutate(
                land_weight = land_cov / 100,
                initial = lost_species + shared_species,
                final = gained_species + shared_species,
                changes = lost_species + gained_species
        )

# Ogni riga deve rappresentare una combinazione cella × gruppo.
stopifnot(
        anyDuplicated(dat[c("CellCode", "taxGroup")]) == 0L,
        all(dat$lost_species >= 0),
        all(dat$shared_species >= 0),
        all(dat$gained_species >= 0)
)

summary_groups <- dat |>
        group_by(taxGroup) |>
        summarise(
                n_cells = n_distinct(CellCode),
                initial_occ = sum(initial),
                final_occ = sum(final),
                lost_occ = sum(lost_species),
                gained_occ = sum(gained_species),
                
                # Totali ponderati per la frazione terrestre della cella
                initial_w = sum(initial * land_weight),
                lost_w = sum(lost_species * land_weight),
                gained_w = sum(gained_species * land_weight),
                .groups = "drop"
        ) |>
        mutate(
                net_occ = gained_occ - lost_occ,
                loss_pct = if_else(
                        initial_occ > 0, 100 * lost_occ / initial_occ, NA_real_
                ),
                net_pct = if_else(
                        initial_occ > 0, 100 * net_occ / initial_occ, NA_real_
                ),
                loss_pct_weighted = if_else(
                        initial_w > 0, 100 * lost_w / initial_w, NA_real_
                ),
                net_pct_weighted = if_else(
                        initial_w > 0,
                        100 * (gained_w - lost_w) / initial_w,
                        NA_real_
                )
        )

summary_groups


# 2. Gruppi dell'analisi principale --------------------------------

groups_h2 <- c("Grasslands", "Forests", "Sclerophyllous scrub")

h2 <- dat |>
        filter(taxGroup %in% groups_h2) |>
        mutate(taxGroup = factor(taxGroup, levels = groups_h2))

# Grasslands è il riferimento.
# La perdita è definita solo se il gruppo era presente inizialmente.
d_loss <- h2 |> filter(initial > 0)

# Il bilancio è definito solo in presenza di guadagni o perdite.
# Include le colonizzazioni di celle con initial = 0.
d_balance <- h2 |> filter(changes > 0)


# 3. Modelli -------------------------------------------------------

m_loss <- glm(
        cbind(lost_species, shared_species) ~
                total_pa_land_cov * taxGroup,
        family = quasibinomial(),
        weights = land_weight,
        data = d_loss
)

m_balance <- glm(
        cbind(gained_species, lost_species) ~
                total_pa_land_cov * taxGroup,
        family = quasibinomial(),
        weights = land_weight,
        data = d_balance
)

V_loss <- vcovCL(
        m_loss, cluster = d_loss$CellCode, type = "HC1"
)

V_balance <- vcovCL(
        m_balance, cluster = d_balance$CellCode, type = "HC1"
)

# Coefficienti con errori standard corretti per cella
coeftest(m_loss, vcov. = V_loss)
coeftest(m_balance, vcov. = V_balance)

# Test complessivo delle interazioni
m_loss_add <- update(
        m_loss, . ~ total_pa_land_cov + taxGroup
)

m_balance_add <- update(
        m_balance, . ~ total_pa_land_cov + taxGroup
)

waldtest(
        m_loss_add, m_loss,
        vcov = V_loss, test = "Chisq"
)

waldtest(
        m_balance_add, m_balance,
        vcov = V_balance, test = "Chisq"
)


# 4. Predizioni con IC coerenti con la correzione per cella ----------

predict_h2 <- function(model, V, data) {
        
        # Predizioni entro il range osservato di ciascun gruppo
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

pred_loss <- predict_h2(m_loss, V_loss, d_loss)
pred_balance <- predict_h2(m_balance, V_balance, d_balance)

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

p_loss <- plot_h2(pred_loss, "Habitat loss probability")

p_balance <- plot_h2(
        pred_balance |>
                mutate(
                        estimate = 2 * estimate - 1,
                        lower = 2 * lower - 1,
                        upper = 2 * upper - 1
                ),
        "Habitat gain–loss balance"
) +
        geom_hline(
                yintercept = 0, linetype = "dashed",
                colour = "grey45"
        ) +
        scale_y_continuous(
                limits = c(-1, 1),
                breaks = seq(-1, 1, 0.5),
                labels = scales::label_number()
        )

p_loss+
p_balance


#### H2b MOSAIC ANALYSIS ####
# 1. Prepare one row per cell ------------------------------------------

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

mosaic_data <- mosaic_long |>
        pivot_wider(
                id_cols = c(CellCode, total_pa_land_cov, land_weight),
                names_from = group,
                values_from = c(lost, gained, initial),
                names_glue = "{group}_{.value}"
        )

# Missing cell–group records are retained as NA, not converted to zero.

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

fit_mosaic <- function(gain_group, loss_group, gain_label, loss_label) {
        
        d <- mosaic_data |>
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
                drop_na() |>
                filter(
                        is.finite(coverage),
                        is.finite(land_weight),
                        land_weight > 0,
                        initial_loss > 0
                )
        
        # Loss must be possible in the predictor group.
        # Initial absence of the response group is allowed.
        
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

# 4. Run all six comparisons -------------------------------------------

mosaic_results <- lapply(seq_len(nrow(comparisons)), function(i) {
        fit_mosaic(
                gain_group = comparisons$gain_group[i],
                loss_group = comparisons$loss_group[i],
                gain_label = comparisons$gain_label[i],
                loss_label = comparisons$loss_label[i]
        )
})

names(mosaic_results) <- paste(
        comparisons$gain_group,
        comparisons$loss_group,
        sep = "_gain_vs_"
)

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


d <- mosaic_results$scrub_gain_vs_grass$data

d |>
        group_by(loss) |>
        summarise(
                n_cells = n(),
                n_gains = sum(gain),
                gain_frequency = mean(gain),
                .groups = "drop"
        )

summary(
        mosaic_results$scrub_gain_vs_grass$model
)$dispersion

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
