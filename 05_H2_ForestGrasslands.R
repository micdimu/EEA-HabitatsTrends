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
#        filter(taxGroup %in% c("Forests", "Sclerophyllous scrub", "Grasslands")) |> 
        mutate(taxGroup = as.factor(taxGroup))


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


