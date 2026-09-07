library(tidyverse)
library(sf)
library(terra)
library(patchwork)
library(ggeffects)
source("Source.R")

#### Load data ####

temporal_beta <- read.csv("processed/temporal_beta_cell_by_cell.csv") |> 
        mutate(net_change = richness_2013_2018 - richness_2007_2012)

grid_pa <- st_read("data/PAs/derived/europe_10km_protected_area_coverage.gpkg")


grid_T <- grid_pa |>
        right_join(temporal_beta, by = c("CellCode" = "cell_id")) 


#### Figure 1 Bivariate map #####

beta_pa_map <- bivariate_map(
        data = grid_T,
        geometry = grid_T,
        x = total_pa_cov,
        y = beta_jaccard,
        x_breaks = c(10, 30, 70),
        y_breaks = c(0.25, 0.50, 0.75),
        xlab = "Protected areas (%)",
        ylab = "Beta diversity (turnover)"
)

beta_pa_map$plot

layout <- c(
        area(t = 0, l = 0, b = 12, r = 12),
        area(t = 2.5, l = 3, b = 3.5, r = 4)
)

mapleg <- beta_pa_map$map + beta_pa_map$legend + plot_layout(design = layout)

mapleg


#### H1a We expected higher temporal beta diversity where PAs coverage is low ####
#### H1b We expected less habitat lost where PAs coverage is high

##### H1a #### 

beta_model <- glm(formula = beta_jaccard ~ log(total_pa_cov  + 1), data = grid_T, family = quasibinomial(link = "logit"), weights = land_cov)

# plot_glm(data = grid_T, model = beta_model, point_size = .1)
summary(beta_model)

ggpredict(beta_model) |> 
        plot()



##### H1b ####

lost_model <- glm(formula = lost_species ~ log(total_pa_cov + 1), data = grid_T, family = poisson, weights = land_cov)

plot_glm(data = grid_T, model = lost_model, point_size = .1)
plot_small_effect(
        model = beta_model,
        data = grid_T,
        x = "n2k_cov",
        y = "beta_jaccard",
        xlab = "Protected-area cover (%)",
        ylab = expression(beta[Jaccard])
)
summary(lost_model)

ggpredict(lost_model) |> 
        plot()




gain_model <- glm(formula = gained_species ~ total_pa_cov, data = grid_T, family = poisson, weights = land_cov)

summary(gain_model)

ggpredict(gain_model) |> 
        plot()

net_model <- glm(formula = net_change ~ total_pa_cov, data = grid_T, family = gaussian, weights = land_cov)

summary(net_model)

ggpredict(net_model) |> 
        plot()





#### 



