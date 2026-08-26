library(tidyverse)
library(sf)
library(terra)
library(patchwork)
library(ggeffects)
source("Source.R")
#### load data ####

temporal_beta_group <- read.csv("processed/temporal_group_cell_by_cell.csv") |> 
        mutate(net_change = richness_2013_2018 - richness_2007_2012) |> 
        drop_na()

grid_eu <- st_read("data/EU_grid/europe_10km.shp")

PAs_cov <- read.csv("data/df_habitat_new.csv") |> 
        select(CellCode, n2k_cov, other_pa_cov, total_pa_cov)

grid_group <- grid_eu |>
        left_join(temporal_beta_group, by = c("CellCode" = "cell_id")) |> 
        left_join(PAs_cov) |> 
        filter(!is.na(total_pa_cov)) |> 
        filter(taxGroup %in% c("Forests", "Sclerophyllous scrub", "Grasslands"))

lostG_model <- glm(formula = lost_species ~ taxGroup , data = grid_group, family = poisson)               

summary(lostG_model)

ggpredict(lostG_model) |> 
        plot()

### testare se l'aumento delle foreste o delle shrub è dovuto a una diminuzione delle praterie