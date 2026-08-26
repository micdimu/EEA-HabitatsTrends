library(tidyverse)
library(sf)
library(terra)
library(patchwork)
library(ggeffects)
source("Source.R")

#### load data ####


habitat_trend <- read.csv("processed/temporal_each_habitat.csv") 

habitat_trends |> 
        filter(netto_pct_2012 < 250) |> 
        ggplot(mapping = aes(x = log(Ncelle_2012), y = loss)) +
  geom_smooth(method = "lm") +
  geom_point() +
  theme_bw() 

habitat_trends
glm(
netto ~ log(Ncelle_2012), 
  data = habitat_trends |> 
          filter(netto_pct_2012 < 250) 
) |> 
        ggeffect() |> 
        plot()+
        geom_point(data = habitat_trends |> 
                     filter(netto_pct_2012 < 250), 
                   mapping = aes(x = Ncelle_2012, y = netto)) +
        theme_bw() +
        scale_x_log10()+
        labs(x = "Habitat area in 2012 (ha)", y = "Loss")
