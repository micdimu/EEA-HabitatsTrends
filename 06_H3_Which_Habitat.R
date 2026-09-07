library(tidyverse)
library(sf)
library(terra)
library(patchwork)
library(ggeffects)
source("Source.R")

#### load data ####


habitat_trend <- read.csv("processed/temporal_each_habitat.csv") 

habitat_trend |> 
        filter(netto_pct_2012 < 250) |> 
        pivot_longer(contains("pct")) |> 
        ggplot(mapping = aes(x = log(Ncelle_2012), y = value)) +
        geom_point() +
        theme_bw() +
        facet_wrap(~name)
        
        
habitat_trend |> 
        filter(netto_pct_2012 < 250) |> 
        pivot_longer(c("netto", "loss", "gain")) |> 
        ggplot(mapping = aes(x = log(Ncelle_2012), y = value)) +
        geom_point() +
        theme_bw() +
        facet_wrap(~name)




habitat_trend
glm(
netto ~ log(Ncelle_2012), 
  data = habitat_trend |> 
          filter(netto_pct_2012 < 250) 
) |> 
        ggeffect() |> 
        plot()+
        geom_point(data = habitat_trend |> 
                     filter(netto_pct_2012 < 250), 
                   mapping = aes(x = Ncelle_2012, y = netto)) +
        theme_bw() +
        scale_x_log10()+
        labs(x = "Habitat area in 2012 (ha)", y = "Loss")
