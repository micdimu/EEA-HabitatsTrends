library(tidyverse)
library(tidylog)
library(vegan)
library(terra)
library(sf)
library(furrr)
library(progressr)


#### load and check data ####

st_layers("data/report_eea/eea_2007-2012/Art17-2007-2012_GPKG/art17_0712_public_r02.gpkg")
st_layers("data/report_eea/eea_2013-2018/Art17-2013-2018_GPKG/art17_2013_2018_public.gpkg")

report2012 <- st_read("data/report_eea/eea_2007-2012/Art17-2007-2012_GPKG/art17_0712_public_r02.gpkg",
              layer = "Art17_habitats_distribution_2007_2012_EU")

report2018 <- st_read("data/report_eea/eea_2013-2018/Art17-2013-2018_GPKG/art17_2013_2018_public.gpkg",
              layer = "Art17_habitats_distribution_2013_2018_EU")

grid_eu <- st_read("data/EU_grid/europe_10km.shp")

# same CRS  

st_crs(grid_eu) == st_crs(report2012)

st_crs(grid_eu) == st_crs(report2018)

#### intersect grid's centroids ####

cent <- grid_eu |>   
        st_centroid()

hab2018_cell <- st_join(
        cent,
        report2018,
        join = st_intersects
        ) |>
        st_drop_geometry()  |>
        transmute(CellCode,
                  habitat = habitatcode,   # cambia nome colonna
                  period = "2013_2018") |>
        drop_na()


hab2012_cell <- st_join(
        cent,
        report2012,
        join = st_intersects
        ) |>
        st_drop_geometry()  |>
        transmute(CellCode,
                  habitat = habitatcodeEU,   # cambia nome colonna
                  period = "2007_2012") |>
        drop_na()

dat <- bind_rows(hab2018_cell, hab2012_cell) |> 
        distinct(CellCode, period, habitat)  |> 
        filter(!is.na(CellCode), !is.na(habitat))

##### check data before run the analysis 

dat |> 
        count(period)

dat |> 
        summarise(n_cells = n_distinct(CellCode),
                  n_species = n_distinct(habitat))

dat |> 
        count(CellCode, period) |> 
        count(period)

#### build the community matrix #####

comm <- dat |> 
        mutate(presence = 1) |> 
        unite(sample_id, CellCode, period, remove = FALSE) |> 
        select(sample_id, habitat, presence) |> 
        pivot_wider(names_from = habitat,
                    values_from = presence,
                    values_fill = 0)

meta <- comm |> 
        select(sample_id)  |> 
        separate(sample_id, into = c("cell_id", "period"), sep = "_(?=20)")

comm_mat <- comm |> 
        select(-sample_id)  |> 
        as.data.frame()

rownames(comm_mat) <- comm$sample_id


#### run the dissimilarity #####

plan(multisession, workers = parallel::detectCores() - 3)
handlers(global = TRUE)
handlers("progress")

cells <- meta %>%
        count(cell_id) |> 
        filter(n == 2) |> 
        pull(cell_id)

# compute temporal trend - NOT TO RUN - Time consuming
# 
# temporal_beta <- with_progress({
#         p <- progressor(along = cells)
#         future_map_dfr(cells, temporal_bgl, comm_mat = comm_mat, meta = meta)
# })

# # if(!dir.exists("processed")){
#         dir.create("processed")
# }
# 
# write.csv(temporal_beta,
#           "processed/temporal_beta_cell_by_cell.csv",
#           row.names = FALSE)

# temporal_beta <- read.csv("processed/temporal_beta_cell_by_cell.csv")

# temporal_beta$beta_jaccard |> hist()

#### temporal trends per each habitat group ####

cc <- report2012 |> 
        as.data.frame() |> 
        select(habitat = habitatcodeEU, taxGroup) |> 
        distinct() |> 
        drop_na()

unique(report2018$habitatcode)[!(unique(report2018$habitatcode) %in% cc$habitat)]
# new habitats defined between the two reports 32A0 e 6540 - not included in the analysis


temporal_beta_group <- with_progress({
        p <- progressor(along = cells)
        future_map_dfr(cells[1:10], temporal_bgl_group, comm_mat = comm_mat, meta = meta, group = cc)
})

write.csv(temporal_beta_group,
          "processed/temporal_group_cell_by_cell.csv",
          row.names = FALSE)

#### compute temporal trends for each habitat ####

habitat_trend <- habitat_temporal_trend(
        comm_mat = comm_mat,
        meta = meta
)

write.csv(habitat_trend,
          "processed/temporal_each_habitat.csv",
          row.names = FALSE)
