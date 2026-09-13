library(tidyverse)
library(tidylog)
library(vegan)
library(terra)
library(sf)
library(furrr)
library(progressr)

source("SourceCodes/Source02_DataPreparation.R")

#### load and check data ####
st_layers("data/report_eea/eea_2007-2012/Art17-2007-2012_GPKG/ART17_2007_2012_public_r01.gpkg")
st_layers("data/report_eea/eea_2007-2012/Art17-2007-2012_GPKG/art17_0712_public_r02.gpkg")


report2012_r01 <- st_read("data/report_eea/eea_2007-2012/Art17-2007-2012_GPKG/ART17_2007_2012_public_r01.gpkg",
                          layer = "Art17_habitats_distribution_2007_2012_EU")

report2012_r02 <- st_read("data/report_eea/eea_2007-2012/Art17-2007-2012_GPKG/art17_0712_public_r02.gpkg",
                          layer = "Art17_habitats_distribution_2007_2012_EU")


grid_eu <- st_read("data/EU_grid/europe_10km.shp")

# same CRS  

st_crs(grid_eu) == st_crs(report2012_r01)
st_crs(grid_eu) == st_crs(report2012_r02)

#### intersect grid's centroids ####

cent <- grid_eu |>   
        st_centroid()

r02_cell <- st_join(
        cent,
        report2012_r02,
        join = st_intersects
) |>
        st_drop_geometry()  |>
        transmute(CellCode,
                  habitat = habitatcodeEU,   # cambia nome colonna
                  period = "r02") |>
        drop_na()


r01_cell <- st_join(
        cent,
        report2012_r01,
        join = st_intersects
) |>
        st_drop_geometry()  |>
        transmute(CellCode,
                  habitat = habitatcode,   # cambia nome colonna
                  period = "r01") |>
        drop_na()


dat <- bind_rows(r01_cell, r02_cell) |> 
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
        separate(sample_id, into = c("cell_id", "period"), sep = "_")

comm_mat <- comm |> 
        select(-sample_id)  |> 
        as.data.frame()

rownames(comm_mat) <- comm$sample_id


#### run the dissimilarity - NOT TO RUN - Load the file below #####

# plan(multisession, workers = parallel::detectCores() - 3)
# handlers(global = TRUE)
# handlers("progress")
# 
# cells <- meta %>%
#         count(cell_id) |> 
#         filter(n == 2) |> 
#         pull(cell_id)
# 
# r_beta <- with_progress({
#         p <- progressor(along = cells)
#         future_map_dfr(cells, temporal_bgl, comm_mat = comm_mat, meta = meta)
# })
# if(!dir.exists("processed")){
#         dir.create("processed")
# }

# write.csv(r_beta,
#    "processed/rcheck.csv",
#    row.names = FALSE)

#### import the dissimilarity ######

r_beta <- read.csv("processed/rcheck.csv")


non_equal <- r_beta |> 
        filter(lost_species > 0)

non_equal

diss_cells <- r02_cell |> 
        filter(CellCode %in% non_equal$cell_id) |> 
        bind_rows(r01_cell |>
                           filter(CellCode %in% non_equal$cell_id)) |> 
        group_by(habitat, CellCode) |> 
        summarise(cc = n()) |> 
        filter(cc == 1)


#### check with last report ####
st_layers("data/report_eea/eea_2013-2018/Art17-2013-2018_GPKG/art17_2013_2018_public.gpkg")

report2018 <- st_read("data/report_eea/eea_2013-2018/Art17-2013-2018_GPKG/art17_2013_2018_public.gpkg",
                      layer = "Art17_habitats_distribution_2013_2018_EU")

st_crs(grid_eu) == st_crs(report2018)

#### intersect grid's centroids ####
r2018 <- st_join(
        cent,
        report2018,
        join = st_intersects
) |>
        st_drop_geometry()  |>
        transmute(CellCode,
                  habitat = habitatcode,   # cambia nome colonna
                  period = "r2018") |>
        drop_na() 

r2018_f <- r2018 |> 
        filter(CellCode %in% non_equal$cell_id)

diss_check_in2018 <- diss_cells |> 
        left_join(r2018_f)

diss_check_in2018 |> 
        count(habitat)


r02_cell |> 
        bind_rows(r01_cell) |> 
        bind_rows(r2018) |> 
        filter(habitat %in% unique(diss_check_in2018$habitat)) |> 
        count(habitat, period) |> 
        pivot_wider(names_from = period,  values_from = n)
