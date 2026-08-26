#### 01. Packages and functions ####

library(tidyverse)
library(sf)
library(terra)

source("Source_PAs.R")


#### 02. Input paths ####

grid_file <- file.path(
        "data",
        "EU_grid",
        "europe_10km.shp"
)

out_dir <- file.path(
        "data",
        "PAs",
        "derived"
)

crop_dir <- file.path(
        out_dir,
        "cropped"
)

wdpa_crop_files <- file.path(
        crop_dir,
        paste0(
                "WDPA_Aug2026_Europe_",
                0:2,
                ".gpkg"
        )
)

wdpa_n2k_files <- file.path(
        crop_dir,
        paste0(
                "WDPA_Aug2026_N2K_",
                0:2,
                ".gpkg"
        )
)

n2k_eea_crop_file <- file.path(
        crop_dir,
        "Natura2000_end2024_Europe.gpkg"
)


# check inputs

input_check <- c(
        grid = file.exists(grid_file),
        wdpa_0 = file.exists(wdpa_crop_files[1]),
        wdpa_1 = file.exists(wdpa_crop_files[2]),
        wdpa_2 = file.exists(wdpa_crop_files[3]),
        wdpa_n2k_0 = file.exists(wdpa_n2k_files[1]),
        wdpa_n2k_1 = file.exists(wdpa_n2k_files[2]),
        wdpa_n2k_2 = file.exists(wdpa_n2k_files[3]),
        n2k_eea = file.exists(n2k_eea_crop_file)
)

if (!all(input_check)) {
        stop(
                "One or more derived PA files are missing. ",
                "Run 01_prepare_PA_vectors.R first."
        )
}


#### 03. Output paths ####

tile_dir <- file.path(
        out_dir,
        "tiles"
)

dir.create(
        tile_dir,
        recursive = TRUE,
        showWarnings = FALSE
)

#### 04. Analysis parameters ####

crs_area <- 3035

grid_resolution <- 10000
raster_resolution <- 100
tile_size <- 200000

id_col <- "CellCode"

#### 05. Load grid and check geometry ####

# The raster approach below assumes that the vector layer contains
# complete 10 x 10 km grid cells.
#
# Coastal cells may later be corrected using a terrestrial land mask,
# but the underlying grid geometries should still be full squares.

grid_10km <- st_read(
        grid_file,
        quiet = TRUE) |>
        st_make_valid() |>
        st_transform(crs_area)

#### Check grid identifier ####

if (!id_col %in% names(grid_10km)) {
        stop(
                "Grid ID column '",
                id_col,
                "' was not found."
        )
}

if (anyDuplicated(grid_10km[[id_col]])) {
        stop(
                "Grid-cell IDs are not unique."
        )
}


grid_areas <- as.numeric(
        st_area(grid_10km))

expected_area <- grid_resolution^2

relative_area_error <- abs(
        grid_areas - expected_area
) / expected_area

cat(
        "\nMedian grid-cell area:",
        median(grid_areas, na.rm = TRUE) / 1e6,
        "km2\n"
)

cat(
        "Maximum relative deviation from 100 km2:",
        max(relative_area_error, na.rm = TRUE),
        "\n"
)

if (max(relative_area_error, na.rm = TRUE) > 0.001) {
        warning(
                "Some grid polygons differ substantially from 100 km2. ",
                "Check whether the grid contains clipped cells before proceeding."
        )
}

# The temporary raster resolution must divide exactly into 10 km

if (grid_resolution %% raster_resolution != 0) {
        stop(
                "The raster resolution must divide exactly into the ",
                "10-km grid resolution."
        )
}

aggregation_factor <- as.integer(
        grid_resolution / raster_resolution
)

cat(
        "Raster cells per grid-cell side:",
        aggregation_factor,
        "\n"
)

cat(
        "Raster cells per 10-km grid cell:",
        aggregation_factor^2,
        "\n"
)

#### 06. Create spatial processing tiles ####

# Tile boundaries are defined relative to the actual origin of
# the 10-km grid, rather than to an arbitrary EPSG:3035 origin.
#
# Because tile_size is an exact multiple of 10 km, tile boundaries
# remain aligned with the original grid.

grid_bbox <- st_bbox(
        grid_10km
)

grid_origin_x <- as.numeric(
        grid_bbox["xmin"]
)

grid_origin_y <- as.numeric(
        grid_bbox["ymin"]
)

if (tile_size %% grid_resolution != 0) {
        stop(
                "tile_size must be an exact multiple of ",
                "the 10-km grid resolution."
        )
}

grid_centroids <- st_centroid(
        grid_10km
)

xy <- st_coordinates(
        grid_centroids
)

grid_10km <- grid_10km |>
        mutate(
                .cx = xy[, 1],
                .cy = xy[, 2],
                
                tile_x = floor(
                        (.cx - grid_origin_x) /
                                tile_size
                ),
                
                tile_y = floor(
                        (.cy - grid_origin_y) /
                                tile_size
                ),
                
                tile_id = paste(
                        tile_x,
                        tile_y,
                        sep = "_"
                )
        )

tiles <- unique(
        grid_10km$tile_id
)

cat(
        "\nNumber of processing tiles:",
        length(tiles),
        "\n"
)


#### 07. Select a representative tile for testing ####

# Use the tile containing the largest number of grid cells,
# rather than simply taking the first tile in the dataset.

tile_test <- grid_10km |>
        st_drop_geometry() |>
        count(
                tile_id,
                sort = TRUE
        ) |>
        slice(1) |>
        pull(
                tile_id
        )

cat(
        "\nTest tile:",
        tile_test,
        "\n"
)


#### 08. Test total WDPA coverage ####

test_wdpa_total <- process_coverage_tile(
        tile_name = tile_test,
        grid = grid_10km,
        vector_files = wdpa_crop_files,
        output_name = "TEST_WDPA_total",
        coverage_name = "total_pa_cov",
        id = id_col,
        raster_res = raster_resolution,
        overwrite = TRUE
)


#### 09. Test Natura 2000 coverage from WDPA ####

test_wdpa_n2k <- process_coverage_tile(
        tile_name = tile_test,
        grid = grid_10km,
        vector_files = wdpa_n2k_files,
        output_name = "TEST_WDPA_N2K",
        coverage_name = "n2k_wdpa_cov",
        id = id_col,
        raster_res = raster_resolution,
        overwrite = TRUE
)


#### 10. Test official EEA Natura 2000 coverage ####

test_eea_n2k <- process_coverage_tile(
        tile_name = tile_test,
        grid = grid_10km,
        vector_files = n2k_eea_crop_file,
        output_name = "TEST_EEA_N2K",
        coverage_name = "n2k_eea_cov",
        id = id_col,
        raster_res = raster_resolution,
        overwrite = TRUE
)

print(
        summary(
                test_eea_n2k$n2k_eea_cov
        )
)


#### 11. Compare test-tile results ####

test_comparison <- test_wdpa_total |>
        left_join(
                test_wdpa_n2k,
                by = id_col
        ) |>
        left_join(
                test_eea_n2k,
                by = id_col
        ) |>
        mutate(
                n2k_difference =
                        n2k_wdpa_cov -
                        n2k_eea_cov,
                
                n2k_absolute_difference =
                        abs(
                                n2k_difference
                        )
        )

print(
        summary(
                test_comparison[
                        c(
                                "total_pa_cov",
                                "n2k_wdpa_cov",
                                "n2k_eea_cov",
                                "n2k_difference"
                        )
                ]
        )
)


#### 12. Validate test tile ####

n_test_cells <- sum(
        grid_10km$tile_id == tile_test
)

test_checks <- c(
        
        correct_rows =
                nrow(test_comparison) == n_test_cells,
        
        unique_cells =
                n_distinct(test_comparison[[id_col]]) == n_test_cells,
        
        no_missing =
                !anyNA(
                        test_comparison[
                                c(
                                        "total_pa_cov",
                                        "n2k_wdpa_cov",
                                        "n2k_eea_cov"
                                )
                        ]
                ),
        
        wdpa_range =
                all(
                        test_comparison$total_pa_cov >= 0 &
                                test_comparison$total_pa_cov <= 100
                ),
        
        wdpa_n2k_range =
                all(
                        test_comparison$n2k_wdpa_cov >= 0 &
                                test_comparison$n2k_wdpa_cov <= 100
                ),
        
        eea_n2k_range =
                all(
                        test_comparison$n2k_eea_cov >= 0 &
                                test_comparison$n2k_eea_cov <= 100
                ),
        
        wdpa_internal_consistency =
                all(
                        test_comparison$n2k_wdpa_cov <=
                                test_comparison$total_pa_cov + 0.05
                )
)

print(test_checks)

if (!all(test_checks)) {
        
        failed_checks <- names(
                test_checks
        )[!test_checks]
        
        stop(
                "\nTest tile validation failed:\n",
                paste(
                        failed_checks,
                        collapse = "\n"
                ),
                "\n\nFull analysis aborted."
        )
}

message(
        "\nTest tile successfully validated.",
        "\nStarting full European analysis..."
)

#### 13. Run the full European analysis ####

# Total WDPA protected-area coverage
wdpa_total_coverage <- map_dfr(
        tiles,
        ~ process_coverage_tile(
                tile_name = .x,
                grid = grid_10km,
                vector_files = wdpa_crop_files,
                output_name = "WDPA_total",
                coverage_name = "total_pa_cov",
                id = id_col,
                raster_res = raster_resolution
        )
)

# Natura 2000 coverage identified within WDPA
wdpa_n2k_coverage <- map_dfr(
        tiles,
        ~ process_coverage_tile(
                tile_name = .x,
                grid = grid_10km,
                vector_files = wdpa_n2k_files,
                output_name = "WDPA_N2K",
                coverage_name = "n2k_wdpa_cov",
                id = id_col,
                raster_res = raster_resolution
        )
)

# Official EEA Natura 2000 coverage
eea_n2k_coverage <- map_dfr(
        tiles,
        ~ process_coverage_tile(
                tile_name = .x,
                grid = grid_10km,
                vector_files = n2k_eea_crop_file,
                output_name = "EEA_N2K",
                coverage_name = "n2k_eea_cov",
                id = id_col,
                raster_res = raster_resolution
        )
)

#### 14. Check complete coverage outputs ####

coverage_checks <- tibble(
        dataset = c(
                "WDPA total",
                "WDPA Natura 2000",
                "EEA Natura 2000"
        ),
        n_rows = c(
                nrow(wdpa_total_coverage),
                nrow(wdpa_n2k_coverage),
                nrow(eea_n2k_coverage)
        ),
        n_unique_cells = c(
                n_distinct(wdpa_total_coverage[[id_col]]),
                n_distinct(wdpa_n2k_coverage[[id_col]]),
                n_distinct(eea_n2k_coverage[[id_col]])
        ),
        n_missing = c(
                sum(is.na(wdpa_total_coverage$total_pa_cov)),
                sum(is.na(wdpa_n2k_coverage$n2k_wdpa_cov)),
                sum(is.na(eea_n2k_coverage$n2k_eea_cov))
        )
)

print(coverage_checks)

if (
        any(coverage_checks$n_rows != nrow(grid_10km)) ||
        any(coverage_checks$n_unique_cells != nrow(grid_10km)) ||
        any(coverage_checks$n_missing > 0)
) {
        stop(
                "Coverage outputs are incomplete or contain duplicated/missing grid cells."
        )
}


#### 15. Join coverage metrics to the original grid ####

grid_pa <- grid_10km |>
        select(
                -starts_with("."),
                -tile_x,
                -tile_y,
                -tile_id
        ) |>
        left_join(
                wdpa_total_coverage,
                by = id_col
        ) |>
        left_join(
                wdpa_n2k_coverage,
                by = id_col
        ) |>
        left_join(
                eea_n2k_coverage,
                by = id_col
        ) |>
        mutate(
                n2k_difference =
                        n2k_wdpa_cov -
                        n2k_eea_cov,
                
                n2k_absolute_difference =
                        abs(
                                n2k_difference
                        )
        )


#### 16. Final quality control ####

cat(
        "\nTotal WDPA coverage:\n"
)

print(
        summary(
                grid_pa$total_pa_cov
        )
)

cat(
        "\nWDPA Natura 2000 coverage:\n"
)

print(
        summary(
                grid_pa$n2k_wdpa_cov
        )
)

cat(
        "\nEEA Natura 2000 coverage:\n"
)

print(
        summary(
                grid_pa$n2k_eea_cov
        )
)

cat(
        "\nDifference between WDPA and EEA Natura 2000:\n"
)

print(
        summary(
                grid_pa$n2k_difference
        )
)

n2k_comparison_summary <- grid_pa |>
        st_drop_geometry() |>
        summarise(
                mean_difference =
                        mean(
                                n2k_difference,
                                na.rm = TRUE
                        ),
                
                mean_absolute_difference =
                        mean(
                                n2k_absolute_difference,
                                na.rm = TRUE
                        ),
                
                median_absolute_difference =
                        median(
                                n2k_absolute_difference,
                                na.rm = TRUE
                        ),
                
                max_absolute_difference =
                        max(
                                n2k_absolute_difference,
                                na.rm = TRUE
                        ),
                
                correlation =
                        cor(
                                n2k_wdpa_cov,
                                n2k_eea_cov,
                                use = "complete.obs"
                        )
        )

print(
        n2k_comparison_summary
)


#### 17. Save final outputs ####

st_write(
        grid_pa,
        file.path(
                out_dir,
                "europe_10km_protected_area_coverage.gpkg"
        ),
        delete_dsn = TRUE,
        quiet = TRUE
)

write.csv(
        grid_pa |>
                st_drop_geometry() |>
                select(
                        all_of(id_col),
                        total_pa_cov,
                        n2k_wdpa_cov,
                        n2k_eea_cov,
                        n2k_difference,
                        n2k_absolute_difference
                ),
        file.path(
                out_dir,
                "europe_10km_protected_area_coverage.csv"
        ),
        row.names = FALSE
)

write.csv(
        n2k_comparison_summary,
        file.path(
                out_dir,
                "N2K_WDPA_vs_EEA_summary.csv"
        ),
        row.names = FALSE
)


