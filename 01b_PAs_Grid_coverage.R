#### 01. Packages and functions ####

library(tidyverse)
library(sf)

source("SourceCodes/Source01_PAs.R")


#### 02. Paths ####

grid_file <- file.path(
        "data",
        "EU_grid",
        "europe_10km.shp"
)

land_file <- file.path(
        "data",
        "EU_boundary",
        "EU_boundary.gpkg"
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

tile_dir <- file.path(
        out_dir,
        "tiles"
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

dir.create(
        tile_dir,
        recursive = TRUE,
        showWarnings = FALSE
)


#### 03. Analysis parameters ####

crs_area <- 3035

grid_resolution <- 10000
raster_resolution <- 100
tile_size <- 200000

id_col <- "CellCode"

coverage_sources <- list(
        total_pa_cov = wdpa_crop_files,
        n2k_wdpa_cov = wdpa_n2k_files,
        n2k_eea_cov = n2k_eea_crop_file
)

coverage_cols <- c(
        names(coverage_sources),
        "land_cov"
)

tolerance <- 0.05


#### 04. Check input files ####

input_files <- c(
        grid_file,
        land_file,
        unlist(
                coverage_sources,
                use.names = FALSE
        )
)

missing_files <- input_files[
        !file.exists(input_files)
]

if (length(missing_files) > 0) {
        stop(
                "Missing input files:\n",
                paste(
                        missing_files,
                        collapse = "\n"
                )
        )
}


#### 05. Load grid and terrestrial boundary ####

grid_10km <- st_read(
        grid_file,
        quiet = TRUE
) |>
        st_make_valid() |>
        st_transform(crs_area)

land_boundary <- st_read(
        land_file,
        quiet = TRUE
) |>
        repair_geometries() |>
        st_transform(crs_area)

# Geometry only: country attributes are not needed
land_boundary <- land_boundary[, 0, drop = FALSE]


#### 06. Check grid ####

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

area_error <- max(
        abs(
                as.numeric(
                        st_area(grid_10km)
                ) -
                        grid_resolution^2
        ) /
                grid_resolution^2,
        na.rm = TRUE
)

if (area_error > 0.001) {
        warning(
                "Some grid polygons differ substantially ",
                "from 100 km2."
        )
}

if (grid_resolution %% raster_resolution != 0) {
        stop(
                "raster_resolution must divide exactly ",
                "into grid_resolution."
        )
}

if (tile_size %% grid_resolution != 0) {
        stop(
                "tile_size must be an exact multiple ",
                "of grid_resolution."
        )
}


#### 07. Create processing tiles ####

grid_bbox <- st_bbox(
        grid_10km
)

xy <- st_coordinates(
        st_centroid(
                grid_10km
        )
)

grid_10km <- grid_10km |>
        mutate(
                tile_x = floor(
                        (
                                xy[, 1] -
                                        as.numeric(
                                                grid_bbox["xmin"]
                                        )
                        ) /
                                tile_size
                ),
                
                tile_y = floor(
                        (
                                xy[, 2] -
                                        as.numeric(
                                                grid_bbox["ymin"]
                                        )
                        ) /
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

message(
        "Processing tiles: ",
        length(tiles)
)


#### 08. Coverage validation ####

validate_coverage <- function(
                x,
                expected_rows) {
        
        values <- as.matrix(
                x[
                        coverage_cols
                ]
        )
        
        c(
                correct_rows =
                        nrow(x) == expected_rows,
                
                unique_cells =
                        n_distinct(
                                x[[id_col]]
                        ) == expected_rows,
                
                no_missing =
                        !anyNA(values),
                
                valid_range =
                        all(
                                values >= 0 &
                                        values <= 100,
                                na.rm = TRUE
                        ),
                
                wdpa_internal_consistency =
                        all(
                                x$n2k_wdpa_cov <=
                                        x$total_pa_cov +
                                        tolerance,
                                na.rm = TRUE
                        ),
                
                terrestrial_consistency =
                        all(
                                x$total_pa_cov <=
                                        x$land_cov +
                                        tolerance &
                                        
                                        x$n2k_wdpa_cov <=
                                        x$land_cov +
                                        tolerance &
                                        
                                        x$n2k_eea_cov <=
                                        x$land_cov +
                                        tolerance,
                                na.rm = TRUE
                        )
        )
}


#### 09. Test representative tile ####

tile_test <- grid_10km |>
        st_drop_geometry() |>
        count(
                tile_id,
                sort = TRUE
        ) |>
        slice(1) |>
        pull(tile_id)

message(
        "Test tile: ",
        tile_test
)

test_coverage <- process_coverage_tile(
        tile_name = tile_test,
        grid = grid_10km,
        coverage_sources = coverage_sources,
        boundary = land_boundary,
        output_name = "TEST_PA_coverage",
        tile_dir = tile_dir,
        id = id_col,
        grid_res = grid_resolution,
        raster_res = raster_resolution,
        overwrite = TRUE
)

print(
        summary(
                test_coverage[
                        coverage_cols
                ]
        )
)

test_checks <- validate_coverage(
        test_coverage,
        sum(
                grid_10km$tile_id ==
                        tile_test
        )
)

print(
        test_checks
)

if (!all(test_checks)) {
        stop(
                "\nTest tile validation failed:\n",
                paste(
                        names(test_checks)[
                                !test_checks
                        ],
                        collapse = "\n"
                ),
                "\n\nFull analysis aborted."
        )
}

message(
        "Test successfully validated. ",
        "Starting full analysis..."
)


#### 10. Calculate European coverage ####

coverage <- map_dfr(
        tiles,
        ~ process_coverage_tile(
                tile_name = .x,
                grid = grid_10km,
                coverage_sources = coverage_sources,
                boundary = land_boundary,
                output_name = "PA_coverage",
                tile_dir = tile_dir,
                id = id_col,
                grid_res = grid_resolution,
                raster_res = raster_resolution
        )
)


#### 11. Validate complete output ####

coverage_checks <- validate_coverage(
        coverage,
        nrow(grid_10km)
)

print(
        coverage_checks
)

if (!all(coverage_checks)) {
        stop(
                "\nFull coverage validation failed:\n",
                paste(
                        names(coverage_checks)[
                                !coverage_checks
                        ],
                        collapse = "\n"
                )
        )
}


#### 12. Join coverage to grid ####

grid_pa <- grid_10km |>
        select(
                -tile_x,
                -tile_y,
                -tile_id
        ) |>
        left_join(
                coverage,
                by = id_col
        ) |>
        mutate(
                # PA coverage relative to terrestrial surface
                total_pa_land_cov =
                        if_else(
                                land_cov > 0,
                                pmin(
                                        100,
                                        100 *
                                                total_pa_cov /
                                                land_cov
                                ),
                                NA_real_
                        ),
                
                n2k_wdpa_land_cov =
                        if_else(
                                land_cov > 0,
                                pmin(
                                        100,
                                        100 *
                                                n2k_wdpa_cov /
                                                land_cov
                                ),
                                NA_real_
                        ),
                
                n2k_eea_land_cov =
                        if_else(
                                land_cov > 0,
                                pmin(
                                        100,
                                        100 *
                                                n2k_eea_cov /
                                                land_cov
                                ),
                                NA_real_
                        ),
                
                n2k_difference =
                        n2k_wdpa_cov -
                        n2k_eea_cov,
                
                n2k_absolute_difference =
                        abs(
                                n2k_difference
                        )
        )


#### 13. Final quality control ####

print(
        summary(
                grid_pa |>
                        st_drop_geometry() |>
                        select(
                                all_of(
                                        coverage_cols
                                ),
                                total_pa_land_cov,
                                n2k_wdpa_land_cov,
                                n2k_eea_land_cov,
                                n2k_difference
                        )
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


#### 14. Save outputs ####

output_base <- file.path(
        out_dir,
        "europe_10km_protected_area_coverage"
)

st_write(
        grid_pa,
        paste0(
                output_base,
                ".gpkg"
        ),
        delete_dsn = TRUE,
        quiet = TRUE
)

write_csv(
        grid_pa |>
                st_drop_geometry(),
        paste0(
                output_base,
                ".csv"
        )
)

write_csv(
        n2k_comparison_summary,
        file.path(
                out_dir,
                "N2K_WDPA_vs_EEA_summary.csv"
        )
)


#### 15. Filter current PA coverage to EU27 before Croatia accession ####

coverage_file <- file.path(
        "data",
        "PAs",
        "derived",
        "europe_10km_protected_area_coverage.gpkg"
)

output_file <- file.path(
        "data",
        "PAs",
        "derived",
        "europe_10km_protected_area_coverage_2012EU27.gpkg"
)


# EU27 Member States during the 2007-2012 reporting period

eu27_2012 <- c(
        "AUT", "BEL", "BGR", "CYP", "CZE", "DNK", "EST",
        "FIN", "FRA", "DEU", "GRC", "HUN", "IRL", "ITA",
        "LVA", "LTU", "LUX", "MLT", "NLD", "POL", "PRT",
        "ROU", "SVK", "SVN", "ESP", "SWE", "GBR"
)

# Read data #
grid_pa <- st_read(
        coverage_file,
        quiet = TRUE
)

land_boundary$ISO_A3

land_boundary_EU27 <- land_boundary |>
        filter(
                ISO_A3 %in% eu27_2012
        ) |>
        st_make_valid() |>
        st_transform(
                st_crs(grid_pa)
        ) |>
        st_union()


# Retain grid cells overlapping EU27 terrestrial territory

grid_pa_eu27 <- st_filter(
        grid_pa,
        land_boundary_EU27,
        .predicate = st_intersects
)


#### 16. Save EU27 before Croatia accession ####

st_write(
        grid_pa_eu27,
        output_file,
        delete_dsn = TRUE,
        quiet = TRUE
)

cat(
        "Original cells:",
        nrow(grid_pa),
        "\nEU27-2012 cells:",
        nrow(grid_pa_eu27),
        "\n"
)
