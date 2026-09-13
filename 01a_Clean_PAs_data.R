#### 01. Packages and functions ####
library(tidyverse)
library(sf)
library(terra)

source("SourceCodes/Source01_PAs.R")

#### 02. Input paths ####

# European 10-km grid
grid_file <- file.path(
        "data",
        "EU_grid",
        "europe_10km.shp"
)



# Protected Planet / WDPA - August 2026
# The dataset is split into three shapefiles
wdpa_files <- file.path(
        "data",
        "PAs",
        "wdpa",
        paste0("WDPA_Aug2026_Public_shp_", 0:2),
        "WDPA_Aug2026_Public_shp-polygons.shp"
)

# Official EEA Natura 2000 dataset
n2k_eea_file <- file.path(
        "data",
        "PAs",
        "n2k",
        "SHP",
        "Natura2000_end2024_epsg3035.shp"
)


#### 03. Output paths ####

out_dir <- file.path(
        "data",
        "PAs",
        "derived"
)

crop_dir <- file.path(
        out_dir,
        "cropped"
)

dir.create(
        out_dir,
        recursive = TRUE,
        showWarnings = FALSE
)

dir.create(
        crop_dir,
        recursive = TRUE,
        showWarnings = FALSE
)

# Lambert Azimuthal Equal-Area projection for Europe
crs_area <- 3035

#### 04. Check input files ####

input_check <- c(
        grid = file.exists(grid_file),
        wdpa_0 = file.exists(wdpa_files[1]),
        wdpa_1 = file.exists(wdpa_files[2]),
        wdpa_2 = file.exists(wdpa_files[3]),
        n2k_eea = file.exists(n2k_eea_file)
)

print(input_check)

if (!all(input_check)) {
        stop("One or more input files were not found.")
}


#### 05. Read and check the European grid ####

grid_10km <- st_read(
        grid_file,
        quiet = TRUE
) |>
        st_make_valid() |>
        st_transform(crs_area)

cat(
        "\nNumber of 10-km grid cells:",
        nrow(grid_10km),
        "\n"
)

cat(
        "Grid CRS:",
        st_crs(grid_10km)$input,
        "\n"
)

wdpa_fields <- c(
        "WDPAID",
        "WDPA_PID",
        "NAME",
        "ORIG_NAME",
        "DESIG",
        "DESIG_ENG",
        "DESIG_TYPE",
        "STATUS",
        "ISO3"
)

#### 06.1 Crop WDPA to study extent ####

wdpa_crop_files <- file.path(
        crop_dir,
        paste0(
                "WDPA_Aug2026_Europe_",
                0:2,
                ".gpkg"
        )
)

walk2(
        wdpa_files,
        wdpa_crop_files,
        ~ crop_vector_to_grid(
                file = .x,
                grid = grid_10km,
                output_file = .y,
                keep_fields = wdpa_fields,
                exclude_proposed = TRUE
        )
)


#### 06.2 Crop the official EEA Natura 2000 dataset ####

n2k_eea_crop_file <- file.path(
        crop_dir,
        "Natura2000_end2024_Europe.gpkg"
)

crop_vector_to_grid(
        file = n2k_eea_file,
        grid = grid_10km,
        output_file = n2k_eea_crop_file,
        keep_fields = character(0),
        exclude_proposed = FALSE
)



#### 07 Mask with the EU Boundaries

land_boundary <- st_read(
        file.path(
                "data",
                "EU_boundary",
                "EU_boundary.gpkg"
        ),
        quiet = TRUE
)


vector_crop_files <- c(
        wdpa_crop_files,
        n2k_eea_crop_file
)


walk(
        vector_crop_files,
        ~ crop_vector_to_boundary(
                input_file = .x,
                boundary = land_boundary,
                output_file = .x
        )
)


#### 08. Extract Natura 2000 sites from WDPA ####

# Natura 2000 sites may occur under different English
# designations in WDPA.
#
# The main categories are:
# - Special Protection Area
# - Site of Community Importance
# - Special Area of Conservation
# - records explicitly containing "Natura 2000"

n2k_pattern <- paste(
        c(
                "Natura[[:space:]]*2000",
                "Special Protection Area",
                "Site of Community Importance",
                "Special Area of Conservation"
        ),
        collapse = "|"
)


wdpa_n2k_files <- file.path(
        crop_dir,
        paste0(
                "WDPA_Aug2026_N2K_",
                0:2,
                ".gpkg"
        )
)

walk2(
        wdpa_crop_files,
        wdpa_n2k_files,
        extract_wdpa_n2k
)


#### 09. Inspect Natura 2000 designations found in WDPA ####

# This table is useful for checking that the Natura 2000
# identification captured the expected WDPA designations.

wdpa_n2k_designations <- map_dfr(
        wdpa_n2k_files,
        function(f) {
                
                x <- st_read(
                        f,
                        quiet = TRUE
                ) |>
                        st_drop_geometry()
                
                available <- intersect(
                        c(
                                "DESIG_ENG",
                                "DESIG_TYPE"
                        ),
                        names(x)
                )
                
                x |>
                        count(
                                across(
                                        all_of(available)
                                ),
                                name = "n"
                        )
        }
) |>
        group_by(
                across(
                        -n
                )
        ) |>
        summarise(
                n = sum(n),
                .groups = "drop"
        ) |>
        arrange(
                desc(n)
        )

print(
        wdpa_n2k_designations,
        n = Inf
)

write.csv(
        wdpa_n2k_designations,
        file.path(
                out_dir,
                "WDPA_Natura2000_designations.csv"
        ),
        row.names = FALSE
)

