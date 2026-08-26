#### 01. Packages ####

library(sf)
library(terra)
library(dplyr)
library(purrr)


#' Repair invalid spatial geometries
#'
#' Attempts to repair invalid geometries in an sf object using
#' `sf::st_make_valid()`. A global repair is attempted first. If this
#' operation fails, invalid features are repaired individually.
#'
#' Features that remain invalid after repair, as well as empty
#' geometries, are removed from the returned object.
#'
#' @param x An `sf` object.
#'
#' @return An `sf` object containing only valid, non-empty geometries.
#'
#' @details
#' The function reports the number of invalid geometries before and
#' after repair. Individual repair is used only as a fallback when
#' `st_make_valid()` fails on the complete object.
#'

repair_geometries <- function(x) {
        
        valid_before <- sf::st_is_valid(x)
        
        n_invalid_before <- sum(
                is.na(valid_before) |
                        !valid_before
        )
        
        message(
                "Invalid geometries before repair: ",
                n_invalid_before
        )
        
        if (n_invalid_before == 0) {
                return(x)
        }
        
        
        #### First attempt: repair the complete object ####
        
        repaired <- tryCatch(
                sf::st_make_valid(x),
                error = function(e) NULL
        )
        
        if (!is.null(repaired)) {
                
                valid_after <- sf::st_is_valid(
                        repaired
                )
                
                message(
                        "Invalid geometries after repair: ",
                        sum(
                                is.na(valid_after) |
                                        !valid_after
                        )
                )
                
                repaired <- repaired[
                        !is.na(valid_after) &
                                valid_after,
                ]
                
                repaired <- repaired[
                        !sf::st_is_empty(repaired),
                ]
                
                return(repaired)
        }
        
        
        #### Fallback: repair invalid features individually ####
        
        message(
                "Global st_make_valid() failed. ",
                "Trying invalid features individually..."
        )
        
        geom <- sf::st_geometry(x)
        
        invalid_idx <- which(
                is.na(valid_before) |
                        !valid_before
        )
        
        repair_success <- logical(
                length(invalid_idx)
        )
        
        for (i in seq_along(invalid_idx)) {
                
                j <- invalid_idx[i]
                
                g <- tryCatch(
                        sf::st_make_valid(
                                geom[j]
                        ),
                        error = function(e) NULL
                )
                
                if (is.null(g)) {
                        next
                }
                
                is_valid_g <- sf::st_is_valid(g)
                
                if (
                        length(is_valid_g) == 1 &&
                        !is.na(is_valid_g) &&
                        is_valid_g
                ) {
                        
                        geom[j] <- g
                        repair_success[i] <- TRUE
                }
        }
        
        sf::st_geometry(x) <- geom
        
        
        #### Final validity check ####
        
        valid_after <- sf::st_is_valid(x)
        
        message(
                "Successfully repaired individually: ",
                sum(repair_success)
        )
        
        message(
                "Geometries discarded: ",
                sum(
                        is.na(valid_after) |
                                !valid_after
                )
        )
        
        x <- x[
                !is.na(valid_after) &
                        valid_after,
        ]
        
        x <- x[
                !sf::st_is_empty(x),
        ]
        
        x
}


#' Safely save an R object as an RDS file
#'
#' Saves an object to a temporary file and renames it only after the
#' write operation has completed successfully. This prevents incomplete
#' or corrupted cache files from being mistaken for valid results after
#' an interrupted analysis.
#'
#' @param x R object to save.
#' @param file Character string giving the destination `.rds` file.
#'
#' @return The file path invisibly.
#'
save_rds_safe <- function(
                x,
                file) {
        
        temp_file <- paste0(
                file,
                ".tmp"
        )
        
        if (file.exists(temp_file)) {
                unlink(temp_file)
        }
        
        saveRDS(
                x,
                temp_file
        )
        
        if (file.exists(file)) {
                unlink(file)
        }
        
        success <- file.rename(
                temp_file,
                file
        )
        
        if (!success) {
                stop(
                        "Could not save RDS file: ",
                        file
                )
        }
        
        invisible(file)
}


#' Crop and clean a polygon dataset to the analysis grid
#'
#' Reads only features overlapping the spatial extent of a supplied
#' grid, repairs invalid geometries, retains polygonal components,
#' optionally removes proposed protected areas, transforms the data
#' to the grid CRS, and crops the result to the grid extent.
#'
#' @param file Character string giving the input vector file.
#' @param grid An `sf` object defining the spatial extent and target CRS.
#' @param output_file Character string giving the output vector file.
#' @param keep_fields Character vector containing attribute fields to
#'   retain. If `NULL`, all fields are retained. An empty character
#'   vector retains geometry only.
#' @param exclude_proposed Logical. If `TRUE`, records with
#'   `STATUS == "PROPOSED"` are removed when the field is available.
#'
#' @return The output file path invisibly.
#'
#' @details
#' A spatial filter is applied in the native CRS of the input dataset
#' before reading the complete object. This substantially reduces
#' memory requirements for large continental datasets such as WDPA.
#'
#' Geometry collections are reduced to their polygon components.
#' Features that cannot be repaired are discarded.
#'
crop_vector_to_grid <- function(
                file,
                grid,
                output_file,
                keep_fields = NULL,
                exclude_proposed = FALSE) {
        
        message(
                "\nReading: ",
                file
        )
        
        
        #### Read source CRS ####
        
        source_proxy <- terra::vect(
                file,
                proxy = TRUE
        )
        
        source_crs <- terra::crs(
                source_proxy
        )
        
        
        #### Build spatial filter ####
        
        grid_filter <- sf::st_bbox(grid) |>
                sf::st_transform(
                        crs = source_crs,
                        densify = 101
                ) |>
                sf::st_as_sfc()
        
        
        #### Read overlapping features ####
        
        x <- sf::st_read(
                file,
                wkt_filter = sf::st_as_text(
                        grid_filter
                ),
                quiet = TRUE
        )
        
        message(
                "Features read: ",
                format(
                        nrow(x),
                        big.mark = ","
                )
        )
        
        message(
                "Original geometry types: ",
                paste(
                        names(
                                table(
                                        sf::st_geometry_type(x)
                                )
                        ),
                        collapse = ", "
                )
        )
        
        
        #### Remove proposed sites ####
        
        if (
                exclude_proposed &&
                "STATUS" %in% names(x)
        ) {
                
                x <- x[
                        is.na(x$STATUS) |
                                toupper(x$STATUS) != "PROPOSED",
                ]
        }
        
        
        #### Repair invalid geometries ####
        
        x <- repair_geometries(x)
        
        
        #### Keep polygonal geometries only ####
        
        geom_type <- as.character(
                sf::st_geometry_type(x)
        )
        
        polygon_direct <- geom_type %in%
                c(
                        "POLYGON",
                        "MULTIPOLYGON"
                )
        
        geometry_collection <-
                geom_type == "GEOMETRYCOLLECTION"
        
        x_polygon <- x[
                polygon_direct,
        ]
        
        if (any(geometry_collection)) {
                
                x_collection <- suppressWarnings(
                        sf::st_collection_extract(
                                x[
                                        geometry_collection,
                                ],
                                "POLYGON"
                        )
                )
                
                x_polygon <- rbind(
                        x_polygon,
                        x_collection
                )
        }
        
        x <- x_polygon
        
        x <- x[
                !sf::st_is_empty(x),
        ]
        
        message(
                "Features after geometry cleaning: ",
                format(
                        nrow(x),
                        big.mark = ","
                )
        )
        
        
        #### Transform to grid CRS ####
        
        x <- sf::st_transform(
                x,
                sf::st_crs(grid)
        )
        
        
        #### Check validity after transformation ####
        
        valid_after_transform <- sf::st_is_valid(x)
        
        n_invalid <- sum(
                is.na(valid_after_transform) |
                        !valid_after_transform
        )
        
        message(
                "Invalid geometries after transformation: ",
                n_invalid
        )
        
        if (n_invalid > 0) {
                x <- repair_geometries(x)
        }
        
        
        #### Crop to study extent ####
        
        x <- suppressWarnings(
                sf::st_crop(
                        x,
                        sf::st_bbox(grid)
                )
        )
        
        x <- x[
                !sf::st_is_empty(x),
        ]
        
        
        #### Retain required attributes ####
        
        if (!is.null(keep_fields)) {
                
                available_fields <- intersect(
                        keep_fields,
                        names(x)
                )
                
                x <- x[
                        ,
                        available_fields,
                        drop = FALSE
                ]
        }
        
        
        #### Final checks ####
        
        message(
                "Final features: ",
                format(
                        nrow(x),
                        big.mark = ","
                )
        )
        
        message(
                "Final geometry types: ",
                paste(
                        names(
                                table(
                                        sf::st_geometry_type(x)
                                )
                        ),
                        collapse = ", "
                )
        )
        
        
        #### Save ####
        
        sf::st_write(
                x,
                output_file,
                delete_dsn = TRUE,
                quiet = TRUE
        )
        
        message(
                "Saved: ",
                output_file
        )
        
        invisible(output_file)
}


#' Clip a vector dataset to a spatial boundary
#'
#' Clips polygon features to a supplied boundary while preserving
#' their original attributes. Boundary polygons are dissolved before
#' intersection so internal administrative borders do not split the
#' input geometries.
#'
#' @param input_file Character string giving the input vector file.
#' @param boundary An `sf` object defining the clipping boundary.
#' @param output_file Character string giving the output vector file.
#'
#' @return The output file path invisibly.
#'
#' @details
#' The boundary is transformed to the CRS of the input dataset and
#' dissolved with `st_union()` before clipping. Features completely
#' outside the boundary are removed, whereas features crossing the
#' boundary retain only their intersecting portion.
#'
crop_vector_to_boundary <- function(
                input_file,
                boundary,
                output_file) {
        
        message(
                "\nClipping to boundary: ",
                input_file
        )
        
        
        #### Read input ####
        
        x <- sf::st_read(
                input_file,
                quiet = TRUE
        )
        
        
        #### Prepare boundary ####
        
        boundary <- boundary |>
                sf::st_make_valid() |>
                sf::st_transform(
                        sf::st_crs(x)
                )
        
        boundary <- sf::st_union(
                boundary
        )
        
        
        #### Retain intersecting features ####
        
        x <- sf::st_filter(
                x,
                boundary,
                .predicate = sf::st_intersects
        )
        
        
        #### Clip geometries to boundary ####
        
        x <- suppressWarnings(
                sf::st_intersection(
                        x,
                        boundary
                )
        )
        
        
        #### Repair and retain polygon geometries ####
        
        x <- repair_geometries(x)
        
        geom_type <- as.character(
                sf::st_geometry_type(x)
        )
        
        x <- x[
                geom_type %in%
                        c(
                                "POLYGON",
                                "MULTIPOLYGON"
                        ),
        ]
        
        x <- x[
                !sf::st_is_empty(x),
        ]
        
        
        #### Save ####
        
        sf::st_write(
                x,
                output_file,
                delete_dsn = TRUE,
                quiet = TRUE
        )
        
        message(
                "Features after boundary clipping: ",
                format(
                        nrow(x),
                        big.mark = ","
                )
        )
        
        message(
                "Saved: ",
                output_file
        )
        
        invisible(output_file)
}




#' Extract Natura 2000 records from a WDPA dataset
#'
#' Identifies Natura 2000 records in a processed WDPA vector layer
#' using the English and original designation fields.
#'
#' @param input_file Character string giving the input WDPA vector file.
#' @param output_file Character string giving the output vector file.
#'
#' @return The output file path invisibly.
#'
#' @details
#' Records are identified using the following designation terms:
#' `"Natura 2000"`, `"Special Protection Area"`,
#' `"Site of Community Importance"`, and
#' `"Special Area of Conservation"`.
#'
#' When `DESIG_TYPE` is available, records are retained when the
#' designation type is `"Regional"` or missing.
#'
extract_wdpa_n2k <- function(
                input_file,
                output_file) {
        
        n2k_pattern <- paste(
                c(
                        "Natura[[:space:]]*2000",
                        "Special Protection Area",
                        "Site of Community Importance",
                        "Special Area of Conservation"
                ),
                collapse = "|"
        )
        
        x <- sf::st_read(
                input_file,
                quiet = TRUE
        )
        
        if (!"DESIG_ENG" %in% names(x)) {
                stop(
                        "DESIG_ENG is missing from: ",
                        input_file
                )
        }
        
        designation_text <- x$DESIG_ENG
        
        if ("DESIG" %in% names(x)) {
                
                designation_text <- paste(
                        designation_text,
                        x$DESIG
                )
        }
        
        is_n2k <- grepl(
                n2k_pattern,
                designation_text,
                ignore.case = TRUE
        )
        
        if ("DESIG_TYPE" %in% names(x)) {
                
                is_n2k <- is_n2k &
                        (
                                is.na(x$DESIG_TYPE) |
                                        toupper(
                                                x$DESIG_TYPE
                                        ) == "REGIONAL"
                        )
        }
        
        n2k <- x[
                is_n2k &
                        !is.na(is_n2k),
        ]
        
        message(
                basename(input_file),
                ": ",
                nrow(n2k),
                " Natura 2000 records"
        )
        
        sf::st_write(
                n2k,
                output_file,
                delete_dsn = TRUE,
                quiet = TRUE
        )
        
        invisible(output_file)
}


#' Calculate polygon coverage for one processing tile
#'
#' Calculates the percentage of each 10-km grid cell covered by one
#' or more polygon datasets. Protected-area polygons are dissolved
#' before rasterisation so overlapping designations are counted only
#' once.
#'
#' @param tile_name Character string identifying the processing tile.
#' @param grid An `sf` object containing the complete analysis grid.
#'   The object must contain a `tile_id` column.
#' @param vector_files Character vector containing one or more polygon
#'   datasets in the same CRS as `grid`.
#' @param output_name Character string used to create the cache
#'   subdirectory.
#' @param coverage_name Character string used as the output coverage
#'   column name.
#' @param tile_dir Character string giving the parent directory used
#'   for cached tile results.
#' @param id Character string identifying the unique grid-cell ID.
#'   Defaults to `"CellCode"`.
#' @param grid_res Numeric grid-cell resolution in metres.
#'   Defaults to `10000`.
#' @param raster_res Numeric temporary raster resolution in metres.
#'   Defaults to `100`.
#' @param overwrite Logical. If `TRUE`, cached results are ignored and
#'   the tile is recalculated.
#'
#' @return A data frame containing one row per grid cell and the
#'   requested percentage-coverage variable.
#'
#' @details
#' Only polygons intersecting the current tile are read from disk.
#' Polygon overlaps are dissolved with `st_union()` before
#' rasterisation.
#'
#' Coverage is rasterised at `raster_res` resolution using fractional
#' pixel coverage and aggregated to the original grid resolution.
#'
#' Tile results are cached as RDS files. Corrupted cache files are
#' automatically discarded and recalculated. New cache files are
#' written through a temporary file to reduce the risk of corruption
#' if the analysis is interrupted.
#'
process_coverage_tile <- function(
                tile_name,
                grid,
                vector_files,
                output_name,
                coverage_name,
                tile_dir,
                id = "CellCode",
                grid_res = 10000,
                raster_res = 100,
                overwrite = FALSE) {
        
        
        #### Initial checks ####
        
        if (!"tile_id" %in% names(grid)) {
                stop(
                        "The grid does not contain a 'tile_id' column."
                )
        }
        
        if (!id %in% names(grid)) {
                stop(
                        "Grid ID column '",
                        id,
                        "' was not found."
                )
        }
        
        missing_files <- vector_files[
                !file.exists(vector_files)
        ]
        
        if (length(missing_files) > 0) {
                stop(
                        "Missing vector files:\n",
                        paste(
                                missing_files,
                                collapse = "\n"
                        )
                )
        }
        
        
        #### Initialise tile output ####
        
        message(
                "\n----------------------------------------"
        )
        
        message(
                output_name,
                " | Tile ",
                tile_name
        )
        
        output_subdir <- file.path(
                tile_dir,
                output_name
        )
        
        dir.create(
                output_subdir,
                recursive = TRUE,
                showWarnings = FALSE
        )
        
        output_file <- file.path(
                output_subdir,
                paste0(
                        tile_name,
                        ".rds"
                )
        )
        
        
        #### Read cached result ####
        
        if (
                file.exists(output_file) &&
                !overwrite
        ) {
                
                cached_result <- tryCatch(
                        readRDS(output_file),
                        error = function(e) NULL
                )
                
                if (!is.null(cached_result)) {
                        
                        message(
                                "Cached result found."
                        )
                        
                        return(
                                cached_result
                        )
                }
                
                message(
                        "Corrupted cache found. ",
                        "Deleting and recomputing tile."
                )
                
                unlink(
                        output_file
                )
        }
        
        
        #### Select grid cells belonging to tile ####
        
        grid_tile <- grid |>
                dplyr::filter(
                        tile_id == tile_name
                )
        
        if (nrow(grid_tile) == 0) {
                stop(
                        "No grid cells found for tile: ",
                        tile_name
                )
        }
        
        tile_bbox <- sf::st_bbox(
                grid_tile
        )
        
        tile_filter <- sf::st_as_sfc(
                tile_bbox
        )
        
        
        #### Read polygons overlapping tile ####
        
        polygon_list <- purrr::map(
                vector_files,
                function(f) {
                        
                        x <- sf::st_read(
                                f,
                                wkt_filter = sf::st_as_text(
                                        tile_filter
                                ),
                                quiet = TRUE
                        )
                        
                        if (nrow(x) == 0) {
                                return(NULL)
                        }
                        
                        x <- suppressWarnings(
                                sf::st_crop(
                                        x,
                                        tile_bbox
                                )
                        )
                        
                        x <- x[
                                !sf::st_is_empty(x),
                        ]
                        
                        if (nrow(x) == 0) {
                                return(NULL)
                        }
                        
                        # Keep geometry only, regardless of the
                        # geometry-column name.
                        x[
                                ,
                                0,
                                drop = FALSE
                        ]
                }
        ) |>
                purrr::compact()
        
        
        #### Return zero if no polygons occur in tile ####
        
        if (length(polygon_list) == 0) {
                
                result <- grid_tile |>
                        sf::st_drop_geometry() |>
                        dplyr::transmute(
                                "{id}" :=
                                        .data[[id]],
                                
                                "{coverage_name}" :=
                                        0
                        )
                
                save_rds_safe(
                        result,
                        output_file
                )
                
                return(result)
        }
        
        
        #### Dissolve overlapping polygons ####
        
        polygons_tile <- dplyr::bind_rows(
                polygon_list
        )
        
        message(
                "Polygon features in tile: ",
                format(
                        nrow(polygons_tile),
                        big.mark = ","
                )
        )
        
        polygon_union <- sf::st_union(
                sf::st_geometry(
                        polygons_tile
                )
        )
        
        polygon_union <- sf::st_sf(
                value = 1,
                geometry = polygon_union,
                crs = sf::st_crs(grid)
        )
        
        
        #### Create temporary raster ####
        
        fact <- grid_res /
                raster_res
        
        if (fact != as.integer(fact)) {
                stop(
                        "Raster resolution does not divide ",
                        "exactly into the grid resolution."
                )
        }
        
        fact <- as.integer(
                fact
        )
        
        template <- terra::rast(
                xmin = as.numeric(
                        tile_bbox["xmin"]
                ),
                xmax = as.numeric(
                        tile_bbox["xmax"]
                ),
                ymin = as.numeric(
                        tile_bbox["ymin"]
                ),
                ymax = as.numeric(
                        tile_bbox["ymax"]
                ),
                resolution = raster_res,
                crs = sf::st_crs(grid)$wkt
        )
        
        if (
                terra::ncol(template) %% fact != 0 ||
                terra::nrow(template) %% fact != 0
        ) {
                stop(
                        "Temporary raster is not aligned ",
                        "with the analysis grid."
                )
        }
        
        
        #### Rasterise polygon coverage ####
        
        coverage_raster <- terra::rasterize(
                terra::vect(
                        polygon_union
                ),
                template,
                field = "value",
                background = 0,
                cover = TRUE
        )
        
        
        #### Aggregate coverage to grid resolution ####
        
        coverage_grid <- terra::aggregate(
                coverage_raster,
                fact = fact,
                fun = "mean",
                na.rm = TRUE
        ) * 100
        
        names(
                coverage_grid
        ) <- coverage_name
        
        
        #### Extract coverage values ####
        
        grid_points <- suppressWarnings(
                sf::st_point_on_surface(
                        grid_tile
                )
        )
        
        extracted <- terra::extract(
                coverage_grid,
                terra::vect(
                        grid_points
                )
        )
        
        result <- grid_tile |>
                sf::st_drop_geometry() |>
                dplyr::transmute(
                        "{id}" :=
                                .data[[id]],
                        
                        "{coverage_name}" :=
                                extracted[[coverage_name]])
        
        
        #### Numerical quality control ####
        
        values <- result[[coverage_name]]
        
        if (
                any(
                        values < -1e-6 |
                        values > 100 + 1e-6,
                        na.rm = TRUE
                )
        ) {
                stop(
                        "Coverage outside the expected 0-100% range."
                )
        }
        
        result[[coverage_name]] <- pmin(
                100,
                pmax(
                        0,
                        values
                )
        )
        
        
        #### Save result ####
        
        save_rds_safe(
                result,
                output_file
        )
        
        
        #### Release large temporary objects ####
        
        rm(
                polygons_tile,
                polygon_union,
                coverage_raster,
                coverage_grid
        )
        
        gc()
        
        result
}