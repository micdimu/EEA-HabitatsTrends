# Functions for script 03: Figure 1 and Figure S1.
# bivariate_map uses explicit package namespaces.
# density_panel requires ggplot2 to be attached (also loaded by tidyverse).

#' Draw a bivariate choropleth map and its legend
#'
#' Classifies two numeric variables and maps their joint classes using a
#' bivariate colour palette. Returns the map and legend separately as well
#' as a combined plot with an inset legend.
#'
#' @param data A data frame or sf object containing the ID and both variables.
#' @param x,y Numeric columns, supplied as bare names or character strings.
#' @param geometry An sf object containing `id`. If NULL, uses the geometry
#'   of `data`, which must then be an sf object.
#' @param id Name of the column used to join attributes and geometry.
#'   Identifiers should uniquely identify features in each input.
#' @param x_breaks,y_breaks Class boundaries: NULL for quantiles, `dim - 1`
#'   internal thresholds, or `dim + 1` complete boundaries.
#' @param dim Number of classes on each axis; defaults to 4.
#' @param pal Name of a palette supported by biscale for the chosen dimension.
#' @param flip_axes Whether to flip the palette axes in the map and legend.
#' @param xlab,ylab Legend axis labels; defaults to the variable names.
#' @param legend_position Named vector with `left`, `bottom`, `right` and
#'   `top` coordinates for the inset legend, relative to the map panel.
#' @param scale_bar Whether to add a scale bar using ggspatial.
#' @param save_path Optional output path for the combined plot. The saving
#'   call passes LZW compression and is intended for TIFF output.
#' @param width,height Dimensions used when saving the combined plot.
#' @param units Units of the saved dimensions; defaults to "mm".
#' @param dpi Resolution of the saved plot; defaults to 300.
#'
#' @return A list containing `plot` (map with inset legend), `map`, `legend`,
#'   `data` (joined sf object with class assignments), and `breaks` (a list
#'   containing the x and y class boundaries).
#'
#' @details
#' Internal helpers `make_breaks()` and `legend_values()` construct class
#' boundaries and legend labels; they are local to this function.
#' Quantile boundaries are computed from finite values. Internal thresholds
#' and quantile boundaries receive infinite outer bounds. Duplicate
#' boundaries cause an error. Classes are right-closed, so a value exactly
#' equal to an internal threshold belongs to the lower class.
#'
#' Unclassified attribute rows are excluded before the left join to geometry;
#' geometries without a class are displayed in white. Legend labels show
#' upper class boundaries, replacing an infinite boundary with the observed
#' extreme among classified values.
#'
#' Requires dplyr, rlang, sf, ggplot2, biscale and patchwork, plus ggspatial
#' when `scale_bar = TRUE`. The function does not transform the input CRS.
bivariate_map <- function(data,
                          x,
                          y,
                          geometry = NULL,
                          id = "CellCode",
                          x_breaks = NULL,
                          y_breaks = NULL,
                          dim = 4,
                          pal = "Brown2",
                          flip_axes = FALSE,
                          xlab = NULL,
                          ylab = NULL,
                          legend_position = c(
                                  left = 0.04,
                                  bottom = 0.68,
                                  right = 0.30,
                                  top = 0.96
                          ),
                          scale_bar = TRUE,
                          save_path = NULL,
                          width = 400,
                          height = 320,
                          units = "mm",
                          dpi = 300) {
        
        # Permette sia x = total_pa_cov sia x = "total_pa_cov"
        x_sym <- rlang::ensym(x)
        y_sym <- rlang::ensym(y)
        
        x_name <- rlang::as_string(x_sym)
        y_name <- rlang::as_string(y_sym)
        
        if (is.null(xlab)) xlab <- x_name
        if (is.null(ylab)) ylab <- y_name
        
        # Se data è già sf, usa direttamente la sua geometria
        if (inherits(data, "sf") && is.null(geometry)) {
                
                geometry_sf <- data |>
                        dplyr::select(dplyr::all_of(id))
                
                attributes_df <- sf::st_drop_geometry(data)
                
        } else {
                
                if (is.null(geometry)) {
                        stop(
                                "`geometry` deve essere fornito quando `data` non è un oggetto sf."
                        )
                }
                
                if (!inherits(geometry, "sf")) {
                        stop("`geometry` deve essere un oggetto sf.")
                }
                
                geometry_sf <- geometry |>
                        dplyr::select(dplyr::all_of(id))
                
                attributes_df <- if (inherits(data, "sf")) {
                        sf::st_drop_geometry(data)
                } else {
                        data
                }
        }
        
        required_columns <- c(id, x_name, y_name)
        
        if (!all(required_columns %in% names(attributes_df))) {
                missing_columns <- setdiff(
                        required_columns,
                        names(attributes_df)
                )
                
                stop(
                        "Colonne mancanti in `data`: ",
                        paste(missing_columns, collapse = ", ")
                )
        }
        
        # Crea i limiti delle classi
        make_breaks <- function(values, breaks, dim, variable_name) {
                
                values <- values[
                        is.finite(values)
                ]
                
                if (length(values) == 0) {
                        stop(
                                "La variabile `",
                                variable_name,
                                "` non contiene valori validi."
                        )
                }
                
                # Default: classi basate sui quantili
                if (is.null(breaks)) {
                        
                        breaks <- stats::quantile(
                                values,
                                probs = seq(0, 1, length.out = dim + 1),
                                na.rm = TRUE,
                                names = FALSE
                        )
                        
                        # Garantisce che siano inclusi anche i valori estremi
                        breaks[1] <- -Inf
                        breaks[length(breaks)] <- Inf
                        
                } else if (length(breaks) == dim - 1) {
                        
                        # Sono state fornite soltanto le soglie interne
                        breaks <- c(-Inf, breaks, Inf)
                        
                } else if (length(breaks) != dim + 1) {
                        
                        stop(
                                "`",
                                variable_name,
                                "_breaks` deve contenere ",
                                dim - 1,
                                " soglie interne oppure ",
                                dim + 1,
                                " estremi completi."
                        )
                }
                
                if (anyDuplicated(breaks)) {
                        stop(
                                "Non è possibile costruire ",
                                dim,
                                " classi per `",
                                variable_name,
                                "`: alcuni limiti coincidono."
                        )
                }
                
                breaks
        }
        
        map_values <- attributes_df |>
                dplyr::transmute(
                        !!id := .data[[id]],
                        .x = as.numeric(!!x_sym),
                        .y = as.numeric(!!y_sym)
                )
        
        x_limits <- make_breaks(
                values = map_values$.x,
                breaks = x_breaks,
                dim = dim,
                variable_name = x_name
        )
        
        y_limits <- make_breaks(
                values = map_values$.y,
                breaks = y_breaks,
                dim = dim,
                variable_name = y_name
        )
        
        classified_data <- map_values |>
                dplyr::mutate(
                        x_class = cut(
                                .x,
                                breaks = x_limits,
                                labels = seq_len(dim),
                                include.lowest = TRUE,
                                right = TRUE
                        ),
                        y_class = cut(
                                .y,
                                breaks = y_limits,
                                labels = seq_len(dim),
                                include.lowest = TRUE,
                                right = TRUE
                        ),
                        bi_class = paste(
                                as.integer(x_class),
                                as.integer(y_class),
                                sep = "-"
                        )
                ) |>
                dplyr::filter(
                        !is.na(x_class),
                        !is.na(y_class)
                )
        
        map_sf <- geometry_sf |>
                dplyr::left_join(
                        classified_data,
                        by = id
                )
        
        # Valori visualizzati nella legenda: limite superiore di ciascuna classe
        legend_values <- function(limits, values) {
                
                labels <- limits[-1]
                
                labels[is.infinite(labels) & labels > 0] <-
                        max(values, na.rm = TRUE)
                
                labels[is.infinite(labels) & labels < 0] <-
                        min(values, na.rm = TRUE)
                
                format(
                        signif(labels, 3),
                        trim = TRUE,
                        scientific = FALSE
                )
        }
        
        legend_breaks <- list(
                bi_x = legend_values(
                        x_limits,
                        classified_data$.x
                ),
                bi_y = legend_values(
                        y_limits,
                        classified_data$.y
                )
        )
        
        map_plot <- ggplot2::ggplot() +
                ggplot2::geom_sf(
                        data = map_sf,
                        ggplot2::aes(fill = bi_class),
                        linewidth = 0,
                        colour = NA,
                        show.legend = FALSE
                ) +
                biscale::bi_scale_fill(
                        pal = pal,
                        dim = dim,
                        flip_axes = flip_axes,
                        na.value = "white"
                ) +
                biscale::bi_theme(
                        bg_color = "#FFFFFF"
                ) +
                ggplot2::theme(
                        panel.background = ggplot2::element_rect(
                                fill = "white",
                                colour = NA
                        ),
                        plot.background = ggplot2::element_rect(
                                fill = "white",
                                colour = NA
                        ),
                        panel.grid = ggplot2::element_blank()
                )
        
        if (scale_bar) {
                map_plot <- map_plot +
                        ggspatial::annotation_scale(
                                location = "br",
                                bar_cols = c("black", "white"),
                                text_family = "sans",
                                line_width = 0.5,
                                height = grid::unit(0.20, "cm"),
                                pad_x = grid::unit(0.25, "cm"),
                                pad_y = grid::unit(0.25, "cm")
                        )
        }
        
        legend_plot <- biscale::bi_legend(
                pal = pal,
                dim = dim,
                flip_axes = flip_axes,
                xlab = xlab,
                ylab = ylab,
                size = 8,
                breaks = legend_breaks
        ) +
                ggplot2::theme(
                        axis.text.x = ggplot2::element_text(
                                angle = 35,
                                hjust = 1
                        )
                )
        
        combined_plot <- map_plot +
                patchwork::inset_element(
                        legend_plot,
                        left = legend_position[["left"]],
                        bottom = legend_position[["bottom"]],
                        right = legend_position[["right"]],
                        top = legend_position[["top"]],
                        align_to = "panel"
                )
        
        if (!is.null(save_path)) {
                ggplot2::ggsave(
                        filename = save_path,
                        plot = combined_plot,
                        width = width,
                        height = height,
                        units = units,
                        dpi = dpi,
                        compression = "lzw"
                )
        }
        
        list(
                plot = combined_plot,
                map = map_plot,
                legend = legend_plot,
                data = map_sf,
                breaks = list(
                        x = x_limits,
                        y = y_limits
                )
        )
}

#' Draw a density panel with the sample median
#'
#' Displays a kernel density estimate, a rug of observations, a dashed median
#' line and a median annotation. Optionally adds a reference line at zero.
#'
#' @param data A data frame containing the variable to plot.
#' @param variable Character string naming a numeric column in `data`.
#' @param x_label Label for the horizontal axis.
#' @param fill Fill colour for the density curve.
#' @param zero_line Whether to add a solid reference line at zero.
#' @param x_limits Currently unused. Retained for compatibility with the
#'   calling script; supplying it does not restrict the plotted range.
#'
#' @return A ggplot object representing one density panel.
#'
#' @details
#' The median ignores missing values and is annotated to two decimal places.
#' Density smoothing uses `adjust = 4`. Missing observations are removed
#' by the density and rug layers. Density estimation requires sufficient
#' finite observations; no additional input checks are performed here.
#' Requires ggplot2 to be attached before the function is called.
density_panel <- function(data,
                          variable,
                          x_label,
                          fill,
                          zero_line = FALSE,
                          x_limits = NULL) {
        
        mediana <- median(data[[variable]], na.rm = TRUE)
        
        p <- ggplot(data, aes(x = .data[[variable]])) +
                geom_density(
                        fill = fill,
                        colour = "black",
                        linewidth = 0.45,
                        alpha = 0.75,
                        adjust = 4,
                        na.rm = TRUE
                ) +
                geom_rug(
                        sides = "b",
                        alpha = 0.20,
                        linewidth = 0.25,
                        na.rm = TRUE
                ) +
                geom_vline(
                        xintercept = mediana,
                        linewidth = 0.55,
                        linetype = "dashed"
                ) +
                annotate(
                        "text",
                        x = Inf,
                        y = Inf,
                        label = paste0("Median = ", round(mediana, 2)),
                        hjust = 1.08,
                        vjust = 1.5,
                        size = 2.8
                ) +
                labs(
                        x = x_label,
                        y = "Density"
                ) +
                theme_classic(base_size = 9) +
                theme(
                        axis.title = element_text(size = 9),
                        axis.text = element_text(size = 8, colour = "black"),
                        axis.line = element_line(linewidth = 0.4),
                        axis.ticks = element_line(linewidth = 0.4),
                        plot.margin = margin(6, 8, 5, 6)
                )
        
        if (zero_line) {
                p <- p +
                        geom_vline(
                                xintercept = 0,
                                linewidth = 0.45,
                                colour = "grey35"
                        )
        }
        
        p
}
