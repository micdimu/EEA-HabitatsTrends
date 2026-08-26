library(tidyverse)
library(biscale)

#### temporal functions #####

temporal_bgl <- function(id, comm_mat, meta) {
        
        p()
        
        rows <- meta$cell_id == id
        m <- comm_mat[rows, , drop = FALSE]
        m <- m[order(meta$period[rows]), , drop = FALSE]
        
        a <- m[1, ]
        b <- m[2, ]
        
        shared <- sum(a == 1 & b == 1)
        lost   <- sum(a == 1 & b == 0)
        gained <- sum(a == 0 & b == 1)
        
        beta_jaccard <- (lost + gained) / (shared + lost + gained)
        
        data.frame(
                cell_id = id,
                beta_jaccard = beta_jaccard,
                richness_2007_2012 = sum(a),
                richness_2013_2018 = sum(b),
                shared_species = shared,
                lost_species = lost,
                gained_species = gained
        )
}

temporal_bgl_group <- function(id, comm_mat, meta, group) {
        
        rows <- meta$cell_id == id
        
        m <- comm_mat[rows, , drop = FALSE]
        m <- m[order(meta$period[rows]), , drop = FALSE]
        
        a <- m[1, ]
        b <- m[2, ]
        
        # Codici habitat come caratteri
        group$habitat <- as.character(group$habitat)
        
        # Calcolo separato per ciascun taxGroup
        risultati <- lapply(unique(group$taxGroup), function(g) {
                
                habitat_g <- group$habitat[group$taxGroup == g]
                
                # Colonne della matrice appartenenti al taxGroup
                cols <- colnames(comm_mat) %in% habitat_g
                
                a_g <- a[cols]
                b_g <- b[cols]
                
                shared <- sum(a_g == 1 & b_g == 1)
                lost   <- sum(a_g == 1 & b_g == 0)
                gained <- sum(a_g == 0 & b_g == 1)
                
                totale <- shared + lost + gained
                
                beta_jaccard <- if (totale == 0) {
                        NA_real_
                } else {
                        (lost + gained) / totale
                }
                
                data.frame(
                        cell_id = id,
                        taxGroup = g,
                        beta_jaccard = beta_jaccard,
                        richness_2007_2012 = sum(a_g),
                        richness_2013_2018 = sum(b_g),
                        shared_species = shared,
                        lost_species = lost,
                        gained_species = gained
                )
        })
        
        do.call(rbind, risultati)
}


#### temporal each habitat ####

habitat_temporal_trend <- function(comm_mat, meta) {
        
        # Controlli essenziali
        if (nrow(comm_mat) != nrow(meta)) {
                stop("comm_mat e meta devono avere lo stesso numero di righe")
        }
        
        periods <- sort(unique(meta$period))
        
        if (length(periods) != 2) {
                stop("meta$period deve contenere esattamente due periodi")
        }
        
        if (anyDuplicated(meta[c("cell_id", "period")])) {
                stop("Ogni combinazione cell_id-period deve comparire una sola volta")
        }
        
        # Celle presenti in almeno uno dei due periodi
        cells <- unique(meta$cell_id)
        
        # Matrici vuote: celle × habitat
        m1 <- matrix(
                0,
                nrow = length(cells),
                ncol = ncol(comm_mat),
                dimnames = list(cells, colnames(comm_mat))
        )
        
        m2 <- m1
        
        rows1 <- meta$period == periods[1]
        rows2 <- meta$period == periods[2]
        
        # Inserisce i dati allineandoli per cell_id
        m1[match(meta$cell_id[rows1], cells), ] <-
                as.matrix(comm_mat[rows1, , drop = FALSE])
        
        m2[match(meta$cell_id[rows2], cells), ] <-
                as.matrix(comm_mat[rows2, , drop = FALSE])
        
        # Trasformazione in presenza/assenza
        m1 <- m1 > 0
        m2 <- m2 > 0
        
        ncelle_2012 <- colSums(m1)
        ncelle_2018 <- colSums(m2)
        
        shared <- colSums(m1 & m2)
        loss   <- colSums(m1 & !m2)
        gain   <- colSums(!m1 & m2)
        
        netto <- gain - loss
        
        data.frame(
                habitat = colnames(comm_mat),
                Ncelle_2012 = ncelle_2012,
                Ncelle_2018 = ncelle_2018,
                shared = shared,
                loss = loss,
                gain = gain,
                netto = netto,
                
                loss_pct_2012 = ifelse(
                        ncelle_2012 == 0,
                        NA_real_,
                        loss / ncelle_2012 * 100
                ),
                
                gain_pct_2012 = ifelse(
                        ncelle_2012 == 0,
                        NA_real_,
                        gain / ncelle_2012 * 100
                ),
                
                netto_pct_2012 = ifelse(
                        ncelle_2012 == 0,
                        NA_real_,
                        netto / ncelle_2012 * 100
                ),
                
                row.names = NULL
        )
}

#### plot glms results ####

plot_glm <- function(model,
                     data,
                     x = NULL,
                     y = NULL,
                     xlab = NULL,
                     ylab = NULL,
                     level = 0.95,
                     n = 200,
                     ylim = NULL,
                     point_alpha = 0.25,
                     point_size = 1.4) {
        
        if (!inherits(model, "glm")) {
                stop("`model` deve essere un oggetto glm.")
        }
        
        # Ricava risposta e predittore dalla formula
        formula_vars <- all.vars(stats::formula(model))
        
        if (is.null(y)) {
                y <- formula_vars[1]
        }
        
        if (is.null(x)) {
                x <- formula_vars[2]
        }
        
        if (!all(c(x, y) %in% names(data))) {
                stop("Le variabili `x` e `y` devono essere presenti in `data`.")
        }
        
        # Sequenza del predittore
        newdata <- data.frame(
                x_seq = seq(
                        min(data[[x]], na.rm = TRUE),
                        max(data[[x]], na.rm = TRUE),
                        length.out = n
                )
        )
        
        names(newdata) <- x
        
        # Predizioni sulla scala del link
        pred <- stats::predict(
                model,
                newdata = newdata,
                type = "link",
                se.fit = TRUE
        )
        
        z <- stats::qnorm(1 - (1 - level) / 2)
        linkinv <- stats::family(model)$linkinv
        
        newdata$.fit <- linkinv(pred$fit)
        newdata$.lower <- linkinv(pred$fit - z * pred$se.fit)
        newdata$.upper <- linkinv(pred$fit + z * pred$se.fit)
        
        p <- ggplot2::ggplot(
                data,
                ggplot2::aes(
                        x = .data[[x]],
                        y = .data[[y]]
                )
        ) +
                ggplot2::geom_point(
                        alpha = point_alpha,
                        size = point_size
                ) +
                ggplot2::geom_ribbon(
                        data = newdata,
                        ggplot2::aes(
                                x = .data[[x]],
                                ymin = .data$.lower,
                                ymax = .data$.upper
                        ),
                        inherit.aes = FALSE,
                        alpha = 0.18
                ) +
                ggplot2::geom_line(
                        data = newdata,
                        ggplot2::aes(
                                x = .data[[x]],
                                y = .data$.fit
                        ),
                        inherit.aes = FALSE,
                        linewidth = 0.8
                ) +
                ggplot2::labs(
                        x = xlab %||% x,
                        y = ylab %||% y
                ) +
                ggplot2::theme_classic(base_size = 11)
        
        if (!is.null(ylim)) {
                p <- p +
                        ggplot2::coord_cartesian(ylim = ylim)
        }
        
        p
}


plot_small_effect <- function(model,
                              data,
                              x,
                              y,
                              xlab = x,
                              ylab = y,
                              n = 200,
                              nsim = 5000,
                              level = 0.95,
                              full_ylim = c(0, 1),
                              seed = 123) {
        
        stopifnot(
                inherits(model, "glm"),
                x %in% names(data),
                y %in% names(data)
        )
        
        set.seed(seed)
        
        # Griglia del predittore
        x_range <- range(data[[x]], na.rm = TRUE)
        
        newdata <- data.frame(
                x_value = seq(
                        x_range[1],
                        x_range[2],
                        length.out = n
                )
        )
        
        names(newdata) <- x
        
        # Matrice del modello
        X <- model.matrix(
                delete.response(terms(model)),
                newdata
        )
        
        beta_hat <- coef(model)
        vcov_hat <- vcov(model)
        linkinv <- family(model)$linkinv
        
        # Simulazioni dalla distribuzione dei coefficienti
        beta_sim <- MASS::mvrnorm(
                n = nsim,
                mu = beta_hat,
                Sigma = vcov_hat
        )
        
        # Predizioni simulate: righe = valori di x, colonne = simulazioni
        eta_sim <- X %*% t(beta_sim)
        mu_sim <- linkinv(eta_sim)
        
        # Predizione centrale
        fit <- as.numeric(
                linkinv(X %*% beta_hat)
        )
        
        alpha <- (1 - level) / 2
        
        lower <- apply(
                mu_sim,
                1,
                quantile,
                probs = alpha
        )
        
        upper <- apply(
                mu_sim,
                1,
                quantile,
                probs = 1 - alpha
        )
        
        # Variazione rispetto al valore minimo osservato di x
        delta_sim <- sweep(
                mu_sim,
                MARGIN = 2,
                STATS = mu_sim[1, ],
                FUN = "-"
        ) * 100
        
        delta_fit <- (fit - fit[1]) * 100
        
        delta_lower <- apply(
                delta_sim,
                1,
                quantile,
                probs = alpha
        )
        
        delta_upper <- apply(
                delta_sim,
                1,
                quantile,
                probs = 1 - alpha
        )
        
        pred_data <- newdata |>
                dplyr::mutate(
                        fit = fit,
                        lower = lower,
                        upper = upper,
                        delta = delta_fit,
                        delta_lower = delta_lower,
                        delta_upper = delta_upper
                )
        
        total_delta <- tail(pred_data$delta, 1)
        
        # Pannello A: dati grezzi e modello sulla scala completa
        p_raw <- ggplot2::ggplot(
                data,
                ggplot2::aes(
                        x = .data[[x]],
                        y = .data[[y]]
                )
        ) +
                ggplot2::geom_point(
                        alpha = 0.18,
                        size = 1
                ) +
                ggplot2::geom_ribbon(
                        data = pred_data,
                        ggplot2::aes(
                                x = .data[[x]],
                                ymin = lower,
                                ymax = upper
                        ),
                        inherit.aes = FALSE,
                        alpha = 0.18
                ) +
                ggplot2::geom_line(
                        data = pred_data,
                        ggplot2::aes(
                                x = .data[[x]],
                                y = fit
                        ),
                        inherit.aes = FALSE,
                        linewidth = 0.8
                ) +
                ggplot2::labs(
                        x = xlab,
                        y = ylab
                ) +
                ggplot2::coord_cartesian(
                        ylim = full_ylim
                ) +
                ggplot2::theme_classic(base_size = 11)
        
        # Pannello B: dimensione dell'effetto
        p_effect <- ggplot2::ggplot(
                pred_data,
                ggplot2::aes(
                        x = .data[[x]],
                        y = delta
                )
        ) +
                ggplot2::geom_hline(
                        yintercept = 0,
                        linewidth = 0.4,
                        linetype = "dashed"
                ) +
                ggplot2::geom_ribbon(
                        ggplot2::aes(
                                ymin = delta_lower,
                                ymax = delta_upper
                        ),
                        alpha = 0.18
                ) +
                ggplot2::geom_line(
                        linewidth = 0.8
                ) +
                ggplot2::geom_rug(
                        data = data,
                        ggplot2::aes(x = .data[[x]]),
                        inherit.aes = FALSE,
                        sides = "b",
                        alpha = 0.12,
                        linewidth = 0.2
                ) +
                ggplot2::annotate(
                        "text",
                        x = Inf,
                        y = Inf,
                        label = sprintf(
                                "Overall predicted change = %.3f p.p.",
                                total_delta
                        ),
                        hjust = 1.05,
                        vjust = 1.4,
                        size = 3.2
                ) +
                ggplot2::labs(
                        x = xlab,
                        y = paste0(
                                "Predicted change in ",
                                ylab,
                                "\n(percentage points)"
                        )
                ) +
                ggplot2::theme_classic(base_size = 11)
        
        p_raw + p_effect +
                patchwork::plot_annotation(
                        tag_levels = "a"
                ) &
                ggplot2::theme(
                        plot.tag = ggplot2::element_text(
                                face = "bold",
                                size = 12
                        )
                )
}

#### geographic distribution #####

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
                          title = NULL,
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
        if (is.null(title)) title <- paste(ylab, "and", xlab)
        
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
                ggplot2::labs(title = title) +
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
