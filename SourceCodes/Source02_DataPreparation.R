# Shared functions for scripts 1 and 2: report checks and data preparation.
# All three functions use base R; no packages are attached by this source.
# Function bodies and output names are unchanged from the original Source.R.
# temporal_bgl requires the progress callback p defined by the calling script.
# temporal_bgl_group does not call a progress callback.
# Rows are compared in the order of meta$period.
# In script 1, the legacy richness columns refer to r01 and r02 of 2007-2012.
# In script 2, they refer to the actual 2007-2012 and 2013-2018 periods.
# habitat_temporal_trend uses all cells in meta and fills missing periods with 0.

#' Calculate habitat change between two records of a grid cell
#'
#' Computes habitat richness, shared habitats, losses, gains and temporal
#' Jaccard dissimilarity for one cell across two periods or report revisions.
#'
#' @param id A single cell identifier matching `meta$cell_id`.
#' @param comm_mat A numeric presence/absence matrix or data frame, with
#'   samples in rows and habitats in columns. Values must be 0 or 1, without
#'   missing values. Rows must correspond to `meta` in the same order.
#' @param meta A data frame containing `cell_id` and `period`. The selected
#'   cell must have exactly two rows, one for each period or revision.
#'
#' @return A one-row data frame containing `cell_id`, `beta_jaccard`,
#'   `richness_2007_2012`, `richness_2013_2018`, `shared_species`,
#'   `lost_species` and `gained_species`. Despite their legacy names,
#'   the species columns contain habitat counts.
#'
#' @details
#' Records are ordered by `meta$period`; losses and gains are measured from
#' the first to the second record. Period labels must sort in the intended
#' comparison order. The richness columns refer to these ordered records:
#' in script 1 they represent r01 and r02 of the same 2007-2012 report;
#' in script 2 they represent the two reporting periods.
#'
#' Jaccard dissimilarity is the number of lost plus gained habitats divided
#' by the number of habitats present in either record. If both records are
#' empty, the result is `NaN`. Input requirements are not checked internally.
#'
#' Calls an externally defined function `p()` once per cell to update
#' progress. This callback must be accessible when the function runs.
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

#' Calculate temporal habitat change by habitat group
#'
#' Computes richness, shared habitats, losses, gains and Jaccard
#' dissimilarity separately for each habitat group within one grid cell.
#'
#' @param id A single cell identifier matching `meta$cell_id`.
#' @param comm_mat A numeric presence/absence matrix or data frame, with
#'   samples in rows and habitat codes as column names. Values must be 0 or
#'   1, without missing values. Rows must match `meta` in the same order.
#' @param meta A data frame containing `cell_id` and `period`. The selected
#'   cell must have exactly two rows, one for each period.
#' @param group A data frame containing `habitat` codes and their `taxGroup`
#'   assignments, without missing values. Habitat codes are converted to
#'   character and matched to the column names of `comm_mat`.
#'
#' @return A data frame with one row per distinct `taxGroup`, containing
#'   `cell_id`, `taxGroup`, `beta_jaccard`, `richness_2007_2012`,
#'   `richness_2013_2018`, `shared_species`, `lost_species` and
#'   `gained_species`. The species columns contain habitat counts.
#'
#' @details
#' Records are ordered by `meta$period`; the two richness columns refer to
#' the first and second records, respectively. Period labels must sort in
#' chronological order. Losses and gains are measured within each group.
#'
#' Habitats absent from `group` are excluded from these calculations.
#' Jaccard dissimilarity is `NA` when no habitat in a group is present in
#' either period, including groups with no matching matrix columns.
#' Other input requirements are not checked internally.
#'
#' This function does not call a progress callback.
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

#' Summarise temporal changes in occupancy for each habitat
#'
#' Counts occupied cells in two periods and distinguishes persistent
#' occurrences, losses and gains, including changes hidden by stable
#' total occupancy.
#'
#' @param comm_mat A numeric matrix or data frame with samples in rows and
#'   habitat codes as column names. Positive values indicate presence;
#'   zero indicates absence. Values should be non-negative and non-missing.
#'   Rows must correspond to `meta` in the same order.
#' @param meta A data frame containing `cell_id` and `period`, with exactly
#'   two distinct periods and no duplicated cell-period combinations.
#'
#' @return A data frame with one row per habitat and columns `habitat`,
#'   `Ncelle_2012`, `Ncelle_2018`, `shared`, `loss`, `gain`, `netto`,
#'   `loss_pct_2012`, `gain_pct_2012` and `netto_pct_2012`.
#'   Net change (`netto`) equals gains minus losses. Percentages use
#'   first-period occupancy as their denominator and are `NA` when it is zero.
#'
#' @details
#' Periods are sorted using `sort()`: the first is assigned to the 2012
#' output columns and the second to the 2018 columns. Labels must therefore
#' sort in the intended chronological order.
#'
#' All cells occurring in either period are included. A missing cell-period
#' record is filled with zeros for every habitat, thus treating that record
#' as absence rather than unknown occupancy. The function does not restrict
#' calculations to cells recorded in both periods.
#'
#' Values are converted to presence/absence using `> 0`. Shared occurrences
#' are cells occupied in both periods; losses and gains are cells occupied
#' only in the first or second period, respectively.
#'
#' Checks equal row counts in `comm_mat` and `meta`, exactly two periods,
#' and unique cell-period combinations. Row alignment and value validity
#' remain the caller's responsibility.
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
