# 06a | Prepare H3 data, fit models and save inputs for 06b.
# Run from the project root.
library(tidyverse)
library(sf)
library(glmmTMB)

#### load data ####

comm_mat <-  read.csv("processed/comm_mat.csv", row.names = 1, check.names = FALSE)

meta <- comm_mat |> 
        rownames_to_column("sample_id") |>
        select(sample_id)  |> 
        separate(sample_id, into = c("cell_id", "period"), sep = "_(?=20)")

grid_pa <- st_read("data/PAs/derived/europe_10km_protected_area_coverage.gpkg")

##### inspect data #####

# Check object classes and dimensions.
class(comm_mat)
dim(comm_mat)

class(meta)
dim(meta)

# Inspect a small subset of the community matrix.
comm_mat[
        seq_len(min(6L, nrow(comm_mat))),
        seq_len(min(8L, ncol(comm_mat))),
        drop = FALSE
]

# Inspect row and column identifiers.
head(rownames(comm_mat))
head(colnames(comm_mat))

# Inspect metadata fields, data types, and first records.
str(meta)
head(meta)


###### Check identifiers and reporting periods ######

# Inspect reporting periods and record counts.
table(meta$period, useNA = "ifany")

# Verify metadata completeness and consistency with matrix row identifiers.
stopifnot(
        nrow(comm_mat) == nrow(meta),
        !anyNA(meta$cell_id),
        !anyNA(meta$period),
        !anyDuplicated(names(comm_mat)),
        setequal(unique(meta$period), c("2007_2012", "2013_2018"))
)

expected_ids <- paste0(
        meta$cell_id, "_",
        meta$period
)

stopifnot(identical(rownames(comm_mat), expected_ids))

# Each cell must have at most one record per reporting period.
stopifnot(!anyDuplicated(meta[c("cell_id", "period")]))

# Habitat records must contain only presence, absence, or missing values.
stopifnot(all(vapply(
        comm_mat,
        function(x) is.numeric(x) && all(is.na(x) | x %in% c(0, 1)),
        logical(1)
)))


##### Match cells between reporting periods ####

# Identify matrix rows belonging to each reporting period.
rows_t0 <- which(meta$period == "2007_2012")
rows_t1 <- which(meta$period == "2013_2018")

cells_t0 <- meta$cell_id[rows_t0]
cells_t1 <- meta$cell_id[rows_t1]

# The input comm_mat.csv must be rebuilt by the preparation script after
# completing missing reporting periods with zero reported habitats.
# This step matches the two records; it does not itself fill missing records.
# Stop on an outdated, incomplete input rather than silently exclude cells.
stopifnot(setequal(cells_t0, cells_t1))

# Retain and align cells represented in both completed periods.
common_cells <- intersect(cells_t0, cells_t1)

stopifnot(length(common_cells) > 0L)

# Report the number of shared and unmatched cells.
c(
        cells_t0 = length(cells_t0),
        cells_t1 = length(cells_t1),
        common_cells = length(common_cells),
        only_t0 = sum(!cells_t0 %in% cells_t1),
        only_t1 = sum(!cells_t1 %in% cells_t0)
)

# Extract both periods in exactly the same cell order.
comm_t0 <- as.matrix(
        comm_mat[rows_t0[match(common_cells, cells_t0)], , drop = FALSE]
)

comm_t1 <- as.matrix(
        comm_mat[rows_t1[match(common_cells, cells_t1)], , drop = FALSE]
)

rownames(comm_t0) <- common_cells
rownames(comm_t1) <- common_cells

# Store binary data as integers to reduce memory usage.
storage.mode(comm_t0) <- "integer"
storage.mode(comm_t1) <- "integer"

# Verify alignment and count missing habitat records.
stopifnot(identical(dimnames(comm_t0), dimnames(comm_t1)))

dim(comm_t0)

c(
        missing_t0 = sum(is.na(comm_t0)),
        missing_t1 = sum(is.na(comm_t1))
)

#### Build the habitat-level loss dataset ####

# Count occupied cells at baseline for each habitat.
occupancy_t0 <- colSums(comm_t0)

# Locate baseline presences without converting the full matrices to long format.
presence_index <- which(comm_t0 == 1L, arr.ind = TRUE)

# Extract follow-up status for the same habitat–cell combinations.
h3_data <- data.frame(
        cell_id = rownames(comm_t0)[presence_index[, 1]],
        habitat = colnames(comm_t0)[presence_index[, 2]],
        loss = 1L - comm_t1[presence_index],
        n_cells_t0 = unname(occupancy_t0[presence_index[, 2]]),
        stringsAsFactors = FALSE
)

# Log-transform baseline occupancy while retaining the original cell counts.
h3_data$log_n_cells_t0 <- log(h3_data$n_cells_t0)

# Verify that every baseline presence has exactly one outcome.
stopifnot(
        nrow(h3_data) == sum(occupancy_t0),
        all(h3_data$loss %in% c(0L, 1L)),
        all(h3_data$n_cells_t0 > 0L),
        !anyDuplicated(h3_data[c("cell_id", "habitat")])
)

# Remove the temporary index to release memory.
rm(presence_index)

##### Inspect habitat loss records ####

head(h3_data)

c(
        observations = nrow(h3_data),
        habitats = length(unique(h3_data$habitat)),
        cells = length(unique(h3_data$cell_id)),
        lost = sum(h3_data$loss == 1L),
        maintained = sum(h3_data$loss == 0L)
)

# Descriptive counts above are recalculated from the current input.

##### Add protected-area coverage ####

# Extract cell-level coverage and coordinates without polygon geometries.
pa_data <- grid_pa |>
        sf::st_drop_geometry() |>
        dplyr::transmute(
                cell_id = CellCode,
                total_pa_cov,
                total_pa_land_cov,
                land_cov,
                x_km = (EofOrigin + 5000) / 1000,
                y_km = (NofOrigin + 5000) / 1000
        )

# Join protection metrics to each habitat–cell occurrence.
h3_data <- h3_data |>
        dplyr::left_join(pa_data, by = "cell_id")

# Inspect h3_data
glimpse(h3_data)

# Inspect coverage values before defining the model dataset.
summary(h3_data[c("total_pa_cov", "total_pa_land_cov", "land_cov")])

##### Prepare model variables ####

# Express protection in units of 10 percentage points.
h3_data$pa_10 <- h3_data$total_pa_cov / 10

# Centre log occupancy using one value per habitat.
# This avoids giving widespread habitats more weight in the centring constant.
occupancy_reference <- mean(
        h3_data$log_n_cells_t0[!duplicated(h3_data$habitat)]
)

h3_data$log_occupancy_c <- h3_data$log_n_cells_t0 -
        occupancy_reference

# Encode grouping variables as factors.
h3_data$habitat <- factor(h3_data$habitat)
h3_data$cell_id <- factor(h3_data$cell_id)

#### Fit models ####

##### Additive model ####

# Allow baseline loss and protection slopes to vary among habitats.
# Include a cell intercept to account for habitats sharing the same cell.
h3_additive <- glmmTMB(
        loss ~ pa_10 + log_occupancy_c +
                (1 + pa_10 | habitat) +
                (1 | cell_id),
        family = binomial(link = "logit"),
        data = h3_data,
        na.action = na.fail
)


##### Interaction model ####

# Test whether the protection slope varies systematically with occupancy.
# Retain the same random-effects structure as in the additive model.
h3_interaction <- glmmTMB(
        loss ~ pa_10 * log_occupancy_c +
                (1 + pa_10 | habitat) +
                (1 | cell_id),
        family = binomial(link = "logit"),
        data = h3_data,
        na.action = na.fail
)


###### Inspect model estimates and convergence ####

summary(h3_interaction)

# Both convergence codes should be 0 and both Hessian checks TRUE.
c(
        additive_convergence = h3_additive$fit$convergence,
        interaction_convergence = h3_interaction$fit$convergence
)

c(
        additive_hessian = h3_additive$sdr$pdHess,
        interaction_hessian = h3_interaction$sdr$pdHess
)


###### Test the protection–occupancy interaction ####

# Compare nested models using a likelihood-ratio test.
# Interpret this comparison only if both models converged successfully.
anova(h3_additive, h3_interaction)

###### model export #####

dir.create("processed/models", recursive = TRUE, showWarnings = FALSE)

# Save fitted models using moderate compression.
saveRDS(
        h3_additive,
        "processed/models/h3_additive.rds",
        compress = "gzip"
)

saveRDS(
        h3_interaction,
        "processed/models/h3_interaction.rds",
        compress = "gzip"
)

# Preserve the centring constant required for subsequent predictions.
saveRDS(
        occupancy_reference,
        "processed/models/h3_occupancy_reference.rds"
)

# Save the exact analysis data used for figures and descriptive tables.
saveRDS(
        h3_data,
        "processed/models/h3_analysis_data.rds",
        compress = "gzip"
)
