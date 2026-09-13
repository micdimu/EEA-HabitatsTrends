library(tidyverse)
library(sf)
library(patchwork)
library(sandwich)
library(lmtest)

source("SourceCodes/Source05_H2.R")

#### Load data ####

# Run from the project root. Keep double absences (0 -> 0): their
# beta_jaccard is undefined (NA), but their gain/loss counts are valid zeros.
# Filter missing values only for variables required by each analysis.
temporal_beta_group <- read.csv("processed/temporal_group_cell_by_cell.csv", check.names = FALSE) |>
        mutate(net_change = richness_2013_2018 - richness_2007_2012)

grid_pa <- st_read("data/PAs/derived/europe_10km_protected_area_coverage.gpkg")


grid_group <- grid_pa |>
        right_join(temporal_beta_group, by = c("CellCode" = "cell_id")) |> 
        mutate(taxGroup = as.factor(taxGroup))

#### H2a: habitat loss and gain-loss balance ####
# Counts refer to habitat types within cells, not habitat area.
# Associations between group changes do not establish direct replacement.


# 1. Prepare data and summarise ALL habitat groups ----------------------

dat <- sf::st_drop_geometry(grid_group) |>
        filter(
                !is.na(CellCode),
                !is.na(taxGroup),
                is.finite(total_pa_land_cov),
                is.finite(land_cov),
                land_cov > 0,
                !is.na(lost_species),
                !is.na(shared_species),
                !is.na(gained_species)
        ) |>
        mutate(
                land_weight = land_cov / 100,
                initial = lost_species + shared_species,
                final = gained_species + shared_species,
                changes = lost_species + gained_species
        )

# Each row must represent a unique cell-by-group combination.
stopifnot(
        anyDuplicated(dat[c("CellCode", "taxGroup")]) == 0L,
        all(dat$lost_species >= 0),
        all(dat$shared_species >= 0),
        all(dat$gained_species >= 0)
)

# Occurrences count habitat-by-cell presences, not cells occupied by a group.
summary_groups <- dat |>
        group_by(taxGroup) |>
        summarise(
                n_cells = n_distinct(CellCode),
                initial_occ = sum(initial),
                final_occ = sum(final),
                lost_occ = sum(lost_species),
                gained_occ = sum(gained_species),
                
                # Habitat-occurrence totals weighted by the terrestrial cell fraction.
                initial_w = sum(initial * land_weight),
                lost_w = sum(lost_species * land_weight),
                gained_w = sum(gained_species * land_weight),
                .groups = "drop"
        ) |>
        mutate(
                net_occ = gained_occ - lost_occ,
                loss_pct = if_else(
                        initial_occ > 0, 100 * lost_occ / initial_occ, NA_real_
                ),
                net_pct = if_else(
                        initial_occ > 0, 100 * net_occ / initial_occ, NA_real_
                ),
                loss_pct_weighted = if_else(
                        initial_w > 0, 100 * lost_w / initial_w, NA_real_
                ),
                net_pct_weighted = if_else(
                        initial_w > 0,
                        100 * (gained_w - lost_w) / initial_w,
                        NA_real_
                )
        )

summary_groups


# 2. Select the three groups for the main analysis ----------------------

groups_h2 <- c("Grasslands", "Forests", "Sclerophyllous scrub")

h2 <- dat |>
        filter(taxGroup %in% groups_h2) |>
        mutate(taxGroup = factor(taxGroup, levels = groups_h2))

# Grasslands is the reference level. Loss is defined only when at least
# one habitat of the group was present initially.
d_loss <- h2 |> filter(initial > 0)

# Gain-loss balance is defined only when at least one change occurred.
# Colonisations from initial absence remain eligible for this analysis.
d_balance <- h2 |> filter(changes > 0)


# 3. Fit grouped quasibinomial models -----------------------------------
# Loss models lost / (lost + shared); balance models gained / (gained + lost).
# The two-column responses supply their denominators automatically.
# land_weight supplies an additional weight based on terrestrial coverage.

m_loss <- glm(
        cbind(lost_species, shared_species) ~
                total_pa_land_cov * taxGroup,
        family = quasibinomial(),
        weights = land_weight,
        data = d_loss
)

m_balance <- glm(
        cbind(gained_species, lost_species) ~
                total_pa_land_cov * taxGroup,
        family = quasibinomial(),
        weights = land_weight,
        data = d_balance
)

V_loss <- vcovCL(
        m_loss, cluster = d_loss$CellCode, type = "HC1"
)

V_balance <- vcovCL(
        m_balance, cluster = d_balance$CellCode, type = "HC1"
)

# Cluster-robust inference accounts for multiple groups in the same cell.
# It does not account for spatial dependence between neighbouring cells.
coeftest(m_loss, vcov. = V_loss)
coeftest(m_balance, vcov. = V_balance)

# Joint Wald tests of the group-by-coverage interactions.
m_loss_add <- update(
        m_loss, . ~ total_pa_land_cov + taxGroup
)

m_balance_add <- update(
        m_balance, . ~ total_pa_land_cov + taxGroup
)

waldtest(
        m_loss_add, m_loss,
        vcov = V_loss, test = "Chisq"
)

waldtest(
        m_balance_add, m_balance,
        vcov = V_balance, test = "Chisq"
)


# 4. Predict with pointwise 95% cell-cluster-robust confidence intervals --



pred_loss <- predict_h2(m_loss, V_loss, d_loss)
pred_balance <- predict_h2(m_balance, V_balance, d_balance)



p_loss <- plot_h2(pred_loss, "Habitat loss probability")

# Transform balance to (gained - lost) / (gained + lost), ranging from -1 to 1.
# This is not net change relative to initial richness.
p_balance <- plot_h2(
        pred_balance |>
                mutate(
                        estimate = 2 * estimate - 1,
                        lower = 2 * lower - 1,
                        upper = 2 * upper - 1
                ),
        "Habitat gain–loss balance"
) +
        geom_hline(
                yintercept = 0, linetype = "dashed",
                colour = "grey45"
        ) +
        scale_y_continuous(
                limits = c(-1, 1),
                breaks = seq(-1, 1, 0.5),
                labels = scales::label_number()
        )

p_loss + p_balance


