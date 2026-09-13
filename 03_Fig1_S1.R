library(tidyverse)
library(patchwork)
library(sf)
library(biscale)
library(ggspatial)

source("SourceCodes/Source03_Fig1_S1.R")

#### Load data ####

temporal_beta <- read.csv("processed/temporal_beta_cell_by_cell.csv") |> 
        mutate(net_change = richness_2013_2018 - richness_2007_2012) 
        

grid_pa <- st_read("data/PAs/derived/europe_10km_protected_area_coverage_2012EU27.gpkg")

grid_T <- grid_pa |>
        filter(land_cov > 0) |>
        right_join(temporal_beta, by = c("CellCode" = "cell_id")) 



temporal_beta$cell_id[(!(temporal_beta$cell_id %in% grid_pa$CellCode))]
all((temporal_beta$cell_id %in% grid_pa$CellCode))
sum(is.na(grid_T$richness_2007_2012))

#### Figure 1 Bivariate map #####

beta_pa_map <- bivariate_map(
        data = grid_T,
        geometry = grid_T,
        x = total_pa_cov,
        y = beta_jaccard,
        x_breaks = c(10, 30, 70),
        y_breaks = c(0.25, 0.50, 0.75),
        xlab = "Protected areas (%)",
        ylab = "Temporal Beta diversity",
        pal = "DkBlue2",
)

beta_pa_map$plot

layout <- c(
        patchwork::area(t = 0, l = 0, b = 12, r = 12),
        patchwork::area(t = 2.5, l = 1.8, b = 3.5, r = 2.8)
)

figure_1 <- beta_pa_map$map + beta_pa_map$legend + plot_layout(design = layout)

figure_1

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

ggsave(
        filename = file.path(
                "output",
                "figures",
                "Figure_1.tiff"
        ),
        plot = figure_1,
        width = 180,
        height = 165,
        units = "mm",
        dpi = 600,
        device = "tiff",
        compression = "lzw",
        bg = "white"
)

##### fig 1 supplementary ####

temporal_beta |> 
        colnames()

plot_specs <- tribble(
        ~variable,             ~x_label,                               ~fill,      ~zero_line, ~x_limits,
        "beta_jaccard",        "Beta diversity (Jaccard)",              "#5B7083",  FALSE,      list(c(0, 1)),
        "richness_2007_2012",  "Habitat richness, 2007–2012",          "#849EAD",  FALSE,      list(NULL),
        "richness_2013_2018",  "Habitat richness, 2013–2018",          "#849EAD",  FALSE,      list(NULL),
        "shared_species",      "Shared Habitat",                       "#9A9A9A",  FALSE,      list(NULL),
        "lost_species",        "Lost Habitat",                         "#B27870",  FALSE,      list(NULL),
        "gained_species",      "Gained Habitat",                       "#6F9278",  FALSE,      list(NULL),
        "net_change",          "Net change in Habitat richness",       "#AA8E61",  TRUE,       list(NULL)
)

plots <- pmap(
        plot_specs,
        \(variable, x_label, fill, zero_line, x_limits) {
                density_panel(
                        data = temporal_beta,
                        variable = variable,
                        x_label = x_label,
                        fill = fill,
                        zero_line = zero_line,
                        x_limits = x_limits
                )
        }
) 



design <- "
AAAAAABBBBBB
CCCCCCDDDDDD
EEEEFFFFGGGG
"

fig_S1 <- wrap_plots( plots[c(2,3,1,7,4,5,6)], ncol = 3 ) +
        plot_layout(design = design)+
        plot_annotation(tag_levels = "a") &
        theme(plot.tag = element_text(face = "bold", size = 11),
              plot.tag.position = c(0.02, 0.98))

fig_S1

ggsave(
        filename = "output/supplementary/Figure_S1.tiff",
        plot = fig_S1,
        device = "tiff",
        width = 260,
        height = 240,
        units = "mm",
        dpi = 600,
        compression = "lzw",
        bg = "white"
)
