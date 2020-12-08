set_environment <- function(){
    library(dplyr)
    library(tidyverse)
    library(tidyr)
    library(ggplot2)
    library(ComplexHeatmap)
}

dh_heatmap <- function(dt, col_names){

    # HM annotations
    col_ha = HeatmapAnnotation(
    condition = col_names,
    simple_anno_size = unit(0.5, "cm"),
    show_legend = FALSE,
    show_annotation_name = FALSE
    )

    row_ha = HeatmapAnnotation(
        Subtype=dt$subtype,
        simple_anno_size = unit(0.5, "cm"),
        which = 'row',
        show_annotation_name = FALSE
    )

    # Get index for cols of interest
    idx <- which(colnames(dt) %in% col_names)

    # Plot heatmap
    print(Heatmap(scale(dt[,idx]), 
        name = "z-score",
        heatmap_legend_param = list(color_bar = "continuous"), 
        show_row_dend = T,
        show_column_dend = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_split = dt$subtype,
        column_split = col_names,
        right_annotation = row_ha,
        top_annotation = col_ha, 
        row_title=NULL))
}



L2FC_barplot <- function(dt, max_metabolites = 20, msubtype) {

    # Color function
    color_levels <- function(dt) {
        colors <- c()
        if (any(dt$state == "DN")) {
          colors <- c(colors, "lightblue")
        }
        if (any(dt$state == "UP")) {
          colors <- c(colors, "#DC143C")
        }
        colors
    }

    # Select top metabolites by subtype
    dt <- dt %>% 
            filter(subtype==msubtype) %>%
            arrange(desc(abs(L2R))) %>%
            dplyr::slice(1:max_metabolites)

    print(ggplot(dt, aes(reorder(metabolite, L2R), L2R)) +
            geom_col(aes(fill = state), width = 0.5, color = "black") +
            scale_size_manual(values = c(0, 1), guide = "none") +
            coord_flip() +
            labs(
                x = 'Metabolite', 
                y = "Log2 fold-change",
                title = paste0("Top ", max_metabolites, 
                    " deregulated metabolites: ", msubtype)
            ) +
            theme_bw() +
            scale_fill_manual(values = color_levels(dt))
    )

}


