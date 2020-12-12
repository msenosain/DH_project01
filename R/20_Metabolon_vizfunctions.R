set_environment <- function(){
    library(dplyr)
    library(tidyverse)
    library(tidyr)
    library(ggplot2)
    library(ComplexHeatmap)
}

dh_heatmap <- function(dt, col_names, scale = TRUE){

    # HM annotations
    col_ha = HeatmapAnnotation(
    condition = col_names,
    simple_anno_size = unit(0.5, "cm"),
    show_legend = FALSE,
    show_annotation_name = FALSE
    )

    row_ha = HeatmapAnnotation(
        Subtype=dt$Super_Pathway,
        simple_anno_size = unit(0.5, "cm"),
        which = 'row',
        show_annotation_name = FALSE
    )

    # Get index for cols of interest
    idx <- which(colnames(dt) %in% col_names)

    # Scaling
    if(scale){
        x <- scale(dt[,idx])
        name_l <- "z-score"
    } else {
        x <- dt[,idx]
        name_l <- "value"
    }
    # Plot heatmap
    print(Heatmap(x, 
        name = name_l,
        heatmap_legend_param = list(color_bar = "continuous"), 
        show_row_dend = T,
        show_column_dend = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_split = dt$Super_Pathway,
        column_split = col_names,
        right_annotation = row_ha,
        top_annotation = col_ha, 
        row_title=NULL))
}



L2FC_barplot <- function(dt, 
    max_items = 20, # metabolites or subpathways
    superpathway, 
    state_col, 
    L2R_col,
    sum_subpathway = FALSE) {

    # Color function
    color_levels <- function(dt, state_col) {
        colors <- c()
        if (any(dt[,state_col] == "DN")) {
          colors <- c(colors, "lightblue")
        }
        if (any(dt[,state_col] == "UP")) {
          colors <- c(colors, "#DC143C")
        }
        colors
    }

    colnames(dt)[which(colnames(dt) == L2R_col)] <- 'L2R'
    colnames(dt)[which(colnames(dt) == state_col)] <- 'state'

    if(sum_subpathway){
        plot_title <- paste0("Top ", max_items,
            " deregulated sub-pathways: ", superpathway)
        x_lab = 'Sub-pathway'
        colnames(dt)[which(colnames(dt) == 'Sub_Pathway')] <- 'item'

    } else {
        plot_title <- paste0("Top ", max_items,
            " deregulated metabolites: ", superpathway)
        x_lab = 'Metabolite'
        colnames(dt)[which(colnames(dt) == 'Metabolite')] <- 'item'
    }
    # Select top metabolites or subpathways per superpathway
    dt <- dt %>% 
        filter(Super_Pathway==superpathway) %>%
        arrange(desc(abs(L2R))) %>%
        dplyr::slice(1:max_items)

    
    print(ggplot(dt, aes(reorder(item, L2R), L2R)) +
            geom_col(aes(fill = state), width = 0.5, color = "black") +
            scale_size_manual(values = c(0, 1), guide = "none") +
            coord_flip() +
            labs(
                x = x_lab, 
                y = "Log2 fold-change",
                title = plot_title
            ) +
            theme_bw() +
            scale_fill_manual(values = color_levels(dt))
    )

}
