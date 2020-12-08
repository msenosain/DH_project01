# Select top n metabolites
n <- 50
top_met <- data_DH[1:n,]

ha = HeatmapAnnotation(
    condition = c('OE', 'WT'),
    simple_anno_size = unit(0.5, "cm"),
    show_legend = FALSE,
    show_annotation_name = FALSE
)

row_ha = HeatmapAnnotation(
    Subtype=top_met$subtype,
    simple_anno_size = unit(0.5, "cm"),
    which = 'row',
    show_annotation_name = FALSE
)

Heatmap(scale(top_met[,3:4]), name = "z-score",
        heatmap_legend_param = list(color_bar = "continuous"), 
        show_row_dend = T,
        show_column_dend = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_split = top_met$subtype,
        column_split = c('OE', 'WT'),
        right_annotation = row_ha,
        top_annotation = ha, 
        row_title=NULL)


max_metabolites <- 50
curated_metabolites <- data_DH %>%
            dplyr::slice(1:max_metabolites)
#curated_pathways['leadingEdge'] <- NULL
ggplot(curated_metabolites, aes(reorder(metabolite, L2R), L2R)) +
    geom_col(aes(fill = state), width = 0.5, color = "black") +
    scale_size_manual(values = c(0, 1), guide = "none") +
    #geom_label(aes(label = round(padj, 4)), size = 3) +
    coord_flip() +
    labs(
        x = 'Metabolite', 
        y = "Log2 fold-change"
        #title = str_c(pathways_title, " pathways: ", condition_name),
        #subtitle = str_c("(Cutoff: p.adj <", cutoff, ")")
    ) +
    theme_bw() 
    #scale_fill_manual(values = color_levels(curated_metabolites))