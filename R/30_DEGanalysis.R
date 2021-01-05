rnaseq_counts <- read.delim("/Users/senosam/Documents/Massion_lab/Others/Dalton/data/RNAseq/BEAS2B/RnaSeq_Dalton4676.count.txt")    
ref = data.frame(cbind(sample_ID = c("OverExpression1", "OverExpression2", "OverExpression3", "WT1", "WT2", "WT3"), 
        Batch = c(1,2,3,1,2,3), 
        condition = c('OE', 'OE', 'OE', 'WT', 'WT', 'WT')))



environment_set <- function(){
    library(EnhancedVolcano)
    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(DESeq2)
    library(sva)
    library(limma)
    library(ComplexHeatmap)
    library(fgsea)
    library(tibble)
    library(gage)
    library(ggplot2)
    library(org.Hs.eg.db)
    library(clusterProfiler)
    library(biomaRt)
    # Check BiocManager::valid()
}

preprocess_rna <- function(rna_all, 
    ref, 
    correct_batch=TRUE, 
    lowvargenesrm = TRUE){

    # Remove low variance genes from counts
    if(lowvargenesrm){
        variances <- apply(rna_all[, 8:ncol(rna_all)], 1, var)
        sdv <- apply(rna_all[, 8:ncol(rna_all)], 1, sd)
        q1 <- quantile(variances, na.rm = T)["25%"]
        idx <- which(is.na(variances) | variances <= q1 | sdv == 0)
        rna_all <- rna_all[-idx, ]
    }
    
    counts_all <- rna_all[, 8:ncol(rna_all)]
    rownames(counts_all) <- rna_all$Feature

    # Normalize
    ref$Batch <- as.factor(ref$Batch)
    dds <- DESeqDataSetFromMatrix(
      countData = counts_all,
      colData = ref,
      design = ~Batch, tidy = F  
    )
    dds <- estimateSizeFactors(dds)
    sizeFactors(dds)
    vsd <- vst(dds)


    pbatch_bf <- plotPCA(vsd, intgroup=c('condition', 'Batch'), returnData = TRUE)
    #pbatch_bf <- plotPCA(vsd, "condition") + labs(fill = "condition") + ggtitle("Batch raw")
    pbatch_bf$condition <- factor(pbatch_bf$condition)
    pbatch_bf$Batch <- factor(pbatch_bf$Batch)
    percentVar <- round(100 * attr(pbatch_bf, "percentVar")) 
    pbatch_bf <- ggplot(pbatch_bf, aes(x = PC1, y = PC2, color = condition, shape = Batch)) + 
      geom_point(size =3, aes(fill=condition)) + 
      scale_shape_manual(values=c(21,22,23)) + 
      xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
      ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
      ggtitle("PCA before batch correction") +
      labs(color='condition', shape='Batch')
    
    # Remove batch effect
    if(correct_batch){
        assay(vsd) <- limma::removeBatchEffect(assay(vsd), 
            batch=vsd$Batch)
        pbatch_af <- plotPCA(vsd, intgroup=c('condition', 'Batch'), returnData = TRUE)
        #pbatch_af <- plotPCA(vsd, "condition") + labs(fill = "condition") + ggtitle("Batch after BE removal")
        pbatch_af$condition <- factor(pbatch_af$condition)
        pbatch_af$Batch <- factor(pbatch_af$Batch)
        percentVar <- round(100 * attr(pbatch_af, "percentVar")) 
        pbatch_af <- ggplot(pbatch_af, aes(x = PC1, y = PC2, color = condition, shape = Batch)) + 
          geom_point(size =3, aes(fill=condition)) + 
          scale_shape_manual(values=c(21,22,23)) + 
          xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
          ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
          ggtitle("PCA after batch correction") +
          labs(color='condition', shape='Batch')
    }
    vsd_mat <- assay(vsd)

    # Save pre-processed data in a list
    dge_preprocessed <- list(ref=ref, rna_all=rna_all, counts_all=counts_all, 
        vsd_mat=vsd_mat, pbatch_bf=pbatch_bf, pbatch_af=pbatch_af)

    dge_preprocessed
}


DE_analysis <- function(ls_preprocessed, 
    RefBased=FALSE,
    NewCondition=FALSE,
    NewCondition_df, # condition label into ref
    cond_nm, # gene or column name
    two_levels=c('WT','OE'), # high low, dead alive...
    reference,
    batch_cor = T
    ){

    # Unlist
    ref=ls_preprocessed$ref
    rna_all=ls_preprocessed$rna_all
    counts_all=ls_preprocessed$counts_all
    vsd_mat=ls_preprocessed$vsd_mat

    message('Unlist done')

    
    if(RefBased){
        meta_data <- ref
        meta_data$sample_ID <- as.character(meta_data$sample_ID)
        names(meta_data)[names(meta_data) == cond_nm] <- 'Condition'
        meta_data <- meta_data %>% filter(Condition %in% two_levels)
    }

    if(NewCondition){
        meta_data <- NewCondition_df
        names(meta_data)[names(meta_data) == cond_nm] <- 'Condition'
        meta_data$Condition <- as.character(meta_data$Condition)
        meta_data <- meta_data %>% filter(Condition %in% two_levels)
    }


    counts_all <- as.matrix(counts_all[,meta_data$sample_ID])
    ref <- ref %>% filter(sample_ID %in% meta_data$sample_ID)
    vsd_mat <- vsd_mat[,meta_data$sample_ID]

    message('Filtering done')

    if(batch_cor) {
        dds <- DESeqDataSetFromMatrix(
        countData = counts_all,
        colData = meta_data,
        design = ~Batch + Condition, tidy = F)
    } else{
        dds <- DESeqDataSetFromMatrix(
        countData = counts_all,
        colData = meta_data,
        design = ~Condition, tidy = F)
    }
    message('Design done')

    dds <- estimateSizeFactors(dds)

    # Generate a transformed matrix with gene symbols
    ens2symbol <- data.frame(cbind(ENSEMBL=as.character(rna_all$Feature), 
        symbol=as.character(rna_all$Feature_gene_name)))
    vsd_mat_sym <- data.frame(vsd_mat) %>% mutate(gene=rownames(.)) %>% as_tibble()
    vsd_mat_sym <- inner_join(vsd_mat_sym, ens2symbol, by=c("gene"="ENSEMBL"))
    rownames(vsd_mat_sym) <- make.names(vsd_mat_sym$symbol, unique=TRUE)

    message('vsd symbols done')

    dds$Condition <- relevel(dds$Condition, ref = reference)
    dds <- DESeq(dds, parallel = F)

    message('DESeq done')

    res <- results(dds, contrast = c('Condition', two_levels))
    res <- lfcShrink(dds, contrast = c('Condition',two_levels),
        res=res, type = 'normal')
    res_df <- data.frame(res) %>% mutate(gene=rownames(.)) %>% as_tibble()

    res_df <- res_df %>% 
      inner_join(., ens2symbol, by=c("gene"="ENSEMBL")) %>%
      data.frame()
    rownames(res_df) <- make.names(res_df$symbol, unique=TRUE)

    message('res symbols done')

    DE_res <- list(dds=dds, vsd_mat=vsd_mat, vsd_mat_sym=vsd_mat_sym, 
        ens2symbol=ens2symbol, res=res, res_df=res_df, 
        ref=ref, meta_data=meta_data)

    message('list done')

    DE_res
}

fgsea_analysis <- function(DE_res){

    # Unlist
    res_df=DE_res$res_df

    ranks <- res_df %>%
      dplyr::select(symbol,stat) %>%
      na.omit() %>% 
      distinct() %>% 
      group_by(symbol) %>% 
      summarize(stat=mean(stat))

    ranks <- deframe(ranks)

    hallmark_path <- '/Users/senosam/Documents/Massion_lab/RNASeq_summary/GSEA/msigdb_v7.1_GMTs/h.all.v7.1.symbols.gmt'
    c1_path <- '/Users/senosam/Documents/Massion_lab/RNASeq_summary/GSEA/msigdb_v7.1_GMTs/c1.all.v7.1.symbols.gmt'
    c2_path <- '/Users/senosam/Documents/Massion_lab/RNASeq_summary/GSEA/msigdb_v7.1_GMTs/c2.all.v7.1.symbols.gmt'
    c3_path <- '/Users/senosam/Documents/Massion_lab/RNASeq_summary/GSEA/msigdb_v7.1_GMTs/c3.all.v7.1.symbols.gmt'
    c4_path <- '/Users/senosam/Documents/Massion_lab/RNASeq_summary/GSEA/msigdb_v7.1_GMTs/c4.all.v7.1.symbols.gmt'
    c5_path <- '/Users/senosam/Documents/Massion_lab/RNASeq_summary/GSEA/msigdb_v7.1_GMTs/c5.all.v7.1.symbols.gmt'
    c6_path <- '/Users/senosam/Documents/Massion_lab/RNASeq_summary/GSEA/msigdb_v7.1_GMTs/c6.all.v7.1.symbols.gmt'
    c7_path <- '/Users/senosam/Documents/Massion_lab/RNASeq_summary/GSEA/msigdb_v7.1_GMTs/c7.all.v7.1.symbols.gmt'
    msig_path <- '/Users/senosam/Documents/Massion_lab/RNASeq_summary/GSEA/msigdb_v7.1_GMTs/msigdb.v7.1.symbols.gmt'

    fgsea_fixed <- function(pthw_path){
        res <- fgsea(pathways=gmtPathways(pthw_path), stats=ranks, nperm=1000) %>%
                    as_tibble %>%
                    arrange(padj, desc(abs(NES)))
        res$state <- ifelse(res$NES > 0, "up", "down")
        res$leadingEdge <- sapply(res$leadingEdge, . %>% {
          str_c(., collapse = " ")})

        res
    }

    res_hm <- fgsea_fixed(pthw_path=hallmark_path)
    res_c1 <- fgsea_fixed(pthw_path=c1_path)
    res_c2 <- fgsea_fixed(pthw_path=c2_path)
    res_c3 <- fgsea_fixed(pthw_path=c3_path)
    res_c4 <- fgsea_fixed(pthw_path=c4_path)
    res_c5 <- fgsea_fixed(pthw_path=c5_path)
    res_c6 <- fgsea_fixed(pthw_path=c6_path)
    res_c7 <- fgsea_fixed(pthw_path=c7_path)
    res_msg <- fgsea_fixed(pthw_path=msig_path)

    fgsea_res <- list(res_hm=res_hm, res_c1=res_c1, res_c2=res_c2, 
        res_c3=res_c3, res_c4=res_c4, res_c5=res_c5, res_c6=res_c6,
        res_c7=res_c7, res_msg=res_msg)

    fgsea_res
}

# plots
heatmap_200 <- function(res_df, vsd_mat, meta_data, pData_rnaseq, n_genes=200,
    padj_cutoff = 0.05, l2fc_cutoff=1.5, ha_custom=NULL, row_km=2){
    res_df <- data.frame(res_df) %>%
        filter(abs(log2FoldChange) > l2fc_cutoff) %>%
        filter(padj < padj_cutoff) %>%
        arrange(padj)

    if(nrow(res_df)>n_genes){
        gene_list <- res_df %>%
            head(n_genes) %>%
            dplyr::select(gene) %>%
            pull()
    }else{
        gene_list <- res_df %>%
            dplyr::select(gene) %>%
            pull()        
    }

    # Use normalized data (vsd_mat) for plot
    filtered_res <- data.frame(vsd_mat) %>%
      filter(gene %in% gene_list) %>%
      as.data.frame()
    rownames(filtered_res) <- make.names(filtered_res$symbol, unique=TRUE)
    filtered_res$gene <- NULL
    filtered_res$symbol <- NULL

    ha = HeatmapAnnotation(

        CANARY = as.factor(pData_rnaseq$CANARY),
        #gender = as.factor(pData_rnaseq$Gender),
        Condition = as.factor(meta_data$Condition),

        simple_anno_size = unit(0.5, "cm")
    )

    if(!is.null(ha_custom)){
        ha <- ha_custom
    }

    print(Heatmap(as.matrix(filtered_res), name = "mat", 
        #column_km = 2, 
      #row_km = row_km,
      column_split =as.factor(meta_data$Condition),
      heatmap_legend_param = list(color_bar = "continuous"), 
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 8), top_annotation = ha))
}

volcano_plot <- function(res_df, gene=NULL, p_title, pCutoff=0.001, FCcutoff=1.5){

    if (!is.null(gene)){
        k <- which(res_df$gene %in% gene)
        res_df <- res_df[-k,]
    }
    print(EnhancedVolcano(res_df,
        lab = rownames(res_df),
        x = 'log2FoldChange',
        y = 'pvalue',
        pCutoff=pCutoff,
        FCcutoff=FCcutoff,
        xlim = c(-5,5),
        pointSize=2,
        labSize=4,
        title = p_title))
}

fgsea_plot <- function(fgsea_res, pathways_title, cutoff = 0.05, 
    max_pathways = 30, condition_name){

        color_levels <- function(fgsea_res) {
            colors <- c()
            if (any(fgsea_res$state == "down")) {
              colors <- c(colors, "lightblue")
            }
            if (any(fgsea_res$state == "up")) {
              colors <- c(colors, "#DC143C")
            }
            colors
        }


        if (!is.null(cutoff)) {
            fgsea_res <- fgsea_res %>% filter(padj < cutoff)
        }

        curated_pathways <- fgsea_res %>%
            dplyr::slice(1:max_pathways)
        curated_pathways['leadingEdge'] <- NULL
        print(ggplot(curated_pathways, aes(reorder(pathway, NES), NES)) +
            geom_col(aes(fill = state), width = 0.5, color = "black") +
            scale_size_manual(values = c(0, 1), guide = "none") +
            geom_label(aes(label = round(padj, 4)), size = 3) +
            coord_flip() +
            labs(
                x = 'Pathway', 
                y = "Normalized Enrichment Score",
                title = str_c(pathways_title, " pathways: ", condition_name),
                subtitle = str_c("(Cutoff: p.adj <", cutoff, ")")
            ) +
            theme_bw() +
            scale_fill_manual(values = color_levels(curated_pathways)))

        curated_pathways
}

