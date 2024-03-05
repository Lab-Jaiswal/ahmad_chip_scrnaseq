library(BiocParallel)
library(BiocSingular)
library(tidySingleCellExperiment)

library(scran)
library(scuttle)
library(igraph)
library(corral)
library(PCAtools)
library(bluster)
library(batchelor)
library(RNOmni)
library(showtext)

library(readxl)
library(magrittr)
library(tidyverse)

iterative_filter <- function(sce_object){
    low_counts <- isOutlier(sce_object$sum, log = TRUE, type = "lower")
    low_features <- isOutlier(sce_object$detected, log = TRUE, type = "lower")
    high_mt <- isOutlier(sce_object$subsets_percent_mt_percent, type = "higher")
    excluded_cells <- low_counts | low_features | high_mt
    n_outlier <- which(excluded_cells) %>% length
    print(n_outlier)
    if(n_outlier >= floor(0.01 * ncol(sce_object))) {
        sce_object_filtered <- sce_object[,!excluded_cells]
        iterative_filter(sce_object_filtered)
    } else {
        return(sce_object)
    }
}

filter_osca <- function(sce_object, gene_annot = NULL, filter_cells = T, iterative_filter = F) {
    # If gene_annot has been passed in, join it to the rowData of the SCE object
    if (!is.null(gene_annot)) {
        rowData(sce_object) %<>% as.data.frame %>% left_join(gene_annot) %>% DataFrame
    }

    # Filter cells that are outliers based on QC metrics
    if(filter_cells == T) {
        # Compute QC metrics
        mito_genes <- which(rowData(sce_object)$SEQNAME == "MT")
        sce_object_qc <- scuttle::addPerCellQC(sce_object, subset = list(percent_mt = mito_genes))
        colData(sce_object_qc)$log_sum <- log2(sce_object_qc$sum)
        colData(sce_object_qc)$log_detected <- log2(sce_object_qc$detected)
        colData(sce_object_qc)$Cell_ID <- str_c(sce_object_qc$Sample, sce_object_qc$Barcode, sep = "_")
        colnames(sce_object_qc) <- sce_object_qc$Cell_ID

        if (iterative_filter == T) {
            sce_object_filtered <- iterative_filter(sce_object)
        } else {
            low_counts <- scuttle::isOutlier(sce_object_qc$sum, log = TRUE, type = "lower")
            low_features <- scuttle::isOutlier(sce_object_qc$detected, log = TRUE, type = "lower")
            high_mt <- scuttle::isOutlier(sce_object_qc$subsets_percent_mt_percent, type = "higher")
            excluded_cells <- low_counts | low_features | high_mt

            sce_object_filtered <- dplyr::filter(sce_object_qc, !excluded_cells)
        }

        # Set rownames to be gene symbols instead of Ensembl IDs, and remove genes with duplicate symbols
        symbol_vector <- rowData(sce_object_filtered)$Symbol
        duplicated_symbols <- duplicated(symbol_vector) | duplicated(symbol_vector, fromLast = TRUE)
        sce_object_filtered <- sce_object_filtered[!duplicated_symbols, ]
        rownames(sce_object_filtered) <- rowData(sce_object_filtered)$Symbol
    } else {
        sce_object_filtered <- sce_object
    }

    # Subset genes to only include ones which meet our variability cutoff
    set.seed(12345)
    sce_object_gene_var <- scran::modelGeneVarByPoisson(sce_object_filtered)
    sce_object_gene_var_top <- scran::getTopHVGs(sce_object_gene_var)
    sce_object_gene_cv2 <- scran::modelGeneCV2(sce_object_filtered)
    sce_object_gene_cv2_top <- scran::getTopHVGs(sce_object_gene_cv2, var.field = "ratio", var.threshold = 1)
    sce_object_gene_shared <- intersect(sce_object_gene_var_top, sce_object_gene_cv2_top)

    # Add annotation for which genes are highly variable
    rowSubset(sce_object_filtered, "hvg") <- sce_object_gene_shared

    # Log normalize counts for downstream analyses
    sce_object_filtered %<>% scuttle::logNormCounts()
    sce_object_filtered
}

iterative_pca <- function(n_pcs, standardized_data, pca_object) {
    var_explained <- pca_object$sdev^2
    if (min(var_explained) > attr(n_pcs, "limit")) {
        n_pcs <- ceiling(1.5 * n_pcs)
        print(str_c("Minimum eigenvalue ",  signif(min(var_explained), digits = 3), " is larger than the threshold ", signif(attr(n_pcs, "limit"), 3), ". Increasing rank to ", n_pcs))
        set.seed(1)
        pca_object <- PCAtools::pca(standardized_data, BSPARAM = BiocSingular::IrlbaParam(), rank = n_pcs)
        iterative_pca(n_pcs, standardized_data, pca_object)
    } else {
        n_pcs_final <- sum(pca_object$sdev^2 > attr(n_pcs, "limit"))
        attr(n_pcs_final, "limit") <- attr(n_pcs, "limit")
        return(list(n_pcs = n_pcs_final, pca_standardized = pca_object))
    }
}

# Use corral to estimate PCA from counts directly.
estimate_corral <- function(sce_object, method = "irlba") {
    # Document this line better
    # Todo: let user choose between inverse normal and scaling and centering
    standardized_data <- SingleCellExperiment::counts(sce_object)[rowSubset(sce_object, "hvg"),] %>% 
        corral::corral_preproc() %>% apply(1, RNOmni::RankNorm) %>% t()

    if (method == "irlba") {
        init_npcs <- dim(standardized_data) %>% min() %>% sqrt() %>% ceiling()
        set.seed(1)
        pca_standardized_init <- PCAtools::pca(standardized_data, BSPARAM = BiocSingular::IrlbaParam(), rank = init_npcs)
        n_pcs <- PCAtools::chooseGavishDonoho(standardized_data, var.explained = pca_standardized_init$sdev^2, noise = 1)
        pca_output <- iterative_pca(n_pcs, standardized_data, pca_standardized_init)

        pca_standardized <- pca_output$pca_standardized
        n_pcs <- pca_output$n_pcs
    } else {
        pca_standardized <- PCAtools::pca(standardized_data, BSPARAM = BiocSingular::ExactParam())
        n_pcs <- PCAtools::chooseGavishDonoho(standardized_data, var.explained = pca_standardized$sdev^2, noise = 1)
    }
    print(str_c("Chose ", n_pcs, " PCs"))

    # Do rank normalization on eigenvectors to make UMAP consistent  
    # TODO: Make this optional
    SingleCellExperiment::reducedDim(sce_object, "corral") <- apply(pca_standardized$rotated[,1:n_pcs], 2, RNOmni::RankNorm)
    sce_object
}

# Estimate UMAP with corral by default
estimate_umap <- function(sce_object, reduced_dim = "corral", random_seed = 1L) {
    set.seed(random_seed)
    sce_object <- scater::runUMAP(sce_object, dimred = reduced_dim)
    sce_object
}

# Infer clusters using Leiden algorith on SNN graph.  k = 10 is default parameter
leiden_clustering <- function(sce_object, reduced_dim = "corral", k = 10, type = "jaccard", resolution = 1.0, random_seed = 1L) {
    blus_param <- bluster::SNNGraphParam(k = k, type = type, cluster.fun = "leiden",
        cluster.args = list(objective_function = "modularity", resolution = resolution))
    set.seed(random_seed)
    colLabels(sce_object) <- scran::clusterCells(sce_object, use.dimred = reduced_dim, BLUSPARAM = blus_param)
    sce_object
}

umap_plot <- function(color_by, prefix = "", sce_object, suffix = "", exprs_values = "logcounts",
                      continuous_scale = TRUE, gene_expression = TRUE,
                      facet_var = NULL, facet_ncol = NULL, add_legend = TRUE, font_name = "Noto Sans",
                      label_var = NULL, scale_name = NULL, file_name = NULL, plot_name = NULL,
                      width = 6, height = 6, random_seed = 1L, ...) {

    if (continuous_scale) {
        continuous_values <- scater::retrieveCellInfo(sce_object, color_by, exprs_values = exprs_values)
        sce_object <- sce_object[,order(continuous_values$value)]
    }

    showtext_auto()
    umap_plot <- scater::plotUMAP(object = sce_object, 
                          color_by = color_by, 
                          other_fields = facet_var,
                          add_legend = add_legend, 
                          exprs_values = exprs_values, ...)

    if (!is.null(scale_name)) {
        scale_name_final <- scale_name
    } else {
        if (gene_expression) {
            scale_name_final <- expression(paste(log[2.0], "(counts)"))
        } else {
            scale_name_final <- color_by
        }
    }

    if (continuous_scale) {
        umap_plot <- umap_plot + scale_color_gradient(low = "lightgrey", high = "blue", name = scale_name_final)
    } else {
        #unique_values <- retrieveCellInfo(sce_object, color_by) %>% extract2("value") %>% unique() %>% length()
        #unique_colors <- grDevices::colors() %>% str_subset("gr[ae]y", negate = TRUE)
        #umap_plot <- umap_plot + scale_color_discrete(unique_colors, name = scale_name_final)
    }

    if (!is.null(plot_name)) {
        umap_plot <- umap_plot + ggtitle(plot_name) + theme(plot.title = element_text())
    } else {
        umap_plot <- umap_plot + ggtitle(color_by)
    }

    if (add_legend) {
        width <- width + 1
    }

    if (!is.null(facet_var) && !is.null(facet_ncol)) {
        facet_formula <- str_c("~ ", facet_var) %>% as.formula()
        umap_plot <- umap_plot + facet_wrap(facet_formula, ncol = facet_ncol) +
            theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
        n_facets <- scater::retrieveCellInfo(sce_object, facet_var) %>% extract2("value") %>% unique %>% length
        facet_rows <- ceiling(n_facets / facet_ncol)
        height <- height * facet_rows
        
        width <- width * facet_ncol
    }

    if (!is.null(font_name)) {
        umap_plot <- umap_plot + theme(text = element_text(family = font_name))
    }

    if (!is.null(file_name)) {
        file_name_final <- str_c(prefix, file_name, suffix)
    } else {
        file_name_final <- str_c(prefix, color_by, suffix)
    }

    ggsave(file_name_final, umap_plot, width = width, height = height, device = cairo_pdf, units = "in")
}

# Plot colData (usually cell type names, cluster labels or QC measures)
coldata_plot <- function(y_var, prefix = "", sce_object, x_var, fill_var = NULL, suffix = "", scale_name = NULL, y_label = NULL, x_label = NULL, font_name = "Noto Sans", file_name = NULL, angled_text = FALSE, celltype = NULL, height = 6, width = 6) {
    # Use plotColData from scater to  make a ggplot object from columns of colData
    coldata_ggplot <- scater::plotColData(object = sce_object,
        y = y_var,
        x = x_var,
        colour_by = fill_var
    )

    # Add angled text to x-axis if axis labels are too long
    if (angled_text) {
        coldata_ggplot <- coldata_ggplot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }

    # Override the y-axis label with the value specified by the user
    if (!is.null(y_label)) {
        coldata_ggplot <- coldata_ggplot + ylab(y_label)
    }

    # Override the x-axis label with the value specified by the user
    if (!is.null(x_label)) {
        if (x_label == FALSE) {
            coldata_ggplot <- coldata_ggplot + theme(axis.title.x = element_blank())
        } else {
            coldata_ggplot <- coldata_ggplot + xlab(y_label)
        }
    }

    # Override scale name with value specified by the user
    if (!is.null(scale_name)) {
        scale_name_final <- scale_name
        coldata_ggplot <- coldata_ggplot + scale_fill_discrete(name = scale_name_final)
    }

    # Add a user specified font name
    if (!is.null(font_name)) {
        coldata_ggplot <- coldata_ggplot + theme_classic(base_family = font_name) +
            theme(text = element_text(family = font_name))
    } else {
        coldata_ggplot <- coldata_ggplot + theme_classic()
    }

    # Add a border around the plot
    coldata_ggplot <- coldata_ggplot + theme(panel.border = element_rect(color = "black", fill = NA))  + scale_color_manual(values = c("black"))

    # Add _violinplot to file name
    if (!is.null(file_name)) {
        file_name_final <- str_c(prefix, file_name, suffix, "_violinplot")
    } else {
        file_name_final <- str_c(prefix, y_var, suffix, "_violinplot")
    }

    # Save plot with ggplot
    ggsave(file_name_final, coldata_ggplot, width = width, height = height, units = "in", device = cairo_pdf)
}

qc_umap_plots <- function(prefix, sce_object, qc_var, ...) {
    umap_plot(qc_var, prefix, sce_object, ...)
    umap_plot(qc_var, prefix, sce_object, facet_var = "Sample", facet_ncol = 4, suffix = "_sample_names", ...)
}

qc_plots <- function(sce_object, prefix, doublets = TRUE, facet_ncol, point_size = 0.3) {
    # Plot clusters
    qc_umap_plots(prefix, sce_object, "label", continuous_scale = FALSE, add_legend = FALSE, file_name = "clusters", text_by = "label", point_size)
    # Plot samples
    umap_plot("Sample", prefix, sce_object, continuous_scale = FALSE, add_legend = TRUE,
              file_name = "sample_names_single", plot_name = "Sample")
    umap_plot("Sample", prefix, sce_object, continuous_scale = FALSE,
              facet_var = "Sample", facet_ncol = 4, 
              file_name = "sample_names", add_legend = FALSE, plot_name = "Sample")

    # Plot QC measures
    qc_umap_plots(prefix, sce_object, "subsets_percent_mt_percent", scale_name = "% MT", file_name = "percent_mt", gene_expression = FALSE)
    qc_umap_plots(prefix, sce_object, "log_sum", gene_expression = FALSE)
    qc_umap_plots(prefix, sce_object, "log_detected", gene_expression = FALSE)

    if (doublets) {
        qc_umap_plots(prefix, sce_object, "doublet_scores", gene_expression = FALSE)
        n_clusters <- unique(sce_object$label) %>% length
        width <- ceiling(n_clusters / 2)
        coldata_plot("doublet_scores", prefix, sce_object, x_var = "label", width = width, y_label = "Doublet scores")
    }
}

add_celltype_labels <- function(sce_object, celltype_labels, new_celltypes = F) {
    sce_metadata <- colData(sce_object) %>% as_tibble
    if (new_celltypes == T) {
        sce_metadata %<>% dplyr::select(-celltype)
    }
    celltype_labels_df <- tibble(label = factor(names(celltype_labels)), celltype = factor(celltype_labels))
    sce_metadata_cell_labels <- left_join(sce_metadata, celltype_labels_df)
    colData(sce_object) <- DataFrame(sce_metadata_cell_labels)
    sce_object
}


add_sample_column <- function(sce_object, sample_name, sample_column, hvg_name, shared_genes) {
    colData(sce_object) %<>% as_tibble() %>% dplyr::mutate("{{sample_column}}" := sample_name) %>% DataFrame()
    rowData(sce_object) %<>% as_tibble() %>% dplyr::select(-{{hvg_name}}) %>% DataFrame()
    SingleCellExperiment::reducedDim(sce_object, "corral") <- NULL
    sce_object[shared_genes,]
}

get_hvgs <- function(sce_object, hvg_column) {
    if (!is_in(hvg_column, colnames(rowData(sce_object)))) {
        stop(str_c("Error: highly variable gene column ", hvg_column, " does not exist in the column names of rowData."))
    }
    rownames(sce_object)[rowSubset(sce_object, hvg_column)]
}

# Merge batches using MNN method from batchelor package
# TODO: Make actual batch correction optional!
merge_batches <- function(sce_list, batch_column = NULL, combine_hvgs = "union", hvg_column = "hvg", BSPARAM = IrlbaParam()) {
    if (combine_hvgs == "union") {
        shared_hvgs <- map(sce_list, get_hvgs, hvg_column) %>% reduce(union)
    } else {
        shared_hvgs <- map(sce_list, get_hvgs, hvg_column) %>% reduce(intersect)
    }
    
    shared_genes <- map(sce_list, rownames) %>% reduce(intersect)
    dimred_ncols <- map(sce_list, reducedDim, "corral") %>% map_int(ncol) %>% max

    if (is.null(batch_column)) {
        sce_list_samples <- imap(sce_list, add_sample_column, "batch", hvg_column, shared_genes)
        batch_column <- "batch"
        rm(sce_list)
        gc()
    } else {
        sce_list_samples <- sce_list
        rm(sce_list)
        gc()
    }

    # Combine SCE objects into a single object with a batch column
    sce_combined <- reduce(sce_list_samples, cbind)
    rowData(sce_combined)[[hvg_column]] <- is_in(rownames(sce_combined), shared_hvgs)

    # Normalize batches to each other before they are combined
    sce_batch_norm <- batchelor::multiBatchNorm(sce_combined,
        batch = colData(sce_combined)[[batch_column]],
        normalize.all = TRUE,
        preserve.single = TRUE,
        subset.row = rowData(sce_combined)[[hvg_column]]
    )
    rm(sce_combined)
    gc()

    set.seed(1L)
    fastmnn_params <- batchelor::FastMnnParam(auto.merge = TRUE,
        d = dimred_ncols,
        cos.norm = FALSE,
        BSPARAM = BSPARAM
    )
    sce_merge <- batchelor::batchCorrect(sce_batch_norm,
        PARAM = fastmnn_params,
        batch = colData(sce_batch_norm)[[batch_column]],
        subset.row = rowData(sce_batch_norm)[[hvg_column]],
        correct.all = TRUE
    )
    counts(sce_merge) <- counts(sce_batch_norm)
    logcounts(sce_merge) <- logcounts(sce_batch_norm)
    colData(sce_merge) <- colData(sce_batch_norm)
    rowData(sce_merge) <- cbind(rowData(sce_batch_norm), rowData(sce_merge))
    sce_merge
}

AGS_sce_macros_batch1 <- read_rds("../clustering/AGS_sce_filtered_macros_batch1.rda")
AGS_sce_macros_batch2 <- read_rds("../clustering/AGS_sce_filtered_macros_batch2.rda")

# Plot marker genes
marker_genes = c(
    "CLDN5", # Endothelial cells
    "CALD1", # Smooth muscle cells

    "CD68", # macrophages
    "CD14", # macrophages/monocytes
    "FCGR3A", # monocytes

    "S100A8", # inflammatory macrophages
    "IL1B", # inflammatory macrophages
    "CD300E", # inflammatory macrophages
    "CSF3R", # inflammatory macrophages

    "APOE", # foam cells
    "TREM2", # foam cells
    "LPL", # foam cells
    "ABCA1", # foam cells
    "FABP5", # foam cells
    "FTL", # foam cells
    "PLTP", # foam cells

    "LYVE1", # TR macrophages
    "MRC1", # TR macrophages
    "FOLR2", # TR macrophages

    "HLA-DRB1", # MHC-hi macrophages
    "HLA-DQB1", # MHC-hi macrophages

    "FCER1A", # classical DCs
    "PLD4", # plasmacytoid DCs
    "XCR1", # XCR1+ DCs

    "FCGR3B", # neutrophils

    "MS4A2", # basophils/neutrophils
    "GATA2", # basophils/neutrophils
    "HPGD", # basophils

    "CD3D", # T cells
    "CD3E", # T cells
    "CD4", # CD4 T cells
    "CD8A", # CD8 T cells
    "TRDC", # Gamma-delta T-cells, NK cells
    "TRGC2", # Gamma-delta T-cells, NK cells
    "FOXP3", # T-regs
    "GZMB", # Effector T-cells
    "GZMH", # Effector T-cells
    "GZMK", # Effector T-cells
    "NKG7", # Effector T-cells
    "TCL1A", # Effector T-cells
    "IL7R", # Naive T-cells

    "NCAM1", # NK cells
    "XCL1", # active NK cells

    "CD19", # B cells
    "CD79A", # B cells

    "SDC1", # Plasma cells

    "MKI67" # Cycling cells
)

AGS_sce_leiden_batch1 <- filter_osca(AGS_sce_macros_batch1, iterative_filter = TRUE) %>%
    estimate_corral() %>%
    estimate_umap() %>%
    leiden_clustering()
qc_plots(AGS_sce_leiden_batch1, "all_cells_batch1/", doublets = F)
catch <- map(marker_genes, umap_plot, "all_cells_batch1/", AGS_sce_leiden_batch1, height = 6, width = 7)

AGS_sce_leiden_batch2 <- filter_osca(AGS_sce_macros_batch2, iterative_filter = TRUE) %>%
    estimate_corral() %>%
    estimate_umap() %>%
    leiden_clustering()
qc_plots(AGS_sce_leiden_batch2, "all_cells_batch2/", doublets = F)
catch <- map(marker_genes, umap_plot, "all_cells_batch2/", AGS_sce_leiden_batch2, height = 6, width = 7)

AGS_sce_macros_batch1_good <- dplyr::filter(AGS_sce_leiden_batch1, !is_in(label, c(4)))
AGS_sce_leiden_batch1_good <- filter_osca(AGS_sce_macros_batch1_good, filter_cells = FALSE) %>%
    estimate_corral() %>%
    estimate_umap() %>%
    leiden_clustering()
qc_plots(AGS_sce_leiden_batch1_good, "all_cells_batch1_good/", doublets = F)
catch <- map(marker_genes, umap_plot, "all_cells_batch1_good/", AGS_sce_leiden_batch1_good, height = 6, width = 7)

AGS_sce_macros_batch2_good <- dplyr::filter(AGS_sce_leiden_batch2, !is_in(label, c(6)))
AGS_sce_leiden_batch2_good <- filter_osca(AGS_sce_macros_batch2_good, filter_cells = FALSE) %>%
    estimate_corral() %>%
    estimate_umap() %>%
    leiden_clustering()
qc_plots(AGS_sce_leiden_batch2_good, "all_cells_batch2_good/", doublets = F)
catch <- map(marker_genes, umap_plot, "all_cells_batch2_good/", AGS_sce_leiden_batch2_good, height = 6, width = 7)

# Create a named list of all separate batches 
AGS_sce_list <- list("1" = AGS_sce_leiden_batch1_good, "2" = AGS_sce_leiden_batch2_good)

AGS_sce_macros_merge <- merge_batches(AGS_sce_list)
AGS_sce_leiden <- estimate_umap(AGS_sce_macros_merge, "corrected") %>%
    leiden_clustering("corrected", resolution = 0.4)
qc_plots(AGS_sce_leiden, "all_cells/", doublets = F)
catch <- map(marker_genes, umap_plot, "all_cells/", AGS_sce_leiden, height = 6, width = 7)

weird_cells <- dplyr::filter(AGS_sce_leiden, is_in(label, c(3,10)))
AGS_sce_leiden_batch2_good2 <- dplyr::filter(AGS_sce_leiden_batch2_good, !is_in(Cell_ID, colData(weird_cells)$Cell_ID))

AGS_sce_list2 <- list("1" = AGS_sce_leiden_batch1_good, "2" = AGS_sce_leiden_batch2_good2)

AGS_sce_macros_merge2 <- merge_batches(AGS_sce_list2)
AGS_sce_leiden2 <- estimate_umap(AGS_sce_macros_merge2, "corrected") %>%
    leiden_clustering("corrected", resolution = 0.4)
qc_plots(AGS_sce_leiden2, "all_cells/", doublets = F)
catch <- map(marker_genes, umap_plot, "all_cells/", AGS_sce_leiden2, height = 6, width = 7)

AGS_markers <- findMarkers(AGS_sce_leiden2, test = "t", direction = "up", pval.type = "some", assay.type = "reconstructed", row.data = rowData(AGS_sce_leiden))
AGS_markers_df <- as.list(AGS_markers) %>% map(as_tibble) %>% bind_rows(.id = "Cluster")
write_tsv(AGS_markers_df, "macro_markers.tsv")

# Note: 4 may be low quality, but does not have clear cell markers
celltype_labels <- c(
    "1" = "Foam cells",
    "2" = "Foam cells",
    "3" = "Classical DCs",
    "4" = "Foam cells",
    "5" = "Classical DCs",
    "6" = "MHC-hi macrophages",
    "7" = "Inflammatory macrophages",
    "8" = "Foam cells",
    "9" = "Foam cells",
    "10" = "Foam cells"
)

AGS_sce_celltypes <- add_celltype_labels(AGS_sce_leiden2, celltype_labels, new_celltypes = TRUE)

AGS_silhouette <- reducedDim(AGS_sce_celltypes, "corrected") %>% 
    approxSilhouette(AGS_sce_celltypes$celltype) %>% as_tibble
colData(AGS_sce_celltypes)$silhouette_width <- AGS_silhouette$width
write_rds(AGS_sce_celltypes, "AGS_sce_celltypes.rda")
