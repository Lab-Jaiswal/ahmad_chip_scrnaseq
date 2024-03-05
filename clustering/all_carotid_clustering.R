library(AnnotationHub)
library(AnnotationDbi)
library(BiocParallel)
library(BiocSingular)

library(DropletUtils)
library(scater)
library(scran)
library(scuttle)
library(scDblFinder)
library(celldex)
library(SingleR)
library(corral)
library(PCAtools)
library(batchelor)
library(bluster)
library(tidySingleCellExperiment)

library(RNOmni)
library(VGAM)
library(lme4)
library(glmmTMB)
library(betareg)
library(RColorBrewer)
library(scales)
library(igraph)
library(meta)
library(showtext)
library(forestplot)
library(openxlsx)

library(rlist)
library(ggrepel)
library(broom)
library(readxl)
library(magrittr)
library(tidyverse)

add_dup_column <- function(gene_rows, grouping_df, dup_column) {
    gene_rows$new_names <- str_c(dup_column, 1:nrow(gene_rows))
    gene_rows
}

# Group gene annotation by Ensembl ID, and create new columns for Ensembl IDs with multiple EntrezIDs
deduplicate_ids <- function(gene_annot, group_column = "GENEID", dup_column = "ENTREZID") {
    gene_annot_nodups <- group_by(gene_annot, !!sym(group_column)) %>% group_modify(add_dup_column, dup_column)
}

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

filter_osca <- function(sce_object, gene_annot = NULL, filter_cells = T, iterative_filter = F, ) {
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

# Compute doublet scores using number of PCs identified above.  k = 50 is default parameter for nearest neighbors
compute_doublets <- function(sce_object, n_neighbors = 50, BSPARAM = IrlbaParam()) {
    n_pcs <- reducedDim(sce_object, "corral") %>% ncol()
    doublet_scores <- scDblFinder::computeDoubletDensity(sce_object, k = n_neighbors, dims = n_pcs, BSPARAM = BSPARAM)
    colData(sce_object)$doublet_scores <- doublet_scores
    sce_object    
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

# Use SingleR to roughly assign cell types from bulk reference data
singler_annotation <- function(sce_object, sce_bulk, label_name) {
    singler_labels <- SingleR::SingleR(counts(sce_object), sce_bulk, labels = sce_bulk$label.fine)
    singler_labels_df <- as_tibble(singler_labels, rownames = "Barcode") 
    colData(sce_object)[[label_name]] <- factor(singler_labels_df$pruned.labels)
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

celltype_bar_plot <- function(prefix, sce_object, width = 10, height = 5, file_name = NULL, relative = FALSE) {
    bar_plot <- ggcells(sce_object, aes(x = Sample, fill = celltype)) 

    if (relative) { 
        bar_plot <- bar_plot + geom_bar(width = 0.7, color = "black", position = position_fill()) 
        file_name <- str_c(prefix, "celltype_barplot_fill")  
    } else {
        bar_plot <- bar_plot + geom_bar(width = 0.7, color = "black", position = position_stack()) 
        file_name <- str_c(prefix, "celltype_barplot_stack")  
    }

    bar_plot <- bar_plot + theme_classic(base_family = "Noto Sans") +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1),
              panel.border = element_rect(color = "black", fill = NA)
        ) +
        scale_fill_brewer(palette = "Paired")

    ggsave(file_name, bar_plot, width = width, height = height, units = "in", device = cairo_pdf)
}


celltype_plots <- function(sce_object, prefix, fill_var = NULL) {
    qc_umap_plots(prefix, sce_object, "celltype", text_by = "celltype", file_name = "celltypes", continuous_scale = FALSE, add_legend = FALSE, plot_name = "", point_size = 0.5, text_size = 3)
    celltype_bar_plot(prefix, sce_object, relative = FALSE)
    celltype_bar_plot(prefix, sce_object, relative = TRUE)

    coldata_plot("doublet_scores", prefix, sce_object, fill_var = fill_var, suffix = "_celltypes", x_var = "celltype", angled_text = TRUE, x_label = "Cell type", y_label = "Doublet scores")
    coldata_plot("subsets_percent_mt_percent", prefix, sce_object, fill_var = fill_var, suffix = "_celltypes", x_var = "celltype", angled_text = T, y_label = "% MT", file_name = "percent_mt")
    coldata_plot("log_sum", prefix, sce_object, fill_var = fill_var, suffix = "_celltypes", x_var = "celltype", angled_text = T, y_label = expression(paste(log[2], "(total counts)")))
    coldata_plot("log_detected", prefix, sce_object, fill_var = fill_var, suffix = "_celltypes", x_var = "celltype", angled_text = T, y_label = expression(paste(log[2], "(detected_genes)")))
    coldata_plot("silhouette_width", prefix, sce_object, fill_var = fill_var, suffix = "_celltypes", x_var = "celltype", angled_text = T, y_label = "Silhouette width")
}

celltype_proportion_boxplot <- function(celltype_name, prefix, sce_object, x_var, suffix = "", scale_name = NULL, 
                             file_name = NULL, y_label = NULL, 
                             angled_text = F, height = 4, width = 12) {
    cell_metadata <- colData(sce_object) %>% as_tibble 
    sample_metadata <- dplyr::select(cell_metadata, Sample, celltype, !!x_var)
    sample_abundances <- group_by(sample_metadata, Sample, celltype) %>%
        summarise(n_cells = n())
    total_cells <- group_by(sample_abundances, Sample) %>% summarise(total_cells = sum(n_cells))
    sample_abundances_df <- left_join(sample_abundances, total_cells) %>% left_join(sample_metadata) %>% dplyr::filter(celltype == celltype_name)
    sample_abundances_df$prop_cells <- sample_abundances_df$n_cells / sample_abundances_df$total_cells

    box_plot <- ggplot(sample_abundances_df, aes(Mutation_1, prop_cells, fill = Mutation_1)) +
        geom_boxplot() +
        theme_classic(base_family = "Noto Sans") +
        theme(axis.title.x = element_blank(),
            legend.position = "none",
            panel.border = element_rect(color = "black", fill = NA)
        )

    if (!is.null(y_label)) {
        box_plot <- box_plot + ylab(y_label)
    }

    if (!is.null(scale_name)) {
        scale_name_final <- scale_name
        box_plot <- box_plot + scale_fill_discrete(name = scale_name_final)
    }

    celltype_format <- str_replace_all(celltype_name, " ", "_") %>% str_replace_all("/", "_")
    if (!is.null(file_name)) {
        file_name_final <- str_c(prefix, file_name, suffix, "_boxplot")
    } else {
        file_name_final <- str_c(prefix, celltype_format, suffix, "_boxplot")
    }
    ggsave(file_name_final, box_plot, width = width, height = height, units = "in", device = cairo_pdf)
}

metadata_scatterplot <- function(celltype_name, prefix, sce_object, x_var, suffix = "", 
                             file_name = NULL, y_label = NULL, x_label = NULL, height = 5, width = 6) {
    cell_metadata <- colData(sce_object) %>% as_tibble 
    sample_metadata <- dplyr::select(cell_metadata, Sample, celltype, !!x_var)
    sample_abundances <- group_by(sample_metadata, Sample, celltype) %>%
        summarise(n_cells = n())
    total_cells <- group_by(sample_abundances, Sample) %>% summarise(total_cells = sum(n_cells))
    sample_abundances_df <- left_join(sample_abundances, total_cells) %>% left_join(sample_metadata) %>% dplyr::filter(celltype == celltype_name)
    sample_abundances_df$prop_cells <- sample_abundances_df$n_cells / sample_abundances_df$total_cells

    scatter_plot <- ggplot(sample_abundances_df, aes_string(x_var, "prop_cells")) +
        geom_point() +
        stat_smooth(method = "lm") +
        theme_classic(base_family = "Noto Sans") +
        theme(panel.border = element_rect(color = "black", fill = NA))

    if (!is.null(y_label)) {
        scatter_plot <- scatter_plot + ylab(y_label)
    }

    if (!is.null(x_label)) {
        scatter_plot <- scatter_plot + xlab(x_label)
    }

    celltype_format <- str_replace_all(celltype_name, " ", "_") %>% str_replace_all("/", "_")
    if (!is.null(file_name)) {
        file_name_final <- str_c(prefix, file_name, suffix, "_scatterplot")
    } else {
        file_name_final <- str_c(prefix, celltype_format, suffix, "_scatterplot")
    }
    ggsave(file_name_final, scatter_plot, width = width, height = height, units = "in", device = cairo_pdf)
}

# TODO: throw an error if cell type is not in name sof top_table_list
volcano_plot <- function(celltype, top_table_list, filename, width = 5, height = 4, logfc_cutoff = NA, cutoff = 0.05, cutoff_column = "adj.P.Val", log_column = "logFC", xlabel = expression(paste(log[2.0], " fold change"))) {
    top_table <- as_tibble(top_table_list[[celltype]])
    top_table_filter <- dplyr::filter(top_table, !is.na(P.Value)) %>% arrange(P.Value)
    top_table_filter$Significant <- top_table_filter[[cutoff_column]] < cutoff
    
    if (!is.na(logfc_cutoff)) {
        sig_logfc <- abs(top_table_filter[[log_column]]) > logfc_cutoff
        top_table_filter$Significant <- top_table_filter$Significant & sig_logfc
    }
    top_table_filter$Log_Pvalue <- -log10(top_table_filter$P.Value)
    top_table_filter$Symbol_new <- NA
    top_table_filter$Symbol_new[1L:10L] <- top_table_filter$Symbol[1L:10L]
    top_table_sig_pos <- dplyr::filter(top_table_filter, Significant & !!ensym(log_column) > 0)
    top_table_sig_neg <- dplyr::filter(top_table_filter, Significant & !!ensym(log_column) < 0)

    p <- ggplot() +
        geom_point(data = top_table_filter, color = "gray", aes_string(x = log_column, y = "Log_Pvalue")) +
        geom_point(data = top_table_sig_pos, color = "blue", aes_string(x = log_column, y = "Log_Pvalue")) +
        geom_point(data = top_table_sig_neg, color = "red", aes_string(x = log_column, y = "Log_Pvalue")) +
        theme_bw(base_family = "Noto Sans") +
        #geom_hline(yintercept = -log10(0.05 / nrow(top_table_filter)), color = "red", linetype = "dashed") +
        geom_text_repel(data = top_table_filter, aes_string(x = log_column, y = "Log_Pvalue"), 
            label = top_table_filter$Symbol_new, 
            color = muted("red"), 
            segment.color = "black", 
            box.padding = 0.5, 
            min.segment.length = 0.0, 
            size = 3.0, 
            nudge_y = 0.25, 
            seed = 12345L, 
            family = "Noto Sans") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background = element_blank(),
              strip.text = element_text(size = 11.0, face = "bold"),
              legend.title = element_blank(),
              panel.border = element_rect(size = 1.0, color = "black"),
              plot.background = element_blank(),
              plot.title = element_text(hjust = 0.5, face = "bold")) +
        xlab(xlabel) + ylab(expression(paste("-", log[10.0], " p-value"))) 
    ggsave(filename, p, width = width, height = height, units = "in", device = cairo_pdf)
}

pseudobulk_counts <- function(sce_object, min_ncells = 10L) {
    celltype_metadata <- colData(sce_object) %>% as_tibble() %>% dplyr::select(Sample, celltype) %>% S4Vectors::DataFrame()
    sce_pseudobulk <- scuttle::aggregateAcrossCells(sce_object, id = celltype_metadata)
    sce_pseudobulk_filter <- sce_pseudobulk[, sce_pseudobulk$ncells >= min_ncells]
    sce_pseudobulk_filter
}

# Call wrapper for pseudobulk differential expression with edgeR
edger_bulk_dge <- function(pseudobulk_object, coefficient) {
    edger_dge <- scran::pseudoBulkDGE(pseudobulk_object,
        label = pseudobulk_object$Celltype,
        design = ~ batch + Mutation_1,
        coef = coefficient,
        condition = pseudobulk_object$Status,
        row.data = rowData(pseudobulk_object),
        method = "edgeR"
    )
    edger_dge
}

# Call wrapper for pseudobulk differential expression with limma
limma_bulk_dge <- function(pseudobulk_object, coefficient) {
    limma_dge <- scran::pseudoBulkDGE(pseudobulk_object,
        label = pseudobulk_object$celltype,
        design = ~ batch + Mutation_1,
        coef = coefficient,
        row.data = rowData(pseudobulk_object),
        method = "voom"
    )
    limma_dge
}

# Use openxlsx to make tables summarizing differential expression analysis
export_de_table <- function(index, workbook_name, full_de_table, wb) {
    # Only keep genes that were actually tested and columns we need
    de_table <- as_tibble(full_de_table) %>%
        dplyr::filter(!is.na(P.Value)) %>% 
        dplyr::select(ID, Symbol, GENEBIOTYPE:GENEIDVERSION, hvg, logFC:B)

    de_table %<>% mutate(across(logFC:B, signif, digits = 3L)) %>% arrange(P.Value) # Only show 3 significant digits
    # Identify pvalue and logFC columns
    pval_cols <- colnames(de_table) %>% str_detect("P.Value") %>% which()
    logfc_cols <- colnames(de_table) %>% str_detect("logFC") %>% which()

    addWorksheet(wb, sheetName = workbook_name) # Add worksheet
    writeDataTable(wb, sheet = index, x = de_table) # Write data to worksheet

    # Make p < 0.05 red, and fill positive logFC with green and negative logFC with red
    sig_style <- createStyle(fontColour = "red")
    conditionalFormatting(wb, index, cols = pval_cols,
                          rows = seq_len(nrow(de_table)), rule = "<0.05", style = sig_style)
    conditionalFormatting(wb, index, cols = logfc_cols,
                          rows = seq_len(nrow(de_table)), style = c("#63BE7B", "white", "red"),
                          type = "colourScale")

    # Set column widths to reasonable values
    setColWidths(wb, index, cols = 1L:3L, widths = "auto")
    setColWidths(wb, index, cols = 4L, widths = 45L)
    setColWidths(wb, index, cols = 5L:ncol(de_table), widths = "auto")

    pageSetup(wb, index, orientation = "landscape", fitToWidth = TRUE) # Set to landscape view
    freezePane(wb, index, firstRow = TRUE) # Freeze pane
}

run_glmmtmb <- function(.x, .y, outcome, model_string, metric_name, glmmtmb_family, odds_ratio = TRUE) {
    model_data <- .x
    model_formula <- str_c(outcome, model_string, sep = " ~ ") %>% as.formula
    vglm_fit <- glmmTMB::glmmTMB(model_formula, family = glmmtmb_family, data = model_data)
    vglm_summary <- summary(vglm_fit)
    vglm_summary_df <- as_tibble(vglm_summary$coefficients$cond, rownames = "term")
    vglm_summary_df$conf.low <- vglm_summary_df$Estimate - (1.96 * vglm_summary_df$`Std. Error`)
    vglm_summary_df$conf.high <- vglm_summary_df$Estimate + (1.96 * vglm_summary_df$`Std. Error`)
    if (odds_ratio == TRUE) {
        vglm_summary_df$OR <- exp(vglm_summary_df$Estimate)
        vglm_summary_df$OR_se <- exp(vglm_summary_df$`Std. Error`)
        vglm_summary_df$OR_conf.low <- exp(vglm_summary_df$conf.low)
        vglm_summary_df$OR_conf.high <- exp(vglm_summary_df$conf.high)
    }
    vglm_summary_df$metric_name <- metric_name
    vglm_summary_df
}

run_vglm <- function(.x, .y, outcome, model_string, metric_name, vgam_family, odds_ratio = TRUE) {
    model_data <- .x
    model_formula <- str_c(outcome, model_string, sep = " ~ ") %>% as.formula
    vglm_fit <- VGAM::vglm(model_formula, family = vgam_family, data = model_data, trace = TRUE)
    vglm_summary <- summary(vglm_fit, HDEtest = FALSE)
    vglm_summary <- as_tibble(vglm_summary@coef3, rownames = "term")
    vglm_summary$conf.low <- vglm_summary$Estimate - (1.96 * vglm_summary$`Std. Error`)
    vglm_summary$conf.high <- vglm_summary$Estimate + (1.96 * vglm_summary$`Std. Error`)
    if (odds_ratio == TRUE) {
        vglm_summary$OR <- exp(vglm_summary$Estimate)
        vglm_summary$OR_se <- exp(vglm_summary$`Std. Error`)
        vglm_summary$OR_conf.low <- exp(vglm_summary$conf.low)
        vglm_summary$OR_conf.high <- exp(vglm_summary$conf.high)
    }
    vglm_summary$metric_name <- metric_name
    vglm_summary
}

run_betabinomial_vglm <- function(...) {
    #if (mixed_effects == TRUE) {
        run_vglm(vgam_family = VGAM::betabinomial(lrho = "rhobitlink"), odds_ratio = TRUE, ...)
    #} else {
        #run_glmmtmb(glmmtmb_family = glmmTMB::betabinomial(), odds_ratio = TRUE)
    #}
}

run_fisherzlink_glm <- function(...) {
    run_vglm(vgam_family = uninormal(lmean = fisherzlink()), odds_ratio = FALSE, ...)
}

run_negbinomial_glm <- function(...) {
    if (mixed_effects == TRUE) {
        run_vglm(vgam_family = negbinomial(), odds_ratio = TRUE, ...)
    } else {
        run_glmmtmb(glmmtmb_family = nbinom1(), odds_ratio = TRUE)
    }
}

run_beta_glm <- function(outcome, model_string, model_data, metric_name, odds_ratio = TRUE) {
    if (min(model_data[[outcome]]) == 0 || max(model_data[[outcome]]) == 1) {
        model_data[[outcome]] <- (model_data[[outcome]] * (length(model_data[[outcome]]) - 1) + 0.5) / length(model_data[[outcome]])
    }
    full_model_string <- str_c(outcome, model_string, sep = " ~ ")
    model_formula <- as.formula(full_model_string)
    glm_fit <- betareg::betareg(full_model_string, model_data)
    glm_summary <- tidy(glm_fit, conf.int = T)
    if (odds_ratio == TRUE) {
        glm_summary$OR <- exp(glm_summary$estimate)
        glm_summary$OR_se <- exp(glm_summary$estimate_se)
        glm_summary$OR_conf.low <- exp(glm_summary$conf.low)
        glm_summary$OR_conf.high <- exp(glm_summary$conf.high)
    }
    glm_summary$metric_name <- metric_name
    glm_summary
}

run_lm <- function(.x, .y, outcome, model_string, metric_name, rank_norm = TRUE, logit = FALSE, log = FALSE, arcsin = FALSE, odds_ratio = FALSE, mixed_effects = FALSE) {
    model_data <- .x
    #transform_args <- c(rank_norm, logit, log, arcsin)
    #if (sum(transform_args) != 1) {j
        #return("Error: one and only one of rank_norm, logit, log, or arcsin must be set to TRUE.")
    #}
    #if ((min(model_data[[outcome]]) < 0 || max(model_data[[outcome]]) > 1) && (logit == TRUE || arcsin == TRUE )) {
        #return("Error: outcome contains values <0 or >1 and logit or arcsin are set to TRUE. Set these values to 0 or 1 or remove them.")
    #}
    #if (min(model_data[[outcome]]) < 0 && log == TRUE) {
        #return("Error: outcome contains values <0 and log is set to TRUE.  Transform these values to be greater than 0 or remove them.")
    #}
    if (rank_norm == TRUE) {
        model_data[[outcome]] <- RNOmni::RankNorm(model_data[[outcome]])
    }
    if (log == TRUE) {
        model_data[[outcome]] <- log2(model_data[[outcome]])
    }
    if (logit == TRUE) {
        if (min(model_data[[outcome]]) == 0 || max(model_data[[outcome]] == 1)) {
            model_data[[outcome]] <- (model_data[[outcome]] * (length(model_data[[outcome]]) - 1) + 0.5) / length(model_data[[outcome]])
        }
        model_data[[outcome]] <- log(model_data[[outcome]] / (1 - model_data[[outcome]]))
    }
    if (arcsin == TRUE) {
        model_data[[outcome]] <- asin(sqrt(model_data[[outcome]]))
    }
    full_model_string <- str_c(outcome, model_string, sep = " ~ ")
    model_formula <- as.formula(full_model_string)
    if (mixed_effects == FALSE) {
        lm_fit <- lm(full_model_string, model_data)
        lm_summary <- tidy(lm_fit, conf.int = T)
        if (odds_ratio == TRUE) {
            if (logit == TRUE || log == TRUE) {
                lm_summary$OR <- exp(lm_summary$estimate)
                lm_summary$OR_se <- exp(lm_summary$se)
                lm_summary$OR_conf.low <- exp(lm_summary$conf.low)
                lm_summary$OR_conf.high <- exp(lm_summary$conf.high)
            } else {
                if (arcsin == TRUE) {
                    lm_summary$OR <- sqrt(sin(lm_summary$estimate))
                    lm_summary$OR_se <- sqrt(sin((lm_summary$se)))
                    lm_summary$OR_conf.low <- sqrt(sin((lm_summary$conf.low)))
                    lm_summary$OR_conf.high <- sqrt(sin((lm_summary$conf.high)))
                }
            }
        }
    } else {
        lm_fit <- lmer(full_model_string, model_data)
        lm_summary <- as.data.frame(summary(lm_fit)$coefficients) %>% rownames_to_column("term")
        if (odds_ratio == TRUE) {
            if (logit == TRUE || log == TRUE) {
                lm_summary$OR <- exp(lm_summary$estimate)
                lm_summary$OR_se <- exp(lm_summary$se)
                #lm_summary$OR_conf.low <- exp(lm_summary$conf.low)
                #lm_summary$OR_conf.high <- exp(lm_summary$conf.high)
            } else {
                if (arcsin == TRUE) {
                    lm_summary$OR <- sqrt(sin(lm_summary$estimate))
                    lm_summary$OR_se <- sqrt(sin((lm_summary$se)))
                    #lm_summary$OR_conf.low <- sqrt(sin((lm_summary$conf.low)))
                    #lm_summary$OR_conf.high <- sqrt(sin((lm_summary$conf.high)))
                }
            }
        }
    }
    lm_summary
}

run_ranknorm_lm <- function(..., mixed_effects = FALSE) {
    run_lm(rank_norm = TRUE, logit = FALSE, log = FALSE, arcsin = FALSE, odds_ratio = FALSE, ...)
}

run_logit_lm <- function(...) {
    run_lm(rank_norm = FALSE, logit = TRUE, log = FALSE, arcsin = FALSE, odds_ratio = TRUE, ...)
}

run_log_lm <- function(...) {
    run_lm(rank_norm = FALSE, logit = FALSE, log = TRUE, arcsin = FALSE, odds_ratio = TRUE, ...)
}

run_arcsin_lm <- function(...) {
    run_lm(rank_norm = FALSE, logit = FALSE, log = FALSE, arcsin = TRUE, odds_ratio = FALSE, ...)
}

combine_detected <- function(metadata_df, sample_df, sce_counts) {
    sample_counts <- counts(sce_counts[, is_in(sce_counts$CellID, metadata_df$CellID)]) %>% as.matrix()
    sample_nonzero <- apply(sample_counts > 0L, 1L, any) %>% which() %>% length()
    tibble(detected = sample_nonzero)
}

avg_silhouette <- function(metadata_df, sample_df) {
    metadata_df$Z_silhouette_width <- VGAM::fisherzlink(metadata_df$silhouette_width)
    mean_silhouette_width <- mean(metadata_df$Z_silhouette_width) %>% VGAM::fisherzlink(inverse = TRUE)
    tibble(mean_silhouette_width = mean_silhouette_width)
}

run_binomial_glm <- function(outcome_val, outcome, model_string, model_data, odds_ratio = TRUE, mixed_effects = FALSE) {
    outcome_lgl <- str_c(outcome, "_lgl")
    model_data[[outcome_lgl]] <- model_data[[outcome]] == outcome_val
    full_model_string <- str_c(outcome_lgl, model_string, sep = " ~ ")
    if (mixed_effects == TRUE) {
        glm_fit <- glmer(full_model_string, model_data, family = "binomial")
        # Why are confidence intervals so slow for glm?
        glm_summary <- as.data.frame(summary(glm_fit)$coefficients) %>% rownames_to_column("term")
        if (odds_ratio == TRUE) {
            glm_summary$OR <- exp(glm_summary$Estimate)
            glm_summary$OR_se <- exp(glm_summary$`Std. Error`)
            #glm_summary$OR_conf.low <- exp(glm_summary$conf.low)
            #glm_summary$OR_conf.high <- exp(glm_summary$conf.high)
        }
    } else {
        glm_fit <- glm(full_model_string, model_data, family = "binomial")
        # Why are confidence intervals so slow for glm?
        glm_summary <- tidy(glm_fit, conf.int = F)
        if (odds_ratio == TRUE) {
            glm_summary$OR <- exp(glm_summary$estimate)
            glm_summary$OR_se <- exp(glm_summary$std.error)
            #glm_summary$OR_conf.low <- exp(glm_summary$conf.low)
            #glm_summary$OR_conf.high <- exp(glm_summary$conf.high)
        }
    }
    glm_summary[[outcome]] <- outcome_val
    glm_summary
}

cell_composition_da <- function(sce_object, model_string, proportion_model = "betabinomial_vglm", meta_contingency = FALSE, mixed_effects = FALSE) {
    #if (!is_in(proportion_model, c("betabinomial_vglm", "beta_glm", "ranknorm_lm", "logit_lm", "arcsin_lm", "fisher_contingency", "chisq_contingency"))) {
        #return("Error: proportion_model must be one of betabinomial_vglm, beta_glm, ranknorm_lm, logit_lm or arcsin_lm")
    #}

    cell_metadata <- colData(sce_object) %>% as_tibble
    model_string_split <- str_remove_all(model_string, " ") %>% str_remove_all("~") %>% str_split("\\+") %>% unlist
    if (proportion_model != "binomial_glm") {
        sample_metadata <- dplyr::select(cell_metadata, Sample, celltype, one_of(model_string_split))
        sample_abundances <- group_by(sample_metadata, Sample, celltype) %>%
            summarise(n_cells = n()) %>% pivot_wider(names_from = "celltype", values_from = "n_cells") %>%
            mutate(across(everything(), replace_na, 0)) %>% 
            pivot_longer(!Sample, names_to = "celltype", values_to = "n_cells")
        total_cells <- group_by(sample_abundances, Sample) %>% summarise(total_cells = sum(n_cells)) 
        sample_metadata_join <- dplyr::select(sample_metadata, -celltype) %>% distinct
        sample_abundances_df <- left_join(sample_abundances, total_cells) %>% left_join(sample_metadata_join)

        if(is_in(proportion_model, c("beta_glm", "ranknorm_lm", "logit_lm", "arcsin_lm"))) {
            sample_abundances_df$prop_cells <- sample_abundances_df$n_cells / sample_abundances_df$total_cells
        }

        switch(proportion_model,
            betabinomial_vglm = {proportion_model_func <- run_betabinomial_vglm},
            beta_glm = {proportion_model_func <- run_beta_glm},
            logit_lm = {proportion_model_func <- run_logit_lm},
            ranknorm_lm = {proportion_model_func <- run_ranknorm_lm},
            arcsin_lm = {proportion_model_func <- run_arcsin_lm}
        )
        otherwise_df <- tibble(term = NA, Estimate = NA, "Std. Error" = NA, "z value" = NA, "Pr(>|z|)" = NA)
        possibly_run_proportion_model <- possibly(proportion_model_func, otherwise = otherwise_df)

        if (proportion_model == "betabinomial_vglm") {
            outcome <- "cbind(n_cells, total_cells - n_cells)"
        } else {
            outcome <- "prop_cells"
        }

        celltype_results <- group_by(sample_abundances_df, celltype) %>% 
            group_modify(proportion_model_func, 
                outcome = outcome,
                model_string = model_string,
                metric_name = "Differential abundance"
            )
    } else {
        if (is_in(proportion_model, c("fisher_test, prop_test"))) {
            #if (length(model_string_split) > 1 ) {
                #return("Error: fisher_contingency and chisq_contingency do not support covariates")
            #}
            #if (!is.factor(cell_metadata[[model_string_split[1]]])) {
                #return("Error: the predictor for fisher_contingency and chisq_contingency must be a factor")
            #}
            #predictor <- model_string_split[1]
            #if (meta_contingency) {
                #control_level <- levels(cell_metadata[[predictor]])[1]
                #metadata_control <- filter(cell_metadata, .data[[predictor]] == control_level)
                #metadata_other <- filter(cell_metadata, .data[[predictor]] != control_level)
                #cell_counts_control <- group_by(metadata_control, Celltype) %>% summarise(n = n())
                #cell_counts_other <- group_by(metadata_other, Celltype, Sample, {{ predictor }}) %>%
                    #group_modify(contingency_test, 
                    #)
            #}
        } else {
            celltype_results <- unique(cell_metadata$celltype) %>%
                map_df(run_binomial_glm,
                    outcome = "celltype",
                    model_string = model_string,
                    model_data = cell_metadata,
                    mixed_effects = mixed_effects
                )
        }
    }
    celltype_results
}

cluster_metrics_model <- function(sce_object, model_string, aggregate = FALSE, proportion_model = "betabinomial_vglm", counts_model = "negbinomial_glm", silhouette_model = "fisherzlink_glm") {
    #if (!is_in(proportion_model, c("betabinomial_vglm", "beta_glm", "ranknorm_lm", "logit_lm", "arcsin_lm"))) {
        #return("Error: proportion_model must be one of betabinomial_vglm, beta_glm, ranknorm_lm, logit_lm or arcsin_lm")
    #}
    #if (!is_in(counts_model, c("negbinomial_glm, log_lm, ranknorm_lm"))) {
        #return("Error: counts_model must be one of negbinomial_glm, log_lm, ranknorm_lm")
    #}
    #if (!is_in(silhouette_model, c("fisherzlink_glm, ranknorm_lm"))) {
        #return("Error: silhouette_model must be one of fisherzlink_lm ranknorm_lm")
    #}
    #if (!is_in(proportion_model, c("beta_glm", "ranknorm_lm", "logit_lm", "arcsin_lm")) && aggregate == FALSE) {
        #return("Error: aggregate must be TRUE when proportion_model is any of beta_glm, ranknorm_lm, logit_lm, arcsin_lm")
    #}

    # If aggregate is set to TRUE, pseudobulk metrics per sample similar to what is done for differential expression.
    # Percent MT reads and total UMIs are directly summed.
    # Detected UMIs are the union of UMIs detected in each cell.
    # Silhouette width is mean silhouette width across cells.

    cell_metadata <- colData(sce_object) %>% as_tibble

    if (aggregate == TRUE) {
        model_string_split <- str_remove_all(model_string, " ") %>% str_split("\\+") %>% unlist
        metadata_reduce <- dplyr::select(cell_metadata, Sample, celltype, one_of(model_string_split)) %>% distinct()
        mt_percent_celltype <- group_by(cell_metadata, Sample, celltype) %>% 
            summarise(subsets_percent_mt_sum = sum(subsets_percent_mt_sum), sum = sum(sum)) %>% 
            left_join(metadata_reduce)
        silhouette_celltype <- group_by(cell_metadata, Sample, celltype) %>% 
            group_modify(avg_silhouette) %>% 
            left_join(metadata_reduce)
        sum_celltype <- group_by(cell_metadata, Sample, celltype) %>% 
            summarise(sum = sum(sum)) %>%
            left_join(metadata_reduce)
        detected_celltype <- group_by(cell_metadata, Sample, celltype) %>% 
            group_modify(combine_detected, sce_object) %>% 
            left_join(metadata_reduce)
        silhouette_outcome <- "mean_silhouette_width"
    } else {
        mt_percent_celltype <- cell_metadata
        silhouette_celltype <- cell_metadata
        sum_celltype <- cell_metadata
        detected_celltype <- cell_metadata
        silhouette_outcome <- "silhouette_width"
    }

    if (proportion_model == "betabinomial_vglm") {
        mt_percent_outcome <- "cbind(subsets_percent_mt_sum, sum - subsets_percent_mt_sum)"
    } else {
        mt_percent_outcome <- "subsets_percent_mt_sum / sum"
    }

    # Assign the modelling function used for each type of metric.
    # Percent MT uses proportion_model
    # Total UMIs and uniquely detected UMIs uses counts_model
    # Silhouette width uses silhouette_model
    switch(proportion_model,
        betabinomial_vglm = {proportion_model_func <- run_betabinomial_vglm},
        beta_glm = {proportion_model_func <- run_beta_vglm},
        logit_lm = {proportion_model_func <- run_logit_lm},
        ranknorm_lm = {proportion_model_func <- run_ranknorm_lm},
        arcsin_lm = {proportion_model_func <- run_arcsin_lm}
    )
    switch(counts_model,
        negbinomial_glm = {counts_model_func <- run_negbinomial_glm},
        log_glm = {counts_model_func <- run_log_lm},
        ranknorm_lm = {counts_model_func <- run_ranknorm_lm}
    )
    switch(silhouette_model,
        fisherzlink_glm = {silhouette_model_func <- run_fisherzlink_glm},
        ranknorm_glm = {silhouette_model_func <- run_ranknorm_glm}
    )

    # Wrap calls to each model in possibly to gracefully capture models that failed to converge.
    # Set a "fallback" dataframe with the otherwise argument to fill in NAs for failed models.
    otherwise_df <- tibble(term = NA, Estimate = NA, "Std. Error" = NA, "z value" = NA, "Pr(>|z|)" = NA)
    possibly_run_proportion_model <- possibly(proportion_model_func, otherwise = otherwise_df)
    possibly_run_counts_model <- possibly(counts_model_func, otherwise = otherwise_df)
    possibly_run_silhouette_model <- possibly(silhouette_model_func, otherwise = otherwise_df)

    # Call wrapped functions for each metric
    mt_percent_summary_df <- group_by(mt_percent_celltype, celltype) %>%
        group_modify(possibly_run_proportion_model,
            outcome = mt_percent_outcome,
            model_string = model_string,
            metric_name = "Proportion reads in MT genes"
        )
    sum_vglm_summary_df <- group_by(sum_celltype, celltype) %>%
        group_modify(possibly_run_counts_model,
            outcome = "sum",
            model_string = model_string,
            metric_name = "Total UMIs"
        )
    detected_vglm_summary_df <- group_by(detected_celltype, celltype) %>%
        group_modify(possibly_run_counts_model,
            outcome = "detected",
            model_string = model_string,
            metric_name = "Uniquely detected UMIs"
        )
    silhouette_vglm_summary_df <- group_by(detected_celltype, celltype) %>%
        group_modify(possibly_run_silhouette_model,
            outcome = silhouette_outcome,
            model_string = model_string,
            metric_name = "Silhouette width"
        )

    # Combine outputs with bind_rows
    metadata_summary <- bind_rows(mt_percent_summary_df, silhouette_vglm_summary_df, sum_vglm_summary_df, detected_vglm_summary_df)
    metadata_summary
}

cluster_metrics <- function(sce_object, model_string, aggregate = FALSE) {
    metadata_df <- colData(sce_object) %>% as_tibble()
    model_string_split <- str_remove_all(model_string, " ") %>% str_split("\\+") %>% unlist
    metrics_long <- dplyr::select(metadata_df, Barcode, Sample, CellID, Celltype, sum, detected, subsets_percent_mt_sum, silhouette_width, one_of(model_string_split))

    metrics_vglm <- group_by(metrics_long, Celltype) %>% group_modify(cluster_metrics_model, model_string, aggregate)
    metrics_vglm
}

contingency_test <- function(metadata_df, celltype_df, control_df) {
    control_celltype_df <- dplyr::filter(control_df, Celltype == celltype_df$Celltype[1])
    total_control <- sum(control_df$n)
    total_celltype_control <- control_celltype_df$n[1]
    total_notcelltype_control <- total_control - total_celltype_control

    total_chip <- metadata_df$total_cells[1]
    total_celltype_chip <- metadata_df$n[1]
    total_notcelltype_chip <- total_chip - total_celltype_chip

    contingency_table <- c(total_celltype_chip, total_notcelltype_chip, total_celltype_control, total_notcelltype_control) %>% matrix(ncol = 2) %>% t
    fisher_result <- fisher.test(contingency_table)
    fisher_result_df <- tibble(
        OR = fisher_result$estimate,
        OR_CI95lo = fisher_result$conf.int[1],
        OR_CI95hi = fisher_result$conf.int[2],
        p_value = fisher_result$p.value
    )
    fisher_result_df$Z <- qnorm(fisher_result_df$p_value/2, lower.tail = FALSE) * sign(log(fisher_result_df$OR))
    fisher_result_df
}

meta_results <- function(grouped_rows, grouping_df, te_column, se_column) {
    meta_output <- meta::metagen(grouped_rows[[te_column]], grouped_rows[[se_column]])
    meta_output_df <- tibble(te_fixed = meta_output$TE.common,
        se_te_fixed = meta_output$seTE.common,
        pval_fixed = meta_output$pval.common,
        te_random = meta_output$TE.random,
        se_te_random = meta_output$seTE.random,
        pval_random = meta_output$pval.random,
    )
    meta_output_df
}

get_sizes <- function(objects, object_names) {
    formatted_sizes <- map_chr(objects, format, units = "auto")
    names(formatted_sizes) <- object_names
    raw_sizes <- reduce(objects, c)
    names(raw_sizes) <- object_names
    raw_sizes_sorted <- sort(raw_sizes)
    formatted_sizes_sorted <- formatted_sizes[match(names(raw_sizes_sorted), names(formatted_sizes))]
    formatted_sizes_sorted
}

sample_dirs_AGS1 <- c(
    AGS10 = "/oak/stanford/groups/sjaiswal/Herra/210310_A00742_0209_BHY7NWDSXY/cellranger_output/AGS10/outs/filtered_feature_bc_matrix",
    AGS11 = "/oak/stanford/groups/sjaiswal/Herra/210310_A00742_0209_BHY7NWDSXY/cellranger_output/AGS11/outs/filtered_feature_bc_matrix",
    AGS12 = "/oak/stanford/groups/sjaiswal/Herra/210310_A00742_0209_BHY7NWDSXY/cellranger_output/AGS12/outs/filtered_feature_bc_matrix",
    AGS13 = "/oak/stanford/groups/sjaiswal/Herra/210310_A00742_0209_BHY7NWDSXY/cellranger_output/AGS13/outs/filtered_feature_bc_matrix",
    AGS14 = "/oak/stanford/groups/sjaiswal/Herra/210310_A00742_0209_BHY7NWDSXY/cellranger_output/AGS14/outs/filtered_feature_bc_matrix",
    AGS15 = "/oak/stanford/groups/sjaiswal/Herra/210310_A00742_0209_BHY7NWDSXY/cellranger_output/AGS15/outs/filtered_feature_bc_matrix",
    AGS16 = "/oak/stanford/groups/sjaiswal/Herra/210310_A00742_0209_BHY7NWDSXY/cellranger_output/AGS16/outs/filtered_feature_bc_matrix",
    AGS18 = "/oak/stanford/groups/sjaiswal/Herra/210310_A00742_0209_BHY7NWDSXY/cellranger_output/AGS18/outs/filtered_feature_bc_matrix",
    AGS19 = "/oak/stanford/groups/sjaiswal/Herra/210310_A00742_0209_BHY7NWDSXY/cellranger_output/AGS19/outs/filtered_feature_bc_matrix",
    AGS20 = "/oak/stanford/groups/sjaiswal/Herra/210310_A00742_0209_BHY7NWDSXY/cellranger_output/AGS20/outs/filtered_feature_bc_matrix"
)

sample_dirs_AGS2 <- c(
    AGS_1B = "/oak/stanford/groups/sjaiswal/Herra/H202SC23041042/H5KHVDSX7/outs/fastq_path/H5KHVDSX7/cellranger_output_1A/1A_Gene_S8/outs/filtered_feature_bc_matrix"
)

sample_dirs_AGS3 <- c(
    AGS_2B = "/oak/stanford/groups/sjaiswal/Herra/H202SC23041042/H5KHVDSX7/outs/fastq_path/H5KHVDSX7/cellranger_output_2A_31A_12A/outs/per_sample_outs/2B/count/sample_filtered_feature_bc_matrix",
    AGS_12B = "/oak/stanford/groups/sjaiswal/Herra/H202SC23041042/H5KHVDSX7/outs/fastq_path/H5KHVDSX7/cellranger_output_2A_31A_12A/outs/per_sample_outs/12B/count/sample_filtered_feature_bc_matrix"
)

sample_dirs_AGS4 <- c(
    AGS_3B = "/oak/stanford/groups/sjaiswal/Herra/H202SC23041042/H5KHVDSX7/outs/fastq_path/H5KHVDSX7/cellranger_output_3A_28A/outs/per_sample_outs/3B/count/sample_filtered_feature_bc_matrix",
    AGS_28B = "/oak/stanford/groups/sjaiswal/Herra/H202SC23041042/H5KHVDSX7/outs/fastq_path/H5KHVDSX7/cellranger_output_3A_28A/outs/per_sample_outs/28B/count/sample_filtered_feature_bc_matrix"
)

sample_dirs_AGS5 <- c(
    AGS_4B = "/oak/stanford/groups/sjaiswal/Herra/H202SC23041042/H5KHVDSX7/outs/fastq_path/H5KHVDSX7/cellranger_output_4A_25A/outs/per_sample_outs/4B/count/sample_filtered_feature_bc_matrix",
    AGS_25B = "/oak/stanford/groups/sjaiswal/Herra/H202SC23041042/H5KHVDSX7/outs/fastq_path/H5KHVDSX7/cellranger_output_4A_25A/outs/per_sample_outs/25B/count/sample_filtered_feature_bc_matrix"
)

sample_dirs_AGS6 <- c(
    AGS_5B = "/oak/stanford/groups/sjaiswal/Herra/H202SC23041042/H5KHVDSX7/outs/fastq_path/H5KHVDSX7/cellranger_output_5A_29A_N/outs/per_sample_outs/5B/count/sample_filtered_feature_bc_matrix",
    AGS_29B = "/oak/stanford/groups/sjaiswal/Herra/H202SC23041042/H5KHVDSX7/outs/fastq_path/H5KHVDSX7/cellranger_output_5A_29A_N/outs/per_sample_outs/29B/count/sample_filtered_feature_bc_matrix"
)

sample_dirs_AGS7 <- c(
    AGS_6B = "/oak/stanford/groups/sjaiswal/Herra/H202SC23041042/H5KHVDSX7/outs/fastq_path/H5KHVDSX7/cellranger_output_7A_6A_30A_16A/outs/per_sample_outs/6B/count/sample_filtered_feature_bc_matrix",
    AGS_7B = "/oak/stanford/groups/sjaiswal/Herra/H202SC23041042/H5KHVDSX7/outs/fastq_path/H5KHVDSX7/cellranger_output_7A_6A_30A_16A/outs/per_sample_outs/7B/count/sample_filtered_feature_bc_matrix",
    AGS_16B = "/oak/stanford/groups/sjaiswal/Herra/H202SC23041042/H5KHVDSX7/outs/fastq_path/H5KHVDSX7/cellranger_output_7A_6A_30A_16A/outs/per_sample_outs/16B/count/sample_filtered_feature_bc_matrix",
    AGS_30B = "/oak/stanford/groups/sjaiswal/Herra/H202SC23041042/H5KHVDSX7/outs/fastq_path/H5KHVDSX7/cellranger_output_7A_6A_30A_16A/outs/per_sample_outs/30B/count/sample_filtered_feature_bc_matrix"
)

sample_dirs_AGS8 <- c(
    AGS_8B = "/oak/stanford/groups/sjaiswal/Herra/H202SC23041042/H5KHVDSX7/outs/fastq_path/H5KHVDSX7/cellranger_output_8A_19A/outs/per_sample_outs/8B/count/sample_filtered_feature_bc_matrix",
    AGS_19B = "/oak/stanford/groups/sjaiswal/Herra/H202SC23041042/H5KHVDSX7/outs/fastq_path/H5KHVDSX7/cellranger_output_8A_19A/outs/per_sample_outs/19B/count/sample_filtered_feature_bc_matrix"
)

sample_dirs_AGS9 <- c(
    AGS_11B = "/oak/stanford/groups/sjaiswal/Herra/H202SC23041042/H5KHVDSX7/outs/fastq_path/H5KHVDSX7/cellranger_output_15A_11A/outs/per_sample_outs/11B/count/sample_filtered_feature_bc_matrix",
    AGS_15B = "/oak/stanford/groups/sjaiswal/Herra/H202SC23041042/H5KHVDSX7/outs/fastq_path/H5KHVDSX7/cellranger_output_15A_11A/outs/per_sample_outs/15B/count/sample_filtered_feature_bc_matrix"
)

sample_dirs_AGS10 <- c(
    AGS_21B = "/oak/stanford/groups/sjaiswal/Herra/H202SC23041042/H5KHVDSX7/outs/fastq_path/H5KHVDSX7/cellranger_output_22A_21A/outs/per_sample_outs/21B/count/sample_filtered_feature_bc_matrix",
    AGS_22B = "/oak/stanford/groups/sjaiswal/Herra/H202SC23041042/H5KHVDSX7/outs/fastq_path/H5KHVDSX7/cellranger_output_22A_21A/outs/per_sample_outs/22B/count/sample_filtered_feature_bc_matrix"
)

sample_dirs_AGS11 <- c(
    AGS_24B = "/oak/stanford/groups/sjaiswal/Herra/H202SC23041042/H5KHVDSX7/outs/fastq_path/H5KHVDSX7/cellranger_output_27A_24A/outs/per_sample_outs/24B/count/sample_filtered_feature_bc_matrix",
    AGS_27B = "/oak/stanford/groups/sjaiswal/Herra/H202SC23041042/H5KHVDSX7/outs/fastq_path/H5KHVDSX7/cellranger_output_27A_24A/outs/per_sample_outs/27B/count/sample_filtered_feature_bc_matrix"
)

# Read data and add gene annotation information
# Batch 1 were samples run on standard 10X libraries for each sample.
# Batch 2 are all samples run together in a hashtag multiplexing library.
 # All batch 1 samples are loaded by this line.
AGS_sce1 <- DropletUtils::read10xCounts(sample_dirs_AGS1)
# Each sample in batch 2 is loaded separately
AGS_sce2 <- DropletUtils::read10xCounts(sample_dirs_AGS2) 
AGS_sce3 <- DropletUtils::read10xCounts(sample_dirs_AGS3)
AGS_sce4 <- DropletUtils::read10xCounts(sample_dirs_AGS4)
AGS_sce5 <- DropletUtils::read10xCounts(sample_dirs_AGS5)
AGS_sce6 <- DropletUtils::read10xCounts(sample_dirs_AGS6)
AGS_sce7 <- DropletUtils::read10xCounts(sample_dirs_AGS7)
AGS_sce8 <- DropletUtils::read10xCounts(sample_dirs_AGS8)
AGS_sce9 <- DropletUtils::read10xCounts(sample_dirs_AGS9)
AGS_sce10 <- DropletUtils::read10xCounts(sample_dirs_AGS10)
AGS_sce11 <- DropletUtils::read10xCounts(sample_dirs_AGS11)

# Identify shared gene names across samples.  This should be very similar as all samples use the same reference.
AGS_batch2_rownames <- list(
    rownames(AGS_sce2),
    rownames(AGS_sce3),
    rownames(AGS_sce4),
    rownames(AGS_sce5),
    rownames(AGS_sce6),
    rownames(AGS_sce7),
    rownames(AGS_sce8),
    rownames(AGS_sce9),
    rownames(AGS_sce10),
    rownames(AGS_sce11)
) %>% reduce(intersect)

# Combine samples for batch 2 into one SingleCellExperiment object because they had to be loaded separately. Only keep the shared genes.
AGS_batch2_sce <- cbind(
    AGS_sce2[AGS_batch2_rownames,],
    AGS_sce3[AGS_batch2_rownames,],
    AGS_sce4[AGS_batch2_rownames,],
    AGS_sce5[AGS_batch2_rownames,],
    AGS_sce6[AGS_batch2_rownames,],
    AGS_sce7[AGS_batch2_rownames,],
    AGS_sce8[AGS_batch2_rownames,],
    AGS_sce9[AGS_batch2_rownames,],
    AGS_sce10[AGS_batch2_rownames,],
    AGS_sce11[AGS_batch2_rownames,]
)

# Load annotation for genes.  Used primarily to get mitochondrial gene names from Ensembl IDs, but also useful for annotating differential expression results later. 
annotation_hub <- AnnotationHub::AnnotationHub()
ens_104 <- annotation_hub[["AH95744"]] # Load Ensembl annotation version 104
annotation_columns <- c( 
    "GENEID",
    "GENEBIOTYPE",
    "DESCRIPTION",
    "SEQNAME",
    "GENESEQSTART",
    "GENESEQEND",
    "SEQSTRAND",
    "CANONICALTRANSCRIPT",
    "GENEIDVERSION",
    "ENTREZID"
)
# Annotate gene metadata with Entrez IDs for downstream analysis and properly handle duplicate Entrez IDs
gene_annot_batch1 <- AnnotationDbi::select(ens_104, keys = rownames(AGS_sce1), keytype = "GENEID", columns = annotation_columns) %>% 
    deduplicate_ids %>% 
    pivot_wider(names_from = "new_names", values_from = "ENTREZID")
colnames(gene_annot_batch1)[1] <- "ID"
gene_annot_batch2 <- AnnotationDbi::select(ens_104, keys = rownames(AGS_batch2_sce), keytype = "GENEID", columns = annotation_columns) %>% 
    deduplicate_ids %>% 
    pivot_wider(names_from = "new_names", values_from = "ENTREZID")
colnames(gene_annot_batch2)[1] <- "ID"

# Add sample phenotype information 
AGS_pheno <- read_csv("../AGS_all_metadata.csv") 
colnames(AGS_pheno)[1] <- "Sample"
colnames(AGS_pheno)[14:16] <- str_c(colnames(AGS_pheno)[14:16], "_gt")

# Join samples phenotype data to colData and replace it in the SCE object
colData(AGS_sce1) %<>% as.data.frame %>% left_join(AGS_pheno) %>% DataFrame
colData(AGS_batch2_sce) %<>% as.data.frame %>% left_join(AGS_pheno) %>% DataFrame

# Do a first pass of filtering, PCA, doublet estimation, UMAP, and clustering on each batch separately
tic("Processing batch 1")
AGS_sce_leiden_batch1 <- filter_osca(AGS_sce1, gene_annot_batch1) %>%
    estimate_corral() %>%
    compute_doublets() %>%
    estimate_umap() %>% 
    leiden_clustering()
toc()

tic("Processing batch 2")
AGS_sce_leiden_batch2 <- filter_osca(AGS_batch2_sce, gene_annot_batch2) %>%
    estimate_corral() %>%
    compute_doublets() %>%
    estimate_umap() %>%
    leiden_clustering()
toc()

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

# First pass of clustering and plotting is to identify non-hematopoietic cell clusters and remove them
catch <- map(marker_genes, umap_plot, "all_cells_batch1/", AGS_sce_leiden_batch1)
catch <- map(marker_genes, umap_plot, "all_cells_batch2/", AGS_sce_leiden_batch2)

qc_plots(AGS_sce_leiden_batch1, "all_cells_batch1/")
qc_plots(AGS_sce_leiden_batch2, "all_cells_batch2/")

# Remove cells which have markers of non-immune cells
AGS_sce_immune_batch1 <- dplyr::filter(AGS_sce_leiden_batch1, !is_in(label, c(19)))
AGS_sce_immune_batch2 <- dplyr::filter(AGS_sce_leiden_batch2, !is_in(label, c(14,22,18)))

# Repeat PCA, UMAP and clustering on data after removing non-hematopoietic cells
AGS_sce_immune_leiden_batch1 <- estimate_corral(AGS_sce_immune_batch1) %>% 
    estimate_umap %>% 
    leiden_clustering
AGS_sce_immune_leiden_batch2 <- estimate_corral(AGS_sce_immune_batch2) %>% 
    estimate_umap %>% 
    leiden_clustering

qc_plots(AGS_sce_immune_leiden_batch1, "immune_cells_batch1/")
qc_plots(AGS_sce_immune_leiden_batch2, "immune_cells_batch2/")

catch <- map(marker_genes, umap_plot, "immune_cells_batch1/", AGS_sce_immune_leiden_batch1)
catch <- map(marker_genes, umap_plot, "immune_cells_batch2/", AGS_sce_immune_leiden_batch2)

# Load reference data and use SingleR to classify cells
monaco <- celldex::MonacoImmuneData()
lm22_se <- read_rds("../../new_code/lm22_se.rda")

AGS_sce_singler_batch1 <- singler_annotation(AGS_sce_immune_leiden_batch1, monaco, "monaco_celltypes", "immune_cells_batch1/monaco_labels") %>% 
    singler_annotation(lm22_se, "lm22_celltypes", "immune_cells_batch1/lm22_labels")
AGS_sce_singler_batch2 <- singler_annotation(AGS_sce_immune_leiden_batch2, monaco, "monaco_celltypes", "immune_cells_batch2/monaco_labels") %>% 
    singler_annotation(lm22_se, "lm22_celltypes", "immune_cells_batch2/lm22_labels")

# Plot SingleR results
umap_plot("monaco_celltypes", "immune_cells_batch1/", AGS_sce_singler_batch1, continuous_scale = FALSE, gene_expression = FALSE, add_legend = FALSE, text_by = "monaco_celltypes", text_size = 3)
umap_plot("monaco_celltypes", "immune_cells_batch2/", AGS_sce_singler_batch2, continuous_scale = FALSE, gene_expression = FALSE, add_legend = FALSE, text_by = "monaco_celltypes", text_size = 3)
umap_plot("lm22_celltypes", "immune_cells_batch1/", AGS_sce_singler_batch1, continuous_scale = FALSE, gene_expression = FALSE, add_legend = FALSE, text_by = "lm22_celltypes", text_size = 3)
umap_plot("lm22_celltypes", "immune_cells_batch2/", AGS_sce_singler_batch2, continuous_scale = FALSE, gene_expression = FALSE, add_legend = FALSE, text_by = "lm22_celltypes", text_size = 3)

AGS_markers_batch1 <- scran::findMarkers(AGS_sce_singler_batch1, test = "t", direction = "up", pval.type = "some", assay.type = "logcounts", row.data = rowData(AGS_sce_singler_batch1))
AGS_markers_batch2 <- scran::findMarkers(AGS_sce_singler_batch2, test = "t", direction = "up", pval.type = "some", assay.type = "logcounts", row.data = rowData(AGS_sce_singler_batch2))

umap_plot("label", "immune_cells_batch1/", AGS_sce_singler_batch1, suffix = "_legend", continuous_scale = FALSE, gene_expression = FALSE, text_by = "label")
umap_plot("label", "immune_cells_batch2/", AGS_sce_singler_batch2, suffix = "_legend", continuous_scale = FALSE, gene_expression = FALSE, text_by = "label")

# Assign rough cell types
celltype_labels_batch1 <- c(
    "1" = "T-cells/NK cells",
    "2" = "T-cells/NK cells",
    "3" = "T-cells/NK cells",
    "4" = "T-cells/NK cells",
    "5" = "XCR1+ DCs",
    "6" = "T-cells/NK cells",
    "7" = "T-cells/NK cells",
    "8" = "Macrophages/cDCs",
    "9" = "Macrophages/cDCs",
    "10" = "T-cells/NK cells",
    "11" = "T-cells/NK cells",
    "12" = "Macrophages/cDCs",
    "13" = "B-cells",
    "14" = "pDCs",
    "15" = "Plasma cells",
    "16" = "T-cells/NK cells",
    "17" = "Basophils",
    "18" = "T-cells/NK cells",
    "19" = "T-cells/NK cells",
    "20" = "T-cells/NK cells"
)

celltype_labels_batch2 <- c(
    "1" = "Macrophages c/DCs",
    "2" = "Macrophages/cDCs",
    "3" = "Macrophages/cDCs",
    "4" = "B-cells",
    "5" = "T-cells/NK cells",
    "6" = "T-cells/NK cells",
    "7" = "T-cells/NK cells",
    "8" = "pDCs",
    "9" = "T-cells/NK cells",
    "10" = "T-cells/NK cells",
    "11" = "T-cells/NK cells",
    "12" = "T-cells/NK cells",
    "13" = "Macrophages/cDCs",
    "14" = "Macrophages/cDCs",
    "15" = "T-cells/NK cells",
    "16" = "T-cells/NK cells",
    "17" = "Macrophages/cDCs"
)

AGS_sce_celltypes_batch1 <- add_celltype_labels(AGS_sce_immune_leiden_batch1, celltype_labels_batch1)
AGS_sce_celltypes_batch2 <- add_celltype_labels(AGS_sce_immune_leiden_batch2, celltype_labels_batch2)

# Separate out macrophages and T-cells
colData(AGS_sce_celltypes_batch1)$CellID <- str_c(colData(AGS_sce_celltypes_batch1)$Sample, colData(AGS_sce_celltypes_batch1)$Barcode, sep = ":")
AGS_sce_filtered_macros_batch1 <- dplyr::filter(AGS_sce_celltypes_batch1, celltype == "Macrophages/cDCs")
write_rds(AGS_sce_filtered_macros_batch1, "AGS_sce_filtered_macros_batch1.rda")
AGS_sce_filtered_tcells_batch1 <- dplyr::filter(AGS_sce_celltypes_batch1, celltype == "T-cells/NK cells")
write_rds(AGS_sce_filtered_tcells_batch1, "AGS_sce_filtered_tcells_batch1.rda")

colData(AGS_sce_celltypes_batch2)$CellID <- str_c(colData(AGS_sce_celltypes_batch2)$Sample, colData(AGS_sce_celltypes_batch2)$Barcode, sep = ":")
AGS_sce_filtered_macros_batch2 <- dplyr::filter(AGS_sce_celltypes_batch2, celltype == "Macrophages/cDCs")
write_rds(AGS_sce_filtered_macros_batch2, "AGS_sce_filtered_macros_batch2.rda")
AGS_sce_filtered_tcells_batch2 <- dplyr::filter(AGS_sce_celltypes_batch2, celltype == "T-cells/NK cells")
write_rds(AGS_sce_filtered_tcells_batch2, "AGS_sce_filtered_tcells_batch2.rda")

# Select cells which are not macrophages or T-cells
AGS_sce_filtered_other_batch1_cells <- !str_detect(colData(AGS_sce_celltypes_batch1)$celltype, fixed("T-cells")) & !str_detect(colData(AGS_sce_celltypes_batch1)$celltype, fixed("Macrophages"))
AGS_sce_filtered_other_batch2_cells <- !str_detect(colData(AGS_sce_celltypes_batch2)$celltype, fixed("T-cells")) & !str_detect(colData(AGS_sce_celltypes_batch2)$celltype, fixed("Macrophages"))

# Load subclustered macrophages and T-cells
AGS_macros_annot <- read_rds("../macrophage_clustering/AGS_sce_celltypes.rda")
AGS_tcells_annot <- read_rds("../tcell_clustering/AGS_sce_celltypes.rda")
subclustered_cells <- union(colData(AGS_macros_annot)$CellID, colData(AGS_tcells_annot)$CellID)

# Only keep cells which were kept in subclustering or were not macrophages or T-cells
AGS_sce_subclustered_batch1 <- dplyr::filter(AGS_sce_celltypes_batch1, AGS_sce_filtered_other_batch1_cells | is_in(CellID, subclustered_cells))
AGS_sce_subclustered_batch2 <- dplyr::filter(AGS_sce_celltypes_batch2, AGS_sce_filtered_other_batch2_cells | is_in(CellID, subclustered_cells))

AGS_sce_list <- list("1" = AGS_sce_subclustered_batch1, "2" = AGS_sce_subclustered_batch2)
AGS_sce_filtered_merge <- merge_batches(AGS_sce_list)
write_rds(AGS_sce_filtered_merge, "AGS_sce_filtered_merge.rda")

AGS_combined_filter <-  estimate_umap(AGS_sce_filtered_merge, "corrected") %>%
    leiden_clustering("corrected")
qc_plots(AGS_combined_filter, "subclustered_cells/", doublets = F)
catch <- map(marker_genes, umap_plot, "subclustered_cells/", AGS_combined_filter)
AGS_markers_combined_filter <- scran::findMarkers(AGS_combined_filter, test = "t", direction = "up", pval.type = "some", assay.type = "logcounts", row.data = rowData(AGS_combined_filter))
AGS_markers_combined_filter_df <- as.list(AGS_markers_combined_filter) %>% map(as_tibble) %>% bind_rows(.id = "Cluster")
write_rds(AGS_markers_combined_filter_df, "AGS_markers_combined_filter_df.rds")

all_samples_celltype_labels <- c(
    "1" = "Active NK cells",
    "2" = "Cytotoxic T-cells",
    "3" = "Macrophages",
    "4" = "Effector memory T-cells",
    "5" = "Macrophages",
    "6" = "Macrophages",
    "7" = "Resting NK cells",
    "8" = "B-cells",
    "9" = "Plasmacytoid dendritic cells",
    "10" = "Plasma cells",
    "11" = "Effector memory T-cells",
    "12" = "Macrophages",
    "13" = "TCL1A+ T-cells",
    "14" = "Basophils",
    "15" = "Unknown",
    "16" = "Effector memory T-cells",
    "17" = "Effector memory T-cells",
    "18" = "T-regs",
    "19" = "Effector memory T-cells",
    "20" = "Macrophages",
    "21" = "Junk"
)

AGS_combined_filter_celltype <- add_celltype_labels(AGS_combined_filter, all_samples_celltype_labels, new_celltypes = TRUE)
AGS_silhouette <- reducedDim(AGS_combined_filter_celltype, "corrected") %>% 
    bluster::approxSilhouette(AGS_combined_filter_celltype$celltype) %>% as_tibble
colData(AGS_combined_filter_celltype)$silhouette_width <- AGS_silhouette$width
celltype_plots(AGS_combined_filter_celltype, "subclustered_cells/")

write_rds(AGS_combined_filter_celltype, "AGS_combined_filter_celltype.rds")

AGS_combined_filter_celltype <- read_rds("./AGS_combined_filter_celltype.rds")

# Cell composition analysis
colData(AGS_combined_filter_celltype)$Mutation_1[AGS_combined_filter_celltype$VAF_1 < 0.02] <- NA
AGS_combined_filter_celltype$VAF_1[AGS_combined_filter_celltype$VAF_1 < 0.02] <- NA
colData(AGS_combined_filter_celltype)$Mutation_1[colData(AGS_combined_filter_celltype)$Mutation_1 == "PPM1D"] <- "DDR"
colData(AGS_combined_filter_celltype)$Mutation_1[colData(AGS_combined_filter_celltype)$Mutation_1 == "TP53"] <- "DDR"
colData(AGS_combined_filter_celltype)$Mutation_1 %<>% replace_na("Control") %>% factor(levels = c("Control", "DNMT3A", "TET2", "ASXL1", "DDR"))
colData(AGS_combined_filter_celltype)$VAF_1 %<>% replace_na(0)

colData(AGS_combined_filter_celltype)$chip <- colData(AGS_combined_filter_celltype)$Mutation_1 != "Control"

colData(AGS_combined_filter_celltype)$chip_ddr <- as.character(colData(AGS_combined_filter_celltype)$Mutation_1)
colData(AGS_combined_filter_celltype)$chip_ddr[!is_in(colData(AGS_combined_filter_celltype)$chip_ddr, c("Control", "DDR"))] <- "CHIP"
colData(AGS_combined_filter_celltype)$chip_ddr %<>% replace_na("Control") %>% factor(levels = c("Control", "CHIP", "DDR"))

AGS_combined_filter_celltype_noASXL1 <- dplyr::filter(AGS_combined_filter_celltype, Mutation_1 != "ASXL1")

cell_composition_mutation <- cell_composition_da(AGS_combined_filter_celltype, "Mutation_1 + batch")
cell_composition_mutation_noASXL1 <- cell_composition_da(AGS_combined_filter_celltype_noASXL1, "Mutation_1 + batch")
cell_composition_mutation_binomial <- cell_composition_da(AGS_combined_filter_celltype, "Mutation_1 + batch + Sample", proportion_model = "binomial_glm")
cell_composition_mutation_binomial_mixed <- cell_composition_da(AGS_combined_filter_celltype, "Mutation_1 + batch + (1|Sample)", proportion_model = "binomial_glm", mixed_effects = TRUE)
cell_composition_chip <- cell_composition_da(AGS_combined_filter_celltype, "chip + batch")
cell_composition_chip_ddr <- cell_composition_da(AGS_combined_filter_celltype, "chip_ddr + batch")

AGS_combined_filter_celltype_chip <- dplyr::filter(AGS_combined_filter_celltype, VAF_1 > 0)
colData(AGS_combined_filter_celltype_chip)$VAF_1 %<>% RankNorm
cell_composition_mutation_vaf <- cell_composition_da(AGS_combined_filter_celltype_chip, "VAF_1 + batch")
cell_composition_mutation_vaf_ranknorm <- cell_composition_da(AGS_combined_filter_celltype_chip, "VAF_1 + batch", proportion_model = "ranknorm_lm")

AGS_combined_filter_celltype_chip_noAGS14 <- dplyr::filter(AGS_combined_filter_celltype, Sample != "AGS14")
cell_composition_mutation_vaf_noAGS14 <- cell_composition_da(AGS_combined_filter_celltype_chip_noAGS14, "VAF_1 + batch")

# Macrophages only
colData(AGS_macros_annot)$Mutation_1[AGS_macros_annot$VAF_1 < 0.02] <- NA
AGS_macros_annot$VAF_1[AGS_macros_annot$VAF_1 < 0.02] <- NA
colData(AGS_macros_annot)$Mutation_1[colData(AGS_macros_annot)$Mutation_1 == "PPM1D"] <- "DDR"
colData(AGS_macros_annot)$Mutation_1[colData(AGS_macros_annot)$Mutation_1 == "TP53"] <- "DDR"
colData(AGS_macros_annot)$Mutation_1 %<>% replace_na("Control") %>% factor(levels = c("Control", "DNMT3A", "TET2", "ASXL1", "DDR"))
colData(AGS_macros_annot)$chip <- colData(AGS_macros_annot)$Mutation_1 != "Control"
colData(AGS_macros_annot)$chip_ddr <- as.character(colData(AGS_macros_annot)$Mutation_1)
colData(AGS_macros_annot)$chip_ddr[!is_in(colData(AGS_macros_annot)$chip_ddr, c("Control", "DDR"))] <- "CHIP"
colData(AGS_macros_annot)$chip_ddr %<>% replace_na("Control") %>% factor(levels = c("Control", "CHIP", "DDR"))
colData(AGS_macros_annot)$VAF_1 %<>% replace_na(0)
#colData(AGS_macros_annot)$VAF_1 %<>% RankNorm

AGS_macros_noAGS14 <- dplyr::filter(AGS_macros_annot, Sample != "AGS14")
AGS_macros_noASXL1 <- dplyr::filter(AGS_macros_annot, Mutation_1 != "ASXL1")
AGS_macros_nosmall <- dplyr::filter(AGS_macros_annot, chip == FALSE | (chip == TRUE & VAF_1 > 0.05))

cell_composition_mutation_macros_noASXL1 <- cell_composition_da(AGS_macros_noASXL1, "Mutation_1 + batch")
cell_composition_mutation_macros <- cell_composition_da(AGS_macros_annot, "Mutation_1 + batch")
cell_composition_mutation_macros_noAGS14 <- cell_composition_da(AGS_macros_annot, "Mutation_1 + batch")
cell_composition_mutation_macros_binomial_mixed <- cell_composition_da(AGS_macros_annot, "Mutation_1 + batch + (1|Sample)", proportion_model = "binomial_glm", mixed_effect = TRUE)
cell_composition_chip_macros <- cell_composition_da(AGS_macros_annot, "chip + batch")
cell_composition_chip_macros_binomial_mixed <- cell_composition_da(AGS_macros_annot, "chip + batch + (1|Sample)", proportion_model = "binomial_glm", mixed_effect = TRUE)
cell_composition_chip_macros_noAGS14 <- cell_composition_da(AGS_macros_noAGS14, "chip + batch")
cell_composition_chip_ddr_macros <- cell_composition_da(AGS_macros_annot, "chip_ddr + batch")
cell_composition_chip_ddr_macros_binomial_mixed <- cell_composition_da(AGS_macros_annot, "chip_ddr + batch + (1|Sample)", proportion_model = "binomial_glm", mixed_effects = TRUE)
cell_composition_chip_ddr_macros_noAGS14 <- cell_composition_da(AGS_macros_noAGS14, "chip_ddr + batch")
cell_composition_chip_ddr_macros_ranknorm <- cell_composition_da(AGS_macros_annot, "chip_ddr + batch", proportion_model = "ranknorm_lm")
cell_composition_chip_ddr_macros_noAGS14_ranknorm <- cell_composition_da(AGS_macros_noAGS14, "chip_ddr + batch", proportion_model = "ranknorm_lm")

cell_composition_chip_macros_nosmall <- cell_composition_da(AGS_macros_nosmall, "chip + batch")

AGS_macros_annot_chip <- dplyr::filter(AGS_macros_annot, VAF_1 > 0)

cell_composition_mutation_macros_vaf <- cell_composition_da(AGS_macros_annot_chip, "VAF_1 + batch")
cell_composition_mutation_macros_vaf_ranknorm <- cell_composition_da(AGS_macros_annot_chip, "VAF_1 + batch", proportion_model = "ranknorm_lm")
cell_composition_mutation_macros_vaf_binomial_mixed <- cell_composition_da(AGS_macros_annot_chip, "VAF_1 + batch + (1|Sample)", proportion_model = "binomial_glm", mixed_effects = TRUE)

AGS_macros_noAGS14_chip <- dplyr::filter(AGS_macros_annot_chip, Sample != "AGS14")
AGS_macros_nosmall_chip <- dplyr::filter(AGS_macros_annot_chip, VAF_1 <= 0.05)
cell_composition_mutation_macros_vaf_noAGS14 <- cell_composition_da(AGS_macros_noAGS14_chip, "VAF_1 + batch")
cell_composition_mutation_macros_vaf_noAGS14_ranknorm <- cell_composition_da(AGS_macros_noAGS14_chip, "VAF_1 + batch", proportion_model = "ranknorm_lm")
cell_composition_mutation_macros_vaf_nosmall <- cell_composition_da(AGS_macros_nosmall_chip, "VAF_1 + batch")
cell_composition_mutation_macros_nolarge_vaf <- cell_composition_da(AGS_macros_no_large, "VAF_1 + batch")

# Cell composition plots
metadata_scatterplot("Active NK cells", "subclustered_cells/", AGS_combined_filter_celltype_chip, x_var = "VAF_1")
metadata_scatterplot("Inflammatory macrophages", "../macrophage_clustering/all_cells/", AGS_macros_annot_chip, x_var = "VAF_1")
metadata_scatterplot("Foam cells", "../macrophage_clustering/all_cells/", AGS_macros_annot_chip, x_var = "VAF_1")
metadata_scatterplot("Inflammatory macrophages", "../macrophage_clustering/all_cells/", suffix = "_no_AGS14", AGS_macros_noAGS14_chip, x_var = "VAF_1")
metadata_scatterplot("Foam cells", "../macrophage_clustering/all_cells/", suffix = "_no_AGS14", AGS_macros_noAGS14_chip, x_var = "VAF_1")
AGS_macros_no_large <- dplyr::filter(AGS_macros_noAGS14_chip, VAF_1 < 0.2)
metadata_scatterplot("Inflammatory macrophages", "../macrophage_clustering/all_cells/", suffix = "_no_large", AGS_macros_no_large, x_var = "VAF_1")
metadata_scatterplot("Foam cells", "../macrophage_clustering/all_cells/", suffix = "_no_large", AGS_macros_no_large, x_var = "VAF_1")

celltype_proportion_boxplot("Inflammatory macrophages", "../macrophage_clustering/all_cells/", AGS_macros_annot_chip, x_var = "Mutation_1")
celltype_bar_plot("../macrophage_clustering/all_cells/", AGS_macros_annot)

# Cell metrics analysis
cell_metrics_mutation <- cluster_metrics_model(AGS_combined_filter_celltype, "Mutation_1 + batch")
cell_metrics_mutation_agg <- cluster_metrics_model(AGS_combined_filter_celltype, "Mutation_1 + batch", aggregate = TRUE)
cell_metrics_mutation_macros_agg <- cluster_metrics_model(AGS_macros_annot, "Mutation_1 + batch", aggregate = TRUE)

object_names <- ls()
object_sizes_format <- map(object_names, get) %>% map(object.size) %>% get_sizes(object_names)

# Forest plot
# Remove cell types that are low abundance
AGS_cell_composition_vglm_mutation <- dplyr::filter(cell_composition_mutation_noASXL1, !is_in(celltype, c("Basophils", "TCL1A+ T-cells", "Unknown", "Junk", "T-regs"))) 
AGS_cell_composition_vglm_mutation_filter <- dplyr::filter(AGS_cell_composition_vglm_mutation, str_detect(term, "Mutation"))
# Create adusted p-value
AGS_cell_composition_vglm_mutation_filter$p_adjust <- p.adjust(AGS_cell_composition_vglm_mutation_filter$`Pr(>|z|)`, method = "fdr")

format_pvalue <- function(mantissa, exponent) {
    mantissa <- formatC(mantissa, format = "f", digits = 1)
    exponent <- as.character(exponent)
    pvalue_format <- bquote(.(mantissa) ~ "x" ~ 10^.(exponent)) %>% as.expression
    pvalue_format
}

# Convert p-value to scientific notation
AGS_cell_composition_vglm_mutation_filter$exponent <- as.numeric(AGS_cell_composition_vglm_mutation_filter$`Pr(>|z|)`) %>% log10 %>% floor
AGS_cell_composition_vglm_mutation_filter$mantissa <- as.numeric(AGS_cell_composition_vglm_mutation_filter$`Pr(>|z|)`) * 10^(-1 * AGS_cell_composition_vglm_mutation_filter$exponent)
formatted_pvalues <- map2(AGS_cell_composition_vglm_mutation_filter$mantissa, AGS_cell_composition_vglm_mutation_filter$exponent, format_pvalue) 
# Only use scientific notation if p < 1x10^-4
formatted_pvalues[AGS_cell_composition_vglm_mutation_filter$exponent > -4] <- formatC(AGS_cell_composition_vglm_mutation_filter$`Pr(>|z|)`[AGS_cell_composition_vglm_mutation_filter$exponent > -4], format = "f", digits = 3) 

# Create formatted confidence interval
AGS_cell_composition_vglm_mutation_filter$estimate_format <- formatC(AGS_cell_composition_vglm_mutation_filter$OR, format = "f", digits = 2)
AGS_cell_composition_vglm_mutation_filter$conf_low_format <- formatC(AGS_cell_composition_vglm_mutation_filter$OR_conf.low, format = "f", digits = 2)
AGS_cell_composition_vglm_mutation_filter$conf_high_format <- formatC(AGS_cell_composition_vglm_mutation_filter$OR_conf.high, format = "f", digits = 2)
AGS_cell_composition_vglm_mutation_filter$conf_int <- str_glue_data(AGS_cell_composition_vglm_mutation_filter, "{estimate_format} ({conf_low_format}-{conf_high_format})")

# Relevel cell type names if necessary and arrange by cell type name and mutation name
AGS_cell_composition_vglm_mutation_plot <- AGS_cell_composition_vglm_mutation_filter
AGS_cell_composition_vglm_mutation_plot$term %<>% str_remove_all("Mutation_1") %>% factor(levels = c("DNMT3A", "TET2", "DDR"))
AGS_cell_composition_vglm_mutation_plot %<>% arrange(celltype, term) # Arranging order determines plotting order

# Create a list of columns for the forest plot
celltypes_final <- list(list.prepend(as.list(as.character(AGS_cell_composition_vglm_mutation_plot$celltype)), expression(bold("Cell type"))),
    list.prepend(as.list(as.character(AGS_cell_composition_vglm_mutation_plot$term)), expression(bold("Mutation"))),
    list.prepend(as.list(as.character(AGS_cell_composition_vglm_mutation_plot$conf_int)), expression(bold("OR (CI95)"))),
    list.prepend(as.list(formatted_pvalues), expression(bold("P-value"))))
# To make it look nice, make the cell type name blank after the first row of that cell type
celltypes_final[[1]][is_in(celltypes_final[[2]], c("TET2", "DNMT3A"))] <- NA

# Create matrix of odds ratio and confidence interval, and put NAs in top row for column names
celltypes_matrix <- ungroup(AGS_cell_composition_vglm_mutation_plot) %>% dplyr::select(OR, OR_conf.low, OR_conf.high) %>% mutate_all(as.numeric) %>% as.matrix 
celltypes_matrix_plot <- rbind(c(NA, NA, NA), celltypes_matrix)

# Create PDF of forest plot
cairo_pdf("cell_composition_forest_plot", width = 11, height = 6)
    forestplot(labeltext = celltypes_final, # Give list of colunms as labeltext
        txt_gp = fpTxtGp(label = gpar(fontfamily = "Noto Sans"), # Set custom font and tick size
                        ticks = gpar(cex = 0.8)),
        hrzl_lines = list("2" = gpar(lwd = 2, col = "black"), # Place horizontal lines between groups of mutations
                         "5" = gpar(lwd = 2, col = "black"),
                         "8" = gpar(lwd = 2, col = "black"),
                         "11" = gpar(lwd = 2, col = "black"),
                         "14" = gpar(lwd = 2, col = "black"),
                         "17" = gpar(lwd = 2, col = "black"),
                         "20" = gpar(lwd = 2, col = "black"),
                         "23" = gpar(lwd = 2, col = "black")),
        xticks = log(c(0.125, 0.25, 0.5, 1.0, 2.0, 4.0, signif(8.0))), # Set customs for odds ratio
        col = fpColors(lines = "black"),  # Set line colors
        ci.vertices = TRUE, ci.vertices.height = 0.2, # Set CI vertice size
        lwd.ci = 2, lwd.xaxis = 2, lwd.zero = 2, graphwidth = unit(150, "points"), # Set line widths
        lineheight = unit(50, "points"), xlog = TRUE, # Set line height and log x-axis
        graph.pos = 3, mean = celltypes_matrix_plot, # Set column number for CI plot and provide matrix as argument
        new_page = FALSE, boxsize = 0.25, zero = 1, align = c("l", "l", "l")) # disable new page, set box size, and left alignment
dev.off()

# Differential expression
# Pseudobulk counts and test for differences between CHIP mutations and control
AGS_pseudobulk <- pseudobulk_counts(AGS_combined_filter_celltype)
AGS_limma_dge_dnmt3a <- limma_bulk_dge(AGS_pseudobulk, "Mutation_1DNMT3A")
AGS_limma_dge_tet2 <- limma_bulk_dge(AGS_pseudobulk, "Mutation_1TET2")
AGS_limma_dge_ddr <- limma_bulk_dge(AGS_pseudobulk, "Mutation_1DDR")

AGS_pseudobulk_macros <- pseudobulk_counts(AGS_macros_annot, min_ncells = 3L)
AGS_limma_dge_dnmt3a_macros <- limma_bulk_dge(AGS_pseudobulk_macros, "Mutation_1DNMT3A")
AGS_limma_dge_tet2_macros <- limma_bulk_dge(AGS_pseudobulk_macros, "Mutation_1TET2")
AGS_limma_dge_ddr_macros <- limma_bulk_dge(AGS_pseudobulk_macros, "Mutation_1DDR")

wb <- openxlsx::createWorkbook()
catch <- pmap(list(1:length(AGS_limma_dge_dnmt3a), names(AGS_limma_dge_dnmt3a), AGS_limma_dge_dnmt3a), export_de_table, wb)
saveWorkbook(wb, "DGE_DNMT3A.xlsx", overwrite = TRUE)
 
wb <- openxlsx::createWorkbook()
catch <- pmap(list(1:length(AGS_limma_dge_tet2), names(AGS_limma_dge_tet2), AGS_limma_dge_tet2), export_de_table, wb)
saveWorkbook(wb, "DGE_TET2.xlsx", overwrite = TRUE)

wb <- openxlsx::createWorkbook()
catch <- pmap(list(1:length(AGS_limma_dge_ddr), names(AGS_limma_dge_ddr), AGS_limma_dge_ddr), export_de_table, wb)
saveWorkbook(wb, "DGE_DDR.xlsx", overwrite = TRUE)

wb <- openxlsx::createWorkbook()
catch <- pmap(list(1:length(AGS_limma_dge_dnmt3a_macros), names(AGS_limma_dge_dnmt3a_macros), AGS_limma_dge_dnmt3a_macros), export_de_table, wb)
saveWorkbook(wb, "DGE_DNMT3A.xlsx", overwrite = TRUE)
 
wb <- openxlsx::createWorkbook()
catch <- pmap(list(1:length(AGS_limma_dge_tet2_macros), names(AGS_limma_dge_tet2_macros), AGS_limma_dge_tet2_macros), export_de_table, wb)
saveWorkbook(wb, "DGE_TET2.xlsx", overwrite = TRUE)

wb <- openxlsx::createWorkbook()
catch <- pmap(list(1:length(AGS_limma_dge_ddr_macros), names(AGS_limma_dge_ddr_macros), AGS_limma_dge_ddr_macros), export_de_table, wb)
saveWorkbook(wb, "DGE_DDR.xlsx", overwrite = TRUE)

volcano_plot("Active NK cells", AGS_limma_dge_tet2, "subclustered_cells/active_nk_tet2")
