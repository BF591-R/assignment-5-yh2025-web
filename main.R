library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
library('fgsea')

#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(counts_csv, metafile_csv, selected_times) {
  
    counts_df <- readr::read_tsv(counts_csv)
    meta_df <- readr::read_csv(metafile_csv, show_col_types = FALSE)
    
    meta_df <- meta_df |>
      dplyr::select(samplename, timepoint) |>
      dplyr::filter(timepoint %in% selected_times)
    
    count_sample_cols <- intersect(meta_df$samplename, colnames(counts_df))
    meta_df <- meta_df |>
      dplyr::filter(samplename %in% count_sample_cols)
    
    counts_df <- counts_df |>
      dplyr::select(1, dplyr::all_of(meta_df$samplename))
    
    gene_ids <- counts_df[[1]]
    count_matrix <- counts_df[, -1] |>
      as.data.frame()
    
    count_matrix <- as.matrix(count_matrix)
    rownames(count_matrix) <- gene_ids
    
    meta_df <- meta_df[match(colnames(count_matrix), meta_df$samplename), ]
    meta_df$timepoint <- factor(meta_df$timepoint)
    meta_df$timepoint <- stats::relevel(meta_df$timepoint, ref = "vP0")
    
    coldata <- S4Vectors::DataFrame(
      samplename = meta_df$samplename,
      timepoint = meta_df$timepoint,
      row.names = meta_df$samplename
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = count_matrix),
      colData = coldata
    )
    
    return(se)
}

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
  
    dds <- DESeq2::DESeqDataSet(se, design = design)
    dds <- DESeq2::DESeq(dds)
    res <- DESeq2::results(dds)
    res_df <- as.data.frame(res)
    return(list(dds = dds, results = res_df))

}

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`. Ensure
#' that the column name for your rownames is called "genes". 
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
 
    labeled_res <- deseq2_res |>
      tibble::rownames_to_column(var = "genes") |>
      tibble::as_tibble() |>
      dplyr::mutate(
        volc_plot_status = dplyr::case_when(
          padj < padj_threshold & log2FoldChange > 0 ~ "UP",
          padj < padj_threshold & log2FoldChange < 0 ~ "DOWN",
          TRUE ~ "NS"
        )
      )
    
    return(labeled_res)
    
}

#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
  
    p <- labeled_results |>
      dplyr::filter(!is.na(pvalue)) |>
      ggplot2::ggplot(ggplot2::aes(x = pvalue)) +
      ggplot2::geom_histogram(bins = 50)
    
    return(p)

}

#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
   
    p <- labeled_results |>
      dplyr::filter(!is.na(padj), padj < padj_threshold) |>
      dplyr::filter(!is.na(log2FoldChange)) |>
      ggplot2::ggplot(ggplot2::aes(x = log2FoldChange)) +
      ggplot2::geom_histogram(bins = 50)
    
    return(p)

}

#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
  
    top_genes <- labeled_results |>
      dplyr::filter(!is.na(padj)) |>
      dplyr::arrange(padj) |>
      dplyr::slice_head(n = num_genes) |>
      dplyr::pull(genes)
    
    norm_counts <- DESeq2::counts(dds_obj, normalized = TRUE)
    norm_counts_df <- as.data.frame(norm_counts)
    norm_counts_df <- tibble::rownames_to_column(norm_counts_df, var = "genes")
    
    plot_df <- norm_counts_df |>
      dplyr::filter(genes %in% top_genes) |>
      tidyr::pivot_longer(
        cols = -genes,
        names_to = "samplename",
        values_to = "normalized_count"
      )
    
    plot_df$genes <- factor(plot_df$genes, levels = top_genes)
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = samplename, y = normalized_count)) +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(~ genes, scales = "free_y")
    
    return(p)
  
}

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
  
    p <- labeled_results |>
      dplyr::filter(!is.na(padj), !is.na(log2FoldChange)) |>
      dplyr::mutate(neg_log10_padj = -log10(padj)) |>
      ggplot2::ggplot(
        ggplot2::aes(
          x = log2FoldChange,
          y = neg_log10_padj,
          color = volc_plot_status
        )
      ) +
      ggplot2::geom_point(alpha = 0.6) +
      ggplot2::labs(
        x = "log2 Fold Change",
        y = "-log10(padj)",
        color = "Status"
      )
    
    return(p)

}

#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  
    id_map <- readr::read_tsv(id2gene_path, col_names = FALSE, show_col_types = FALSE)
    colnames(id_map) <- c("genes", "symbol")
    
    merged_df <- labeled_results |>
      dplyr::left_join(id_map, by = "genes") |>
      dplyr::filter(!is.na(symbol), !is.na(log2FoldChange)) |>
      dplyr::distinct(symbol, .keep_all = TRUE) |>
      dplyr::arrange(dplyr::desc(log2FoldChange))
    
    ranked_vector <- merged_df$log2FoldChange
    names(ranked_vector) <- merged_df$symbol
    
    return(ranked_vector)
  }

#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  
  cat("=== DEBUG START run_fgsea ===\n")
  
  pathways <- fgsea::gmtPathways(gmt_file_path)
  
  cat("original pathway sizes:\n")
  print(summary(lengths(pathways)))
  
  pathways <- lapply(pathways, function(x) intersect(x, names(rnk_list)))
  
  cat("after intersect sizes:\n")
  print(summary(lengths(pathways)))
  
  pathways <- pathways[lengths(pathways) >= min_size & lengths(pathways) <= max_size]
  
  cat("after filtering pathways:\n")
  print(summary(lengths(pathways)))
  
  fgsea_res <- fgsea::fgsea(
    pathways = pathways,
    stats = rnk_list
  )
  
  fgsea_res <- as.data.frame(fgsea_res)
  
  cat("raw fgsea size:\n")
  print(summary(fgsea_res$size))
  cat("raw min/max:\n")
  print(min(fgsea_res$size, na.rm=TRUE))
  print(max(fgsea_res$size, na.rm=TRUE))
  
  fgsea_res <- fgsea_res[!is.na(fgsea_res$size), , drop = FALSE]
  fgsea_res <- fgsea_res[fgsea_res$size >= min_size & fgsea_res$size <= max_size, , drop = FALSE]
  
  cat("final size:\n")
  print(summary(fgsea_res$size))
  cat("final min/max:\n")
  print(min(fgsea_res$size))
  print(max(fgsea_res$size))
  
  cat("=== DEBUG END ===\n")
  
  tibble::as_tibble(fgsea_res)
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){
   
    top_pos <- fgsea_results |>
      dplyr::arrange(desc(NES)) |>
      dplyr::slice_head(n = num_paths)
    
    top_neg <- fgsea_results |>
      dplyr::arrange(NES) |>
      dplyr::slice_head(n = num_paths)
    
    plot_df <- dplyr::bind_rows(top_pos, top_neg)
    
    plot_df$pathway <- factor(plot_df$pathway, levels = plot_df$pathway)
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = pathway, y = NES, fill = NES)) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::labs(x = "Pathway", y = "NES")
    
    return(p)
}