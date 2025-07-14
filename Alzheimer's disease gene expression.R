count_matrix <- read.delim("C:/Users/lenovo/Downloads/GSE153873_summary_count.star.txt.gz",
                           +                                                      sep = "\t", header = TRUE, row.names = 1)
> 
  > sample_names <- colnames(count_matrix)
  >  condition <- ifelse(grepl("Old$", sample_names), "Old",
                         +                                            ifelse(grepl("AD$", sample_names), "AD",
                                                                             +                                                                              ifelse(grepl("Young$", sample_names), "Young", NA)))
  > 
    > metadata <- data.frame(     sample = sample_names,
                                  +                             condition = factor(condition, levels = c("Old", "AD", "Young"))
                                  +                         )
    > 
      > rownames(metadata) <- sample_names  # Required for DESeq2
      > if (!requireNamespace("DESeq2", quietly = TRUE)) {
        +          install.packages("BiocManager")
        +          BiocManager::install("DESeq2")
        +      }
      > 
        > library(DESeq2)
      Loading required package: S4Vectors
      Loading required package: stats4
      Loading required package: BiocGenerics
      Loading required package: generics
      
      Attaching package: ‘generics’
      
      The following objects are masked from ‘package:base’:
        
        as.difftime, as.factor, as.ordered, intersect, is.element, setdiff, setequal,
      union
      
      
      Attaching package: ‘BiocGenerics’
      
      The following objects are masked from ‘package:stats’:
        
        IQR, mad, sd, var, xtabs
      
      The following objects are masked from ‘package:base’:
        
        anyDuplicated, aperm, append, as.data.frame, basename, cbind, colnames, dirname,
      do.call, duplicated, eval, evalq, Filter, Find, get, grep, grepl, is.unsorted,
      lapply, Map, mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
      Position, rank, rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
      unsplit, which.max, which.min
      
      
      Attaching package: ‘S4Vectors’
      
      The following object is masked from ‘package:utils’:
        
        findMatches
      
      The following objects are masked from ‘package:base’:
        
        expand.grid, I, unname
      
      Loading required package: IRanges
      
      Attaching package: ‘IRanges’
      
      The following object is masked from ‘package:grDevices’:
        
        windows
      
      Loading required package: GenomicRanges
      Loading required package: GenomeInfoDb
      Loading required package: SummarizedExperiment
      Loading required package: MatrixGenerics
      Loading required package: matrixStats
      
      Attaching package: ‘MatrixGenerics’
      
      The following objects are masked from ‘package:matrixStats’:
        
        colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse, colCounts,
      colCummaxs, colCummins, colCumprods, colCumsums, colDiffs, colIQRDiffs, colIQRs,
      colLogSumExps, colMadDiffs, colMads, colMaxs, colMeans2, colMedians, colMins,
      colOrderStats, colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
      colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads, colWeightedMeans,
      colWeightedMedians, colWeightedSds, colWeightedVars, rowAlls, rowAnyNAs,
      rowAnys, rowAvgsPerColSet, rowCollapse, rowCounts, rowCummaxs, rowCummins,
      rowCumprods, rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
      rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins, rowOrderStats,
      rowProds, rowQuantiles, rowRanges, rowRanks, rowSdDiffs, rowSds, rowSums2,
      rowTabulates, rowVarDiffs, rowVars, rowWeightedMads, rowWeightedMeans,
      rowWeightedMedians, rowWeightedSds, rowWeightedVars
      
      Loading required package: Biobase
      Welcome to Bioconductor
      
      Vignettes contain introductory material; view with 'browseVignettes()'. To cite
      Bioconductor, see 'citation("Biobase")', and for packages 'citation("pkgname")'.
      
      
      Attaching package: ‘Biobase’
      
      The following object is masked from ‘package:MatrixGenerics’:
        
        rowMedians
      
      The following objects are masked from ‘package:matrixStats’:
        
        anyMissing, rowMedians
      > dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                      +                                                              colData = metadata,
                                      +                                                              design = ~ condition)
      > 
        > dds <- dds[rowSums(counts(dds)) > 10, ]
      > dds <- DESeq(dds)
      estimating size factors
      estimating dispersions
      gene-wise dispersion estimates
      mean-dispersion relationship
      final dispersion estimates
      fitting model and testing
      -- replacing outliers and refitting for 202 genes
      -- DESeq argument 'minReplicatesForReplace' = 7 
      -- original counts are preserved in counts(dds)
      estimating dispersions
      fitting model and testing
      > res_AD_vs_Old <- results(dds, contrast = c("condition", "AD", "Old"))
      > res_AD_vs_Old <- results(dds, contrast = c("condition", "AD", "Old"))
      > res_AD_vs_Old <- results(dds, contrast = c("condition", "AD", "Old"))
      > res_AD_vs_Old <- results(dds, contrast = c("condition", "AD", "Old"))
      >  res_ordered <- res_AD_vs_Old[order(res_AD_vs_Old$pvalue), ]
      >  res_ordered <- na.omit(res_ordered)  # Remove rows with NA p-values
      > 
        > head(res_ordered, 5)
      log2 fold change (MLE): condition AD vs Old 
      Wald test p-value: condition AD vs Old 
      DataFrame with 5 rows and 6 columns
      baseMean log2FoldChange     lfcSE      stat      pvalue        padj
      <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
        DUSP6         312.9708      -0.908530  0.152604  -5.95352 2.62434e-09 3.00648e-05
      CRH            15.1804      -2.600392  0.442563  -5.87576 4.20912e-09 3.00648e-05
      LOC101926975    9.0819      -2.092138  0.356166  -5.87405 4.25265e-09 3.00648e-05
      KLF15         834.8493       1.181748  0.213764   5.52828 3.23391e-08 1.41878e-04
      C2orf16        75.9637      -0.896143  0.162275  -5.52236 3.34477e-08 1.41878e-04
      >   ggplot(crh_df, aes(x = condition, y = expression, fill = condition)) +
        +               geom_boxplot(outlier.shape = NA, alpha = 0.6) +
        +               geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
        +         labs(
          +                           title = "CRH Expression across Conditions",
          +                           x = "Condition",
          +                           y = "Normalized Expression (VST)"
          +                      ) +
        +               theme_minimal() +
        +         scale_fill_manual(values = c("Old" = "#FDB863", "AD" = "#B2ABD2", "Young" = "#5E3C99"))
      Error in ggplot(crh_df, aes(x = condition, y = expression, fill = condition)) : 
        could not find function "ggplot"
      
      > library(ggplot2)
      > library(reshape2)
      > 
        > top_genes <- c("DUSP6", "CRH", "LOC101926975", "KLF15", "C2orf16")
      > expr_data <- assay(vsd)[top_genes, ]  # rows = genes, cols = samples
      Error in h(simpleError(msg, call)) : 
        error in evaluating the argument 'x' in selecting a method for function 'assay': object 'vsd' not found
      
      > aes(x = condition, y = expression, fill = condition)
      Aesthetic mapping: 
        * `x`    -> `condition`
      * `y`    -> `expression`
      * `fill` -> `condition`
      > expr_df <- as.data.frame(t(expr_data))  # transpose: rows = samples, cols = genes
      Error in h(simpleError(msg, call)) : 
        error in evaluating the argument 'x' in selecting a method for function 'as.data.frame': error in evaluating the argument 'x' in selecting a method for function 't': object 'expr_data' not found
      
      >  top_genes <- c("DUSP6", "CRH", "LOC101926975", "KLF15", "C2orf16")
      > expr_data <- assay(vsd)[top_genes, ]  # rows = genes, columns = samples
      Error in h(simpleError(msg, call)) : 
        error in evaluating the argument 'x' in selecting a method for function 'assay': object 'vsd' not found
      
      > vsd <- vst(dds)  # dds must be your DESeq2 object (already created with DESeq())
      >  library(DESeq2)
      > vsd <- vst(dds)  # 'dds' should be your DESeq2 object
      >  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = metadata, design = ~ condition)
      >     dds <- DESeq(dds)
      estimating size factors
      estimating dispersions
      gene-wise dispersion estimates
      mean-dispersion relationship
      final dispersion estimates
      fitting model and testing
      -- replacing outliers and refitting for 202 genes
      -- DESeq argument 'minReplicatesForReplace' = 7 
      -- original counts are preserved in counts(dds)
      estimating dispersions
      fitting model and testing
      >     
        > library(DESeq2)
      > dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                      +                                                              colData = metadata,
                                      +                                                             design = ~ condition)
      >  dds <- DESeq(dds)
      estimating size factors
      estimating dispersions
      gene-wise dispersion estimates
      mean-dispersion relationship
      final dispersion estimates
      fitting model and testing
      -- replacing outliers and refitting for 202 genes
      -- DESeq argument 'minReplicatesForReplace' = 7 
      -- original counts are preserved in counts(dds)
      estimating dispersions
      fitting model and testing
      > 
        > 
        > top_genes <- c("DUSP6", "CRH", "LOC101926975", "KLF15", "C2orf16")
      > expr_data <- assay(vsd)[top_genes, ]  # rows = genes, columns = samples
      > expr_df <- as.data.frame(t(expr_data))
      >  expr_df$condition <- colData(vsd)$condition
      > expr_long <- melt(expr_df, id.vars = "condition",
                          +                                        variable.name = "gene",
                          +                                   value.name = "expression")
      > 
        > ggplot(expr_long, aes(x = condition, y = expression, fill = condition)) +
        +          geom_boxplot(outlier.shape = NA, alpha = 0.6) +
        +          geom_jitter(width = 0.2, size = 1.8, alpha = 0.8) +
        +          facet_wrap(~ gene, scales = "free_y", ncol = 2) +
        +       labs(
          +                 title = "Top 5 Differentially Expressed Genes (AD vs Old)",
          +                 x = "Condition",
          +                 y = "Normalized Expression (VST)"
          +              ) +
        +          theme_minimal(base_size = 13) +
        +          scale_fill_manual(values = c("Old" = "#FDB863", "AD" = "#B2ABD2", "Young" = "#5E3C99")) 
      >     theme(
        +                  strip.text = element_text(face = "bold", size = 13),
        +                 axis.text.x = element_text(angle = 30, hjust = 1)
        +            )
      List of 2
      $ axis.text.x:List of 11
      ..$ family       : NULL
      ..$ face         : NULL
      ..$ colour       : NULL
      ..$ size         : NULL
      ..$ hjust        : num 1
      ..$ vjust        : NULL
      ..$ angle        : num 30
      ..$ lineheight   : NULL
      ..$ margin       : NULL
      ..$ debug        : NULL
      ..$ inherit.blank: logi FALSE
      ..- attr(*, "class")= chr [1:2] "element_text" "element"
      $ strip.text :List of 11
      ..$ family       : NULL
      ..$ face         : chr "bold"
      ..$ colour       : NULL
      ..$ size         : num 13
      ..$ hjust        : NULL
      ..$ vjust        : NULL
      ..$ angle        : NULL
      ..$ lineheight   : NULL
      ..$ margin       : NULL
      ..$ debug        : NULL
      ..$ inherit.blank: logi FALSE
      ..- attr(*, "class")= chr [1:2] "element_text" "element"
      - attr(*, "class")= chr [1:2] "theme" "gg"
      - attr(*, "complete")= logi FALSE
      - attr(*, "validate")= logi TRUE
      > 
        > library(ggpubr)
      > 
        > # Define which groups to compare
        > my_comparisons <- list(c("AD", "Old"))
      > 
        > # Create the plot with significance markers
        > ggboxplot(expr_long, x = "condition", y = "expression", fill = "condition",
                    +           palette = c("Old" = "#FDB863", "AD" = "#B2ABD2", "Young" = "#5E3C99"),
                    +           facet.by = "gene", short.panel.labs = TRUE,
                    +           add = "jitter", width = 0.6) +
        +     stat_compare_means(comparisons = my_comparisons,
                                 +                        method = "wilcox.test",  # you can change to "t.test" if appropriate
                                 +                        label = "p.signif",       # shows stars (*, **, ***)
                                 +                        hide.ns = TRUE) +
        +     labs(
          +         title = "CRH and KLF15 Expression with Significance (AD vs Old)",
          +         x = "Condition",
          +         y = "Normalized Expression (VST)"
          +     ) +
        +     theme_minimal(base_size = 13)
      Warning message:
        No shared levels found between `names(values)` of the manual scale and the data's colour
values. 
> 
> # Define genes of interest
> genes_of_interest <- c("CRH", "KLF15")
> 
> # Extract expression data
> expr_data <- assay(vsd)[genes_of_interest, ]
> 
> # Transpose for plotting
> expr_df <- as.data.frame(t(expr_data))
> expr_df$condition <- colData(vsd)$condition  # add condition column
> 
> # Melt to long format
> expr_long <- melt(expr_df, id.vars = "condition",
+                   variable.name = "gene",
+                   value.name = "expression")
> 
> # Define which conditions to compare
> my_comparisons <- list(c("AD", "Old"))
> 
> # Generate the plot
> ggboxplot(expr_long, 
+           x = "condition", 
+           y = "expression", 
+           fill = "condition",
+           palette = c("Old" = "#FDB863", "AD" = "#B2ABD2", "Young" = "#5E3C99"),
+           facet.by = "gene", 
+           short.panel.labs = TRUE,
+           add = "jitter", 
+           width = 0.6) +
+     stat_compare_means(comparisons = my_comparisons,
+                        method = "wilcox.test",
+                        label = "p.signif",
+                        hide.ns = FALSE) +
+     labs(
+         title = "Expression of CRH and KLF15 with Statistical Significance",
+         x = "Condition",
+         y = "Normalized Expression (VST)"
+     ) +
+     theme_minimal(base_size = 13)
Warning message:
No shared levels found between `names(values)` of the manual scale and the data's colour
      values. 
      > 
        > # Plot using ggpubr
        > ggboxplot(expr_long, x = "condition", y = "expression", fill = "condition",
                    +           palette = c("Old" = "#FDB863", "AD" = "#B2ABD2", "Young" = "#5E3C99"),
                    +           facet.by = "gene", short.panel.labs = TRUE, width = 0.6, add = "jitter") +
        +     stat_compare_means(comparisons = list(c("AD", "Old")),
                                 +                        method = "wilcox.test", 
                                 +                        label = "p.signif", 
                                 +                        hide.ns = TRUE) +  # show stars (*, **, ***)
        +     labs(
          +         title = "Expression of CRH and KLF15 with Significance Markers",
          +         x = "Condition",
          +         y = "Normalized Expression (VST)"
          +     ) +
        +     theme_minimal(base_size = 13)
      Warning message:
        No shared levels found between `names(values)` of the manual scale and the data's colour
values. 
> 
> # Load libraries
> library(ggplot2)
> library(reshape2)
> 
> # Define your two genes of interest
> genes_of_interest <- c("CRH", "KLF15")
> 
> # Extract expression data from VST-transformed object
> expr_data <- assay(vsd)[genes_of_interest, ]
> 
> # Transpose and convert to data frame
> expr_df <- as.data.frame(t(expr_data))  # rows = samples, columns = genes
> 
> # Add condition info
> expr_df$condition <- colData(vsd)$condition
> 
> # Convert to long format
> expr_long <- melt(expr_df, id.vars = "condition", 
+                   variable.name = "gene", 
+                   value.name = "expression")
> 
> # Plot boxplots
> ggplot(expr_long, aes(x = condition, y = expression, fill = condition)) +
+     geom_boxplot(outlier.shape = NA, alpha = 0.6) +
+     geom_jitter(width = 0.2, size = 1.8, alpha = 0.8) +
+     facet_wrap(~ gene, scales = "free_y", ncol = 2) +
+     labs(
+         title = "Expression of CRH and KLF15 across AD, Old, and Young Samples",
+         x = "Condition",
+         y = "Normalized Expression (VST)"
+     ) +
+     theme_minimal(base_size = 13) +
+     scale_fill_manual(values = c("Old" = "#FDB863", "AD" = "#B2ABD2", "Young" = "#5E3C99")) +
+     theme(
+         strip.text = element_text(face = "bold", size = 13),
+         axis.text.x = element_text(angle = 30, hjust = 1)
+     )
> 
> 