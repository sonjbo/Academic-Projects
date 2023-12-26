library(shiny)
library(colourpicker)
library(bslib)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(stringr)
library(tools)
library(shinyWidgets)
library(purrr)
library(gplots)
library(RColorBrewer)


read_meta <- function(file) {
  meta <- read.csv(file = file, header = TRUE)
  meta_t <- as.data.frame(t(meta))
  colnames(meta_t) <- meta_t[1, ]
  meta_t <- meta_t[-1, , drop = FALSE]
  
  cols_to_drop <- c('Sample_status', 'Sample_submission_date', 'Sample_last_update_date')
  meta_t <- meta_t[, !colnames(meta_t) %in% cols_to_drop, drop = FALSE]
  
  numeric_cols <- c('Sample_channel_count', 'pmi', 'age_of_death', 'rin', 'mrna_seq_reads')
  meta_t[, numeric_cols] <- lapply(meta_t[, numeric_cols], as.numeric)
  #meta_t %>% rename(gene = Sample_geo_accession)
  return(meta_t)
}

sum_table_gen <- function(meta) {
  summary_df <- data.frame(stringsAsFactors = FALSE)
  col_types <- sapply(meta, class)
  
  for (i in seq_along(col_types)) {
    col_name <- names(col_types)[i]
    col_type <- ifelse(col_types[i] == 'character', 'Factor', 'Double')
    
    # Extract mean and standard deviation for numeric columns
    if (is.numeric(meta[[col_name]])) {
      col_mean <- mean(meta[[col_name]], na.rm = TRUE)
      col_sd <- sd(meta[[col_name]], na.rm = TRUE)
      col_values <- paste(round(col_mean, 2), " (+/- ", round(col_sd, 2), ")", sep = "")
    } else {
      col_values <- toString(unique(meta[[col_name]]))
    }
    
    col_name_cleaned <- toTitleCase(str_replace_all(col_name, "_", " "))
    col_values_cleaned <- toTitleCase(col_values)
    
    summary_row <- data.frame(
      'Column Name' = col_name_cleaned,
      'Type' = col_type,
      'Mean (sd) or Distinct Values' = col_values_cleaned,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    summary_df <- rbind(summary_df, summary_row)
  }
  return(summary_df)
}


#density plot
dens_plot <- function(metadata, x, fill){
  ggplot(data = metadata, aes_string(x = x, fill = fill)) +
    geom_density(alpha = 0.4)+
    labs(x = toTitleCase(str_replace_all(x, "_", " ")),
         y = toTitleCase(str_replace_all(fill, "_", " ")),
         fill = str_to_title(str_replace_all(fill, "_", " ")))
}

#histogram plot
hist_plot <- function(metadata, x, fill){
  ggplot(data = metadata, aes_string(x = x, fill = fill)) +
    geom_histogram(alpha = 0.4, color = "black")+
    labs(x = toTitleCase(str_replace_all(x, "_", " ")),
         y = toTitleCase(str_replace_all(fill, "_", " ")),
         fill = str_to_title(str_replace_all(fill, "_", " ")))
}

#violin plot 
violin_plot <- function(metadata, x, y){
  ggplot(data = metadata, aes_string(x = y, y = x, fill = y)) +
    geom_violin(alpha = 0.4, color = "black")+
    labs(x = toTitleCase(str_replace_all(y, "_", " ")),
         y = toTitleCase(str_replace_all(x, "_", " ")),
         fill = str_to_title(str_replace_all(y, "_", " ")))
}


#working with counts data

read_counts <- function(file){
  counts <- read.csv(file = file, header = TRUE)
  return (counts)
}


summary_counts <- function(counts, var_threshold, non_zero) {
  #filter genes based on variance threshold
  genes_var_filtered <- counts[apply(counts, 1, function(row) var(row, na.rm = TRUE) >= var_threshold), ]
  
  #filter genes based on non-zero samples
  genes_zero_filtered <- counts[rowSums(counts != 0) >= non_zero, ]
  
  #total genes filtered
  #genes_filtered <- counts[apply(counts, 1, function(row) var(row, na.rm = TRUE) >= var_threshold) &
  #rowSums(counts != 0) >= non_zero, ]
  
  genes_filtered <- intersect(genes_var_filtered, genes_zero_filtered)
  
  # Calculate summary values
  total_genes <- nrow(counts)
  total_samples <- ncol(counts[-1])
  genes_after_var_filter <- nrow(genes_var_filtered)
  genes_after_all_filters <- nrow(genes_filtered)
  genes_after_zero_filters <- nrow(genes_zero_filtered)
  
  # Create summary dataframe
  summary_df <- data.frame(
    "Total genes" = total_genes,
    "Total samples" = total_samples,
    "Number of genes passing variance filter"  = genes_after_var_filter,
    "% of genes passing variance filter" = genes_after_var_filter / total_genes *100,
    "Number of genes passing sample filter" = genes_after_zero_filters,
    "% of genes passing sample filter" = genes_after_zero_filters / total_genes *100,
    check.names = FALSE
  )
  
  return(summary_df)
}

plot_pca <- function(data, meta, PC1, PC2) {
  rownames(data) <- data$gene
  data <- data[-1, drop = FALSE]
  pca_results <- prcomp(scale(t(data[])))
  variance <- pca_results$sdev^2
  var_explained <- variance / sum(variance) * 100
  #pca_results <- prcomp(scale(t(data[], center = TRUE, scale = TRUE)))
  pca_sub <- data.frame(sample = colnames(data), PC1 = pca_results$x[, PC1], PC2 = pca_results$x[, PC2])
  #meta$sample <- pca_sub$sample
  merged <- cbind(pca_sub, meta)
  scatter_pca <- ggplot(merged, aes(x = PC1, y = PC2, color = diagnosis)) +
    geom_point() +
    labs(x = paste0("PC", PC1, "-", round(var_explained[PC1], 0), "% Variance"),
         y = paste0("PC", PC2, "-", round(var_explained[PC2], 0), "% Variance"))
  
  return(scatter_pca)
}

plot_variance_vs_mean <- function(data, var_threshold){
  
  genes_var_filtered <- data[apply(data, 1, function(row) var(row, na.rm = TRUE) >= var_threshold), ]
  
  gene_filt <- genes_var_filtered %>%
    pivot_longer(cols = -1, names_to = "sample", values_to = 'values')%>%
    group_by(gene)%>%
    summarize(mean_count = mean(values, na.rm = TRUE), variance = var(values))%>%
    arrange(mean_count)
  
  gene_filt <- gene_filt %>%
    mutate(mean_rank = rank(mean_count, ties.method = "first"))
  
  gene_stats <- data %>%
    pivot_longer(cols = -1, names_to = "sample", values_to = 'values')%>%
    group_by(gene)%>%
    summarize(mean_count = mean(values, na.rm = TRUE), variance = var(values))%>%
    arrange(mean_count)
  
  gene_stats <- gene_stats %>%
    mutate(mean_rank = rank(mean_count, ties.method = "first"))
  
  last_gene_index <- tail(gene_filt$gene, 1)
  
  last_gene_rank <- gene_stats %>%
    filter(gene == last_gene_index)
  
  end_rank <- as.integer(last_gene_rank$mean_rank)
  rows <- as.integer(nrow(gene_filt))
  start_rank <- as.integer(end_rank - rows)
  
  gene_filt <- gene_filt %>%
    mutate(mean_rank = rank(mean_count, ties.method = "first") + start_rank - 1)
  
  gene_stats$filtered <- "FALSE"
  gene_filt$filtered <- "TRUE"
  combined_data <- bind_rows(gene_stats, gene_filt)
  
  p <- ggplot(combined_data, aes(x = mean_rank, y = variance, color = filtered)) +
    geom_point() +
    xlab("Rank(Mean)") +
    ylab("Variance") +
    scale_color_manual(values = c("FALSE" = "#979AED", "TRUE" = "#170549")) +
    scale_y_continuous(trans = "log10")
  #}
  return(p)
}

plot_zero_vs_mean <- function(data, non_zero){
  
  genes_var_filtered <- data[rowSums(data != 0) >= non_zero, ]
  
  gene_filt <- genes_var_filtered %>%
    pivot_longer(cols = -1, names_to = "sample", values_to = 'values')%>%
    group_by(gene) %>%
    summarize(mean_count = mean(values, na.rm = TRUE),
              num_zeros = sum(values == 0)) %>%
    arrange(mean_count)
  
  gene_filt <- gene_filt %>%
    mutate(mean_rank = rank(mean_count, ties.method = "first"))
  
  gene_stats <- data %>%
    pivot_longer(cols = -1, names_to = "sample", values_to = 'values')%>%
    group_by(gene) %>%
    summarize(mean_count = mean(values, na.rm = TRUE),
              num_zeros = sum(values == 0)) %>%
    arrange(mean_count)
  
  gene_stats <- gene_stats %>%
    mutate(mean_rank = rank(mean_count, ties.method = "first"))
  
  last_gene_index <- tail(gene_filt$gene, 1)
  
  last_gene_rank <- gene_stats %>%
    filter(gene == last_gene_index)
  
  end_rank <- as.integer(last_gene_rank$mean_rank)
  rows <- as.integer(nrow(gene_filt))
  start_rank <- as.integer(end_rank - rows)
  
  gene_filt <- gene_filt %>%
    mutate(mean_rank = rank(mean_count, ties.method = "first") + start_rank - 1)
  
  gene_stats$filtered <- "FALSE"
  gene_filt$filtered <- "TRUE"
  combined_data <- bind_rows(gene_stats, gene_filt)
  
  p <- ggplot(combined_data, aes(x = mean_rank, y = num_zeros, color = filtered)) +
    geom_point() +
    xlab("Rank(Mean)") +
    ylab("Number of Zeros") +
    scale_color_manual(values = c("FALSE" = "#979AED", "TRUE" = "#170549")) 
  #scale_y_continuous(trans = "log10")
  #}
  return(p)
}

heatmapplot <- function(data, var_threshold, non_zero) {
  suppressWarnings({
    #filter genes based on variance threshold
    genes_var_filtered <- data[apply(data, 1, function(row) var(row, na.rm = TRUE) >= var_threshold), ]
    #filter genes based on non-zero samples
    genes_filtered <- genes_var_filtered[rowSums(genes_var_filtered != 0) >= non_zero, ]
    
    marker_matrix <- as.matrix(genes_filtered[-1])
    rownames(marker_matrix) <- genes_filtered$gene
    log_marker_matrix <- log2(marker_matrix+1)  
    
    p <- heatmap.2(log_marker_matrix,
                   scale = 'column',
                   col="bluered",
                   key = TRUE,    
                   key.title = 'Expression Level',
                   trace="none")  
    
    return(p)
  })
}

read_deseq <- function(deseq){
  results <- read.csv(deseq)
  return (results)
}

label_res <- function(deseq, padj_threshold) {
  labeled_deseq <- data_frame(deseq)
  deseq$volc_plot_status <- NA 
  deseq$volc_plot_status[deseq$padj >= padj_threshold] <- 'NS'
  deseq$volc_plot_status[deseq$padj < padj_threshold & deseq$log2FoldChange < 0] <- 'DOWN'
  deseq$volc_plot_status[deseq$padj < padj_threshold & deseq$log2FoldChange > 0] <- 'UP'
  return(deseq)
}

plot_vals <- function(labeled_res, X) {
  labeled_res[[X]] <- as.numeric(labeled_res[[X]])
  
  p <- ggplot(labeled_res, aes(x = !!as.name(X))) +
    geom_histogram(binwidth = 0.02, fill = "olivedrab3", color = 'black') +
    labs(title = "Histogram of selected values obtained from DE analysis",
         x = X,
         y = "Count"
    )
  return (p)
}


plot_log2fc <- function(labeled_results, padj_threshold){
  filtered_results <- labeled_results %>% filter(padj < padj_threshold)
  p <- ggplot(filtered_results, aes(x = log2FoldChange)) +
    geom_histogram(binwidth = 0.1, fill = "royalblue1", color = 'black') +
    labs(title = "Histogram of Log2FoldChange for DE genes",
         x = "Log2FoldChange Values",
         y = "Count"
    )
  return (p)
}

scatter_norm_counts <- function(labeled_results, counts, num_genes, num_samples){
  top_genes <- head(labeled_results[order(labeled_results$padj), ], num_genes)
  gene_names <- top_genes$gene
  norm_counts <- subset(counts, counts$gene %in% gene_names)
  row.names(norm_counts) <- norm_counts$gene
  sample_subset <- sample(colnames(norm_counts), num_samples)
  
  norm_counts_subset <- norm_counts[, sample_subset, drop = FALSE]
  
  df_long <- norm_counts_subset %>%
    rownames_to_column(var = 'gene') %>%
    pivot_longer(cols = -gene, names_to = 'Samples', values_to = 'count') %>%
    mutate(gene = as.character(gene)) 
  
  p <- ggplot(df_long, aes(x=gene, y = log10(count), color = Samples)) +
    geom_point() +
    labs(title = "Plot of log10(normalized counts) for top ten DE Genes",
         x = "Genes",
         y = "log10(norm_counts")+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90)) +
    coord_cartesian(ylim = c(0, 1.0))
  
  return(p)
}


plot_volcano <- function(labeled_results, color1, color2, padj_threshold){
  deseq <- data_frame(labeled_results)
  deseq$volc_plot_status <- NA 
  deseq$volc_plot_status[deseq$padj >= padj_threshold] <- 'NS'
  deseq$volc_plot_status[deseq$padj < padj_threshold & deseq$log2FoldChange < 0] <- 'DOWN'
  deseq$volc_plot_status[deseq$padj < padj_threshold & deseq$log2FoldChange > 0] <- 'UP'
  
  label <- ifelse(padj_threshold < -150, 'padj < 1 * 10^-150', paste0('padj < 1 * 10^', padj_threshold))
  
  labeled_results <- deseq%>%
    drop_na(volc_plot_status)
  p <- ggplot(labeled_results, aes(x = log2FoldChange, y = -log10(padj))) + 
    geom_point(aes(color = ifelse(padj < 10^padj_threshold, "Significant", "Not Significant")), size = 1) +
    scale_colour_manual(name = label, values = c("Significant" = color1, "Not Significant" = color2, "grey" = 'grey')) +
    scale_x_continuous(expand = c(0, 0.05)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  return (p)
}

make_ranked_log2fc <- function(labeled_results) {
  labeled_results <- labeled_results %>%
    drop_na()
  rnk_vec <- deframe(labeled_results[c("symbol", "log2FoldChange")])
  rnk_vec <- rnk_vec[order(-rnk_vec)]
  return(rnk_vec)
}

run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  fgsea_file <- fgsea::gmtPathways(gmt_file_path)
  fgsea_res <- fgsea(fgsea_file, rnk_list, minSize = min_size, maxSize=max_size)
  fgsea_res <- as_tibble(fgsea_res)
  return(fgsea_res)
}

read_fgsea <- function(fgsea){
  fgsea <- read.csv(fgsea)
  return (fgsea)
}

top_pathways <- function(fgsea_results, padj_filter, num_paths){
  
  pos <- filter(fgsea_results, padj<10^padj_filter)%>%
    arrange(padj) %>%
    slice_head(n=num_paths)
  #neg <- filter(fgsea_results, <0)%>%
  #  arrange((NES))%>%
  #slice_head(n=num_paths)
  
  #fgsea_top <- bind_rows(pos, neg)
  pos$log_padj <- log10(pos$padj)
  
  p <- ggplot(pos, aes(x = fct_reorder(pathway, -log_padj), y = log_padj, fill = padj > padj_filter)) +
    geom_bar(stat = 'identity') +
    scale_fill_manual(values = c('TRUE' = 'red', 'FALSE' = 'blue')) +
    theme(axis.text = element_text(size=4)) +
    theme(axis.title = element_text(size=6)) +
    theme(plot.title = element_text(size = 10)) +
    ggtitle('Fgsea Results for Hallmark MSigDB Genes') +
    ylab('Normalized Enrichment Score (NES)') +
    xlab('') +
    guides(fill = 'none') +
    coord_flip() 
  return(p)
}

filtered_fgsea <- function(fgsea, padj_filter, NES_filter) {
  fgsea <- filter(fgsea, padj > 10^padj_filter) %>%
    arrange(padj)
  
  if (NES_filter == "positive") {
    fgsea <- filter(fgsea, NES > 0)
  } else if (NES_filter == "negative") {
    fgsea <- filter(fgsea, NES < 0)
  } else if (NES_filter == "all") {
    fgsea <- fgsea
  }
  
  return(fgsea)
}

padjNES <- function(fgsea, padj_threshold) {
  # Filter data based on padj threshold
  filtered_fgsea <- fgsea %>%
    filter(padj <= 10^padj_threshold)
  
  # Create scatter plot
  plot <- ggplot(fgsea, aes(x = NES, y = -log10(padj))) +
    geom_point(aes(color = ifelse(padj <= 10^(padj_threshold), "TRUE", "FALSE")), size = 1) +
    scale_color_manual(values = c("TRUE" = "#CF0071", "FALSE" = "grey"), 
                       name = paste("padj <", 10^padj_threshold),
                       labels = c("TRUE", "FALSE")) +
    labs(title = "Scatter Plot of -log10(padj) vs NES",
         x = "NES",
         y = "-log10(padj)")
  
  return(plot)
}


ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "yeti"),
  titlePanel("Post-mortem Huntington’s disease prefrontal cortex compared with neurologically healthy controls"),
  tabsetPanel(
    tabPanel('Sample',
             sidebarLayout(
               sidebarPanel(
                 fileInput("file", "Sample file", accept = c('.csv', '.tsv')),
                 actionButton("loadData", "Load Data", class = "btn-info")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel('Sample Summary', tableOutput('summary')),
                   tabPanel('Sample Table', dataTableOutput('table')),
                   tabPanel('Distribution Plots',
                            selectInput("plotType", "Select Plot Type", choices = c("Density Plot", "Histogram", "Violin Plot"), selected = "Density Plot"),
                            radioButtons("button1", "Select column to plot",
                                         choices = c("tissue", "age_of_death", "diagnosis", "pmi", "rin", "mrna_seq_reads"), selected = "age_of_death"),
                            radioButtons("button2", "Select column to group by",
                                         choices = c("tissue", "age_of_death", "diagnosis", "pmi", "rin", "mrna_seq_reads"), selected = "diagnosis"),
                            plotOutput("selected_plot")
                   )
                 )
               )
             )
    ),
    tabPanel('Counts',
             sidebarLayout(
               sidebarPanel(
                 fileInput("file2", "Counts file"),
                 actionButton("loadCounts", "Load Data", class = "btn-info"),
                 sliderInput(inputId = 'variance', min = 0, 
                             max = 100,
                             label = 'Select the variance threshold', value = 80, step = 1),
                 sliderInput(inputId = 'non_zero', min = 0, max = 69,
                             label = 'Select the number of non-zero samples', value = 50, step = 1)
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel('Filtered Counts Table', 
                            tableOutput('sum_counts')),
                   tabPanel('Diagnostic Scatter Plots',
                            div(
                              plotOutput('varmeans', width = "800px", height = "300px"),
                              plotOutput('zeromeans', width = "800px", height = "300px"))),
                   tabPanel('PCA', 
                            radioButtons("pc1", "Choose first principal component",
                                         choices = 1:20, selected = 1,  inline = TRUE),
                            radioButtons("pc2", "Choose second principal component",
                                         choices = 1:20, selected = 2,  inline = TRUE),
                            plotOutput('pca')),
                   tabPanel('Heatmap',
                            sliderInput(inputId = 'variance2', min = 0, 
                                        max = 100,
                                        label = 'Select the variance threshold', value = 80, step = 1),
                            sliderInput(inputId = 'non_zero2', min = 0, max = 69,
                                        label = 'Select the number of non-zero samples', value = 50, step = 1),
                            plotOutput('heatmap', width = "600px", height = "600px"))
                 )
               )
             )
    ),
    tabPanel('Differential Expression',
             sidebarLayout(
               sidebarPanel(
                 fileInput("file3", "DE file"),
                 actionButton('loadDE', 'Load Data', class = "btn-info")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel('DE Results Table', dataTableOutput('deseq')),
                   tabPanel('Histograms',
                            radioButtons("button3", "Choose column for the x-axis",
                                         choices = c('log2FoldChange','lfcSE', 'stat', 'pvalue', 'padj'), 
                                         selected = 'pvalue'),
                            plotOutput('vals')),
                   tabPanel('Log2FoldChange', 
                             sliderInput(inputId = 'padj', min = 0, max = 1, 
                                         label = "Select the magnitude of the p adjusted threshold:", value = 0.1, step = 0.01),
                             plotOutput('log2fc')),
                   tabPanel('Scatter Plot',
                            sliderInput(inputId = 'genes', min = 0, max = 100,
                                        label = 'Select the number of top genes to plot:', value = 10, step = 1),
                            sliderInput(inputId = 'samples', min = 0, max = 69,
                                        label = 'Select the number of random samples to plot:', value = 10, step = 1),
                            plotOutput('scatter')),
                            #add option for excluding samples with 0 counts?
                   tabPanel('Volcano Plot',
                            colourInput("base_color", "Base point color", "#C05746"),
                            colourInput("highlight_color", "Highlight point color", "#ADC698"),
                            sliderInput(inputId = 'padj_filter', min = -30, max = 0, 
                                        label = "Select the magnitude of the p adjusted coloring:", value = -10, step = 1),
                            plotOutput('volcano')
                            )
                 )
               )
             )
    ),
    tabPanel('Gene Set Enrichment Analysis',
             sidebarLayout(
               sidebarPanel(
                 fileInput("file4", "GSEA file"),
                 actionButton("loadGSEA", 'Load Fgsea Results', class = "btn-info")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel('Top Pathways',
                            sliderInput(inputId = 'pathways', min = -50, max = 0, 
                                        label = "Select the p-adjusted value to plot top results by", value = -10, step = 1),
                            sliderInput(inputId = 'numpaths', min = 1, max = 100,
                                        label = 'Select the maximum number of pathways to plot', value = 10, step = 1),
                            plotOutput('pathways')),
                   tabPanel('Pathway Table',
                            sliderInput(inputId = 'padjtable', min = -50, max = 0,
                                        label = "Select the p-adjusted value to filter results by", value = -10, step = 1),
                            radioButtons('NESbutton', 'Choose NES pathways to display',
                                         choices = c('all', 'positive', 'negative'),
                                         selected = 'all'),
                            dataTableOutput('toppadj'),
                            downloadButton('downloadTopPadj', 'Download Results', class = 'btn-success')),
                   tabPanel('Scatter Plot',
                            sliderInput(inputId = 'padjscatter', min = -50, max = 0,
                                        label = 'Select p-adjusted value for scatter plot', value = -10, step = 1),
                            plotOutput('padjvsNES'))
                 )
                 )
               )
             )
     )
  )

#Define server logic
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 30 * 1024^2)
  
  sample_data <- reactiveVal(NULL)
  counts_data <- reactiveVal(NULL)
  deseq_data <- reactiveVal(NULL)
  de_choices <- reactiveVal(NULL)
  fgsea_data <-reactiveVal(NULL)

  observeEvent(input$loadData, {
    sample_data(read_meta(input$file$datapath))
  })
  
  observeEvent(input$loadCounts, {
    counts_data(read_counts(input$file2$datapath))
  })
  
  observeEvent(input$loadDE, {
    deseq_data(read_deseq(input$file3$datapath))
    # Update choices for radioButtons
    de_choices(names(deseq_data()))
  })
  
  observeEvent(input$loadGSEA, {
    # Assuming you have the necessary inputs for these functions
    fgsea_data(read_fgsea(input$file4$datapath))  # Adjust accordingly
  })
  
  output$summary <- renderTable({
    sum_table_gen(sample_data())
  })
  
  output$table <- renderDataTable({
    sample_data()
  })
  
  output$selected_plot <- renderPlot({
    if (!is.null(sample_data())) {
      if (input$plotType == "Density Plot") {
        dens_plot(sample_data(), input$button1, input$button2)
      } else if (input$plotType == "Histogram") {
        hist_plot(sample_data(), input$button1, input$button2)
      } else if (input$plotType == "Violin Plot") {
        violin_plot(sample_data(), input$button1, input$button2)
      }
    }
  })
  
  output$sum_counts <- renderTable({
    summary_counts(counts_data(), input$variance, input$non_zero)
  })
  
  output$varmeans <- renderPlot({
    plot_variance_vs_mean(counts_data(), input$variance)
  })
  
  output$zeromeans <- renderPlot({
    plot_zero_vs_mean(counts_data(), input$non_zero)
  })
  
  output$pca <- renderPlot({

    pc1 <- as.numeric(input$pc1)
    pc2 <- as.numeric(input$pc2)

    plot_pca(counts_data(), sample_data(), pc1, pc2)
  })
  
  output$heatmap <- renderPlot({
    heatmapplot(counts_data(), input$variance2, input$non_zero2)
  })
  
  output$deseq <- renderDataTable({
    deseq_data()
  })
  
  #observe({
    # Set choices for radioButtons in 'DE' tab panel
  #  updateRadioButtons(session, "button3", choices = de_choices(), selected = 'pvalue')
  #})
  
  output$vals <- renderPlot({
    plot_vals(deseq_data(), input$button3)
  })
  
  output$log2fc <- renderPlot({
    plot_log2fc(deseq_data(), input$padj)
  })
  
  output$scatter <- renderPlot({
    # Check if both deseq_data and counts_data are available
    scatter_norm_counts(deseq_data(), counts_data(), input$genes, input$samples)
  })
  
  
  output$volcano <- renderPlot({
    plot_volcano(deseq_data(), input$base_color, input$highlight_color, input$padj_filter)
  })
  
  output$pathways <- renderPlot ({
    top_pathways(fgsea_data(), input$pathways, input$numpaths)
  })
  
  output$toppadj <- renderDataTable({
    
    #print (input$padjtable)
    #print (input$NESbutton)
    
    #NES_filter <- input$NESbutton
    
    filtered_fgsea(fgsea_data(), input$padjtable, input$NESbutton)
  })
  
  output$downloadTopPadj <- downloadHandler(
    filename = function() {
      paste("toppadj_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(filtered_fgsea(fgsea_data(), input$padjtable, input$NESbutton), file)
    }
  )

  output$padjvsNES <- renderPlot({
    padjNES(fgsea_data(), input$padjscatter)
  })

}
shinyApp(ui = ui, server = server)
