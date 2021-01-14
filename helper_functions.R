library(tidyverse)
library(viridis)
library(colorspace)
library(rlang)
library(MatrixCorrelation)
library(rhdf5)
library(glue)
library(parallel)
library(data.table)
library(gridExtra)
library(SummarizedExperiment)
library(RColorBrewer)
library(wesanderson)
library(pheatmap)
library(glue)
library(ggalluvial)
library(cowplot)
library(readxl)

# Lookups --------------------------------------------------------

# Create lookups 

qmdiab_path <- "data/QMDiab.xls"
schiz_path  <- "data/Schizophrenia.xls"
aml_path    <- "data/AML.xls"

qm_lookup <- 
  read_xls(qmdiab_path, "Metabolite Annotations")


aml_lookup <- 
  read_xls(aml_path, "Metabolite Annotations")


schizo_lookup <- 
  read_xls(schiz_path, "Metabolite Annotations")


# 03_assess_model_performance ---------------------------------------------

get_bootstrap_scores <- function(recon_vae, recon_pca, ref_data, n_times = 100, n_cores = 5) {
  
  # converting to data.table makes all subsequent calculations faster
  recon_vae <- data.table(recon_vae)
  recon_pca <- data.table(recon_pca)
  ref_data  <- data.table(ref_data)
  
  1:n_times %>% 
    mclapply(function(bs_idx) {
      set.seed(bs_idx)
      if (bs_idx %% 100 == 0) print(glue("Working on index {bs_idx}"))
      sample_idxs <- nrow(ref_data) %>% sample(replace = TRUE)
      
      ref_rand <- ref_data[sample_idxs, ]
      vae_rand <- recon_vae[sample_idxs, ]
      pca_rand <- recon_pca[sample_idxs, ]
      
      vae_cor <- cor(vae_rand) 
      pca_cor <- cor(pca_rand) 
      ref_cor <- cor(ref_rand)
      
      # Calculate MSE between the reconstructed correlation matrix and reference
      # data. Calculate this on the upper triangular matrix
      vae_cor_ut <- vae_cor[upper.tri(vae_cor, diag = FALSE)]
      pca_cor_ut <- pca_cor[upper.tri(pca_cor, diag = FALSE)]
      ref_cor_ut <- ref_cor[upper.tri(ref_cor, diag = FALSE)]
      
      vae_cormat_mse <- (vae_cor_ut - ref_cor_ut)^2 %>% mean()
      pca_cormat_mse <- (pca_cor_ut - ref_cor_ut)^2 %>% mean()
      
      vae_mse <- (vae_rand - ref_rand)^2 %>% unlist() %>% mean()
      pca_mse <- (pca_rand - ref_rand)^2 %>% unlist() %>% mean()
      
      tibble(model = c("VAE", "PCA"),
             mse_cormat = c(vae_cormat_mse, pca_cormat_mse),
             mse        = c(vae_mse, pca_mse)) %>% 
        dplyr::mutate(seed = bs_idx)
    }, 
    mc.preschedule = TRUE, 
    mc.cores = 5, 
    mc.cleanup = TRUE) %>% 
    bind_rows()
  
}

get_mse_plot <- function(data_name, title) {
  
  all_metrics_processed %>% 
    filter(d_dim == 18,
           label == "MSE",
           data == data_name) %>% 
    ggplot(aes(x = "", y = value, color = model)) + 
    geom_boxplot() +
    labs(y = "MSE", color = "Model") +
    # theme(axis.title = element_blank()) +
    facet_grid(cols = vars(label),
               rows = vars(data),
               scales = "free") +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.text.x = element_blank(),
          strip.text.y = element_blank(),
          panel.spacing = unit(1, "lines")) +
    scale_color_manual(values=wes_palette(n=2, name="BottleRocket1")) +
    ggtitle(title)
  
}


# 05_intepret_latent_space.R ----------------------------------------------

# Count number of pathways, required to create visualization 
# gaps in the SAGE value heatmaps

aml_mets <- 
  read_xls("data/AML.xls", "Metabolite Data", col_types = "numeric") %>% 
  colnames()

subpw_counts <- 
  aml_lookup %>% 
  filter(COMP_IDstr %in% aml_mets) %>% 
  dplyr::count(pathway_name = SUB_PATHWAY)

superpw_counts <- 
  aml_lookup %>% 
  filter(COMP_IDstr %in% aml_mets) %>% 
  dplyr::count(pathway_name = SUPER_PATHWAY) %>% 
  arrange()

superpw_sub_counts <- 
  aml_lookup %>% 
  filter(COMP_IDstr %in% aml_mets) %>% 
  dplyr::distinct(SUB_PATHWAY, SUPER_PATHWAY) %>% 
  dplyr::count(pathway_name = SUPER_PATHWAY) %>% 
  arrange()


plot_superpw <- function(df) {
  
  plot.new()
  df %>% 
    pheatmap(
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      # cellwidth = 18, cellheight = 16,
      angle_col = 0,
      color = viridis(n = 100),
      border_color = NA)
  
}

plot_subpw <- function(df, cap) {
  
  breaks_list = seq(0, cap, by = 0.01)
  
  plot.new()
  df %>% 
    pheatmap(
      breaks = breaks_list,
      annotation_row = sage_sub_pw,
      annotation_colors = anno_colors,
      annotation_names_row = FALSE,
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      gaps_row = cumsum(superpw_sub_counts$n),
      # cellwidth = 18, cellheight = 16,
      angle_col = 0,
      color = viridis(n = length(breaks_list)),
      border_color = NA
    )
  
}

plot_met <- function(df) {
  
  plot.new()
  df %>% 
    pheatmap(
      # fontsize = 5,
      show_rownames = FALSE,
      annotation_row = sage_met_pw,
      annotation_colors = anno_colors,
      annotation_names_row = FALSE,
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      gaps_row = cumsum(superpw_counts$n),
      # cellwidth = 18, cellheight = 2,
      angle_col = 0,
      # color = magma(100)
      color = viridis(100)
    )
}


get_alluvial_plots <- function(dim_plot, met_df, subpw_df, superpw_df, n_vis) {
  
  subpw_sage <- 
    subpw_df %>% 
    dplyr::select(superpw = SUPER_PATHWAY, 
                  subpw = pathway_name, 
                  subpw_value := dim_plot) %>% 
    slice_max(n = n_vis, order_by = subpw_value) 
  
  superpw_sage <- 
    superpw_df %>% 
    filter(dim == dim_plot) %>% 
    dplyr::select(superpw = pathway_name, superpw_value = sage_value) 
  
  pw_sage <- left_join(superpw_sage, subpw_sage)
  
  met_sage <- 
    met_df %>% 
    dplyr::select(met = BIOCHEMICAL, 
                  subpw = SUB_PATHWAY, 
                  superpw = SUPER_PATHWAY,
                  met_value := dim_plot) %>% 
    slice_max(n = n_vis, order_by = met_value) 
  
  ordered_sages <- 
    pw_sage %>% 
    inner_join(met_sage) %>% 
    mutate_at(c("superpw", "subpw", "met"), function(col) str_wrap(col, width = str_wrap_width)) %>%
    dplyr::transmute(superpw = fct_reorder(superpw, superpw_value),
                     subpw = if_else(subpw == "Unknown", "Unknown ", subpw),
                     subpw   = fct_reorder(subpw,   subpw_value),
                     met     = fct_reorder(met,     met_value)) %>% 
    mutate_all(fct_rev)
  
  
  
  ordered_sages %>% 
    ggplot(aes(axis1 = superpw, axis2 = subpw, axis3 = met, 
               fill = superpw)) +
    scale_fill_manual(values = alluvial_color %>% unlist()) +
    geom_alluvium(width = 1/2) +
    geom_stratum(color = NA, width = 1/2, alpha = 0.65) +
    ggfittext::geom_fit_text(stat = "stratum", aes(label = after_stat(stratum)), width = 1/2, min.size = 1, contrast = TRUE) +
    theme_minimal() +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank()) +
    xlim(0.75, 3.25) 
  
}



# 07_associate_dimensions_with_diseases ----------------------------------------------

process_qmdiab <- function(df) {
  df %>% 
    dplyr::select(-c("PC1", "PC2", "PC3", "SOMAPC1", "SOMAPC2", "SOMAPC3")) %>%
    gather(latent_dim, value, -c(SampleId, SEX, AGE, DIAB, BMI)) %>%
    mutate(DIAB = if_else(DIAB == 0, "group_1", "group_2")) %>% 
    dplyr::select(latent_dim, group = DIAB, value)
}

process_aml <- function(df, regex_) {
  df %>% 
    dplyr::select(IDENTITY, sex, age, response, matches(regex_)) %>% 
    gather(latent_dim, value, -c(IDENTITY, sex, age, response)) %>%
    mutate(response = if_else(response == 1, "group_1", "group_2")) %>% 
    dplyr::select(latent_dim, group = response, value)
}

process_schizo <- function(df, regex_) {
  df %>% 
    dplyr::select(SAMPLE_ID, GENDER, AGE, Group, matches(regex_)) %>%
    gather(latent_dim, value, -c(SAMPLE_ID, GENDER, AGE, Group)) %>%
    mutate(Group = if_else(Group == "case", "group_1", "group_2")) %>% 
    dplyr::select(latent_dim, group = Group, value) 
}




get_latent_scores <- function(df) {
  
  df_auc <- 
    df %>% 
    group_by(latent_dim, model) %>% 
    dplyr::summarise(auc = pROC::auc(group, value) %>% as.numeric())
  
  df %>%
    group_by(model, latent_dim, group) %>% 
    nest() %>% 
    spread(key = group, value = data) %>% 
    mutate(
      pval_nominal = map2(group_1, group_2, ~{t.test(.x$value, .y$value)$p.value}),
      group_1 = map(group_1, nrow),
      group_2 = map(group_2, nrow)
    ) %>% 
    unnest(cols = c(group_1, group_2, pval_nominal)) %>% 
    ungroup() %>% 
    group_by(model) %>%
    mutate(pval_fdr = p.adjust(pval_nominal, method = "fdr")) %>% 
    inner_join(df_auc, by = c("latent_dim", "model")) %>%
    dplyr::arrange(pval_nominal) %>% 
    ungroup()
}


get_rank_plot <- function(df, title) {
  
  df %>% 
    dplyr::arrange(pval_nominal) %>% 
    mutate(rank = row_number()) %>% 
    ggplot() + 
    geom_point(aes(x = rank, y = -log10(pval_nominal), color = model), size = 1.3) +
    ggtitle(title) +
    scale_color_manual(values=wes_palette(n=2, name="BottleRocket1")) +
    theme_minimal() +
    theme(axis.text=element_text(size = 9)) +
    labs(x = "Rank", y = expression(-log[10](p-value)), color = "Model") +
    theme(
      legend.position = c(0.9, 0.6), # c(0,0) bottom left, c(1,1) top-right.
      legend.background = element_rect(fill = "white", colour = NA)
    )
  
} 

add_small_legend <- function(myPlot, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

get_association_boxplot <- function(df_values, df_scores) {
  
  dim_order <- 
    df_scores %>% 
    .$latent_dim
  
  dim_label <- 
    df_scores %>%
    # dplyr::mutate(labels = glue::glue("dim_{latent_dim}\npval = {signif(pval_nominal, 2)}\nAUC = {round(auc, 2)}")) %>%
    dplyr::mutate(labels = glue::glue("Dim. {latent_dim}\np = {signif(pval_nominal, 2)}")) %>%
    .$labels
  
  df_values %>% 
    # group_by(model)
    dplyr::mutate(latent_dim = factor(latent_dim, levels = dim_order, labels = dim_label)) %>% 
    ggplot(aes(x=model, y=value, fill=group)) + 
    geom_boxplot() +
    labs(y = "Latent value") +
    facet_grid(cols = vars(latent_dim), 
               scale = "free") +
    
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_text(size = 12),
          axis.text.y=element_text(size = 9),
          axis.title.y=element_text(size = 12),
          axis.ticks.x=element_blank(),
          strip.text.x = element_text(size = 3),
          # plot.title = element_text(size=8),
          legend.margin = margin(c(0,0,0,0)),
          panel.spacing = unit(1, "lines")) +
    labs(fill = "Group") 
  
}

get_clinical_var_heatmaps <- function(model_name) {
  qm_enc <- 
    glue("{encode_path}//QMDiab_{model_name}_encoding.csv") %>% 
    read_csv() %>% 
    dplyr::select(-c("PC1", "PC2", "PC3", "SOMAPC1", "SOMAPC2", "SOMAPC3")) %>%
    gather(latent_dim, value, -c(SampleId, SEX, AGE, DIAB, BMI)) %>% 
    mutate(DIAB = if_else(DIAB == 0, FALSE, TRUE),
           SEX = if_else(SEX == 0, "Female", "Male"))
  
  qm_continuous <- 
    qm_enc %>% 
    gather("clinical_variable", "clin_value", -c(SampleId, SEX, DIAB, latent_dim, value))
  
  qm_clin <- read_xls(qmdiab_path, "Clinical Parameters")
    
  clin_vars <- 
    qm_clin %>% 
    right_join(qm_enc, by = "SampleId") %>% 
    bind_rows(qm_continuous)
  
  qm_cor <- 
    clin_vars %>% 
    filter(!is.na(clin_value)) %>% 
    group_by(latent_dim, clinical_variable) %>% 
    dplyr::summarise(pval_nominal = cor.test(clin_value, value, method = "pearson")$p.value,
                     corval = round(cor(clin_value, value, method = "pearson"), 2)) %>%  
    dplyr::mutate(pval_score = -log10(pval_nominal)*sign(corval)) %>%
    ungroup() %>% 
    pivot_wider(id_col = clinical_variable, names_from = latent_dim, values_from = pval_score) %>% 
    column_to_rownames("clinical_variable") 
  
  plot.new()
  break_rng = seq(0, 25, by = 1)
  qm_cor %>% 
    abs() %>% 
    pheatmap(
      # cellwidth = 12, cellheight = 10,
      breaks = break_rng,
      angle_col = 0,
      cutree_cols = 3,
      cutree_rows = 3,
      color = viridis(n = length(break_rng)),
    )
}

get_aml_pheno_associations <- function(df, pheno_info) {
  
  pheno_info
  pheno_info[!pheno_info %in% c("t_8_16_", "del_9q_", "tri_21_", "tri_11_", "t69", "EVI1pos", "IDH2_R172K", "RUNX1", "ASXL1_Methylation", "FLT3_TKD", "PTEN", "TP53")] %>% 
    lapply(function(ph) {
      print(ph)
      df_clean <- 
        df %>% 
        dplyr::select(SAMPLE_NAME, value, latent_dim, !!sym(ph), model) %>% 
        filter(!is.na(!!sym(ph)))
      
      df_clean$group <- str_c("group_", group_indices(df_clean, !!sym(ph)))
      
      df_clean %>%
        dplyr::select(-SAMPLE_NAME, -!!sym(ph)) %>% 
        get_latent_scores() %>% 
        dplyr::mutate(gene = ph)
    }) %>% 
    bind_rows()
}

get_aml_mut_heatmap <- function(df){
  
  plot.new()
  df %>% 
    mutate(pval_score = -log10(pval_nominal)) %>% 
    dplyr::select(latent_dim, pval_score, gene) %>% 
    pivot_wider(names_from = gene, values_from = pval_score) %>% 
    column_to_rownames("latent_dim") %>% 
    t() %>% 
    as.data.frame() %>% 
    pheatmap(
      # cellwidth = 12, cellheight = 10,
      angle_col = 0,
      cutree_cols = 3,
      cutree_rows = 3,
      color = viridis(n = 10)
    )
  
}

