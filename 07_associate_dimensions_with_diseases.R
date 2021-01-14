# Associate latent dimensions to disease groups
source("helper_functions.R")


# Load encodings ----------------------

encode_path <- "results/encodings//"

qm_encoding <- 
  bind_rows(
    glue("{encode_path}//QMDiab_VAE_encoding.csv") %>% 
      read_csv() %>% 
      process_qmdiab() %>% 
      dplyr::mutate(model = "VAE"),
    glue("{encode_path}//QMDiab_PCA_encoding.csv") %>% 
      read_csv() %>% 
      process_qmdiab() %>% 
      dplyr::mutate(model = "PCA")
    
  )
qm_scores <- qm_encoding %>% get_latent_scores()

aml_encoding <- 
  bind_rows(
    glue("{encode_path}//AML_VAE_encoding.csv") %>% 
      read_csv() %>% 
      process_aml("^[0-9]+") %>% 
      dplyr::mutate(model = "VAE"),
    glue("{encode_path}//AML_PCA_encoding.csv") %>% 
      read_csv() %>% 
      process_aml("^[0-9]+") %>% 
      dplyr::mutate(model = "PCA")
  )
aml_scores <- aml_encoding %>% get_latent_scores()

schizo_encoding <- 
  bind_rows(
    glue("{encode_path}//Schizo_VAE_encoding.csv") %>% 
      read_csv() %>% 
      process_schizo("^[0-9]+") %>% 
      dplyr::mutate(model = "VAE"),
    glue("{encode_path}//Schizo_PCA_encoding.csv") %>% 
      read_csv() %>% 
      process_schizo("^[0-9]+") %>% 
      dplyr::mutate(model = "PCA")
  )
schizo_scores <- schizo_encoding %>% get_latent_scores() 



# Univariate analysis -----------------------------------

# Load raw data
qmdiab_path <- "data/QMDiab.xls"
schiz_path  <- "data/Schizophrenia.xls"
aml_path    <- "data/AML.xls"

# Overlap data
qm_met_score <- 
  read_xls(qmdiab_path, "Sample Annotations") %>% 
  bind_cols(read_xls(qmdiab_path, "Metabolite Data", col_types ="numeric")) %>% 
  process_qmdiab() %>% 
  mutate(model = "none") %>% 
  get_latent_scores() %>% 
  left_join(qm_lookup %>% select(COMP_IDstr, BIOCHEMICAL), by = c("latent_dim" = "COMP_IDstr")) %>% 
  dplyr::select(BIOCHEMICAL, pval_nominal, pval_fdr, auc)

aml_met_score <- 
  read_xls(aml_path, "Sample Annotations") %>% 
  bind_cols(read_xls(aml_path, "Metabolite Data", col_types ="numeric")) %>% 
  process_aml("M[0-9]+") %>% 
  mutate(model = "none") %>% 
  get_latent_scores() %>% 
  left_join(aml_lookup %>% select(COMP_IDstr, BIOCHEMICAL), by = c("latent_dim" = "COMP_IDstr")) %>% 
  dplyr::select(BIOCHEMICAL, pval_nominal, pval_fdr, auc)

schizo_met_score <-
  read_xls(schiz_path, "Sample Annotations") %>% 
  bind_cols(read_xls(schiz_path, "Metabolite Data", col_types ="numeric")) %>% 
  drop_na() %>% 
  process_schizo("M[0-9]+") %>% 
  mutate(model = "none") %>% 
  get_latent_scores() %>% 
  left_join(schizo_lookup %>% select(COMP_IDstr, BIOCHEMICAL), by = c("latent_dim" = "COMP_IDstr")) %>% 
  dplyr::select(BIOCHEMICAL, pval_nominal, pval_fdr, auc)



# Create pvalue rank plots for all three datasets ------------------------------------------------

qm_rankplot     <- qm_scores     %>% get_rank_plot("QMDiab")
aml_rankplot    <- aml_scores    %>% get_rank_plot("AML")
schizo_rankplot <- schizo_scores %>% get_rank_plot("Schizo")



# Create association boxplots for all three datasets ------------------------------------------------

# create color palette
dim_pal <- brewer.pal(6, name = "Dark2")
clin_groups <- c("Diabetic", "Non-diabetic",
                 "Schizo.", "Non-schizo.",
                 "Full response", "No response")

brewer_pal <- brewer.pal(length(clin_groups), name = "Set1")
dim_color <- 
  clin_groups %>% 
  sapply(function(grp) brewer_pal[which(clin_groups == grp)],
         USE.NAMES = TRUE,
         simplify = FALSE)

dim_color$`No response` <- "#666666"


# QMDiab
qm_boxplot <- 
  get_association_boxplot(
    qm_encoding %>% 
      filter(model == "VAE" & latent_dim == "9" |
               model == "PCA" & latent_dim == "16") %>% 
      mutate(group = ifelse(group == "group_1", "Non-diabetic", "Diabetic"),
             group = fct_relevel(group, "Diabetic", "Non-diabetic")),
    qm_scores %>% 
      filter(model == "VAE" & latent_dim == "9" |
               model == "PCA" & latent_dim == "16")
  ) 
qm_boxplot <- 
  qm_boxplot %>% 
  add_small_legend(myPlot = ., textSize = 8, spaceLegend = 0.75) +
  scale_fill_manual(values=dim_color)

# Schizo
schizo_boxplot <- 
  get_association_boxplot(
    schizo_encoding %>% 
      mutate(model = fct_relevel(model, "VAE")) %>% 
      filter(model == "VAE" & latent_dim == "11" |
               model == "PCA" & latent_dim == "15") %>% 
      mutate(group = ifelse(group == "group_1", "Schizo.", "Non-schizo."),
             group = fct_relevel(group, "Schizo.", "Non-schizo.")) %>% 
      drop_na(),
    schizo_scores %>% 
      mutate(model = fct_relevel(model, "VAE")) %>% 
      filter(model == "VAE" & latent_dim == "11" |
               model == "PCA" & latent_dim == "15") %>% 
      drop_na()
  )

schizo_boxplot <- 
  schizo_boxplot %>% 
  add_small_legend(myPlot = ., textSize = 8, spaceLegend = 0.75) +
  scale_fill_manual(values=dim_color)


# AML
aml_boxplot <- 
  get_association_boxplot(
    aml_encoding %>% 
      mutate(model = fct_relevel(model, "VAE")) %>% 
      filter(model == "VAE" & latent_dim == "15" |
               model == "PCA" & latent_dim == "10") %>% 
      mutate(group = ifelse(group == "group_1", "Full response", "No response"),
             group = fct_relevel(group, "Full response", "No response")),
    aml_scores %>% 
      mutate(model = fct_relevel(model, "VAE")) %>% 
      filter(model == "VAE" & latent_dim == "15" |
               model == "PCA" & latent_dim == "10")
  ) %>% 
  add_small_legend()
aml_boxplot <- 
  aml_boxplot %>% 
  add_small_legend(myPlot = ., textSize = 8, spaceLegend = 0.75) +
  scale_fill_manual(values=dim_color)



# QMDiab: clinical variable associations -----------------------------

vae_qm_cor_heatmap <- get_clinical_var_heatmaps("VAE")
pca_qm_cor_heatmap <- get_clinical_var_heatmaps("PCA")



# AML: mutation scores and plots -----------------------------------

pheno_data <- 
  read_xls(aml_path, "Mutation Data")

phenos <- 
  pheno_data %>%
  pivot_longer(cols = -SAMPLE_NAME) %>% 
  filter(!is.na(value)) %>% 
  dplyr::count(name, value) %>% 
  # filter for mutations with at least 10 samples
  filter(n >= 10) %>% 
  dplyr::select(-n) %>% 
  dplyr::count(name) %>% 
  filter(n == 2) %>% 
  .$name


aml_vae_pheno <- 
  glue("{encode_path}//AML_VAE_encoding.csv") %>% 
  read_csv() %>% 
  dplyr::select(SAMPLE_NAME, as.character(1:18)) %>% 
  inner_join(pheno_data) %>%
  pivot_longer(cols = matches("^[0-9]+"), names_to = "latent_dim") %>% 
  dplyr::mutate(model = "VAE")
aml_vae_pheno_score <- 
  aml_vae_pheno %>% 
  get_aml_pheno_associations(phenos)

aml_pca_pheno <- 
  glue("{encode_path}//AML_PCA_encoding.csv") %>% 
  read_csv() %>% 
  dplyr::select(SAMPLE_NAME, as.character(1:18)) %>% 
  inner_join(pheno_data) %>%
  pivot_longer(cols = matches("^[0-9]+"), names_to = "latent_dim") %>% 
  dplyr::mutate(model = "PCA") 
aml_pca_pheno_score <- 
  aml_pca_pheno %>% 
  get_aml_pheno_associations(phenos)


aml_vae_pheno_heatmap <- aml_vae_pheno_score %>% get_aml_mut_heatmap()
aml_pca_pheno_heatmap <- aml_pca_pheno_score %>% get_aml_mut_heatmap()


# Plot top 2 associations, i.e. NPM1 and IDH
aml_npm1_boxplot <- 
  get_association_boxplot(
    bind_rows(aml_vae_pheno, aml_pca_pheno) %>% 
      dplyr::select(latent_dim, group = "NPM1", value, model) %>% 
      filter(!is.na(group)) %>% 
      filter(model == "VAE" & latent_dim == "8" |
               model == "PCA" & latent_dim == "8") %>% 
      mutate(model = fct_relevel(model, "VAE", "PCA"),
             latent_dim = if_else(model == "PCA", "8 ", "8")),
    bind_rows(aml_vae_pheno_score, aml_pca_pheno_score) %>% 
      filter(latent_dim == 8,
             gene == "NPM1") %>% 
      mutate(latent_dim = if_else(model == "PCA", "8 ", "8"))
  ) +
  ggtitle("NPM1") +
  theme(strip.text.x = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))

aml_idh_boxplot <- 
  get_association_boxplot(
    bind_rows(aml_vae_pheno, aml_pca_pheno) %>% 
      dplyr::select(latent_dim, group = "IDH", value, model) %>% 
      filter(!is.na(group)) %>% 
      filter(model == "VAE" & latent_dim == "15" |
               model == "PCA" & latent_dim == "8") %>% 
      mutate(model = fct_relevel(model, "VAE", "PCA")),
    bind_rows(aml_vae_pheno_score, aml_pca_pheno_score) %>% 
      filter(model == "VAE" & latent_dim == "15" |
               model == "PCA" & latent_dim == "8") %>% 
      filter(gene == "IDH")
  ) +
  ggtitle("IDH") +
  theme(strip.text.x = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))



# Show plots -------------------------------------------------------

# Latent dimension rank plots
qm_rankplot
schizo_rankplot
aml_rankplot

# Patient group box plots of highest scoring dimensions
qm_boxplot
schizo_boxplot
aml_boxplot

# Type 2 diabetes clinical variable heatmaps
plot.new()
vae_qm_cor_heatmap

plot.new()
pca_qm_cor_heatmap

# AML mutation heatmaps
plot.new()
aml_vae_pheno_heatmap

plot.new()
aml_pca_pheno_heatmap

# AML, NPM1 and IDH box plots
aml_npm1_boxplot
aml_idh_boxplot



