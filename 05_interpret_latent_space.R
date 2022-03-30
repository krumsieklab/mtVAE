source("helper_functions.R")


# Create color palettes ----------------------------------------------------

superpws <- aml_lookup$SUPER_PATHWAY %>% unique()
brewer_pal <- brewer.pal(length(superpws), name = "Set1")
brewer_pal[[1]] <- "#FFBE0B" # Lipid superpathway
brewer_pal[[6]] <- "#EFC3E6" # Nucleotide superpathway
brewer_pal[[length(brewer_pal)]] <- "#adb5bd" # For the Unknown superpathway

superpws_color <- 
  superpws %>% 
  sapply(function(pw) brewer_pal[which(superpws == pw)],
         USE.NAMES = TRUE,
         simplify = FALSE)

anno_colors <- 
  superpws_color %>% 
  # Add opacity alpha of 0.5 (following ggplot convention).
  # This is equivalent to the hex value 80:
  # https://gist.github.com/lopspower/03fb1cc0ac9f32ef38f4
  lapply(function(val) str_c(val, "FF")) %>%
  unlist() %>%
  .[order(names(.))] %>%
  list(SUPER_PATHWAY = .)

# Alluvial plot palette
str_wrap_width <- 20
alluvial_superpws <- 
  superpws %>% 
  str_wrap(width = str_wrap_width)

alluvial_color <- 
  alluvial_superpws %>% 
  sapply(function(pw) brewer_pal[which(alluvial_superpws == pw)],
         USE.NAMES = TRUE,
         simplify = FALSE)

alluvial_color$Lipid <- "#FFBE0B"
alluvial_color$Unknown <- "#adb5bd"



# Load VAE SAGE scores --------------------------------------------------

sage_scores <- 
  list.files("results/sage_values/VAE", full.names = TRUE, pattern = ".csv") %>% 
  lapply(function(fn) read_csv(fn) %>% dplyr::mutate(file = fn)) %>% 
  bind_rows() %>% 
  dplyr::mutate(dim = str_extract(file, "dim_[0-9]+"),
                sage_type = str_extract(file, "met|sub|super")) 

# Metabolites
sage_met_df <- 
  sage_scores %>% 
  filter(sage_type == "met") %>% 
  dplyr::select(-c(pathway_name, file, sage_type)) %>% 
  distinct() 

sage_met_wide <- 
  sage_met_df %>% 
  dplyr::select(-sage_value_sd) %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  left_join(aml_lookup, by = c("metabolite_id" = "COMP_IDstr")) %>% 
  # for visualization purposes, sort by metabolite with highest
  # absolute mean values
  dplyr::mutate(met_mean = rowMeans(abs(select(., starts_with("dim"))))) %>% 
  group_by(SUPER_PATHWAY) %>% 
  dplyr::arrange(met_mean) %>% 
  dplyr::select(-met_mean) %>% 
  ungroup() %>% 
  dplyr::arrange(SUPER_PATHWAY) %>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

# used for heatmap splitting
sage_met_pw <- 
  sage_met_wide %>% 
  dplyr::transmute(BIOCHEMICAL, 
                   SUPER_PATHWAY = factor(SUPER_PATHWAY)) %>% 
  column_to_rownames("BIOCHEMICAL")

sage_met <- 
  sage_met_wide %>% 
  dplyr::select(-c(metabolite_id, SUB_PATHWAY, SUPER_PATHWAY)) %>% 
  column_to_rownames("BIOCHEMICAL")



# Sub-pathways
sage_subpw_df <- 
  sage_scores %>% 
  filter(sage_type == "sub") %>% 
  left_join(aml_lookup %>% distinct(SUB_PATHWAY, SUPER_PATHWAY), 
            by = c("pathway_name" = "SUB_PATHWAY")) 


sage_subpw_wide <- 
  sage_subpw_df %>% 
  dplyr::select(-c(metabolite_id, file, sage_type, sage_value_sd)) %>% 
  distinct() %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  dplyr::mutate(pw_mean = rowMeans(abs(select(., starts_with("dim"))))) %>% 
  group_by(SUPER_PATHWAY) %>% 
  dplyr::arrange(pw_mean) %>% 
  dplyr::select(-pw_mean) %>% 
  ungroup() %>% 
  dplyr::arrange(SUPER_PATHWAY) %>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

# used for heatmap splitting
sage_sub_pw <- 
  sage_subpw_wide %>% 
  dplyr::transmute(pathway_name, 
                   SUPER_PATHWAY = factor(SUPER_PATHWAY)) %>% 
  column_to_rownames("pathway_name")

sage_subpw <- 
  sage_subpw_wide %>% 
  dplyr::select(-SUPER_PATHWAY) %>% 
  column_to_rownames("pathway_name")


# Super-pathways

sage_superpw_df <- 
  sage_scores %>% 
  filter(sage_type == "super")

sage_superpw <- 
  sage_superpw_df %>% 
  dplyr::select(-c(metabolite_id, file, sage_type, sage_value_sd)) %>% 
  distinct() %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  column_to_rownames("pathway_name")%>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())
  


# Load PCA SAGE scores -----------------------------------------------------

pca_sage_scores <- 
  list.files("results/sage_values/PCA", full.names = TRUE) %>% 
  lapply(function(fn) read_csv(fn) %>% dplyr::mutate(file = fn)) %>% 
  bind_rows() %>% 
  dplyr::mutate(dim = str_extract(file, "dim_[0-9]+"),
                sage_type = str_extract(file, "met|sub|super")) 

# metabolites
pca_met_sage_df <- 
  pca_sage_scores %>% 
  filter(sage_type == "met") %>% 
  dplyr::select(-c(pathway_name, file, sage_type)) %>% 
  distinct() 

pca_met_sage_wide <- 
  pca_met_sage_df %>% 
  dplyr::select(-sage_value_sd) %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  left_join(aml_lookup, by = c("metabolite_id" = "COMP_IDstr")) %>% 
  # for visualization purposes, sort by metabolite with highest
  # absolute mean values
  dplyr::mutate(met_mean = rowMeans(abs(select(., starts_with("dim"))))) %>% 
  group_by(SUPER_PATHWAY) %>% 
  dplyr::arrange(met_mean) %>% 
  dplyr::select(-met_mean) %>% 
  ungroup() %>% 
  dplyr::arrange(SUPER_PATHWAY) %>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

# Sub-pathways
pca_sage_subpw_df <- 
  pca_sage_scores %>% 
  filter(sage_type == "sub") %>% 
  left_join(aml_lookup %>% distinct(SUB_PATHWAY, SUPER_PATHWAY), 
            by = c("pathway_name" = "SUB_PATHWAY")) 


pca_sage_subpw_wide <- 
  pca_sage_subpw_df %>% 
  dplyr::select(-c(metabolite_id, file, sage_type, sage_value_sd)) %>% 
  distinct() %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  dplyr::mutate(pw_mean = rowMeans(abs(select(., starts_with("dim"))))) %>% 
  group_by(SUPER_PATHWAY) %>% 
  dplyr::arrange(pw_mean) %>% 
  dplyr::select(-pw_mean) %>% 
  ungroup() %>% 
  dplyr::arrange(SUPER_PATHWAY) %>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

pca_sage_subpw <- 
  pca_sage_subpw_wide %>% 
  dplyr::select(-SUPER_PATHWAY) %>% 
  column_to_rownames("pathway_name")
  

# Super-pathways
pca_sage_superpw_df <- 
  pca_sage_scores %>% 
  filter(sage_type == "super")

pca_sage_superpw <- 
  pca_sage_superpw_df %>% 
  dplyr::select(-c(metabolite_id, file, sage_type, sage_value_sd)) %>% 
  distinct() %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  column_to_rownames("pathway_name")%>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

# Load cosine KPCA SAGE scores -----------------------------------------------------

cosine_sage_scores <- 
  list.files("results/sage_values/cosine", full.names = TRUE) %>% 
  lapply(function(fn) read_csv(fn) %>% dplyr::mutate(file = fn)) %>% 
  bind_rows() %>% 
  dplyr::mutate(dim = str_extract(file, "dim_[0-9]+"),
                sage_type = str_extract(file, "met|sub|super")) 

# metabolites
cosine_met_sage_df <- 
  cosine_sage_scores %>% 
  filter(sage_type == "met") %>% 
  dplyr::select(-c(pathway_name, file, sage_type)) %>% 
  distinct() 

cosine_met_sage_wide <- 
  cosine_met_sage_df %>% 
  dplyr::select(-sage_value_sd) %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  left_join(aml_lookup, by = c("metabolite_id" = "COMP_IDstr")) %>% 
  # for visualization purposes, sort by metabolite with highest
  # absolute mean values
  dplyr::mutate(met_mean = rowMeans(abs(select(., starts_with("dim"))))) %>% 
  group_by(SUPER_PATHWAY) %>% 
  dplyr::arrange(met_mean) %>% 
  dplyr::select(-met_mean) %>% 
  ungroup() %>% 
  dplyr::arrange(SUPER_PATHWAY) %>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

# Sub-pathways
cosine_sage_subpw_df <- 
  cosine_sage_scores %>% 
  filter(sage_type == "sub") %>% 
  left_join(aml_lookup %>% distinct(SUB_PATHWAY, SUPER_PATHWAY), 
            by = c("pathway_name" = "SUB_PATHWAY")) 


cosine_sage_subpw_wide <- 
  cosine_sage_subpw_df %>% 
  dplyr::select(-c(metabolite_id, file, sage_type, sage_value_sd)) %>% 
  distinct() %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  dplyr::mutate(pw_mean = rowMeans(abs(select(., starts_with("dim"))))) %>% 
  group_by(SUPER_PATHWAY) %>% 
  dplyr::arrange(pw_mean) %>% 
  dplyr::select(-pw_mean) %>% 
  ungroup() %>% 
  dplyr::arrange(SUPER_PATHWAY) %>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

cosine_sage_subpw <- 
  cosine_sage_subpw_wide %>% 
  dplyr::select(-SUPER_PATHWAY) %>% 
  column_to_rownames("pathway_name")


# Super-pathways
cosine_sage_superpw_df <- 
  cosine_sage_scores %>% 
  filter(sage_type == "super")

cosine_sage_superpw <- 
  cosine_sage_superpw_df %>% 
  dplyr::select(-c(metabolite_id, file, sage_type, sage_value_sd)) %>% 
  distinct() %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  column_to_rownames("pathway_name")%>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

# Load sigmoid SAGE scores -----------------------------------------------------

sigmoid_sage_scores <- 
  list.files("results/sage_values/sigmoid", full.names = TRUE) %>% 
  lapply(function(fn) read_csv(fn) %>% dplyr::mutate(file = fn)) %>% 
  bind_rows() %>% 
  dplyr::mutate(dim = str_extract(file, "dim_[0-9]+"),
                sage_type = str_extract(file, "met|sub|super")) 

# metabolites
sigmoid_met_sage_df <- 
  sigmoid_sage_scores %>% 
  filter(sage_type == "met") %>% 
  dplyr::select(-c(pathway_name, file, sage_type)) %>% 
  distinct() 

sigmoid_met_sage_wide <- 
  sigmoid_met_sage_df %>% 
  dplyr::select(-sage_value_sd) %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  left_join(aml_lookup, by = c("metabolite_id" = "COMP_IDstr")) %>% 
  # for visualization purposes, sort by metabolite with highest
  # absolute mean values
  dplyr::mutate(met_mean = rowMeans(abs(select(., starts_with("dim"))))) %>% 
  group_by(SUPER_PATHWAY) %>% 
  dplyr::arrange(met_mean) %>% 
  dplyr::select(-met_mean) %>% 
  ungroup() %>% 
  dplyr::arrange(SUPER_PATHWAY) %>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

# Sub-pathways
sigmoid_sage_subpw_df <- 
  sigmoid_sage_scores %>% 
  filter(sage_type == "sub") %>% 
  left_join(aml_lookup %>% distinct(SUB_PATHWAY, SUPER_PATHWAY), 
            by = c("pathway_name" = "SUB_PATHWAY")) 


sigmoid_sage_subpw_wide <- 
  sigmoid_sage_subpw_df %>% 
  dplyr::select(-c(metabolite_id, file, sage_type, sage_value_sd)) %>% 
  distinct() %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  dplyr::mutate(pw_mean = rowMeans(abs(select(., starts_with("dim"))))) %>% 
  group_by(SUPER_PATHWAY) %>% 
  dplyr::arrange(pw_mean) %>% 
  dplyr::select(-pw_mean) %>% 
  ungroup() %>% 
  dplyr::arrange(SUPER_PATHWAY) %>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

sigmoid_sage_subpw <- 
  sigmoid_sage_subpw_wide %>% 
  dplyr::select(-SUPER_PATHWAY) %>% 
  column_to_rownames("pathway_name")


# Super-pathways
sigmoid_sage_superpw_df <- 
  sigmoid_sage_scores %>% 
  filter(sage_type == "super")

sigmoid_sage_superpw <- 
  sigmoid_sage_superpw_df %>% 
  dplyr::select(-c(metabolite_id, file, sage_type, sage_value_sd)) %>% 
  distinct() %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  column_to_rownames("pathway_name")%>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())


# Load rbf SAGE scores -----------------------------------------------------

rbf_sage_scores <- 
  list.files("results/sage_values/rbf", full.names = TRUE) %>% 
  lapply(function(fn) read_csv(fn) %>% dplyr::mutate(file = fn)) %>% 
  bind_rows() %>% 
  dplyr::mutate(dim = str_extract(file, "dim_[0-9]+"),
                sage_type = str_extract(file, "met|sub|super")) 

# metabolites
rbf_met_sage_df <- 
  rbf_sage_scores %>% 
  filter(sage_type == "met") %>% 
  dplyr::select(-c(pathway_name, file, sage_type)) %>% 
  distinct() 

rbf_met_sage_wide <- 
  rbf_met_sage_df %>% 
  dplyr::select(-sage_value_sd) %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  left_join(aml_lookup, by = c("metabolite_id" = "COMP_IDstr")) %>% 
  # for visualization purposes, sort by metabolite with highest
  # absolute mean values
  dplyr::mutate(met_mean = rowMeans(abs(select(., starts_with("dim"))))) %>% 
  group_by(SUPER_PATHWAY) %>% 
  dplyr::arrange(met_mean) %>% 
  dplyr::select(-met_mean) %>% 
  ungroup() %>% 
  dplyr::arrange(SUPER_PATHWAY) %>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

# Sub-pathways
rbf_sage_subpw_df <- 
  rbf_sage_scores %>% 
  filter(sage_type == "sub") %>% 
  left_join(aml_lookup %>% distinct(SUB_PATHWAY, SUPER_PATHWAY), 
            by = c("pathway_name" = "SUB_PATHWAY")) 


rbf_sage_subpw_wide <- 
  rbf_sage_subpw_df %>% 
  dplyr::select(-c(metabolite_id, file, sage_type, sage_value_sd)) %>% 
  distinct() %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  dplyr::mutate(pw_mean = rowMeans(abs(select(., starts_with("dim"))))) %>% 
  group_by(SUPER_PATHWAY) %>% 
  dplyr::arrange(pw_mean) %>% 
  dplyr::select(-pw_mean) %>% 
  ungroup() %>% 
  dplyr::arrange(SUPER_PATHWAY) %>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

rbf_sage_subpw <- 
  rbf_sage_subpw_wide %>% 
  dplyr::select(-SUPER_PATHWAY) %>% 
  column_to_rownames("pathway_name")


# Super-pathways
rbf_sage_superpw_df <- 
  rbf_sage_scores %>% 
  filter(sage_type == "super")

rbf_sage_superpw <- 
  rbf_sage_superpw_df %>% 
  dplyr::select(-c(metabolite_id, file, sage_type, sage_value_sd)) %>% 
  distinct() %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  column_to_rownames("pathway_name")%>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

# Load poly SAGE scores -----------------------------------------------------

poly_sage_scores <- 
  list.files("results/sage_values/poly", full.names = TRUE) %>% 
  lapply(function(fn) read_csv(fn) %>% dplyr::mutate(file = fn)) %>% 
  bind_rows() %>% 
  dplyr::mutate(dim = str_extract(file, "dim_[0-9]+"),
                sage_type = str_extract(file, "met|sub|super")) 

# metabolites
poly_met_sage_df <- 
  poly_sage_scores %>% 
  filter(sage_type == "met") %>% 
  dplyr::select(-c(pathway_name, file, sage_type)) %>% 
  distinct() 

poly_met_sage_wide <- 
  poly_met_sage_df %>% 
  dplyr::select(-sage_value_sd) %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  left_join(aml_lookup, by = c("metabolite_id" = "COMP_IDstr")) %>% 
  # for visualization purposes, sort by metabolite with highest
  # absolute mean values
  dplyr::mutate(met_mean = rowMeans(abs(select(., starts_with("dim"))))) %>% 
  group_by(SUPER_PATHWAY) %>% 
  dplyr::arrange(met_mean) %>% 
  dplyr::select(-met_mean) %>% 
  ungroup() %>% 
  dplyr::arrange(SUPER_PATHWAY) %>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

# Sub-pathways
poly_sage_subpw_df <- 
  poly_sage_scores %>% 
  filter(sage_type == "sub") %>% 
  left_join(aml_lookup %>% distinct(SUB_PATHWAY, SUPER_PATHWAY), 
            by = c("pathway_name" = "SUB_PATHWAY")) 


poly_sage_subpw_wide <- 
  poly_sage_subpw_df %>% 
  dplyr::select(-c(metabolite_id, file, sage_type, sage_value_sd)) %>% 
  distinct() %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  dplyr::mutate(pw_mean = rowMeans(abs(select(., starts_with("dim"))))) %>% 
  group_by(SUPER_PATHWAY) %>% 
  dplyr::arrange(pw_mean) %>% 
  dplyr::select(-pw_mean) %>% 
  ungroup() %>% 
  dplyr::arrange(SUPER_PATHWAY) %>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())

poly_sage_subpw <- 
  poly_sage_subpw_wide %>% 
  dplyr::select(-SUPER_PATHWAY) %>% 
  column_to_rownames("pathway_name")


# Super-pathways
poly_sage_superpw_df <- 
  poly_sage_scores %>% 
  filter(sage_type == "super")

poly_sage_superpw <- 
  poly_sage_superpw_df %>% 
  dplyr::select(-c(metabolite_id, file, sage_type, sage_value_sd)) %>% 
  distinct() %>% 
  pivot_wider(names_from = dim, values_from = sage_value) %>% 
  column_to_rownames("pathway_name")%>% 
  relocate(lapply(1:18, function(dm) str_c("dim_",dm)) %>% unlist())



# Create heatmaps -----------------------------------------------------


# Reassign dimension names, instead of the dim_xx syntax
colnames(sage_subpw) <- 1:18
colnames(sage_superpw) <- 1:18
colnames(sage_met) <- 1:18

# VAE metabolite SAGE dimension-normalized
vae_met_dim_normed <- 
  sage_met %>% 
  abs() %>% 
  scale(center = FALSE) %>% 
  plot_met()

# VAE metabolite SAGE metabolite-normalized
vae_met_met_normed <- 
  sage_met %>% 
  abs() %>% 
  t() %>% 
  as.data.frame() %>% 
  scale(center = FALSE) %>% 
  t() %>% 
  as.data.frame() %>% 
  plot_met()


# VAE sub-pathway SAGE dimension-normalized
vae_subpw_dim_normed <- 
  sage_subpw %>%
  scale(center = FALSE) %>%
  plot_subpw(., max(.))

# VAE sub-pathway SAGE subpathway-normalized
vae_subpw_pw_normed <- 
  sage_subpw %>% 
  t() %>% 
  as.data.frame() %>% 
  scale(center = FALSE) %>% 
  t() %>% 
  as.data.frame() %>% 
  plot_subpw(., max(.))


# VAE super-pathway SAGE dimension-normalized
vae_superpw_dim_normed <- 
  sage_superpw %>% 
  scale(center = FALSE) %>% 
  plot_superpw()

# VAE super-pathway SAGE subpathway-normalized
vae_superpw_pw_normed <- 
  sage_superpw %>% 
  t() %>% 
  as.data.frame() %>% 
  scale(center = FALSE) %>% 
  t() %>% 
  as.data.frame() %>% 
  plot_superpw()


colnames(pca_sage_subpw) <- 1:18

# PCA sub-pathway SAGE dimension-normalized
pca_subpw_dim_normed <- 
  pca_sage_subpw %>%
  scale(center = FALSE) %>%
  plot_subpw(., max(.))

# PCA sub-pathway SAGE subpathway-normalized
pca_subpw_pw_normed <- 
  pca_sage_subpw %>% 
  t() %>% 
  as.data.frame() %>% 
  scale(center = FALSE) %>% 
  t() %>% 
  as.data.frame() %>% 
  plot_subpw(., max(.))

colnames(cosine_sage_subpw) <- 1:18

# cosine sub-pathway SAGE dimension-normalized
cosine_subpw_dim_normed <- 
  cosine_sage_subpw %>%
  scale(center = FALSE) %>%
  plot_subpw(., max(.))

# cosine sub-pathway SAGE subpathway-normalized
cosine_subpw_pw_normed <- 
  cosine_sage_subpw %>% 
  t() %>% 
  as.data.frame() %>% 
  scale(center = FALSE) %>% 
  t() %>% 
  as.data.frame() %>% 
  plot_subpw(., max(.))

colnames(sigmoid_sage_subpw) <- 1:18

# sigmoid sub-pathway SAGE dimension-normalized
sigmoid_subpw_dim_normed <- 
  sigmoid_sage_subpw %>%
  scale(center = FALSE) %>%
  plot_subpw(., max(.))

# sigmoid sub-pathway SAGE subpathway-normalized
sigmoid_subpw_pw_normed <- 
  sigmoid_sage_subpw %>% 
  t() %>% 
  as.data.frame() %>% 
  scale(center = FALSE) %>% 
  t() %>% 
  as.data.frame() %>% 
  plot_subpw(., max(.))

colnames(rbf_sage_subpw) <- 1:18

# rbf sub-pathway SAGE dimension-normalized
rbf_subpw_dim_normed <- 
  rbf_sage_subpw %>%
  scale(center = FALSE) %>%
  plot_subpw(., max(.))

# rbf sub-pathway SAGE subpathway-normalized
rbf_subpw_pw_normed <- 
  rbf_sage_subpw %>% 
  t() %>% 
  as.data.frame() %>% 
  scale(center = FALSE) %>% 
  t() %>% 
  as.data.frame() %>% 
  plot_subpw(., max(.))

colnames(poly_sage_subpw) <- 1:18

# poly sub-pathway SAGE dimension-normalized
poly_subpw_dim_normed <- 
  poly_sage_subpw %>%
  scale(center = FALSE) %>%
  plot_subpw(., max(.))

# poly sub-pathway SAGE subpathway-normalized
poly_subpw_pw_normed <- 
  poly_sage_subpw %>% 
  t() %>% 
  as.data.frame() %>% 
  scale(center = FALSE) %>% 
  t() %>% 
  as.data.frame() %>% 
  plot_subpw(., max(.))

# Create alluvial plots for specific dimensions --------------------------------------


vae_alluvial <- 
  c("dim_12", "dim_11") %>% 
  sapply(function(dim) {
    get_alluvial_plots(dim, 
                       sage_met_wide, 
                       sage_subpw_wide,
                       sage_superpw_df,
                       n_vis = 10)
  },
  USE.NAMES = TRUE,
  simplify = FALSE)

pca_alluvial <- 
  c("dim_16", "dim_15", "dim_10") %>% 
  sapply(function(dim) {
    get_alluvial_plots(dim, 
                       pca_met_sage_wide, 
                       pca_sage_subpw_wide,
                       pca_sage_superpw_df,
                       n_vis = 10)
  },
  USE.NAMES = TRUE,
  simplify = FALSE)



# Display plots -------------------------------------------------------

# VAE
plot.new()
vae_met_dim_normed
plot.new()
vae_met_met_normed

plot.new()
vae_subpw_dim_normed
plot.new()
vae_subpw_pw_normed

plot.new()
vae_superpw_dim_normed
plot.new()
vae_superpw_pw_normed

# PCA
plot.new()
pca_subpw_dim_normed
plot.new()
pca_subpw_pw_normed

# cosine
plot.new()
cosine_subpw_dim_normed
plot.new()
cosine_subpw_pw_normed

#sigmoid
plot.new()
sigmoid_subpw_dim_normed
plot.new()
sigmoid_subpw_pw_normed

#rbf
plot.new()
rbf_subpw_dim_normed
plot.new()
rbf_subpw_pw_normed

#poly
plot.new()
poly_subpw_dim_normed
plot.new()
poly_subpw_pw_normed

# Alluvial plots
vae_alluvial
pca_alluvial
