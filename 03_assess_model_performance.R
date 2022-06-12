source("helper_functions.R")


# Set bootstrap hyperparameters -------------------------------------------

# number of bootstrap samples
n_boot <- 1000

# number of cores for bootstrapping
n_cores <- 5


# Load data ---------------------------------------------------------------

twins_path  <- "data/TwinsUK.xls"

# Load from files
twins_train_data <- read_xls(twins_path,  "Training Set", col_types = "numeric")
twins_test_data  <- read_xls(twins_path,  "Testing Set", col_types = "numeric")

# Load precalculated reconstructions
recon_path <- "results/reconstructions"


# Get reconstruction scores -------------------------------------

# delete precomputed results/reconstruction_scores.csv file to run this part
# NOTE: deleting this file will create a version that potentially does NOT
# contain varying d reconstruction scores. If you want to recreate the
# figure with correlation matrix MSE, please make sure you do
# run 02_reconstruct_data.ipynb in a loop with 
# d = [5, 10, 15, 18, 20, 30, 40, 60, 80, 100, 120, 160, 200]


if(!file.exists("results/reconstruction_scores.csv")) {
  
  twins_recon_files <-
    list.files(recon_path,
               pattern = "Twins",
               recursive = TRUE,
               full.names = TRUE)
  
  d_dims <- str_extract(twins_recon_files, "_d_[0-9]*") %>% unique()
  
  twins_scores <-
    d_dims %>%
    lapply(function(dim_idx) {
      
      print(glue("Working on {dim_idx}"))
      train_recon_files <-
        twins_recon_files %>%
        str_subset(glue("{dim_idx}.csv")) %>%
        str_subset("train")
      
      test_recon_files <-
        twins_recon_files %>%
        str_subset(glue("{dim_idx}.csv")) %>%
        str_subset("test")
      
      # Note: Need to use _PCA because just "PCA" will match the KPCA files as well as the target "PCA" file

      train_scores <-
        get_bootstrap_scores(train_recon_files %>% str_subset("VAE") %>% read_csv(),
                             train_recon_files %>% str_subset("_PCA") %>% read_csv(),
                             twins_train_data,
                             n_times = n_boot,
                             n_cores = n_cores) %>%
        dplyr::mutate(data = "Twins Train",
                      d_dim = str_extract(dim_idx, "[0-9]+") %>% as.numeric())
      
      test_scores <-
        get_bootstrap_scores(test_recon_files %>% str_subset("VAE") %>% read_csv(),
                             test_recon_files %>% str_subset("_PCA") %>% read_csv(),
                             twins_test_data,
                             n_times = n_boot,
                             n_cores = n_cores) %>%
        dplyr::mutate(data = "Twins Test",
                      d_dim = str_extract(dim_idx, "[0-9]+") %>% as.numeric())
      
      bind_rows(train_scores, test_scores)
      
    }) %>%
    bind_rows()
  
  twins_scores %>% write_csv("results/reconstruction_scores.csv")
  
}


# Create plot object, different MSEs for d=18 --------------------------------------

all_metrics <- read_csv("results/reconstruction_scores.csv")

# plot score boxplots

all_metrics_processed <- 
  all_metrics %>% 
  # NOTE: I'm only looking at Twins reconstruction here
  filter(str_detect(data, "Twins")) %>% 
  pivot_longer(cols = c("mse_cormat", 
                        "mse")) %>% 
  dplyr::mutate(
    label = case_when(
      name == "mse_cormat" ~ "MSE (CorMat)",
      name == "mse" ~ "MSE"
    ),
    data = fct_relevel(data, "Twins Train", "Twins Test" )
  )


# Get MSE plots

mse_plot_train <- get_mse_plot("Twins Train", "Training set", all_metrics_processed, 2)
mse_plot_test <- get_mse_plot("Twins Test", "Test set",all_metrics_processed, 2)

# Get correlation matrix MSE (CM-MSE) plots
mse_cormat_plot <- 
  all_metrics_processed %>% 
  filter(label == "MSE (CorMat)",
         d_dim == 18) %>% 
  ggplot(aes(x = "", y = value, color = model)) + 
  geom_boxplot() +
  labs(y = "Correlation Matrix MSE") +
  theme(axis.title = element_blank()) +
  facet_grid(cols = vars(label),
             rows = vars(data),
             scales = "free") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        panel.spacing = unit(1, "lines")) +
  scale_color_manual(values=wes_palette(n=2, name="BottleRocket1"))


# Create plot object, correlation matrix MSE over d --------------------------------------

# NOTE: This only work if you have run 02_reconstruct_data.ipynb
# in a loop with d = [5, 10, 15, 18, 20, 30, 40, 60, 80, 100, 120, 160, 200]

mse_cormat_range_plot <- 
  all_metrics_processed %>% 
  filter(name == "mse_cormat") %>% 
  filter(label == "MSE (CorMat)") %>% 
  dplyr::select(model, seed, data, d_dim, name, value, label) %>% 
  group_by(model, data, d_dim, label) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            min = mean-sd,
            max = mean+sd) %>%
  ggplot(aes(x = d_dim, y = mean, color = model, ymin=min, ymax=max)) +
  geom_line() +
  geom_errorbar(width = 2, color = "black") +
  geom_point(size = 0.5) +
  facet_grid(rows = vars(data),
             cols = vars(label)) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(2, "lines"),
        strip.text.y = element_blank()) +
  scale_color_manual(values=wes_palette(n=2, name="BottleRocket1"))


# display plot objects -------------------------------------------------------

# MSE
mse_plot_train
mse_plot_test

# CM-MSE
mse_cormat_plot
mse_cormat_range_plot


# KPCA extended figure plots

if(!file.exists("results/reconstruction_scores_kpca.csv")) {
  
  twins_recon_files <-
    list.files(recon_path,
               pattern = "Twins",
               recursive = TRUE,
               full.names = TRUE)
  
  d_dims <- str_extract(twins_recon_files, "_d_[0-9]*") %>% unique()
  
  
  kpca_twins_scores <-
    d_dims %>%
    lapply(function(dim_idx) {
      
      print(glue("Working on {dim_idx}"))
      train_recon_files <-
        twins_recon_files %>%
        str_subset(glue("{dim_idx}.csv")) %>%
        str_subset("train")
      
      test_recon_files <-
        twins_recon_files %>%
        str_subset(glue("{dim_idx}.csv")) %>%
        str_subset("test")
      
      train_scores <-
        get_bootstrap_scores_kpca(train_recon_files %>% str_subset("poly_KPC") %>% read_csv(),
                             train_recon_files %>% str_subset("cosine_KPC") %>% read_csv(),
                             train_recon_files %>% str_subset("sigmoid_KPC") %>% read_csv(),
                             train_recon_files %>% str_subset("rbf_KPC") %>% read_csv(),
                             twins_train_data,
                             n_times = n_boot,
                             n_cores = n_cores) %>%
        dplyr::mutate(data = "Twins Train",
                      d_dim = str_extract(dim_idx, "[0-9]+") %>% as.numeric())
      
      test_scores <-
        get_bootstrap_scores_kpca(test_recon_files %>% str_subset("poly_KPC") %>% read_csv(),
                             test_recon_files %>% str_subset("cosine_KPC") %>% read_csv(),
                             test_recon_files %>% str_subset("sigmoid_KPC") %>% read_csv(),
                             test_recon_files %>% str_subset("rbf_KPC") %>% read_csv(),
                             twins_test_data,
                             n_times = n_boot,
                             n_cores = n_cores) %>%
        dplyr::mutate(data = "Twins Test",
                      d_dim = str_extract(dim_idx, "[0-9]+") %>% as.numeric())
      
      bind_rows(train_scores, test_scores)
      
    }) %>%
    bind_rows()
  
  kpca_twins_scores %>% write_csv("results/reconstruction_scores_kcpa.csv")
}

kpca_metrics <- read_csv("results/reconstruction_scores_kcpa.csv")

kpca_metrics_processed <- 
  kpca_metrics %>% 
  # NOTE: I'm only looking at Twins reconstruction here
  filter(str_detect(data, "Twins")) %>% 
  pivot_longer(cols = c("mse_cormat", 
                        "mse")) %>% 
  dplyr::mutate(
    label = case_when(
      name == "mse_cormat" ~ "MSE (CorMat)",
      name == "mse" ~ "MSE"
    ),
    data = fct_relevel(data, "Twins Train", "Twins Test" )
  )
  
# Get KPCA MSE plots

kpca_mse_plot_train <- get_mse_plot("Twins Train", "Training set", kpca_metrics_processed, 4)
kpca_mse_plot_test <- get_mse_plot("Twins Test", "Test set",kpca_metrics_processed, 4)

# Get correlation matrix MSE (CM-MSE) plots
kpca_mse_cormat_plot <- 
  kpca_metrics_processed %>% 
  filter(label == "MSE (CorMat)",
         d_dim == 18) %>% 
  ggplot(aes(x = "", y = value, color = model)) + 
  geom_boxplot() +
  labs(y = "Correlation Matrix MSE") +
  theme(axis.title = element_blank()) +
  facet_grid(cols = vars(label),
             rows = vars(data),
             scales = "free") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        panel.spacing = unit(1, "lines")) +
  scale_color_manual(values=wes_palette(4, name = "Zissou1", type = "continuous"))


# Create plot object, correlation matrix MSE over d --------------------------------------

# NOTE: This only work if you have run 02_reconstruct_data.ipynb
# in a loop with d = [5, 10, 15, 18, 20, 30, 40, 60, 80, 100, 120, 160, 200]

kpca_mse_cormat_range_plot_train <- 
  kpca_metrics_processed %>% 
  filter(name == "mse_cormat") %>% 
  filter(label == "MSE (CorMat)") %>% 
  filter(data=="Twins Train") %>% 
  dplyr::select(model, seed, data, d_dim, name, value, label) %>% 
  group_by(model, data, d_dim, label) %>%
  dplyr::summarise(mean = mean(value),
            sd = sd(value),
            min = mean-sd,
            max = mean+sd) %>%
  ggplot(aes(x = d_dim, y = log10(mean), color = model, ymin=log10(min), ymax=log10(max))) +
  geom_line() +
  geom_errorbar(width = 2, color = "black") +
  geom_point(size = 0.5) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(2, "lines"),
        strip.text.y = element_blank()) +
  scale_color_manual(values=wes_palette(4, name = "Zissou1", type = "continuous"))

kpca_mse_cormat_range_plot_test <- 
  kpca_metrics_processed %>% 
  filter(name == "mse_cormat") %>% 
  filter(label == "MSE (CorMat)") %>% 
  filter(data=="Twins Test") %>% 
  dplyr::select(model, seed, data, d_dim, name, value, label) %>% 
  group_by(model, data, d_dim, label) %>%
  dplyr::summarise(mean = mean(value),
            sd = sd(value),
            min = mean-sd,
            max = mean+sd) %>%
  ggplot(aes(x = d_dim, y =mean, color = model, ymin=min, ymax=max)) +
  geom_line() +
  geom_errorbar(width = 2, color = "black") +
  geom_point(size = 0.5) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(2, "lines"),
        strip.text.y = element_blank()) +
  scale_color_manual(values=wes_palette(4, name = "Zissou1", type = "continuous"))


# display plot objects -------------------------------------------------------

# KPCA MSE
kpca_mse_plot_train
kpca_mse_plot_test

# KPCA CM-MSE
kpca_mse_cormat_plot
kpca_mse_cormat_range_plot_train
kpca_mse_cormat_range_plot_test

