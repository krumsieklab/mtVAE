source("helper_functions.R")

# Get distribution of clinical parameter p-values at finer granularity

# Load encodings ----------------------

vae_reps_path <- "results/replications"
pca_path <- 'results/encodings'

vae_reps_files <-
  list.files(vae_reps_path,
             recursive = TRUE,
             full.names = TRUE)

d_dims <- str_extract(vae_reps_files, "\\d{2}[_]") %>% 
  unique() %>% 
  strsplit("_") %>% 
  unlist()

models<- c("PCA", "RBF KPCA", "Sigmoid KPCA", "Polynomial KPCA", "Cosine KPCA")


# Load PCA top performing dimensions --------------------------------------------------
# NOTE: This code will only run if you have created encodings 
# for the PCA and KPCA models at d = [10,13,15,16,17,18,20]

top_performers_pca<-
  d_dims %>% lapply(function(dim){
    qm_encoding <- 
      bind_rows(
        glue("{pca_path}//QMDiab_PCA_encoding_d{dim}.csv") %>% 
          read_csv() %>% 
          process_qmdiab() %>% 
          dplyr::mutate(model = "PCA"),
        glue("{pca_path}//QMDiab_KPCA_cosine_encoding_d{dim}.csv") %>% 
          read_csv() %>% 
          process_qmdiab() %>% 
          dplyr::mutate(model = "Cosine KPCA"),
        glue("{pca_path}//QMDiab_KPCA_sigmoid_encoding_d{dim}.csv") %>% 
          read_csv() %>% 
          process_qmdiab() %>% 
          dplyr::mutate(model = "Sigmoid KPCA"),
        glue("{pca_path}//QMDiab_KPCA_rbf_encoding_d{dim}.csv") %>% 
          read_csv() %>% 
          process_qmdiab() %>% 
          dplyr::mutate(model = "RBF KPCA"),
        glue("{pca_path}//QMDiab_KPCA_poly_encoding_d{dim}.csv") %>% 
          read_csv() %>% 
          process_qmdiab() %>% 
          dplyr::mutate(model = "Polynomial KPCA")
        
      )
    qm_scores <- qm_encoding %>% get_latent_scores()
    qm_mins <- models %>% lapply(function(mod){log10(min(qm_scores %>% 
                                                           filter(model==mod) %>% .$pval_nominal))})
    qm_df <- tibble(model = models,
                    dim = rep(dim, times = length(models)),
                    disease = rep("qm", times= length(models)),
                    mins = qm_mins %>% unlist() )
    
    aml_encoding <- 
      bind_rows(
        glue("{pca_path}//AML_PCA_encoding_d{dim}.csv") %>% 
          read_csv() %>% 
          process_aml("^[0-9]+") %>% 
          dplyr::mutate(model = "PCA"),
        glue("{pca_path}//AML_KPCA_cosine_encoding_d{dim}.csv") %>% 
          read_csv() %>% 
          process_aml("^[0-9]+") %>% 
          dplyr::mutate(model = "Cosine KPCA"),
        glue("{pca_path}//AML_KPCA_sigmoid_encoding_d{dim}.csv") %>% 
          read_csv() %>% 
          process_aml("^[0-9]+") %>% 
          dplyr::mutate(model = "Sigmoid KPCA"),
        glue("{pca_path}//AML_KPCA_RBF_encoding_d{dim}.csv") %>% 
          read_csv() %>% 
          process_aml("^[0-9]+") %>% 
          dplyr::mutate(model = "RBF KPCA"),
        glue("{pca_path}//AML_KPCA_poly_encoding_d{dim}.csv") %>% 
          read_csv() %>% 
          process_aml("^[0-9]+") %>% 
          dplyr::mutate(model = "Polynomial KPCA")
      )
    aml_scores <- aml_encoding %>% get_latent_scores()
    aml_mins <- models %>% lapply(function(mod){log10(min(aml_scores %>% 
                                                            filter(model==mod) %>% .$pval_nominal))})
    aml_df <- tibble(model = models, 
                     dim = rep(dim, times = length(models)),
                     disease = rep("aml", times= length(models)),
                     mins = aml_mins %>% unlist() )
    
    schizo_encoding <- 
      bind_rows(
        glue("{pca_path}//Schizo_PCA_encoding_d{dim}.csv") %>% 
          read_csv() %>% 
          process_schizo("^[0-9]+") %>% 
          dplyr::mutate(model = "PCA"),
        glue("{pca_path}//Schizo_KPCA_cosine_encoding_d{dim}.csv") %>% 
          read_csv() %>% 
          process_schizo("^[0-9]+") %>% 
          dplyr::mutate(model = "Cosine KPCA"),
        glue("{pca_path}//Schizo_KPCA_sigmoid_encoding_d{dim}.csv") %>% 
          read_csv() %>% 
          process_schizo("^[0-9]+") %>% 
          dplyr::mutate(model = "Sigmoid KPCA"),
        glue("{pca_path}//Schizo_KPCA_rbf_encoding_d{dim}.csv") %>% 
          read_csv() %>% 
          process_schizo("^[0-9]+") %>% 
          dplyr::mutate(model = "RBF KPCA"),
        glue("{pca_path}//Schizo_KPCA_poly_encoding_d{dim}.csv") %>% 
          read_csv() %>% 
          process_schizo("^[0-9]+") %>% 
          dplyr::mutate(model = "Polynomial KPCA")
      )
    schizo_scores <- schizo_encoding %>% get_latent_scores() 
    schizo_mins <- models %>% lapply(function(mod){log10(min(schizo_scores %>% 
                                                               filter(model==mod) %>% .$pval_nominal))})
    schizo_df <-tibble(model = models,
                       dim = rep(dim, times = length(models)),
                       disease = rep("schizo", times= length(models)),
                       mins = schizo_mins %>% unlist() )
    
    bind_rows(qm_df, aml_df, schizo_df)
  }) %>% bind_rows()

# Load VAE top performing dimensions --------------------------------------------------
# NOTE: This code will only run if you have trained a VAE model
# at d = [10,13,15,16,17,18,20] over multiple replications

reps <-  str_extract(vae_reps_files, "[0-9]{1,3}.csv") %>% unique()

top_performers_vae<-
  d_dims %>% lapply(function(dim){
    reps %>% lapply(function(rep){
    qm_encoding <- 
      glue("{vae_reps_path}/QMDiab_VAE_encoding_{dim}_{rep}") %>% 
      read_csv() %>% 
      process_qmdiab() %>% 
      dplyr::mutate(model = "VAE")
    qm_scores <- qm_encoding %>% get_latent_scores()
    qm_min <- min(qm_scores %>% .$pval_nominal)
    
    aml_encoding <- 
      glue("{vae_reps_path}/AML_VAE_encoding_{dim}_{rep}") %>% 
      read_csv() %>% 
      process_aml("^[0-9]+") %>% 
      dplyr::mutate(model = "VAE")
    aml_scores <- aml_encoding %>% get_latent_scores()
    aml_min <- min(aml_scores %>% .$pval_nominal)
    
    schizo_encoding <- 
      glue("{vae_reps_path}/Schizo_VAE_encoding_{dim}_{rep}") %>% 
      read_csv() %>% 
      process_schizo("^[0-9]+") %>% 
      dplyr::mutate(model = "VAE")
    schizo_scores <- schizo_encoding %>% get_latent_scores()
    schizo_min <- min(schizo_scores %>% .$pval_nominal)
    rep_df <- tibble(rep = rep(rep,3),
                     model = rep("VAE", 3),
                     dim = rep(dim, 3),
                     disease = c("qm", "aml", "schizo"),
                     vae_mins = c(qm_min, aml_min, schizo_min))
    
  }) %>% bind_rows()
    }) %>% bind_rows()


# Calculate error bars and plot --------------------------------------------------

df<- top_performers_vae %>% group_by(dim,disease) %>%
  dplyr::summarise(mean = mean(log10(vae_mins)),
            sd = sd(log10(vae_mins)),
            min = mean-sd,
            max = mean+sd) 
  
  ggplot(df,aes(x = dim, y =mean))+
    geom_line(aes(group=1))+geom_errorbar(aes(ymin=min, ymax = max),width = 2, color = "black")+
    geom_line(data=top_performers_pca,aes(x=dim, y=mins, group = model,color=model))+facet_wrap(~disease,scale="free")+
      theme_minimal() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.x = element_blank(),
            panel.spacing = unit(2, "lines"),
            strip.text.y = element_blank()) +
      scale_color_manual(values=brewer.pal(6, name = "Dark2"))
    