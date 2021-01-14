# mtVAE

The repository contains scripts to replicate findings in the paper **Gomari et al., "_Variational autoencoders learn universal latent representations of metabolomics data"_**

<br>

### Requirments:

- [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html)
- [R version >= 4.0.0](https://www.r-project.org/)
- [RStudio](https://rstudio.com/products/rstudio/)

<br>

### 1. git repository

Clone local copy of git repository

`git clone https://github.com/krumsieklab/mtVAE`

(or use a git GUI client of your choice)

<br>

### 2. Environment setup

##### Setup **python** environment

In a terminal, switch to the directory of the local git repository.

`conda env create --force --file environment.yml`
<br>
<br>
##### **R** setup

1. Open `mtvae.Rproj`
2. Run `R_setup.R`

Note: this has been tested under **R** version 4.0.0 and **RStudio** version 1.3.1073
<br>
<br>
<br>

### 3. Activating the conda environment to access jupyter notebooks
`conda activate mtvae_env`

`jupyter notebook`

<br>
<br>

### 4. Instructions for running the scripts

1. Place datasets into [data/](https://gitlab.com/krumsieklab/parviz/mtvae_dev/-/tree/master/data).
2. With access to TwinsUK, Type 2 diabetes (T2D), schizophrenia, acute myeloid leukemia (AML) data, run scripts in increasing order of the file prefixes, starting from **01**_train_VAE.ipynb. 

* **Notes:** 
    * Pre-trained models from **01**_train_VAE.ipynb can be found under [models/](https://gitlab.com/krumsieklab/parviz/mtvae_dev/-/tree/master/models)
    * All R scripts should be run from within RStudio
    * **00**_optimize_VAE_hyperparameters.ipynb can be skipped and should only be used as a guide to select hyperparameters.



<br>

| Name | Description |
| ------ | ------ |
| **00**_optimize_VAE_hyperparameters.ipynb | Optimize for VAE hyperparameters using TwinsUK train data. (Runtime: 1h15m on a MacBook pro)|
| **01**_train_VAE.ipynb | Train VAE model on TwinsUK data and calculate evaluation metrics. Note that this requires access to TwinsUK, which should be requested separately from https://twinsuk.ac.uk/. |
| **02**_reconstruct_data.ipynb | Generate TwinsUK data reconstructions using trained VAE and PCA models. Used for model performance assessments. |
| **03**_assess_model_performance.R | Compute mean squared error (MSE) and correlation matrix MSE (CM-MSE) for VAE and PCA. This includes the calculation of MSE and CM-MSE for varying latent space dimensionality _d_. |
| **04**_calculate_SAGE_values_VAE.ipynb | Calculate VAE SAGE values using TwinsUK test data. This script should be parallelized, due to its long runtime. Pre-computed VAE SAGE values can be found under [results/sage_values](https://gitlab.com/krumsieklab/parviz/mtvae_dev/-/tree/master/results/sage_values). (Runtime: if all instances are parallelized ~7.5h) |
| **04**_calculate_SAGE_values_PCA.ipynb | Calculate PCA SAGE values using TwinsUK test data. Pre-computed PCA SAGE values can be found under [results/sage_values](https://gitlab.com/krumsieklab/parviz/mtvae_dev/-/tree/master/results/sage_values). (Runtime: if all instances are parallelized ~1.5h) |
| **05**_interpret_latent_space.R| Create SAGE value heatmaps and alluvial plots for VAE and PCA. |
| **06**_encode_data.ipynb | Generate type 2 diabetes (T2D), schizophrenia, and acute myeloid leukemia (AML) data encodings using VAE and PCA models. |
| **07**_associate_dimensions_with_diseases.R| Associate VAE and PCA encodings with patient groups from T2D, schizophrenia, and AML data. This includes T2D clinical variables (e.g. HbA1c %) and AML mutations. |


<br>

#### Other files

| Name | Description |
| ------ | ------ |
| models.py | Contains VAE and PCA model classes. |
| metric_functions.py | Functions used for model assessment in python can be found here.  |
| helper_functions.R | R functions that are required for the calculation of evaluation results and the construction of plots can be found here. |

