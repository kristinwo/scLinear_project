library(renv)
# Set up a python environment. I will use the renv package to create a new local environment. In case you want to use your default environment and r library, the first two commands can be scipped.
renv::init("./")
renv::use_python()

install.packages("devtools")
install.packages("reticulate")

# Get sure that the correct python environment is used by reticulate
reticulate::py_config()

#> python:         /path/to/renv/python/virtualenvs/renv-python-3.10/bin/python
#> libpython:      /usr/lib/python3.10/config-3.10-x86_64-linux-gnu/libpython3.10.so
#> pythonhome:     /path/to/renv/python/virtualenvs/renv-python-3.10:/path/to/renv/python/virtualenvs/renv-python-3.10
#> version:        3.10.12 (main, Nov 20 2023, 15:14:05) [GCC 11.4.0]
#> numpy:           [NOT FOUND]
#> 
#> NOTE: Python version was forced by RETICULATE_PYTHON


# Install all python dependencies (The version specified for scikit-learn & numpy should be yoused, if you want to use the pretrained models)
py_modules_list <- c("numpy<1.26.0","joblib","scikit-learn==1.2.0","anndata","warnings","scanpy")
for (m in py_modules_list){
  if(!reticulate::py_module_available(m)){reticulate::py_install(m)}
}
if(!reticulate::py_module_available("pytorch_lightning")){reticulate::py_install("pytorch_lightning", pip = TRUE)} # install with pip


# Test if all python dependencies are available
py_modules_list_available <- c("numpy","joblib","sklearn","anndata","warnings","torch","scanpy","os","scipy","typing", "pytorch_lightning")
py_modules_list_not_available <- c()
for (m in py_modules_list_available){
  if(!reticulate::py_module_available(m)){py_modules_list_not_available <- c(py_modules_list_not_available,m)}
}
if(is.null(py_modules_list_not_available)){ "All python modules are available"}else{
  print(paste0("The following python modules are not available: ", paste0(py_modules_list_not_available, collapse = ", ")))
  print("Try installing them with: reticulate::py_install(missing module name)")
}

#> [1] "All python modules are available"

# if the installation of scLinear fails, try restarting the R session, to reload the reticulate environment

# Install scLinear

if (!require("devtools", quietly = TRUE)){install.packages("devtools")}

devtools::install_github("DanHanh/scLinear")
