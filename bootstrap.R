#/// JUMPn Dependencies Installation \\\#
#################################################################
### Install BiocManager for Bioconductor package installation ###
#################################################################
if (!requireNamespace("BiocManager", quietly = TRUE))         ###
  install.packages("BiocManager",repos="http://cran.us.r-project.org")              ###
#################################################################
#################################################################

################################################
#### Install BiocManager-dependent packages ####
################################################
BiocManager::install("ReactomePA", dependencies=TRUE)          ####
BiocManager::install("org.Hs.eg.db", dependencies=TRUE)        ####
BiocManager::install("AnnotationDbi", dependencies=TRUE)       ####
BiocManager::install("annotate", dependencies=TRUE)            ####
BiocManager::install("WGCNA", dependencies=TRUE)               ####
BiocManager::install("msigdbr", dependencies=TRUE)             ####
################################################
################################################

#######################################################################
###### Install CRAN-dependent packages via vector simultaneously ######
#######################################################################
pkgs <- c('shiny','gplots', 'ggplot2','DT','tidyverse','Rcpp',                        ######
          'plyr','readr','writexl','readxl',                                          ######
          'rvest','httr','igraph','visNetwork',                                       ######
          'stringr','shinyjs','ggnewscale','org.Mm.eg.db',                            ######
          'enrichplot', 'shinyBS','shinydashboard','shinyWidgets')                    ######
install.packages(pkgs, repos="http://cran.us.r-project.org", dependencies=TRUE)       ######
#######################################################################
#######################################################################
