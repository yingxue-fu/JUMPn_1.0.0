#!/bin/bash
show_help() {
    cat <<EOF
This is the JUMP-n pipeline bootstrapping script.  Execute
with no arguments for a standard installation.  An installation of
conda or minicoda is a prerequisite.
Bootstrapping will create a conda environment in this directory for
use with JUMP-n.

EOF
}

# xcode-select --install

printf "creating conda environment $PWD/JUMPn\n"

conda create -p $PWD/JUMPn -y \
  -c bioconda \
  -c conda-forge \
  -c defaults \
  r-base=4.0.0 \
  imagemagick -y \
  r-pdftools -y \
. $(conda env list | tr -d '*' | grep -E '^base' | awk '{print $2;}')/etc/profile.d/conda.sh
conda activate $PWD/JUMPn
# conda install -c conda-forge r=4.0.0 -y
# conda install -c conda-forge imagemagick -y
# conda install -c conda-forge r-pdftools -y

Rscript bootstrap.R

#cd execution
#R -e "shiny::runApp()"

#give the link and just copy paste the link to use JUMPn