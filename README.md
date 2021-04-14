# JUMPn_1.0.0

#############
### JUMPn ###
#############


##############
## Overview ##
##############

JUMPn is a bioinformatic software developed as a part of the JUMP Proteomic Pipeline for the purpose of undertaking co-expression, protein-protein interaction, network, and pathway analysis of proteomic/genomic data sets 

JUMPn implements the widely utilized WGCNA algorithm to identify cohorts of genes/proteins which exhibit similar expression/abundance/quantitative profiles across multiple replicates (the similarity statistic employed here is the Pearson Correlation). 

Once co-expression clusters are identified, the program overlays all of the known protein-protein interactions over the cluster to create a network. The program identifies first level interaction modules of this PPI network by searching for disconnected components. 

After identifying first level modules, if any first level module exceeds the minimum threshold specified by the user (the min_TOM_module_size) a second round of modularization will be undertaken. The second level modules are identified by topological overlap mattrix (TOM) analysis, which identifies dense regions of interactions as displayed on an adjacency matrix. 

After identifying all of the modules, JUMPn constructs visual network displays of each cluster and  each module.  

After network construction,JUMPn searches for enriched pathways at both the co-expression cluster level as well as the modular level. Modularization provably yields more significantly enriched pathways, as evidenced by p-value and Benjiamini and Hochberg FDR. 

##########################
### JUMPn Installation ###
##########################

1. JUMPn insallation requires installation of conda on the user's local machine. If conda has not been previously installed, we recommend using Anaconda or Miniconda 


    a. Anaconda: https://docs.anaconda.com/anaconda/install/ 
    
    b. Miniconda https://docs.conda.io/en/latest/miniconda.html
    
    
2. Download the JUMPn source code from GitHub 
3. Double click to unzip the downloaded file JUMPn_Shiny.zip; a new folder    named JUMPn_Shiny will be created.
4.	Open command line terminal. For example, in MacOS:


    a.	Press two buttons: command button + space button
    
    b.	Type terminal
    
    c.	Press button return
    
    
5.	Create Conda environment for JUMPn. 


    a.	Change the current directory to the JUMPn_Shiny folder: on the             terminal, type cd /path/to/JUMPn_shiny (note: you need to know the         absolute path to the JUMPn_Shiny folder)
    
    b.	Create and activate an empty Conda environment. On the terminal, type the following commands:
    
              conda create -p $PWD/JUMPn -y
              
              conda activate $PWD/JUMPn
              
              
6.	Install JUMPn dependencies. 


    a.	Install R: on the terminal, type:
    
              conda install -c conda-forge r=4.0.0 -y
              
    b.	Install dependency packages: on the terminal, type: 
    
              Rscript bootstrap.R
              
              
              
7.	Launch JUMPn on web browser.


    a.	Change the current directory to the execution folder: on the terminal, type 
    
                 cd execution
                 
    b.	Launch JUMPn: on the terminal, type: 
    
                R -e "shiny::runApp()"
                
    c.	Once the above is executed, the terminal screen will show up Listening on http://127.0.0.1:XXXX (here XXXX indicates 4 random numbers). Copy and paste               http://127.0.0.1:XXXX onto the web browser, on which JUMPn welcome page will show up
    
    
    
 

