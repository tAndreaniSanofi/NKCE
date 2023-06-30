# NKCE
Notebook and computational analysis workstream for **NK cell engager multi-modal data analysis (scRNA-seq+CITE-seq)**. Here, we report on three main tools that have been applied to perfor trajectory and pseudo-time inference, to estimate the RNA velocity and the “regulon”, which depict the cellular differentiation and molecular processes underlying NKCE effector function.

## Trajectory and pseudo-time inference
We propose to use [Monocle3](http://cole-trapnell-lab.github.io/monocle3/) for trajectory and pseudo-time inference. Monocle3 is an analysis toolkit for single-cell RNA-Seq experiments. It uses an algorithm to learn the sequence of gene expression changes each cell must go through as part of a dynamic biological process. Once it has learned the overall "trajectory" of gene expression changes, Monocle3 can place each cell at its proper position in the trajectory.

### Run the analysis
To run the analysis, intall the libraries 
```{r}
devtools::install('monocle3')
BiocManager::install(c('org.Hs.eg.db', 'org.Mm.eg.db'))
devtools::install_github("cole-trapnell-lab/garnett", ref="monocle3")
```
and run the script *monocle3.R* on your data. You will need to update the file with the path to the scRNA-seq count matrix and, cell and genes metadata matrices.

## RNA velocity
The mRNA velocity is estimated based on [La Manno et al 2018](https://www.nature.com/articles/s41586-018-0414-6) on spliced and un-spliced reads for a group of genes that can explain the differentiation of a cell type to another cell type such as CD56dim to CIML. Thus, this will allow to **model the quantifed reads and conduct trajectories analysis**.  

### Basic information
The analysis run by activating the [AIDA-ODS](https://kb-am1.sanofi.com/display/MP/AIDA-ODS) anconda image where the tools to run the analysis have been established. For a full list of packages please check [Sanofi confluence](https://kb-am1.sanofi.com/display/MP/AIDA-ODS) page.

### Analysis starts with exon and introns quantification on Magellan. Initialize an instance with r6i.4xlarge capabilites with an ubuntu terminal and a SLURM cluster with 9 nodes.
```
sbatch NKCE/Tools/Intron_Exon_Quantification.sh
```

### Loom files are then used for trajectory analysis and plot
```
Rscript NKCE/Script/RNA_Velocity.r sample.loom marker_genes.txt plot.pdf
```

## Regulon activity inference
Regulon activity inference is performed using the python implementation of [SCENIC](https://www.nature.com/articles/nmeth.4463). The tool allows to simultaneously reconstruct gene regulatory networks and identify stable cell states from single-cell RNA-seq data. 

### Run the analysis
To run the analysis you will need to install the pySCENIC package on your machine. We reccomend to install it on a conda environment 
```
conda create -y -n pyscenic python=3.10
conda activate pyscenic
pip install pyscenic
```
You will need to download motif rankings files and motif annotation file from https://resources.aertslab.org/cistarget/databases/. Update the *pySCENIC.py* file with your working directory information as well as the directory to your input data and the desired output folders, and run it on python.


