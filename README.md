# NKCE
Notebook and computational analysis workstream for **NK cell engager multi-modal data analysis (scRNA-seq+CITE-seq)**. The mRNA velocity is estimated based on [La Manno et al 2018](https://www.nature.com/articles/s41586-018-0414-6) on spliced and un-spliced reads for a group of genes that can explain whether the cell is in a differentiation state. Thus, this will allow to **model the quantifed reads and conduct trajectories analysis**.  

### Basic information
The analysis run by activating the [AIDA-ODS](https://kb-am1.sanofi.com/display/MP/AIDA-ODS) anconda image where the tools to run the analysis have been established. For a full list of packages please check [Sanofi confluence](https://kb-am1.sanofi.com/display/MP/AIDA-ODS) page.

### Analysis starts with exon and introns quantification on Magellan. Initialize an instance with r6i.4xlarge capabilites with an ubuntu terminal.
```
sbatch NKCE/Tools/Intron_Exon_Quantification.sh
```

### Loom files are then used for trajectory analysis and plot
```
Rscript NKCE/Script/Velocyto.R sample.loom marker_genes.txt plot.pdf
```

