# NKCE
Notebook and computational analysis workstream for **NK cell engager multi-modal data analysis (scRNA-seq+CITE-seq)**. This repository contains command lines bioinformatics software run in parallel for the quantification of spliced and un-spliced reads on bam files for the RNA velocity analysis [La Manno et al 2018](https://www.nature.com/articles/s41586-018-0414-6). Furthermore, it containts notebook to **model the quantifed reads and conduct trajectories analysis**. 

### Basic information
The analysis run by activating the [AIDA-ODS](https://kb-am1.sanofi.com/display/MP/AIDA-ODS) anconda image where the tools to run the analysis have been established. For a full list of packages please check [Sanofi confluence](https://kb-am1.sanofi.com/display/MP/AIDA-ODS) page.

### Analysis starts with exon and introns quantification on Magellan. Initialize an instance with r6i.4xlarge capabilites with an ubuntu terminal.
```
sbatch NKCE/Tools/Intron_Exon_Quantification.sh
```

### Loom are then used for trajectory analysis
```
Rscript NKCE/Script/Velocyto.R
```

