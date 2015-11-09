# plasmaMiRNATissues

This repository contains scripts and data for constructing the figure 6 for the paper:

Qin, Yidan and Yao, Jun and **Wu, Douglas C.** and Nottingham, Ryan M. and Mohr, Sabine and Hunicke-Smith, Scott and Lambowitz, Alan M. High-throughput sequencing of human plasma RNA by using thermostable group II intron reverse transcriptases. *RNA*. 2015

The directory tree are as follow:
├── README.md  
├── miRNACounts  
│   ├── D.counts  
│   ├── D_heatmap_miRNA.pdf  
│   ├── MCBZD.counts  
│   └── MCBZD_heatmap_miRNA.pdf  
├── miRNAcount.R  
└── tissuesDB  
    ├── bloodComponentsMiRNA.csv  
    ├── journal.pone.0041561.s005.csv  
    └── miRNATissue.csv  

miRNACounts/\*.counts files are raw count data.
tissuesDB contains: 

* bloodComponentsMiRNA.csv: a table contains miRNA that expressed in [blood cells](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0041561)
* journal.pone.0041561.s005.csv: a more extensive table that contains microRNA that expressed in [blood cells](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0041561)
* miRNATissue.csv: a table containing qRT-PCR-measured expressions of miRNA across [tissues](http://www.ncbi.nlm.nih.gov/pubmed/17604727)
