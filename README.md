# Visual macro-synteny with VINGO 

                                              Shingo Miyauchi 3Aug21
# Descriptions

Welcome to our visual-omics platform! Visually Integrated Numerous Genres of Omics (VINGO) is a set of custom R scripts that combine and visualise output from our visual integrative omics platform (a.k.a ShingoTools) including SHIN+GO, TINGO, PRINGO, and SynGO. It produces multi-layered plots integrating macrosynteny, secretome, transposable elements, gene clusters for secondary metabolites, and gene expression. 

*Please note that VINGO itself does not produce bioinformatic results. It is a visualisation tool for complex omics data. 

**NOTE 1: Synteny detection**
Identifying syntenic blocks among the species was performed with R package DECIPHER (http://www2.decipher.codes/; see the parameters in Looney et al., 2021). This process is excluded as it is computationally intensive. A high performance computing cluster may be required for calculations. Setting up R environments on the computing cluster is out of scope in this demo.  

**NOTE 2: Genomic features**
We performed prediction of secretome and identification of TEs using in-house pipelines at INRAE Nancy, France (Pellegrin et al. 2015; Morin et al., 2019). Then, clean combined output files were made with PRINGO and TINGO pipelines (Miyauchi et al., 2020). 

# Requirement for R and packages
R Studio,
R 3.6.3

dplyr 1.0.2
karyoploteR 1.12.4
scales 1.1.0

# How to run this demo
1) Download INPUT folder
2) Start up R Studio by clicking an R script
3) Make sure your working directory is in INPUT folder so that R recognises input files. 
4) Read instructions in the scripts and execute the code (line by line).  
5) Figures are generated as output in the same folder. 

# Input files - VINGO_INPUT.zip

1) Gene coordinates with genes, TEs in assembly scaffold 1 to 10 from SynGO

SynGO_Scaff1to10_GenomeFeature.csv

2) Scaffold size calculated from SynGO

SynGO_Scaffolds_Size.csv

3) Synteny identified from SynGO

SynGO_synteny_locations.csv

4) Scaffold of interest

Glocon1_compared_with_7fungi.txt

5) Genes coding for CAZymes in genomes from the CAZy team (www.cazy.org)

CAZymes_CAZy_Total.csv
 
# Output files (for busy people)

VINGO_OUTPUT.zip

# Description of output - Kirisame (drizzle) plots

The script generates a genome-wide plot showing genome assembly scaffolds, genes, transposable elements, syntenic blocks. You can get an eagle view of macrosynteny with some species. The figure made here was used in Looney et al., (2021).

**NOTE 3:** Raw figures generated may not be aesthetic enough and it may require adjustments and beautification with Adobe Illustrator.

# References
1. Looney B, Miyauchi S, Morin E, Drula E, Courty PE, Kohler A, Lindquist E, Kuo A, LaButti K, Pangilinan J, et al. 2021. Evolutionary priming and transition to the ectomycorrhizal habit in an iconic lineage of mushroom-forming fungi: is preadaptation a requirement? bioRxiv: 2021.02.23.432530.
2. Miyauchi S, Kiss E, Kuo A, Drula E, Kohler A, Sánchez-García M, Morin E, Andreopoulos B, Barry KW, Bonito G, et al. 2020. Large-scale genome sequencing of mycorrhizal fungi provides insights into the early evolution of symbiotic traits. Nature Communications 11: 1–17.
3. Morin E, Miyauchi S, San Clemente H, Chen EC, Pelin A, de la Providencia I, Ndikumana S, Beaudet D, Hainaut M, Drula E, et al. 2019. Comparative genomics of Rhizophagus irregularis, R. cerebriforme, R. diaphanus and Gigaspora rosea highlights specific genetic features in Glomeromycotina. New Phytologist 222: 1584–1598.
4. Pellegrin C, Morin E, Martin FM, Veneault-Fourrey C. 2015. Comparative analysis of secretomes from ectomycorrhizal fungi with an emphasis on small-secreted proteins. Frontiers in Microbiology 6.
