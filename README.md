# sky-islands-analysis

This repo describes the code use for the analysis used in the paper ####

Apologize in advance if the code is not polished, but I will try polsih it little by little with time. It should still run anyways and you should be able to reproduce the analyses from the beginning (FASTQ files) to the end figures and stats. 


## Bioinformatics
<img width="3240" height="2167" alt="block1" src="https://github.com/user-attachments/assets/29370f6b-721f-4da7-a692-461948ce24a4" />


The first part of the code includes part of the bioinformatics analysis. 

1) It starts from the PB_nextflow_run.slurm which runs the nextflow pipeline as provided by pacbio in their github repo (https://github.com/pacificbiosciences/HiFi-16S-workflow).
This code must be run on an HPC cluster, as it is very computationally demanding. Please note that the taxonomic and phylogeny tree were not derived by this pipeline, but were made with custom methods:

2) The code tax_ass.R assign taxonomy using RDP classifier trained on GTDB database. To make it less computationally intesive, the code divides the ASV lists in many sublists. 

3) ASV_table_generation.R changes ASVs into ASV names (ASV1 ASV2 etc) and prepares the ASV table for further analysis.

4) tree_construction_no_fastest.slurm is the code used for alignment, gap removal and tree construction.

5) tree_nodes_division.R organizes the tree in different nodes, roots the tree using Archaea and calculates which nodes to use for codiversification analysis (select for nodes in the 1/4 of the tree closer to the tips)

## Miscellaneous codes
<img width="2864" height="1443" alt="Picture1" src="https://github.com/user-attachments/assets/c59f8122-9711-4d14-8481-1a70a0fcd717" />
These codes calculate various things needed for further analysis

1) calculate NDVi values. This code download and average NDVI values
2) tree_pruning_mammals.R This code download and trims host tree
3) this code assign traits to our ASVs



## Alpha/beta Statistics
<img width="3190" height="2374" alt="Picture1" src="https://github.com/user-attachments/assets/c90f15b6-ce9d-4a57-a504-9300622c07ab" />

The following up code is code where more common alpha and beta diversit analysis are performed. 
There are four codes. 
1) alpha_beta_diversity.R calculates basic alpha and beta diversity stats and visualizations
2) boxplots_same_species_same_mountain.R boxplots of visualizations of differences in strain sharing rates between same species and same mountains.
3) strain_sharing_and_co_phylogeny.R strain sharing between species visualization 
4) modelling of ASV species occupancy

## BRMS model
<img width="3396" height="1851" alt="Picture1" src="https://github.com/user-attachments/assets/879f0e69-679f-4777-8995-93f734a36df8" />

This is the code to run BRMS model. 
1) The first step is to generate the metadata
2) then run the models in the HPC. Bewareyou will also need probably a slurm script, very simple for running the three models in the HPC/
3) THe last code is to analyze the model results 

## Codiversification

The last block is codiversification. 
<img width="2970" height="2617" alt="Picture1" src="https://github.com/user-attachments/assets/4ad657c3-b966-4dab-96a0-8b577bc7cee2" />

01_mantel_computation.R — Mantel test (host vs. microbial phylogenetic distance) for every node subtree; the heavy HPC step.
02_node_filtering_taxonomy.R — removes redundant nodes, FDR-corrects p-values, annotates taxonomy, flags dominant ASVs.
03_tree_visualization.R — cophylogeny plot and circular cladogram colored by Mantel's r.
04_trait_association_summary.R — codiversification by family/phylum, linked to GTDB life-history traits.
05_mantel_by_mountain.R — re-runs Mantel tests within each mountain range separately.
06_trait_permutation_tests.R — permutation tests for trait enrichment among codiversifying taxa.
07_codiversification_strength_abundance.R — categorizes ASVs by codiversification strength; compares abundance and occupancy.
08_genus_trait_species_occupancy.R — models host-species occupancy as a function of genus-level microbial traits and codiversification signal.








