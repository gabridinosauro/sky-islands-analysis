# sky-islands-analysis

## Bioinformatics
<img width="3240" height="2167" alt="block1" src="https://github.com/user-attachments/assets/29370f6b-721f-4da7-a692-461948ce24a4" />


The first part of the code includes part of the bioinformatics analysis. 

It starts from the PB_nextflow_run.slurm which runs the nextflow pipeline as provided by pacbio in their github repo (https://github.com/pacificbiosciences/HiFi-16S-workflow).

This code must be run on an HPC cluster, as it is very computationally demanding.

Please note that the taxonomic and phylogeny tree were not derived by this pipeline, but were made with custom methods:

The code tax_ass.R assign taxonomy using RDP classifier trained on GTDB database. To make it less computationally intesive, the code divides the ASV lists in many sublists. 

ASV_table_generation.R changes ASVs into ASV names (ASV1 ASV2 etc) and prepares the ASV table for further analysis.

tree_construction_no_fastest.slurm is the code used for alignment, gap removal and tree construction.

tree_nodes_division.R organizes the tree in different nodes, roots the tree using Archaea and calculates which nodes to use for codiversification analysis (distal 1/4 of the tree) 






