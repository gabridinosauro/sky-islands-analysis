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

## Alpha/beta Statistics

The following up code is code where more common alpha and beta diversit analysis are performed. Followed by a code for plotting strain sharing and a then the codes for modelling strain with BRMS.

## Codiversification

The last block is codiversification. 









