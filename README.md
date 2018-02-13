# Psuedo-LD calculations


> This repository contains a set of scripts that calculate the correlation of genotypes across the genome.  There are two primary programs LD.cpp and LD_chrom.cpp.  The former calculates the average correlations for a range of distances across the genome.  The latter calculates correlations between sites of a given distance in windows across chromosomes.  The other scripts are wrappers that implement this on a cluster computing system.  It is currently designed for a specific cluster, but can be generalized somewhat easily.  

> It is important to note that the genotypic correlations do not strictly measure linkage disequilibrium, although they will be proportional to traditional measures of LD.  

> The basic input is a table of genotypes that is produced by GATK's VariantToTable.  The LD_formatter.py script attempts to parse this table by population and create the two input files (per population) used as input for the main programs (LD.cpp and LD_chrom.cpp).  The 'Wrapper' files do the automation of job submission to the cluster system at the John Innes Centre.

