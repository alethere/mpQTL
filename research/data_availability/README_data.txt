The following folder contains:
- Genotypes: a folder containing genotypes from all NAM populations expressed in three ways:
	- In NAMi_j_biallelic.txt: SNP dosage (one column per individual, one row per marker, each cell
	  is the allelic dosage, from 0 to 4).
	- In NAMi_j_IBD.txt: founder/ancestral/IBD alleles (four columns per individual, one for each 
	  homologue of an individual, one row per marker, each cell containing a numeric allele indicating
	  a founder allele).
	- In NAMi_j_haplotype.txt: haplotype alleles (four columns per individual, one for each homologue
	  of each individual, one row per haploblock, each cell containing a string of 6 letters indicating
	  the haplotype.
- Phenotype.txt: a table containing the simulated phenotypes for the individuals of each population. Each
  column corresponds to a population. The first row corresponds to the first individual of each population, 
  the second row to the second individual, etc. 
- Potato.map: the genetic map corresponding to biallelic and IBD models.
- haplotype.map: the genetic map corresponding to haplotypes.