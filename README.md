# GWAS_hmp_to_plot
#require:

## 1 load function
source("https://raw.githubusercontent.com/suzkami/GWAS_hmp_to_plot/main/GWAS_RtoP_su.R")  
source("https://raw.githubusercontent.com/suzkami/GWAS_hmp_to_plot/main/function_GWAS.R")  

## 2 input:  
(A)Genotype file (hapmap): df.hmp.txt  
(B)Phenotype file (same accession ID as hapmap file) df.txt  
(C)LD calculate by Plink or TASSEL; for candidate region narrow down.  

#example
GWAS(hmpfile = "./6xhmp50362**.hmp.txt**",  
     LDfile = "./GRDA6X50362**.ld**",   
     **phenofile** = "./Test.txt",   
     type = "MLM", ptype_filte = "./ptype.txt")  


