GWAS(hmpfile = "./6xhmp50362.hmp.txt",
     LDfile = "./GRDA6X50362.ld", 
     phenofile = "./Mira_GWAS.txt", 
     type = "MLM", ptype_filte = "./ptype.txt")

setwd("../50362/")

source("https://raw.githubusercontent.com/suzkami/GWAS_hmp_to_plot/main/GWAS_RtoP_su.R")
source("https://raw.githubusercontent.com/suzkami/GWAS_hmp_to_plot/main/function_GWAS.R")