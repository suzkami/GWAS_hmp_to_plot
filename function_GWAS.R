#Function GWAS to plot----

#Need modify:

#2. function selection
#3. violin plot
#4. output list of sigsnps

GWAS <- function(
    input_home = getwd(),
    output_home = getwd(),
    hmpfile = "NA",  
    LDfile = "NA",
    LDGinput = 4, #if LD = "NA", LDG is needed for candidate region selection
    type = "EMMA",
    ptype_filte = "NA",
    phenofile = "NA",
    PvsG = T
){
  #setting env
  requirment(input_home = input_home, hmpfile = hmpfile, LDfile = LDfile, output_home = output_home, 
             ptype_filte = ptype_filte, phenofile = phenofile, PvsG = PvsG)
  
  #transfor input to plot_file
  print(type)
  plot_file <- readinput(type)
  print("plot_file prepared")
  
  #draw CMplot for each
  setwd(paste0(output_home, "/"))
  Man_QQ_plot(plot_file)
  setwd(paste0(input_home, "/"))
  print("Manhhatun plot done")
  print("QQ plot done")
  
  #Calculate LD decay and save in LDA,LDB,LDG & LDD (if 6x)
  LDG <- LDdecay(n)
  print("LD decay done")
  
  #circulation for:
  
  #prepare for siglist output
  siglist <- hmp[1,]
  siglist[1,] <- names(siglist)
  names(siglist)[6:11] <- c("Allele1", "Allele2", "Allele3", "AlleleT", "p", "traits")
  
  #LD heatmap
  for (cirn in 4:ncol(plot_file)) {
    plot_file_cir <<- plot_file[,c(1,2,3,cirn)]
    if (nrow(plot_file_cir[plot_file_cir[, 4] < 0.05/nrow(plot_file_cir), ]) > 0 ){
      # 1. Significant snp selection --> summary for output [rs, Chr, POS, Allele, AlleleA, AlleleB, AlleleC, P, traits]
      sigsnptemp <- sigsnp(plot_file_cir) #<<- for test; change to <- when finished
      siglist <- rbind(siglist, sigsnptemp)
      print(paste0("significant snp collected: ", names(plot_file_cir)[4],"_done"))
      
      # 2. Get candidate region (Sigsnp +- LDG/ABD)
      MTArange <- MTA_range(sigsnptemp, LDG)
      print(paste0("MTA range selection: ", names(plot_file_cir)[4],"_done"))
      
      # 3. Make LD heatmap require input file
      # 4. heatmap
      LDheatmap(MTArange = MTArange)
      print(paste0("Heatmap draw: ", names(plot_file_cir)[4],"_done"))
      
      # 5. pheno vs geno
      if (PvsG == T) {GenovsPheno(MTArange, pheno, ptype)}
      print(paste0("PvsG plot: ", names(plot_file_cir)[4],"_done"))

    #end condition----
      if (cirn == ncol(plot_file) ) {
        print(paste0("Job done, totall ", cirn-3, " traits"))
      }
      
    }else{
      cirn = cirn + 1
      if (cirn > ncol(plot_file) ) {
       print(paste0("Job done, totall ", cirn-4, " traits"))
      }
    }
  }

  #output snp list----
  names(siglist)[6:11] <- c("Allele1", "Allele2", "Allele3", "AlleleT", "pvalue", "traits")
  write.table(siglist, "siglist.txt", row.names = F, col.names = F, sep = "\t", quote = F)
  
  #Phenotype vs best SNP allele of each peak
  #1 get phenotype exist in siglist
  #2 judege discrete or countine  
  #3 need group infor or not
  #4 prepare for violin or bar plot
   
}
