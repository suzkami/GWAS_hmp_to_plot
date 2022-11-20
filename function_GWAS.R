#Function GWAS to plot----

#Need modify:
#1. output in output_home
#2. function selection
#3. violin plot

GWAS <- function(
    input_home = getwd(),
    output_home = getwd(),
    hmp = "NA",  
    LD = "NA",
    LDGinput = 4, #if LD = "NA", LDG is needed for candidate region selection
    n = 180,      #number of individuals in genotype data set
    type = "EMMA"
){
  #setting env
  requirment(input_home, hmp, LD, output_home)
  
  #transfor input to plot_file
  print(type)
  plot_file <- readinput(type)
  print("plot_file prepared")
  
  #draw CMplot for each
  Man_QQ_plot(plot_file)
  print("Manhhatun plot done")
  print("QQ plot done")
  
  #Calculate LD deacy and save in LDA,LDB,LDG & LDD (if 6x)
  LDG <- LDdecay(n)
  print("LD decay done")
  
  #circulation for:
  for (cirn in 4:ncol(plot_file)) {
    plot_file_cir <<- plot_file[,c(1,2,3,cirn)]
    if (nrow(plot_file_cir[plot_file_cir[, 4] < 0.05/nrow(plot_file_cir), ]) > 0 ){
      # 1. Significant snp selection --> summary for output [rs, Chr, POS, Allele, AlleleA, AlleleB, AlleleC, P, traits]
      sigsnptemp <- sigsnp(plot_file_cir) #<<- for test; change to <- when finished
      print(paste0("significant snp collected: ", names(plot_file_cir)[4],"_done"))
      
      # 2. Get candidate region (Sigsnp +- LDG/ABD)
      MTArange <- MTA_range(sigsnptemp, LDG)
      print(paste0("MTA range selection: ", names(plot_file_cir)[4],"_done"))
      
      # 3. Make LD heatmap require input file
      # 4. heatmap
      LDheatmap(MTArange = MTArange)
      print(paste0("Heatmap draw: ", names(plot_file_cir)[4],"_done"))

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
}
