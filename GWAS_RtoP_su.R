#GWAS_FINAL_CODE_BY_SU
#2022_11_16

#function A: setting required
#All the results, hmp and LD file should in one same results document

requirment <- function(
    input_home = getwd(), 
    hmpfile = "NA", 
    LDfile = "NA", 
    output_home = getwd(), 
    LDGinput = 4,
    ptype_filte = "NA",
    phenofile = "NA",
    PvsG = T){
    
  #require library----
    library(tidyverse)
    library(snp.plotter)
    source("https://raw.githubusercontent.com/YinLiLin/CMplot/master/R/CMplot.r")
    
    #require input----
    input_home <<- input_home
    #hmpfile is must needed, can transform from vcf by tassel
    
    hmp <<-  read.table(paste0(input_home, "/", hmpfile), header = T, fill = T, check.names = F)
    
    if (LDfile != "NA"){
      LD <<- read.table(paste0(input_home, "/", LDfile), header = T, fill = T)
      }
   
    if (PvsG == T) {
     ifelse(phenofile == "NA", print("please input the name of phenotype or set PvsG to F"), 
            pheno <<- read.table(paste0(input_home, "/", phenofile), header = T, check.names = F))
      ifelse(ptype_filte == "NA", print("please input the name of ptype_filte or set PvsG to F"),
            ptype <<- read.table(paste0(input_home, "/", ptype_filte), header = T, check.names = F))
    }

    #require output load
    output_home <<- output_home
  }

#read_input----

readinput <- function(type = "EMMA"){
  if (type == "EMMA") {
    #from_road_list EMMAX
    flist <- list.files(path = input_home, pattern = ".ps")
    plot_file <- hmp[,c(1,3,4)]
    for (i in 1:length(flist)) {
      assign(str_sub(flist[i], 1, -4), read.table(paste0(input_home, "/", flist[i]), header = F, fill = T))
      plot_file <- cbind(plot_file, get(str_sub(flist[i], 1, -4))[,4])
      names(plot_file)[i+3] = str_sub(flist[i], 1, -4)
    }
  }else{
    #from_tassel_file
    flist <- list.files(path = input_home, pattern = type)
    plot_file <- hmp[,c(1,3,4)]
    for (i in 1:length(flist)){
      assign(str_sub(flist[i], 1, -5), read.table(paste0(input_home, "/", flist[i]), header = T, fill = T))
      if ("None" == get(str_sub(flist[i], 1, -5))[1,2]){
        assign(str_sub(flist[i], 1, -5), filter(get(str_sub(flist[i], 1, -5)), !(Chr == "0")))
        }
      
      n <- (nrow(get(str_sub(flist[i], 1, -5))))/nrow(hmp)
      
      if (n == 1) {
        plot_file[, ncol(plot_file) + 1] <- get(str_sub(flist[i], 1, -5))$'p'
        names(plot_file)[ncol(plot_file)] = get(str_sub(flist[i], 1, -5))[1,1]
      }
      else{
        nametemp <- get(str_sub(flist[i], 1, -5))[!duplicated(get(str_sub(flist[i], 1, -5))[,1]),1]
        for (ntraits in 1:n) {
          plot_file[, ncol(plot_file) + 1] <- get(str_sub(flist[i], 1, -5))$'p'[((ntraits-1)*nrow(hmp)+1):(nrow(hmp)*(ntraits))]
          names(plot_file)[ncol(plot_file)] = nametemp[ntraits]
        }
      }
    }
  }
  names(plot_file)[1:3] <- c("Marker", "Chr", "Pos")
  plot_file <- arrange(plot_file, Pos)
  plot_file <- arrange(plot_file, Chr)
  return(plot_file)
}


#From Input to plot_file----
#input plot_file; hmp
#----CMplot
Man_QQ_plot <- function(plot_file){
  CMplot(plot_file, plot.type = "m", multracks = F, LOG10 = T,
         threshold=c(0.05/nrow(plot_file),0.01/nrow(plot_file)), pch = 20,
         threshold.lty=c(1,2), threshold.lwd=c(2,2), threshold.col=c("black","black"), 
         amplify = T,
         #判断为几倍体，需要几个颜色
         if (length(plot_file[!duplicated(plot_file$"Chr"),2]) == 14){
           col = c("orange","blue")
           }else{
             col = c("orange", "blue", "darkgreen")
             },
         file="jpg",file.name="",dpi=300, width = 20, height = 6,
         file.output=TRUE, signal.col = c("red","red"), signal.cex = 1,
         main = " ",
         #highlight = Grasdisnp, highlight.cex= 0.8, highlight.pch = 16, highlight.col = "black"
  )
  
  CMplot(plot_file, plot.type="q",col="black",threshold = 0.05/nrow(plot_file),
         threshold.col="red", signal.pch=3, signal.cex=1,
         signal.col="deeppink",conf.int=T,conf.int.col="orange")
}

#LDdecay; return LDA, LDB, LDD, LDALL ----
LDdecay <- function(n){
  if (typeof(LD) != "character") {

    LDtemp <- LD
  
    LDtemp[,ncol(LDtemp)+1] <- LDtemp$BP_B - LDtemp$BP_A
    names(LDtemp)[ncol(LDtemp)] <- "Dist_bp"
    LDtemp[,ncol(LDtemp)+1] <- LDtemp$Dist_bp/1000000
    names(LDtemp)[ncol(LDtemp)] <- "Dist_mb"
    
    distance <- LDtemp$Dist_mb
    LDtemp.data <- LDtemp$R2
    n <- ncol(hmp) - 11 #number of individuals in genotype data set, for LDdecay calculate
    
    HW.st<-c(C=0)
    HW.nonlinear<-nls(LDtemp.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=1000))
    tt<-summary(HW.nonlinear)
    new.rho<-tt$parameters[1]
    fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))
    
    LD_output <- mutate(LDtemp, decay = fpoints)
    LD_output$decay <- abs(LD_output$decay -0.2)
    
    LDG0.2 <- arrange(LD_output, decay)$Dist_mb[1]
    return(LDG0.2)
  }
  else{
    LDG <- LDGinput
    return(LDG)
  }
}


#Significant SNPs ----
Allele_number <- function(sigsnp_hmp){
  Allele1 <- rowSums(sigsnp_hmp[,12:ncol(sigsnp_hmp)] == str_sub(sigsnp_hmp$alleles, 1, 1))
  Allele2 <- rowSums(sigsnp_hmp[,12:ncol(sigsnp_hmp)] == str_sub(sigsnp_hmp$alleles, 3, 3))
  Allele3 <- rowSums(sigsnp_hmp[,12:ncol(sigsnp_hmp)] == str_sub(sigsnp_hmp$alleles, 5, 5))
  AlleleT <- Allele1+Allele2+Allele3
  Allele_num <- data.frame("Allele1" = Allele1, "Allele2" = Allele2, "Allele3" = Allele3, "AlleleT" = AlleleT)
  return(Allele_num)
}
sigsnp <- function(plot_file){
    sigsnp_temp <- plot_file[plot_file[, 4] < 0.05/nrow(plot_file), ]
    sigsnp_hmp <- hmp[hmp$rs %in% sigsnp_temp$Marker, ]
    allele_infor <- Allele_number(sigsnp_hmp)
    sigsnp_hmp[,6:9] <- allele_infor
    sigsnp_hmp[,10] <- sigsnp_temp[,4]
    sigsnp_hmp[,11] <- names(plot_file)[4]
    names(sigsnp_hmp)[6:11] <- c("Allele1", "Allele2", "Allele3", "AlleleT", "p", "traits")
    return(sigsnp_hmp)
}

#select candidate region within LD decay block ----
MTA_range <- function(sigsnp_hmp, LDG){
  #先分染色体，再分染色体块
  chr <- sigsnp_hmp$chrom[!duplicated(sigsnp_hmp$chrom)]
  traitname <- names(plot_file_cir)[4]
  
  for (i in 1:length(chr)) {
    snp_chr <- sigsnp_hmp[sigsnp_hmp$chrom == chr[i],]
    index1 <- 1 #times of change
    index2 <- 1 #times of no change
    if (nrow(snp_chr) == 1) {
      assign(paste0(traitname, chr[i], "MTA_", index1), snp_chr[1,])
    }else if (snp_chr$pos[nrow(snp_chr)] - snp_chr$pos[1] < LDG*1000000){
      assign(paste0(traitname, chr[i], "MTA_", index1), snp_chr)
    }else{
      for (n in 2:(nrow(snp_chr))+1){
        if (n <= nrow(snp_chr)) {
        
          pos_dis <- snp_chr$pos[n] - snp_chr$pos[n-1]
          if (pos_dis > LDG*1000000) {
            assign(paste0(traitname, chr[i], "MTA_", index1), snp_chr[(n-index2):(n-1)])
            index1 <- index1 + 1
            index2 <- 1
          }else{
            index2 <- index2 + 1 
          }
          
        }else{
          pos_dis <- snp_chr$pos[n-1] - snp_chr$pos[n-2]
          if (pos_dis > LDG*1000000) {
            assign(paste0(traitname, chr[i], "MTA_", index1), snp_chr[(n-index2):(n-1)])
          }else{
            assign(paste0(traitname, chr[i], "MTA_", index1), snp_chr[(n-index2-1):(n-1)])
          }
      }
    }
  }
  
  MTA_1_trait <- ls(pattern = traitname)
  list_temp <- list()
  
  for (i in 1:length(MTA_1_trait)) {
    snp <- get(MTA_1_trait[i])[1:10]
    MTA <- filter(hmp, chrom == snp[1,3], pos >= snp[1,4] - LDG*10e5, pos <=  snp[1,4] + LDG*10e5)
    
    for (n in 1:nrow(snp)) {
      loop <- filter(hmp, chrom == snp[n,3], pos >= snp[n,4] - LDG*10e5, pos <=  snp[n,4] + LDG*10e5)
    }
    
    MTA <- rbind(MTA, loop)
    MTA <- MTA[!duplicated(MTA[,1]),]
    MTA$panelLSID <- plot_file_cir[plot_file_cir$Marker %in% MTA$rs, 4]
    MTA$QCcode <- names(plot_file_cir)[4]
    MTA[,6:9] <- Allele_number(MTA)
  
    list_temp[[i]] <- MTA
    names(list_temp)[i] <- MTA_1_trait[i]
  }
  return(list_temp)

}

#heatmap snp_plotter----
LDheatmap <- function(MTArange){
  #text file need for snp plotter
  color <- read.table(paste0(input_home, "/", "palette.txt"), header = F)
  config_fomat <- read.table(paste0(input_home, "/", "config.txt"), header = F, sep = "\n")
  config_in <- data.frame(label = c("SNP.FILE", "COLOR.LIST", "GENOTYPE.FILE", "IMAGE.TITLE", "IMAGE.NAME"), connect = "=", value = "NA")
  
  for (i in 1:length(MTArange)) {
    dir.create(paste0(output_home, "/", names(MTArange)[i]))
    setwd(paste0(output_home, "/", names(MTArange)[i]))
    write.table(color, "palette.txt", row.names = F, col.names = F, sep = "\n", quote = F)
    
    chr_color <- ifelse(str_sub(MTArange[[i]][1,3], 2, 2) == "A", "orange", ifelse(str_sub(MTArange[[i]][1,3], 2, 2) == "B", "blue", "darkgreen"))
    config_in$value <- c(paste0(names(MTArange)[i], "_ss.txt"), chr_color, paste0(names(MTArange)[i], "_geno.txt"), paste0("\" \""), paste0(names(MTArange)[i]))
    config_out <- config_fomat
    config_out[c(1,3,4,5,6),] <- str_c(config_in[,1], config_in[,2], config_in[,3])
    write.table(config_out, "config.txt", row.names = F, col.names = F, sep = "\n", quote = F)
    
    ss <- MTArange[[i]][,c(1,2,4,10)]
    names(ss) <- c("ASSOC", "SNP.NAME",	"LOC",	"SS.PVAL")
    ss$ASSOC <- "+"
    ss$SNP.NAME <- ifelse(ss$SS.PVAL < 0.05/nrow(plot_file_cir), "SIGSNP", "*")
    ss$LOC <- as.integer(ss$LOC)
    ss$SS.PVAL <- as.numeric(ss$SS.PVAL)
    write.table(ss, paste0(names(MTArange)[i], "_ss.txt"), row.names = F, col.names = T, sep = "\t", quote = F)
    
    geno <- data.frame()
    geno[1:(ncol(MTArange[[i]])-11),1:6] <- "0"
    alleleinfor <- MTArange[[i]]
    numericaallele <- data.frame()
    n <- 1
    for (a in 1:nrow(alleleinfor)) {
      for (j in 12:ncol(alleleinfor)) {
        if (alleleinfor[a,j] == str_sub(alleleinfor[a,2],1,1)) {
          numericaallele[j-11, n:(n+1)] = "1"
        }
        else if (alleleinfor[a,j] == str_sub(alleleinfor[a,2],3,3)) {
          numericaallele[j-11, n:(n+1)] = "2"
        }
        else if (alleleinfor[a,j] == "N") {
          numericaallele[j-11, n:(n+1)] = "0"
        }
        else{
          numericaallele[j-11, n] = "1"
          numericaallele[j-11, (n+1)] = "2"
        }
      }
      n <- n + 2
    }
    geno <- cbind(geno, numericaallele)
    write.table(geno, paste0(names(MTArange)[i], "_geno.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
    
    #for each block, LD heatmap   
    snp.plotter(config.file = "config.txt")
    setwd(paste0(output_home))
    print(paste0("LD heatmap for ", names(MTArange)[i], " done"))
  }
  setwd(paste0(input_home))
}


#violion plot / barplot----
#phenotype_file + siglist
#color depend on alleles

#Phenotype vs best SNP allele of each peak
#1 get phenotype exist in siglist
#2 judege discrete or countine
#3 need group infor or not
#4 prepare for violin or bar plot

GenovsPheno <- function(MTArange, pheno, ptype){
  #Thics functiuon can run slef by using the output filed from function_GWAS: 
    #1. siglist for all the sigsnps
    #2. phenotype data used for GWAS
    #3. phenotype data type (eg. discrete or quantity)
  
#----plot required functions, should in preparation step  
  #function inside GenovsPheno, for generate SNPs lsit
  Allele_plot <- function(geno_temp){
    allele <- data.frame(str_split(geno_temp$alleles, "/"))
    allele <- allele[order(allele[,1]),]
    if (geno_temp$assayLSID != (ncol(geno_temp) - 11)) { allele = c(allele, "N") }
  }
  
  #function for plots
  savesvg <- function(plot, name){
    svg(paste0(name, ".svg"))
    print(plot)
    dev.off()
  }
  
  #violin
  violin_plots <- function(plot_temp = plot_temp, allele_color = allele_color, SNPs = SNPs, geno_temp = geno_temp){
    range <- max(plot_temp$Phenotype, na.rm = T) - min(plot_temp$Phenotype, na.rm = T)  
    min_value <- min(plot_temp$Phenotype, na.rm = T) - range*0.15
    max_value <- max(plot_temp$Phenotype, na.rm = T) + range*0.2
      
    p <- ggplot(plot_temp, aes(x = Genotype, y = Phenotype, fill = Genotype))+
            geom_violin(width = 1)+
            scale_x_discrete(limits = SNPs)+
            geom_jitter(height = 0, width = 0.1, alpha=0.8)+
            geom_boxplot(width=0.1, color="red", alpha=0.2, size = 1)+
            scale_fill_manual(values = allele_color)+
            scale_y_continuous(limits = c(min_value, max_value),
                               breaks = ceiling(seq(min_value, max_value, 
                               by = round((max_value - min_value)/5, 2))))+
            theme_bw()+
            labs(x="Aelles", y="Phenotype", title = paste0(geno_temp$QCcode, "_", geno_temp$chrom))+
            theme(legend.position ="none", 
                  panel.background = element_rect(fill = "white", colour = "grey50"),
                  axis.text = element_text(colour = "black", size = rel(2)),
                  axis.title.x = element_text(size = 20),
                  axis.title.y = element_text(size = 20),
                  plot.title = element_text(size = 20))
    
    savesvg(plot = p, name = paste0("violinplot", geno_temp$QCcode, "_", geno_temp$chrom, geno_temp$pos))
  }
  
  #bar plot
  bar_plots <- function(plot_temp = plot_temp, trait_color = trait_color, geno_temp = geno_temp){
    temp <- plot_temp 
    temp <- temp %>% group_by(Genotype, Phenotype) %>% summarise(n = n()) %>% as.data.frame()
    temp_non <- filter(temp, Genotype != 'N')
    temp_n <- filter(temp, Genotype == 'N')
    temp_non <- temp_non[order(temp_non[,1]),]
    
    temp <- rbind(temp_non, temp_n)
    
    p <- ggplot(temp, aes(x = Genotype, y = n, fill = factor(Phenotype)))+
      geom_bar(stat="identity", position = position_stack(), width = 0.6, )+
      scale_fill_manual(values = trait_color)+
      theme_minimal() +
      theme_classic() +
      theme_bw()+
      labs(y = "Numbers",title = " ") +
      theme(legend.position ="bottom", legend.box = "horizontal",
            legend.title = element_blank(),
            legend.text = element_text(size = 18),
            panel.background = element_rect(fill = "white", colour = "grey50"),
            axis.text = element_text(colour = "black", size = rel(1.5)),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 18),
            plot.title = element_text(size = 20),
            text = element_text(family = "Arial"))
      
    savesvg(plot = p, name = paste0("barpolot", geno_temp$QCcode, "_", geno_temp$chrom, geno_temp$pos))
  }
#----  

  #main body
  #input MTArange, pheno, ptype; output plot_temp, allele_color, SNPs, geno_temp
  slist <- MTArange[[1]][1,]
  
  names(ptype) <- c("ID", "type")

  for (i in 1:length(MTArange)){
    ptemp <- MTArange[[i]]
    slist[i,] <- arrange(ptemp, panelLSID)[1,]
  }

  for (i in 1:nrow(slist)) {
    pheno_temp <- pheno[,c(1, which(names(pheno) == slist$QCcode[i]))]
    geno_temp <- slist[i,]
    geno_pheno_temp <- data.frame(t(geno_temp))
    geno_pheno_temp[,2] <- rownames(geno_pheno_temp)
    names(geno_pheno_temp) <- c("Allele","<Trait>")

    plot_temp <- merge(pheno_temp, geno_pheno_temp, by = "<Trait>")
    names(plot_temp) <- c("ID", "Phenotype", "Genotype")
    
    SNPs <- Allele_plot(geno_temp)
    allele_color <- c("A" = "#2121D9", "G" = "#3399FF", "N" = "gray", "C" = "#A945FF", "T" = "#FE2E9A", "+" = "#FF9326", "-" = "#04B404")
    trait_color <- c("#2121D9", "#A945FF", "#FE2E9A", "#3399FF", "#FF9326", "gray")
    #allele_color <- c("A" = "#FF3333", "G" = "#FF7F00", "N" = "gray", "C" = "#4DAF4A", "T" = "#3399FF", "+" = "#FF7F00", "-" = "#4DAF4A")
    #trait_color <- c("#00CC00", "#FF0000", "#FF8000", "#3399FF", "#CCOOCC", "gray")
    #discrete or continue
    if (ptype[geno_temp$QCcode == ptype$ID, 2] != "discrete") {
      violin_plots(plot_temp, allele_color, SNPs, geno_temp)
    }else{
      bar_plots(plot_temp = plot_temp, trait_color = trait_color, geno_temp = geno_temp)
    }
  }
}


















