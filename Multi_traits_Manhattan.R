test_gwas <- read.table("./Resuls/HEADING/GWAS/TTR_KYT_HEADING_6xCC_MLM.txt", header = T, fill = T)
hmp <- read.table("./Genotype/GRDA_re150_36262_0.2_0.05_0.2_pav.hmp.txt", header = T)
plot_file <- hmp[,c(1,3,4)]
p_temp <- test_gwas[,c(1, 3, 4, 7, 15)]
p_temp <- filter(p_temp, !(is.nan(Pos)))
n <- nrow(p_temp)/nrow(plot_file)

if (n == 1) {
    plot_file[, ncol(plot_file) + 1] <- p_temp$'p'
    names(plot_file)[ncol(plot_file)] = p_temp[1,1]
  }else{
    nametemp <- p_temp[!duplicated(p_temp[,1]),1]
    for (ntraits in 1:n) {
      plot_file[, ncol(plot_file) + 1] <- p_temp$'p'[((ntraits-1)*nrow(hmp)+1):(nrow(hmp)*(ntraits))]
      names(plot_file)[ncol(plot_file)] = nametemp[ntraits]
    }
  }

colorlist <- c("#2121D9","#FF9326","#9999FF","#04B404","#0089B2","#A945FF","#B26314","#610B5E","#FE2E9A","#FFFB23","#BFF217")
#for plot
Manhattan_local <-  function(plot_file = plot_file, nametemp = nametemp, outname = nametemp[1], color = colorlist){
  #prepare the parameters
  #plot_pos, for x-aex
  #plot_file, for y and point
  if (ncol(plot_file) - 3 != length(nametemp)) {
    print("nametemp unequal to the number of your trait, please check!")
  }
  
  if (ncol(plot_file) > length(color)) {
    print("color less than the number of your trait, please check!")
  }
  
  plot_pos <- plot_file %>% group_by(chrom) %>% summarise(max(pos))
  plot_file$pos <- plot_file$pos*10/sum(plot_pos$`max(pos)`)
  plot_pos <- mutate(plot_pos, "addpos" = 0.1)
  plot_pos$`max(pos)` <- plot_pos$`max(pos)`*10/sum(plot_pos$`max(pos)`)
  
  for (i in 2:nrow(plot_pos)) {plot_pos$addpos[i] = plot_pos$addpos[i-1] + plot_pos$`max(pos)`[i-1] + 0.1}
  
  chrlist <- plot_file[!duplicated(p_temp[,2]),2]
  for (i in 1:length(chrlist)) {
    plot_file[plot_file$chrom == chrlist[i], 3] <- plot_file[plot_file$chrom == chrlist[i], 3] + as.numeric(plot_pos[plot_pos$chrom == chrlist[i],3],7)
  }
  
  plot_pos <- mutate(plot_pos, "right" = max(plot_file$pos))
  plot_pos$right[1:(nrow(plot_pos) - 1)] <- plot_pos$addpos[2:nrow(plot_pos)] - 0.1
  
  bf = -log10(0.05/nrow(plot_file))
  y_top <- ceiling(-log10(min(plot_file[,4:ncol(plot_file)])))
  
  #draw plots
  svg(str_c(outname, ".svg"), width = 18, height = 8)
  par(bg = "white", mar = c(4.1, 3.1, 3.1, 4.1))
  
  plot(1, xlab = " ", ylab = "", xaxt = "n",
       yaxt = "n", main = "", bty = "n", type = "n", 
       ylim = c(0, y_top), xlim = c(0, max(plot_file$pos)))
  
  #y axis
  rect(xleft = 0, ybottom = -0.05, xright = max(plot_file$pos) +0.5, ytop = 0, border = NA, col = "black")
  rect(xleft = (plot_pos$addpos+plot_pos$right)/2 - 0.01, ybottom = 0, xright = (plot_pos$addpos+plot_pos$right)/2 + 0.01, ytop = 0.05, border = NA, col = "black")
  mtext(text = chrlist, font = 2, side = 1, at = (plot_pos$addpos+plot_pos$right)/2, las = 1, line = 0,col ="black")
  
  #x axis
  rect(xleft = 0, ybottom = 0, xright = 0.02, ytop = y_top +0.5, border = NA, col = "black")
  rect(xleft = -0.05, ybottom = c(1:8)-0.01, xright = 0.01, ytop = c(1:8)+0.01, border = NA, col = "black")
  mtext(text = 1:8, font = 2, side = 2, at = 1:8, las = 1, line = -1,col ="black")
  
  segments(x0 = plot_pos$addpos - 0.05, y0 = 0, y1 = y_top, lty = 2, col = "gray")
  segments(x0 = 0, y0 = 5, x1 = max(plot_file$pos), lty = 2, col = "black")
  segments(x0 = 0, y0 = bf, x1 = max(plot_file$pos), lty = 1, col = "black")
  
  for (i in 1:length(nametemp)) {
    point_temp <- plot_file[,c(1:3, i+3)]
    point_temp <- filter(point_temp, -log10(point_temp[,4]) > 5)
    points(x = point_temp$pos, y = -log10(point_temp[,4]), pch = 19, cex = 0.8, col = color[i])
  }
  
  lab <- as.data.frame(cbind(nametemp, color[1:length(nametemp)]))
  legend(x = 1, y = 1, legend = lab$nametemp, col = lab$V2, pch = 19, ncol = 3, cex = 1.2)
  
  dev.off()
}

Manhattan_local(plot_file, nametemp, "HEADING_36262_MLM", color = colorlist)
