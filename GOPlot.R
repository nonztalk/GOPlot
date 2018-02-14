
GOPlot <- function(file, plot.BP, plot.CC, plot.MF, method, width = 0.7, font.size = 10) {
  # GOPlot visualizes GO enrichment analysis results based on ggplot2
  # Input dataset should contain columns as：Category, Term, Count, PValue, FDR
  #  
  # Parameters：
  #   file：data file
  #   plot.BP：Number of terms in biological process
  #   plot.CC：Number of terms in cellular components
  #   plot.MF：Number of terms in molecular function
  #   method：Standard of significance, "FDR" or "pvalue"
  #   width：Control bar width, default 0.7
  #   font.size：Control font size of y-axis texts, default 10
  
  require(ggplot2)
  
  # Data manipulation
  GO.data <- read.table(file = file, sep = "\t", header = T, stringsAsFactors = F)
  
  GO.data$Type <- NA
  GO.data$Type[grep("BP", GO.data$Category)] <- "Biological Process"
  GO.data$Type[grep("MF", GO.data$Category)] <- "Molecular Function"
  GO.data$Type[grep("CC", GO.data$Category)] <- "Cellular Component"
  
  GO.data$Term <- gsub("^GO.*~", "", GO.data$Term)
  
  if (! is.numeric(plot.BP) | ! is.numeric(plot.CC) | ! is.numeric(plot.MF)) {
    stop("Arguments plot.BP, plot.CC, plot.MF should all be set as numeric")
  }

  if (plot.BP > nrow(GO.data[which(GO.data$Type == "Biological Process"), ])) {
    stop("Argument plot.BP is incorrectly set higher than the record number of Biological Process")
  }
  if (plot.CC > nrow(GO.data[which(GO.data$Type == "Cellular Component"), ])) {
    stop("Argument plot.CC is incorrectly set higher than the record number of Cellular Component")
  }
  if (plot.MF > nrow(GO.data[which(GO.data$Type == "Molecular Function"), ])) {
    stop("Argument plot.CC is incorrectly set higher than the record number of Molecular Function")
  }
  
  if (method != "FDR" & method != "pvalue") {
    stop("Method can only be set as FDR or pvalue")
  }
  
  if (method == "FDR") {

    GO.data <- GO.data[order(GO.data$Type, GO.data$FDR), ]
    GO.data.BP <- GO.data[which(GO.data$Type == "Biological Process"), ][1:plot.BP, ]
    GO.data.CC <- GO.data[which(GO.data$Type == "Cellular Component"), ][1:plot.CC, ]
    GO.data.MF <- GO.data[which(GO.data$Type == "Molecular Function"), ][1:plot.MF, ]
    
    if (! all(GO.data.BP$FDR <= 0.05)) {
      warning("Some chosen records of Biological Process are insignificant (FDR > 0.05) and removed")
    }
    if (! all(GO.data.CC$FDR <= 0.05)) {
      warning("Some chosen records of Cellular Component are insignificant (FDR > 0.05) and removed")
    }
    if (! all(GO.data.MF$FDR <= 0.05)) {
      warning("Some chosen records of Molecular Function are insignificant (FDR > 0.05) and removed")
    }
    
    GO.data.plot <- rbind(GO.data.BP[GO.data.BP$FDR <= 0.05, ], 
                          GO.data.CC[GO.data.CC$FDR <= 0.05, ],
                          as.data.frame(matrix(c(NA, "NA2", rep(NA, ncol(GO.data) - 4), 1, "Cellular Component"), nrow = 1, ncol = ncol(GO.data), dimnames = list("", colnames(GO.data)))),
                          GO.data.MF[GO.data.MF$FDR <= 0.05, ],
                          as.data.frame(matrix(c(NA, "NA1", rep(NA, ncol(GO.data) - 4), 1, "Molecular Function"), nrow = 1, ncol = ncol(GO.data), dimnames = list("", colnames(GO.data)))))
    
    if (all(is.na(GO.data.plot$Category))) {
      stop("There is no record significant under FDR judgement")
    }
    
    GO.data.plot$FDR <- as.numeric(GO.data.plot$FDR)
    GO.data.plot <- GO.data.plot[order(GO.data.plot$Type, -GO.data.plot$FDR), ]
    GO.data.plot$Term <- factor(GO.data.plot$Term, levels = GO.data.plot$Term)
    
    p <- ggplot(GO.data.plot, aes(Term, -log(FDR)))
    p <- p + geom_bar(aes(fill = Type), stat = "identity", width = width) + coord_flip() + theme_bw()
    p <- p + labs(x = "", y = "-ln(FDR)") + scale_fill_hue("") + theme(axis.text.y = element_text(size = font.size))

    labels = c(paste(GO.data.plot$Term[which(GO.data.plot$Type == "Biological Process")][1:plot.BP], " ", "(", GO.data.plot$Count[which(GO.data.plot$Type == "Biological Process")][1:plot.BP], ")", sep = ""),
               "",
               paste(GO.data.plot$Term[which(GO.data.plot$Type == "Cellular Component")][2:(plot.CC + 1)], " ", "(", GO.data.plot$Count[which(GO.data.plot$Type == "Cellular Component")][2:(plot.CC + 1)], ")", sep = ""),
               "",
               paste(GO.data.plot$Term[which(GO.data.plot$Type == "Molecular Function")][2:(plot.MF + 1)], " ", "(", GO.data.plot$Count[which(GO.data.plot$Type == "Molecular Function")][2:(plot.MF + 1)], ")", sep = ""))
    p <- p + scale_x_discrete(breaks = GO.data.plot$Term, labels = labels)
    p <- p + guides(fill = guide_legend(reverse = T))
    
  }
  
  if (method == "pvalue") {
    
    GO.data <- GO.data[order(GO.data$Type, GO.data$PValue), ]
    GO.data.BP <- GO.data[which(GO.data$Type == "Biological Process"), ][1:plot.BP, ]
    GO.data.CC <- GO.data[which(GO.data$Type == "Cellular Component"), ][1:plot.CC, ]
    GO.data.MF <- GO.data[which(GO.data$Type == "Molecular Function"), ][1:plot.MF, ]
    
    if (! all(GO.data.BP$PValue <= 0.05)) {
      warning("Some chosen records of Biological Process are insignificant (PValue > 0.05) and removed")
    }
    if (! all(GO.data.CC$PValue <= 0.05)) {
      warning("Some chosen records of Cellular Component are insignificant (PValue > 0.05) and removed")
    }
    if (! all(GO.data.MF$PValue <= 0.05)) {
      warning("Some chosen records of Molecular Function are insignificant (PValue > 0.05) and removed")
    }
    
    GO.data.plot <- rbind(GO.data.BP[GO.data.BP$PValue <= 0.05, ], 
                          GO.data.CC[GO.data.CC$PValue <= 0.05, ],
                          as.data.frame(matrix(c(NA, "NA2", NA, NA, 1, rep(NA, ncol(GO.data) - 6), "Cellular Component"), nrow = 1, ncol = ncol(GO.data), dimnames = list("", colnames(GO.data)))),
                          GO.data.MF[GO.data.MF$PValue <= 0.05, ],
                          as.data.frame(matrix(c(NA, "NA1", NA, NA, 1, rep(NA, ncol(GO.data) - 6), "Molecular Function"), nrow = 1, ncol = ncol(GO.data), dimnames = list("", colnames(GO.data)))))
    
    if (all(is.na(GO.data.plot$Category))) {
      stop("There is no record significant under p-value judgement")
    }
    
    GO.data.plot$PValue <- as.numeric(GO.data.plot$PValue)
    GO.data.plot <- GO.data.plot[order(GO.data.plot$Type, -GO.data.plot$PValue), ]
    GO.data.plot$Term <- factor(GO.data.plot$Term, levels = GO.data.plot$Term)
    
    p <- ggplot(GO.data.plot, aes(Term, -log(PValue)))
    p <- p + geom_bar(aes(fill = Type), stat = "identity", width = width) + coord_flip() + theme_bw()
    p <- p + labs(x = "", y = "-ln(p-Value)") + scale_fill_hue("") + theme(axis.text.y = element_text(size = font.size))
    labels = c(paste(GO.data.plot$Term[which(GO.data.plot$Type == "Biological Process")][1:plot.BP], " ", "(", GO.data.plot$Count[which(GO.data.plot$Type == "Biological Process")][1:plot.BP], ")", sep = ""),
               "",
               paste(GO.data.plot$Term[which(GO.data.plot$Type == "Cellular Component")][2:(plot.CC + 1)], " ", "(", GO.data.plot$Count[which(GO.data.plot$Type == "Cellular Component")][2:(plot.CC + 1)], ")", sep = ""),
               "",
               paste(GO.data.plot$Term[which(GO.data.plot$Type == "Molecular Function")][2:(plot.MF + 1)], " ", "(", GO.data.plot$Count[which(GO.data.plot$Type == "Molecular Function")][2:(plot.MF + 1)], ")", sep = ""))
    p <- p + scale_x_discrete(breaks = GO.data.plot$Term, labels = labels)
    p <- p + guides(fill = guide_legend(reverse = T))
    
  }
  
  return(p)

}



  
  
  
  
  

