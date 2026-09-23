setwd("~/Documents/2025/Nocturnal_fish/paml_FreeRatio")
dNdS<-read.table("dNdS_pvalue.txt", header = T)
dNdS$FDR<-p.adjust(
  dNdS$Pvalue,
  method = "BH"
)

write.table(dNdS,
            "wilcox_results_FDR.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
summary(dNdS$Pvalue)
