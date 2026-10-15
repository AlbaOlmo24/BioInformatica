library(edgeR)
library(ggplot2)
library(pheatmap)


ruta_matriz <- path.expand("~/proyecto/trabajoR/matriz_conteos.txt")

counts_raw <- read.table(ruta_matriz, header = TRUE, sep = "\t", skip = 1, row.names = 1, check.names = FALSE)

counts_matrix <- counts_raw[, 6:ncol(counts_raw)]

colnames(counts_matrix) <- gsub(".bam$", "", basename(colnames(counts_matrix)))

cat("Se han detectado", ncol(counts_matrix), "muestras en la matriz.\n")
cat("Muestras leídas correctamente:\n")
for (muestra in colnames(counts_matrix)) {
  cat(" - ", muestra, "\n")
}

grupos <- factor(c("CYA", "CYA", "CYA", "UVA", "UVA", "UVA"))
grupos <- relevel(grupos, ref = "CYA")








y <- DGEList(counts = counts_matrix, group = grupos)

keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]









y <- calcNormFactors(y, method="TMM")

counts_normalizados <- cpm(y)
colnames(counts_normalizados) <- colnames(counts_matrix)

ruta_matriz_norm <- path.expand("~/proyecto/trabajoR/matriz_normalizada.txt")
write.table(counts_normalizados, file = ruta_matriz_norm, sep = "\t", quote = FALSE, col.names = NA)








pdf(file = "~/proyecto/trabajoR/grafico_MDS.pdf", width = 7, height = 7)
colores_mds <- c("blue", "red")[grupos]
plotMDS(y, col = colores_mds, pch = 16, main = "MDS Plot - Aspergillus flavus")
legend("topright", legend = levels(grupos), col = c("blue", "red"), pch = 16)
dev.off()

design <- model.matrix(~grupos)
rownames(design) <- colnames(counts_matrix)

y <- estimateDisp(y, design) 









fit <- glmQLFit(y, design) 

qlf <- glmTreat(fit, coef = 2, lfc = 2)

tabla_completa <- topTags(qlf, n = Inf)$table
write.csv(tabla_completa, file = "~/proyecto/trabajoR/resultados_expresion_diferencial_treat.csv")

genes_filtrados <- tabla_completa[tabla_completa$FDR <= 0.05, ]
ruta_salida_degs <- path.expand("~/proyecto/trabajoR/DEGs_FDR005_Treat.csv")
write.csv(genes_filtrados, file = ruta_salida_degs, row.names = TRUE)








pdf(file = "~/proyecto/trabajoR/grafico_Smear.pdf", width = 7, height = 7)

de_genes_tags <- decideTests(qlf, p.value = 0.05, lfc = 2)

plotSmear(qlf, de.tags = rownames(qlf)[de_genes_tags != 0], 
          col = "black", cex = 0.5,
          xlab = "Abundancia Promedio de Expresión (logCPM)",
          ylab = "Magnitud del Cambio (logFC)",
          main = "Smear Plot - Aspergillus flavus (UVA vs CYA)")

abline(h = c(-2, 2), col = "blue", lty = 2, lwd = 1.5)
dev.off()













pdf(file = "~/proyecto/trabajoR/grafico_Volcano.pdf", width = 7, height = 7)

puntos_colores <- rep("black", nrow(tabla_completa))
puntos_colores[tabla_completa$logFC >= 2 & tabla_completa$FDR <= 0.05] <- "red"    
puntos_colores[tabla_completa$logFC <= -2 & tabla_completa$FDR <= 0.05] <- "green" 

plot(tabla_completa$logFC, -log10(tabla_completa$FDR),
     col = puntos_colores, pch = 20, cex = 0.6,
     xlab = "Magnitud del Cambio (log2 Fold Change)", 
     ylab = "Significación Estadística (-log10 FDR)",
     main = "Volcano Plot - Aspergillus flavus (UVA vs CYA)")

abline(h = -log10(0.05), col = "darkgray", lty = 2)
abline(v = c(-2, 2), col = "darkgray", lty = 2)
dev.off()














if(nrow(genes_filtrados) >= 2) {
  cpm_degs <- counts_normalizados[rownames(genes_filtrados), ]
  
  pdf(file = "~/proyecto/trabajoR/grafico_Heatmap.pdf", width = 8, height = 10)
  pheatmap(cpm_degs, 
           cluster_cols = FALSE,         
           cluster_rows = TRUE,           
           scale = "row",                
           show_rownames = (nrow(genes_filtrados) <= 50), 
           main = "Heatmap DEGs - Aspergillus flavus (Z-score)",
           color = colorRampPalette(c("green", "black", "red"))(50)) 
  dev.off()
} else {
  cat("IMPORTANTE!: No se han detectado tantos DEGs necesarios como para construir el Mapa de Calor.\n")
}


cat("Número de genes diferencialmente expresados (DEGs) encontrados:", nrow(genes_filtrados), "\n")







degs_sobreexpresados <- genes_filtrados[genes_filtrados$logFC > 0, ]
ruta_sobre <- path.expand("~/proyecto/trabajoR/DEGs_sobreexpresados_UVA.csv")
write.csv(degs_sobreexpresados, file = ruta_sobre, row.names = TRUE)


degs_infraexpresados <- genes_filtrados[genes_filtrados$logFC < 0, ]
ruta_infra <- path.expand("~/proyecto/trabajoR/DEGs_infraexpresados_UVA.csv")
write.csv(degs_infraexpresados, file = ruta_infra, row.names = TRUE)
