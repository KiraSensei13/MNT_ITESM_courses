options(stringsAsFactors=F)
library(limma)
setwd("/Users/emartinez/Research/TecAgo2018/Clases/Bioinformatica/GSE_Networks/Ejercicio/")

#Functions

#Función para hacer la prueba hipergeométrica
p.overlap.gene_set <- function(total_genes, geneset_genes, significant_genes, lower.tail=FALSE) {
	total_genes_in_set <- length(which(total_genes %in% geneset_genes))
	total_genes_not_in_set <- length(total_genes) - total_genes_in_set
	significant_genes_in_set <- length(which(significant_genes %in% geneset_genes))
	c(phyper(significant_genes_in_set, total_genes_in_set, total_genes_not_in_set, length(significant_genes), lower.tail=lower.tail), significant_genes_in_set)
}

#Función para hacer la prueba hipergeométrica en una lista de gene sets
gene_set_stats <- function(total_genes, gene_set_list, significant_genes){
	gs_stats <- c()
	for(i in 1:length(gene_set_list)){
		aux_gs <- gene_set_list[[i]]
		aux_gs <- aux_gs[which(aux_gs %in% total_genes)]
		gs_stats <- rbind(gs_stats, p.overlap.gene_set(total_genes, aux_gs, significant_genes))
	}
	rownames(gs_stats) <- names(gene_set_list)
	gs_stats <- cbind(gs_stats, p.adjust(gs_stats[,1]))
	colnames(gs_stats) <- c("p.value", "Num", "FDR")
	gs_stats <- gs_stats[-which(gs_stats[,"Num"] < 2),]
	gs_stats <- gs_stats[,c(1,3,2)]
	gs_stats <- gs_stats[order(gs_stats[,1], -gs_stats[,3]),]
	gs_stats
}

#Función para correr limma
limma4DS_fdr <- function(matrix_e, classes_e, classes_names=NA){
	ans_matrix <- t(apply(matrix_e, 1, function(x) c(mean(x[which(classes_e == 1)], na.rm=T), mean(x[which(classes_e == 0)], na.rm=T))))
	f <- factor(classes_e)
	design <- model.matrix(~0+f)
	colnames(design) <- c("Resistant", "Sensitive")	
	fit <- lmFit(matrix_e, design)
	contrast.matrix <- makeContrasts(Sensitive=Sensitive-Resistant, levels=design)
	fit2 <- contrasts.fit(fit, contrast.matrix)
	fit2 <- eBayes(fit2)
	limma_mat <- toptable(fit2, number=nrow(matrix_e))
	ans_matrix <- ans_matrix[match(rownames(limma_mat), rownames(ans_matrix)),]
	ans_matrix <- cbind(ans_matrix, limma_mat[,c(1,3,4)])
	colnames(ans_matrix) <- c("Sensitive","Resistant", "FC", "p.value", "q.value")
	if(!any(is.na(classes_names)))
		colnames(ans_matrix)[c(1,2)] <- classes_names
	ans_matrix
}

#Cargar datos 3 objetos: lgg_datos (datos de expresión), lgg_subtypes (clases de las muestras) y onco_gs (lista de gene sets)
load("TCGA_LGG_data_GSEA.RData")

#Limma: comparar respecto la clase IDHmut-non-codel
aux_classes <- ifelse(lgg_subtypes == "IDHmut-non-codel", 1, 0)
lgg_limma <- limma4DS_fdr(lgg_datos, aux_classes, c("IDHmut-non-codel", "IDHmut-codel"))

#Obtener genes up y down regulated
lgg_up_genes <- rownames(subset(lgg_limma, FC > 0 & q.value < 1e-06))
lgg_dw_genes <- rownames(subset(lgg_limma, FC < 0 & q.value < 1e-06))

#Enrichment analysis usando hipergeométrica
lgg_up_gs <- gene_set_stats(rownames(lgg_datos), onco_gs, lgg_up_genes)
lgg_dw_gs <- gene_set_stats(rownames(lgg_datos), onco_gs, lgg_dw_genes)

#GSEA
#Ordenar genes por signo FC * -log10(p-value)
lgg_limma$Cor <- sign(lgg_limma$FC) * -log10(lgg_limma$p.value)
lgg_cor <- data.frame(Gene=rownames(lgg_limma), Cor=lgg_limma$Cor)
lgg_cor <- lgg_cor[order(lgg_cor$Cor, decreasing=T),]
#Guardar en archivo de texto
write.table(lgg_cor, sep="\t", quote=F, row.names=F, col.names=F, file="LGG_subtypes.rnk")

#Java
rnk_files <- dir(pattern="*.rnk")
rnk_names <- gsub(pattern=".rnk", replacement="", rnk_files)
gsea_jar <- "gsea2-2.2.4.jar"
#Crear directorio
aux_dir <- "GSEAscores/"
aux_gmt <- "c6.all.v7.1.symbols.gmt"
nplot <- 30

for(i in 1:length(rnk_names)){
	java_command <- paste("java -cp", gsea_jar, "-Xmx3072m xtools.gsea.GseaPreranked")
	java_command <- paste(java_command, "-gmx", aux_gmt, "-collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk")
	java_command <- paste(java_command, rnk_files[i], "-scoring_scheme weighted -rpt_label", rnk_names[i], "-include_only_symbols true")
	java_command <- paste(java_command, "-make_sets true -plot_top_x", nplot, "-rnd_seed timestamp -set_max 500 -set_min 5 -zip_report false -out")
	java_command <- paste(java_command,  aux_dir, "-gui false")
	system(java_command)
}
