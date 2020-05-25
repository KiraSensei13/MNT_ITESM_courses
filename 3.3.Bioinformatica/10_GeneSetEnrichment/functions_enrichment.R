setwd("./")

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
