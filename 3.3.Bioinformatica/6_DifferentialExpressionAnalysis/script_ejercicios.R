options(stringsAsFactors=F)
library(DESeq2)
library(limma)

#Función para hacer log2 de datos continuos
glog <- function(y, lambda=1, func=log2) {
	return(func(y+sqrt(y^2+lambda)))
}

#Función para correr LIMMA
#matrix_e: matrix de expresión
#classes_e: vector de grupos como 1 y 0, la comparación se hace Grupo1 - Grupo 0
#classes_names: nombre de clases alternativos, primer elemento es para el grupo 1 y segundo elemento es para el grupo 0
limma4DS_fdr <- function(matrix_e, classes_e, classes_names=NA){
	ans_matrix <- t(apply(matrix_e, 1, function(x) c(mean(x[which(classes_e == 1)], na.rm=T), mean(x[which(classes_e == 0)], na.rm=T))))
	f <- factor(classes_e)
	design <- model.matrix(~0+f)
	colnames(design) <- c("A", "B")
	fit <- lmFit(matrix_e, design)
	contrast.matrix <- makeContrasts(B=B-A, levels=design)
	fit2 <- contrasts.fit(fit, contrast.matrix)
	fit2 <- eBayes(fit2)
	limma_mat <- toptable(fit2, number=nrow(matrix_e))
	ans_matrix <- ans_matrix[match(rownames(limma_mat), rownames(ans_matrix)),]
	ans_matrix <- cbind(ans_matrix, limma_mat[,c(1,3,4)])
	colnames(ans_matrix) <- c("B","A", "FC", "p.value", "q.value")
	if(!any(is.na(classes_names)))
		colnames(ans_matrix)[c(1,2)] <- classes_names
	ans_matrix
}

#Función para correr DESeq2
#matrix_c: matrix de conteos
#classes_c: vector de grupos como 1 y 0, la comparación se hace Grupo1 - Grupo 0
#classes_names: nombre de clases alternativos, primer elemento es para el grupo 1 y segundo elemento es para el grupo 0
DESeq_func <- function(matrix_c, classes_c, classes_names=NA){
	classes_1 <- which(classes_c == 1)
	classes_0 <- which(classes_c == 0)
	aux_data <- matrix_c[,c(classes_1,classes_0)]
	aux_desc <- data.frame(condition=c(rep("B",length(classes_1)),rep("A",length(classes_0))), type=rep("paired-end",c(length(classes_1)+length(classes_0))))
	aux_dds <- DESeqDataSetFromMatrix(countData = aux_data, colData = aux_desc, design = ~condition)
	aux_dds <- DESeq(aux_dds)
	aux_results <- as.data.frame(results(aux_dds))
	aux_results <- aux_results[order(aux_results$pvalue),]
	aux_results
}

setwd("/Users/emartinez/Research/TecAgo2018/Clases/Bioinformatica/Differential_Exp/")
load("Data_ejercicios_ejemplo.RData")

#Filtrar genes
#Conteos
#Para cada gen, contar el número de muestras con mayor a 5 conteos
count_mat_filter <- apply(count_mat, 1, function(x) length(which(x >= 5)))
#Filtrar genes con menor a 2 muestras con más de 5 conteos
count_mat <- count_mat[which(count_mat_filter >= 2),]

#FPKM
fpkm_mat <- glog(fpkm_mat)
#Para cada gen, contar el número de muestras con expresión mayor a 2
fpkm_mat_filter <- apply(fpkm_mat, 1, function(x) length(which(x >= 2)))
#Filtrar genes con menor a 2 muestras con expresión mayor a 2
fpkm_mat <- fpkm_mat[which(fpkm_mat_filter >= 2),]

#Conteos
aux_classes <- rep(1, times=ncol(count_mat))
aux_classes[grep(pattern="Control", x=colnames(count_mat))] <- 0
count_results <- DESeq_func(matrix_c=count_mat, classes_c=aux_classes)

#FPKM
aux_classes <- rep(1, times=ncol(fpkm_mat))
aux_classes[grep(pattern="Control", x=colnames(fpkm_mat))] <- 0
fpkm_results <- limma4DS_fdr(matrix_e=fpkm_mat, classes_e=aux_classes, classes_names=c("Treatment", "Control"))
