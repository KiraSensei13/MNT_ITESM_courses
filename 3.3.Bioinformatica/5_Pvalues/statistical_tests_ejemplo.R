source("/Users/emartinez/Research/TecAgo2018/Clases/Bioinformatica/Statisitcal_Tests/statistical_tests_lib.R")

#Ejemplos
mat_example <- matrix(rnorm(100, mean=7), nrow=10)
plot.densities(mat_example)
mat_example_qn <- quantile.normalization(mat_example)
plot.densities(mat_example_qn)

de.test(x=mat_example_qn, classes=c(rep(1,5),rep(0,5)), test=c("ttest", "kolmogorov", "wilcoxon"))

##########################################################################################

#Siggenes para correr SAM

#Liga al tutorial
#https://bioconductor.org/packages/release/bioc/vignettes/siggenes/inst/doc/siggenes.pdf

#Instalar el paquete
source("http://bioconductor.org/biocLite.R")
biocLite("siggenes")

#Ejemplo
library(siggenes)
data(golub)

#Correr SAM
golub_sam <- sam(golub, golub.cl, method = "d.stat", gene.names = golub.gnames[,2])

#Vector con p-values
length(golub_sam@p.value)

#Mostrar resultados de SAM con cierto valor de delta
golub_sam_sum <- summary(golub_sam, 1.9)

#Acceder a matrix de Genes significativos 
dim(golub_sam_sum@mat.sig)