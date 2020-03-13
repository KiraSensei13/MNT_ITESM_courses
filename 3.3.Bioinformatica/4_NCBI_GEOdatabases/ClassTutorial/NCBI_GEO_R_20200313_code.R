#Slide 4
source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")
biocLite("limma")
biocLite("Biobase")
biocLite("affy")

#Slide 8
source("GEOpatch.R")
gset <- getGEO ("GSE46268" , GSEMatrix =TRUE, destdir="directorio por elegir")

#Slide 9
gset
class(gset)
names(gset)
gset <- gset[[1]]

#Slide 11
class(gset)
length(gset)
slotNames(gset)

class(pData(phenoData(gset)))
class(gset@phenoData@data)
dim(pData(phenoData(gset)))
dim(gset@phenoData@data)

#Slide 12
colnames(pData(phenoData(gset)))
pData(phenoData(gset))[, c(12:13)]

#Slide 13
dim(gset@featureData@data)
gset@featureData@data$ENTREZ_GENE_ID[1:10]
gset@featureData@data$Gene Symbol

fvarLabels(gset)
fvarLabels(gset) <- make.names(fvarLabels(gset))
fvarLabels(gset)
gset@featureData@data$Gene.Symbol[1:10]

#Slide 14
pData(phenoData(gset))[13]
class(pData(phenoData(gset))[13])
pData(phenoData(gset))[13][,1]
class(pData(phenoData(gset))[13][,1])
as.character(pData(phenoData(gset))[13][,1])
gsub(pattern="treatment: ", replacement="", as.character(pData(phenoData(gset))[13][,1]))

#Slide 15
dim(exprs(gset))
exprs(gset)[1:10,]
ex <- exprs(gset)

par(mfrow=c(1,2))
plot(density(ex))
plot(density(log2(ex+1)))

ex_log2 <- log2(ex+1)

#Slide 16
x <- rnorm(10)
t.test(x)
x <- rnorm(10, mean=10, sd=1)
t.test(x)

#Slide 17
x <- rnorm(10, mean=10, sd=1)
x
y <- rnorm(10)
y
t.test(x,y)

#Slide 19
xy_ttest <- t.test(x,y)
xy_ttest
xy_ttest$statistic
xy_ttest$p.value

#Slide 18
compare_groups <- gsub(pattern="treatment: ", replacement="", as.character(pData(phenoData(gset))[13][,1]))
compare_groups

vd_ctrl_ttest <- apply(ex_log2, 1, function(x) {
aux <- t.test(x[which(compare_groups == "1,25a-dihydroxyvitamin D3")], x[which(compare_groups == "control")]);
aux$p.value
})
length(vd_ctrl_ttest)
vd_ctrl_ttest[1:10]
sort(vd_ctrl_ttest)[1:10]
order(vd_ctrl_ttest)[1:10]
cbind(compare_groups, ex_log2[35356,])

#Slide 20
cbind(compare_groups, ex_log2[35356,])

t.test(ex_log2[35356,which(compare_groups == "control")], ex_log2[35356,which(compare_groups == "all-trans retinoic acid")])

#Slide 21
x
y
wilcox.test(x)
wilcox.test(x, y)

xy_wilcoxon <- wilcox.test(x, y)
xy_wilcoxon$statistic
xy_wilcoxon$p.value

#Slide 22
vd_ctrl_wilcoxon <- apply(ex_log2, 1, function(x) {
aux <- wilcox.test(x[which(compare_groups == "1,25a-dihydroxyvitamin D3")], x[which(compare_groups == "control")]);
aux$p.value
})
length(vd_ctrl_wilcoxon)
vd_ctrl_wilcoxon[1:10]
sort(vd_ctrl_wilcoxon)[1:10]
order(vd_ctrl_wilcoxon)[1:10]
cbind(compare_groups, ex_log2[3,])

sort(vd_ctrl_ttest)[1:5]
order(vd_ctrl_ttest)[1:5]
vd_ctrl_wilcoxon[order(vd_ctrl_ttest)[1:5]]

#Slide 23
plot(log10(vd_ctrl_ttest), log10(vd_ctrl_wilcoxon), xlab="t-test", ylab="Wilcoxon")

#Slide 24
library(psych)
data(iris)

cor(iris$Petal.Width, iris$Petal.Length)
cor.test(iris$Petal.Width, iris$Petal.Length)
cor.test(iris$Petal.Width, iris$Petal.Length)$estimate
cor.test(iris$Petal.Width, iris$Petal.Length)$p.value

plot(iris$Petal.Width, iris$Petal.Length)

#Slide 25
plot(rank(iris$Petal.Width), rank(iris$Petal.Length))
cor(iris$Petal.Width, iris$Petal.Length, method="spearman")
cor.test(iris$Petal.Width, iris$Petal.Length, method="spearman")

#Slide 26
pairs.panels(x=iris[,-5], method="pearson", ellipses=F)

#Slide 27
cor_35356 <- cor(ex_log2[35356,], t(ex_log2))
names(cor_35356) <- rownames(ex_log2)
sort(cor_35356)[1:10]
order(cor_35356)[1:10]
sort(cor_35356, decreasing=T)[1:10]
order(cor_35356, decreasing=T)[1:10]

plot(ex_log2[35356,], ex_log2[19000,])
plot(ex_log2[35356,], ex_log2[35013,])

#Slide 28
cbind(compare_groups, ex_log2[35356,], ex_log2[35013,])
cbind(compare_groups, ex_log2[35356,], ex_log2[19000,])

#Slide 29
vd_ctrl_ttest_top_1500 <- ex_log2[order(vd_ctrl_ttest)[1:1500],]
dim(vd_ctrl_ttest_top_1500)
vd_ctrl_ttest_dist <- dist(t(vd_ctrl_ttest_top_1500))
vd_ctrl_ttest_dist

#Slide 30
vd_ctrl_ttest_hh <- hclust(vd_ctrl_ttest_dist, method="average")
plot(vd_ctrl_ttest_hh, labels=compare_groups)

#Slide 31
vd_ctrl_ttest_pca <- prcomp(t(vd_ctrl_ttest_top_1500))
summary(vd_ctrl_ttest_pca)

summary(vd_ctrl_ttest_pca)$importance
summary(vd_ctrl_ttest_pca)$importance[1:2,]

#Slide 32
vd_ctrl_ttest_pca$x
predict(vd_ctrl_ttest_pca)

#Slide 33
plot(vd_ctrl_ttest_pca$x[,1], vd_ctrl_ttest_pca$x[,2], pch=19, col=as.factor(compare_groups))
as.factor(compare_groups)

#Slide 34
library(ComplexHeatmap)
library(circlize)

heatmap(vd_ctrl_ttest_top_1500)
heatmap(vd_ctrl_ttest_top_1500, margins=c(8,5))

#Slide 35
heatmap(vd_ctrl_ttest_top_1500, margins=c(8,5), labRow="")
cbind(compare_groups, colnames(vd_ctrl_ttest_top_1500))
