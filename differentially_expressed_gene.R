################################################################################
#Usage: screen differentially expressed genes
library("affy") #use affy package to do background correcting & Normalizing & Calculating expression
library("limma") #use limma package to find differentially expressed genes
data=ReadAffy(celfile.path="~/../../..") #read .CEL files into data
eset=rma(data) #use rma to normalize data
matrix=as.matrix(eset) #convert S4 file to matrix
write.csv(matrix,"~/../../xxx.csv") #write expression data
design=model.matrix(~ 0+factor(c(1,1,1,1,0,0,0,0)))
colnames(design)=c("group1","group2")
fit=lmFit(eset,design)
contrast.matrix=makeContrasts(group2-group1,levels=design)
fit2=contrasts.fit(fit,contrast.matrix)
fit2=eBayes(fit2)
result=topTable(fit2,number=60000,adjust.method="BH",p.value=1,lfc=0)


