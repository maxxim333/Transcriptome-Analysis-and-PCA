#Set directory to file location
#setwd("C:/Users/Maksym/Desktop/UPM/Genomic Data Analysis/Tema 5 Practice 2")
mydata <- read.csv("table.csv", header = TRUE, sep = ",", dec = ".", row.names = 1) 


##1.- Calculate relationships among transcriptomes

#Pick only the columns you are interest in. They correspond to WT, mutant (shrJ0571) 
#and mutants with each complementation gene.

mydata2 <- mydata[,1:8]
mydata2

#Another way to do it. Needs to specify the row and column names
mydata3 <- cbind(mydata$J0571, mydata$shrJ0571, mydata$shrJ0571.BLJ, mydata$shrJ0571.JKD, mydata$shrJ0571.MGP, mydata$shrJ0571.NUC, mydata$shrJ0571.IME, mydata$shrJ0571.SCR)
colnames(mydata3) <- c("shr","WT","shr+BLJ","shr+JKD","shr+MGP","shr+NUC","shr+IME","shr+SCR")
rownames(mydata3) <-rownames(mydata)


# Perform the Principal Components Analysis
PCA<-princomp(mydata2,cor=FALSE,scores=TRUE)
PCA[1:7]

plot(PCA$loadings, main = "Ground Tissue Sorted Cells") # Vas a tener WT en un grafico cartesiano y el shr en otro sitio. Los que estan entre ellos
#van a ser los que presentan alguna complementarieda
text (PCA$loadings, labels =rownames(PCA$loadings), cex = 0.6,pos=c(4,2,2,2,2,2,2,2))

## 2.- Create transcriptomes to calibrate mutant complementation
#Now we plot intermediate transcriptomes. We are interested in
#transcriptome that recovers 25, 50 and 75% of the WT transcriptome
#We make vectors and multiply shr transcriptome by 0.25, 05 and 0.75.
#These will be artificial transcriptomes but we will be able to see where 
#it falls in PCA graph.

transformed25 = transform (mydata2, threeforths = 0.25 * mydata2$shrJ0571 + 0.75 * mydata2$J0571, half = 0.5 * mydata2$shrJ0571 + 0.5 * mydata2$J0571, quarter = 0.75 * mydata2$shrJ0571 + 0.25 * mydata2$J0571)

PCA2<-princomp(transformed25,cor=FALSE,scores=TRUE)
plot(PCA2$loadings, main = "Ground Tissue Sorted Cells") # Vas a tener WT en un grafico cartesiano y el shr en otro sitio. Los que estan entre ellos
text (PCA2$loadings, labels =rownames(PCA2$loadings), cex = 0.5)

#Here we see where transcriptomes are positioned. For example, "quarter" is located 1/4th
#of the way between mutant and WT.The ones located further from J0571 show an
#overcomplementation (gain of function)

## 3.- Add the SCR domain transcriptome
#What does it mean to have "overcomplementation"?
#We add SCRdomain into the plot. It is a Transc. Factor of stem cells
#J0571 (WT) is located between mutant J0571 and SCRdomain


mydata5 = transform(transformed25, SCRdomain = mydata$SCRdomain)

PCA3<-princomp(mydata5,cor=FALSE,scores=TRUE)
plot(PCA3$loadings, main = "Ground Tissue Sorted Cells") # Vas a tener WT en un grafico cartesiano y el shr en otro sitio. Los que estan entre ellos
text (PCA3$loadings, labels =rownames(PCA3$loadings),cex = 0.5)


##4.-Extract genes contributing the most to sample separation

#First we need to select the most important genes... 
#we take genes with lowest projection values in components 1 and 2 (based on scores)
#Para encontrar los genes mas importantes, cogemos del PCA los scores

#For heatman, RCOlorBrewer is used and function ColorRampPalette
#to define color gradient.
#Defined colors=c("blue","black","red) from less to more expressed.
#We will have (256) tonalities
#Como sacar los genes? Hay que usar la funcion sort(gen) y coge por ejemplo 100 genes. Saca los scores mas grandes para comp1 y mas pequenos a la comp2

#We pich 20 genes with highest scores of one Component and 20 of the other.
H1<-head(sort(PCA3$scores[,"Comp.1"]),decreasing=TRUE,20)
H2<-head(sort(PCA3$scores[,"Comp.2"]),decreasing=TRUE,20)
topVarGenesH <- as.matrix(c(H1,H2))
topVarGenesHExp <- as.matrix(mydata5[rownames(topVarGenesH),]) #Coger los dos vectores en forma de matriz
library(RColorBrewer)
newcolors<-colorRampPalette(colors=c("blue","black","red"))(256)
heatmap(topVarGenesHExp,scale="row",col =newcolors)

#this can also be done with highest proyection values in components 1 and 2
T1<-tail(sort(PCA3$scores[,"Comp.1"]),decreasing=TRUE,20)
T2<-tail(sort(PCA3$scores[,"Comp.2"]),decreasing=TRUE,20)
topVarGenesT <- as.matrix(c(T1,T2))

data4 = transform(transformed25, SCRdomain = mydata$SCRdomain) #Esto lo he puesto yo


topVarGenesTExp <- as.matrix(data4[rownames(topVarGenesT),])
library(RColorBrewer)
newcolors<-colorRampPalette(colors=c("blue","black","red"))(256)
heatmap(topVarGenesTExp,scale="row",col =newcolors)

#The most expressed are the genes of SCRdomain.
#Summary: We can identify the genes that affect globaltranscriptome the most
#using PCA