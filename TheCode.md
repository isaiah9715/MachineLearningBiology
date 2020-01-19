# Let's take a look at the data
```{r}
#Just reading in the data 
NSFX = read.csv("C:/Users/Isaiah1/Documents/Research/WallaceLab/NSFX_DataSets/thedat.csv")
head(NSFX)
dim(NSFX)
```
```

# We need to clean up the data in order to construct some heatmaps  

# Let's start out by cleaning it up
```{r}
#Take out the Gene Symbol, GSE, and Numeric Staging
NSFX_1 = NSFX[,-c(1,3)]

#Setting up a df of stages
stages = (NSFX_1[,1])

#This will now be a matrix of numeric values only
NSFX_1 = t(NSFX_1[,-1])

#Set the columns to the stages. Columns are stages and rows are genes
colnames(NSFX_1) = stages

#What we did was transform the row id's to the correct/corresponding stages
x = as.matrix(t(NSFX_1))

#Transform the stage variable to a factor variable 
stages = as.factor(stages)

#This method that is commented out will let you change the factor levels with maipulation of the numeric vector c()
stages = factor(stages, levels(stages)[c(1,3,4,2,5,6,7)])
list_stages = list(stages)
```  
  
## Aggregating the data
```{r}
#We will group the data by stage
heat = aggregate(x, by=list_stages, FUN = mean, na.rm=TRUE)

#transposing the aggregated data
heat2 = t(heat)
#giving the columns stage names 
colnames(heat2) = heat2[1,]

#getting rid of the extra row
heat2 = heat2[-1,]
#change the heat2 to numeric
mode(heat2) = 'numeric'
```
  
```{r include=FALSE}
#Libraries
library(glmnet)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggdendro)
library(reshape2)
library(grid)
library(gplots)
library(edgeR)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(ggpmisc) 

#Set a seed
set.seed(747)

#Formatiting a blank heatmap 
eaxis<-list(
  showticklabels=FALSE,
  showgrid=FALSE,
  zeroline=FALSE
)
```
  
## Further data preperation 
```{r include=FALSE}
gene.matrix = heat2
stage.matrix = t(heat2)


#We need to make heat00 a data.frame that looks like heat2
#for melt() we use later. The only difference is we have the
#rownames as an actual variable in the data.frame
genes2<-rownames(heat2)
genes2<-as.data.frame(genes2)
heat00<-cbind(genes2,heat2)
heat00<-as.data.frame(heat00)
attach(heat00)

#First make heat2 a data.frame so we can melt it down. Give it a long format which melt will do
#Melting the data into long format is necessary to create a heatmap
heat.long<-melt(heat00, id.vars="genes2")

#This commented out line of code is just rearranging the variales
#heat.long$variable = factor(heat.long$variable, levels(heat.long$variable)[c(1,3,4,2,5,6,7)])

#Taking out missing values
heat.long<-na.omit(heat.long)

```
  
## Actual Heatmap construction (unclustered).
```{r}
#Plotting the unclustered heatmap - no dendograms no order
unclustered = ggplot(heat.long, aes(x=variable, y=genes2)) +
  geom_tile(aes(fill=value)) + xlab("Stages") + ylab("Genes") +
  scale_fill_gradientn(colours=c("blue3","yellow","red3"),
                       values=c(0.0, 0.6, 1.0),
                       guide="colorbar")+
  theme(axis.text.y = element_text(size=8), axis.text.x = element_text(size=8), 
        axis.title.x=element_text(size=10), axis.title.y=element_text(size=10), legend.position='bottom', 
        legend.direction = 'horizontal', legend.text=element_text(size=8), legend.title=element_text(size=10))+
  labs(x="Stages", fill="Value") + ggtitle("Unclustered Heatmap: NSFX Genes")
print(unclustered)
```
  
  
## Construction of the dendograms; required for clustering.  
```{r}
#Gene Dendogram
heat.dendro<-as.dendrogram(hclust(d=dist(x=gene.matrix)))
dendro.plot1<-ggdendrogram(data=heat.dendro, rotate=TRUE)
dendro.plot1<-dendro.plot1 + theme(axis.text.y = element_text(size = 0), axis.text.x=element_text(size=0))
plot(dendro.plot1, xlab=eaxis, ylab=eaxis)

#Stage Dendogram
heat.dendro2<-as.dendrogram(hclust(d=dist(x=stage.matrix)))
dendro.plot2<-ggdendrogram(data=heat.dendro2, rotate=FALSE, labels=FALSE)
dendro.plot2<-dendro.plot2 + theme(axis.text.y = element_text(size = 0), axis.text.x = element_text(size=0))
plot(dendro.plot2, xlab=eaxis, ylab=eaxis)

#Ordering the dendograms 
heat.order<-order.dendrogram(heat.dendro)
heat.order2<-order.dendrogram(heat.dendro2)
heat.long$genes<-factor(x=heat.long$genes, levels=rownames(heat00)[heat.order], ordered=TRUE)
heat.long$variable<-factor(x=heat.long$variable, levels=colnames(heat00[,-1])[heat.order2], ordered=TRUE)
heat.long<-na.omit(heat.long)
```
  
## Construction of the clustered gene/stage heatmaps.  
```{r}
#plot heat map
Clustered<-ggplot(data=heat.long, aes(x=variable, y=genes)) + xlab("Stages") + ylab("Genes") + 
  geom_tile(aes(fill=value)) +
  scale_fill_gradientn(colours=c("blue3","yellow","red3"),
                       values=c(0.0, 0.6, 1.0),
                       guide="colorbar")+
  theme(axis.text.y = element_text(size=8), axis.text.x = element_text(size=8), 
        axis.title.x=element_text(size=10), axis.title.y=element_text(size=10), legend.position='bottom', 
        legend.direction = 'horizontal', legend.text=element_text(size=8), legend.title=element_text(size=10))+
  labs(x="Stages", fill="Value")
print(Clustered)
''' 
