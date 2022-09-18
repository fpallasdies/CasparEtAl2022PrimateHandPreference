#write.csv(sumtab, "sumtab.csv", row.names = FALSE)
sumtab <- read.csv("sumtab.csv", header = TRUE, sep = ';', stringsAsFactors=TRUE)
#litPopulation <- read.csv("LiteraturedataPopLevel.csv", header = TRUE, sep = '\t')

#litPopulation$NameInTree <- NA

#for (i in 1:length(litPopulation$Species)){
#  litPopulation$NameInTree[i] <- paste(litPopulation$Genus[i], litPopulation$Species[i], sep="_")
#} 






predictor <- read.csv("predictortab.csv", header = TRUE, sep = ',', stringsAsFactors=TRUE, fileEncoding="UTF-8-BOM")
nameintree <- predictor$Name.in.tree

predictor <- read.csv("predictortab.csv", header = TRUE, sep = ',', stringsAsFactors=TRUE, row.names = 9, fileEncoding="UTF-8-BOM")
predictor$Name.in.tree <- nameintree





#Phylo Data

library(ape)
library(nlme)
library(geiger)
library(phytools)

arbol <- read.nexus("consensusTree_10kTrees_laterality_strength.nex")
obj <-name.check(arbol,predictor, data.names = predictor$Name.in.tree)



data_with_names <- predictor$HI
names(data_with_names) <- predictor$Name.in.tree



arbol.dir<- read.nexus("consensusTree_10kTrees_laterality_direction.nex")

predictor.dir<- predictor
obj <-name.check(arbol.dir,predictor.dir, data.names = predictor.dir$Name.in.tree)
obj
predictor.dir <- subset(predictor.dir, predictor.dir$Name.in.tree %in% obj$data_not_tree == FALSE )
name.check(arbol.dir,predictor.dir, data.names = predictor.dir$Name.in.tree)


obj <-name.check(arbol,predictor.dir, data.names = predictor.dir$Name.in.tree)
obj
arbol.dir<-drop.tip(arbol, obj$tree_not_data)
name.check(arbol.dir,predictor.dir, data.names = predictor.dir$Name.in.tree)

bm<-corBrownian(1, arbol, form= ~Name.in.tree)
bmdir<-corBrownian(1, arbol.dir, form= ~Name.in.tree)




predictor.nohumans <- subset(predictor, Name.in.tree != "Homo_sapiens")
obj <-name.check(arbol,predictor.nohumans, data.names = predictor.nohumans$Name.in.tree)
obj
arbol.nohumans<-drop.tip(arbol, obj$tree_not_data)
name.check(arbol.nohumans,predictor.nohumans, data.names = predictor.nohumans$Name.in.tree)
bmnohumans<-corBrownian(1, arbol.nohumans, form= ~Name.in.tree)


predictor.dir.nohumans <- subset(predictor.dir, Name.in.tree != "Homo_sapiens")
obj <-name.check(arbol.dir,predictor.dir.nohumans, data.names = predictor.dir.nohumans$Name.in.tree)
obj
arbol.dir.nohumans<-drop.tip(arbol.dir, obj$tree_not_data)
name.check(arbol.dir.nohumans,predictor.dir.nohumans, data.names = predictor.dir.nohumans$Name.in.tree)
bmdirnohumans<-corBrownian(1, arbol.dir.nohumans, form= ~Name.in.tree)



tab <- read.csv("detailedtab.csv", header = TRUE, sep = ',', stringsAsFactors=TRUE)

