
library(nlme)
library(MuMIn)
library(MASS)
library(lme4)

library(ggplot2)
library(modelsummary)

source("Preprocessing.r")

############### 
### ChiSquared Tests

#Creating a combined table of literature Data and our own
HIs<-sumtab$HI
Clade<- sumtab$Clade
TreeName<-sumtab$Treename
characterization<- sumtab$Category
Sex <- sumtab$Sex
Age <- sumtab$Age
Genus<-as.character(sumtab$Genus)
for (i in 1:length(characterization)) {
  if (characterization[i] == "left-handed") {
    characterization[i]<- "left"
  }
  if (characterization[i] == "right-handed") {
    characterization[i]<- "right"
  } 
}

compiledData <- data.frame(HIs, Clade, TreeName, characterization, Genus, Sex, Age)


compiledData$characterization <- as.factor(compiledData$characterization)

levels(compiledData$characterization)

compiledData$TreeName <- as.factor(compiledData$TreeName)
compiledData$Genus <- as.factor(compiledData$Genus)






# Chisquared and t-tests
newworld <-subset(compiledData, compiledData$Clade == "Platyrrhini")
oldworld <-subset(compiledData, compiledData$Clade == "Cercopithecoidea")
homi <-subset(compiledData, compiledData$Clade == "Hominoidea")
homiWOSap <-subset(compiledData, compiledData$Clade == "Hominoidea" & compiledData$Genus != "Homo")


summary(homiWOSap$characterization)
#187          302          413 
# without humans
#0.2073171

summary(oldworld$characterization)
#123           181           170 
#0.259493


summary(newworld$characterization)
#45          129          109 
#0.1731449


pval <- vector()
ttestpval <- vector()
names <- vector()
clade <- vector()
ambifreq <- vector()
HImean <- vector()
HIabs <- vector()

for (i in 1:length(levels(compiledData$TreeName))) {
  print(i)
  tempsub <- subset(compiledData, compiledData$TreeName == levels(compiledData$TreeName)[i] )
  names <- append(names, levels(compiledData$TreeName)[i])
  clade <- append(clade, tempsub$Clade[1])
  
  
  freq <-summary(tempsub$characterization)
  

  
  HImean <- append(HImean, mean(tempsub$HI))
  HIabs <- append(HIabs, mean(abs(tempsub$HI)))
  
  
  if (tempsub$Clade[1] == "Platyrrhini") {
    res <- chisq.test(freq, p=c(0.173145, 0.4134275,0.4134275))
  }
  if (tempsub$Clade[1] == "Cercopithecoidea") {
    #res <- chisq.test(freq, p=c(0.2414634, 0.3792683,0.3792683))
    res <- chisq.test(freq, p=c(0.259493, 0.3702535,0.3702535))
  }
  if (tempsub$Clade[1] == "Hominoidea") {
    res <- chisq.test(freq, p=c(0.2073172, 0.3963414,0.3963414))
  }
  
  pval <- append(pval, res$p.value)
  
  if (!is.na(tempsub$HIs[1])) {
    x <- t.test(tempsub$HIs)
    ttestpval<- append(ttestpval, x$p.value )
    
  }
  else{
    ttestpval<- append(ttestpval, NA )
  }
  



}
chiresults <- data.frame(names,pval, ttestpval, clade, HImean, HIabs)


hist(pval)

hist(ttestpval)


ks.test(pval,"punif",0,1)
ks.test(ttestpval,"punif",0,1)



#Chisquared tests on genus level
pval <- vector()
ttestpval <- vector()
names <- vector()
clade <- vector()
ambifreq <- vector()
HImean <- vector()
HIabs <- vector()

for (i in 1:length(levels(compiledData$Genus))) {
  print(i)
  tempsub <- subset(compiledData, compiledData$Genus == levels(compiledData$Genus)[i] )
  names <- append(names, levels(compiledData$Genus)[i])
  clade <- append(clade, tempsub$Clade[1])
  
  
  freq <-summary(tempsub$characterization)
  
  
  
  HImean <- append(HImean, mean(tempsub$HI))
  HIabs <- append(HIabs, mean(abs(tempsub$HI)))
  
  
  if (tempsub$Clade[1] == "Platyrrhini") {
    res <- chisq.test(freq, p=c(0.173145, 0.4134275,0.4134275)) #(0.159574, 0.420213,0.420213)
  }
  if (tempsub$Clade[1] == "Cercopithecoidea") {
    #res <- chisq.test(freq, p=c(0.2414634, 0.3792683,0.3792683))
    res <- chisq.test(freq, p=c(0.259493, 0.3702535,0.3702535))    #(0.219114, 0.390443,0.390443)
  }
  if (tempsub$Clade[1] == "Hominoidea") {
    res <- chisq.test(freq, p=c(0.2073172, 0.3963414,0.3963414))        #(0.20041, 0.399795,0.399795)
  }
  
  pval <- append(pval, res$p.value)
  
  if (!is.na(tempsub$HIs[1])) {
    x <- t.test(tempsub$HIs)
    ttestpval<- append(ttestpval, x$p.value )
    
  }
  else{
    ttestpval<- append(ttestpval, NA )
  }
  
  
  
  
}
chiresults.genus <- data.frame(names,pval, ttestpval, clade, HImean, HIabs)


######################




############
#Phylogenetic signal
HIlist <- predictor.dir$HI
names(HIlist) <- rownames(predictor.dir)
phylosig(arbol.dir, HIlist, method = "lambda", test = TRUE)


HIabslist <- predictor$HIabs
names(HIabslist) <- rownames(predictor)

phylosig(arbol, HIabslist, method = "lambda", test = TRUE)





#################
###Phylogenetic generalized least squares





glsControl(maxIter = 100, msMaxIter = 100)

#############
#Lateralization strength
model.full<-gls(HIabs~log(female.endocranial.volume)+Ecology+Habitual.tool.use.in.wild.populations, data=predictor,correlation=corPagel(0.885359,arbol, form=~Name.in.tree, fixed=TRUE), method ="ML")


dfirst<-dredge(model.full)
dsec<-get.models(dfirst, subset=TRUE)
dAbs<-model.avg(dsec, revised.var=TRUE)
summary(dAbs)


modelplot(dAbs) 


###### Direction

model.fulldir<-gls(HI~log(female.endocranial.volume)+Ecology+Habitual.tool.use.in.wild.populations, data=predictor.dir, correlation=corPagel(0,arbol.dir, form=~Name.in.tree, fixed=TRUE), method ="ML")


dfirst.dir<-dredge(model.fulldir)
dsec.dir<-get.models(dfirst.dir, subset=TRUE)
ddir<-model.avg(dsec.dir, revised.var=TRUE)
summary(ddir)


modelplot(ddir)





###### without humans
model.fulldirnh<-gls(HI~log(female.endocranial.volume)+Ecology+Habitual.tool.use.in.wild.populations, data=predictor.dir.nohumans, correlation=corPagel(0,arbol.dir.nohumans,  form=~Name.in.tree, fixed=TRUE), method ="ML")

dfirst.dir<-dredge(model.fulldirnh)
dsec.dir<-get.models(dfirst.dir, subset=TRUE)
ddir.nohumans<-model.avg(dsec.dir, revised.var=TRUE)
summary(ddir.nohumans)

modelplot(ddir.nohumans)


##########
#GLM on insertions

tabalt<- subset(tab, is.na(tab$Handed.in.with..)== FALSE)
tabalt$HIInsert <- (tabalt$HIInsert +1) /2

firstmodel<-glmer(HIInsert~ Handed.in.with..+(1 |Individual), data=tabalt, family=binomial )#, na.action="na.omit")

summary(firstmodel)

##########
#Ancestral character estimates

#Direction
fastAnc(arbol.dir, HIlist, vars= TRUE, CI= TRUE)

#Strength
fastAnc(arbol, HIabslist, vars= TRUE, CI= TRUE)

##############################
#############################
###############################
# Test normality:
library(olsrr)
library(car)
shapiro.test(model.full$residuals)
shapiro.test(model.fulldir$residuals)
shapiro.test(model.fulldirnh$residuals)


vif(model.full)
vif(model.fulldir)
vif(model.fulldirnh)




#brms modelling of effects of sex and age cohorts on individual-level laterality
library(brms)
A <- ape::vcv.phylo(arbol) 
sumalt<- subset(sumtab, is.na(sumtab$Sex)== FALSE)
sumalt<- subset(sumalt, is.na(sumalt$Age)== FALSE)

indimodel<-brm(HI~Sex+Age+ (1|gr(Treename, cov = A )), data=sumalt, data2 = list(A = A))
summary(indimodel)
#plot(indimodel)
#pp_check(indimodel, type = "ecdf_overlay")



indimodel.abs<-brm(AbsHI~Sex+Age+ (1|gr(Treename, cov = A )), data=sumalt, data2 = list(A = A))
summary(indimodel.abs)
#plot(indimodel.abs)
#pp_check(indimodel.abs, type = "ecdf_overlay")





sumplat<- subset(sumalt, sumalt$Clade == "Platyrrhini")

indimodel.plat<-brm(HI~Sex+Age+ (1|gr(Treename, cov = A )), data=sumplat, data2 = list(A = A),iter = 12000)
summary(indimodel.plat)
#plot(indimodel.plat)
#pp_check(indimodel.plat, type = "ecdf_overlay")



indimodel.abs.plat<-brm(AbsHI~Sex+Age+ (1|gr(Treename, cov = A )), data=sumplat, data2 = list(A = A),iter = 5000)
summary(indimodel.abs.plat)
#plot(indimodel.abs.plat)
#pp_check(indimodel.abs.plat, type = "ecdf_overlay")



sumcer<- subset(sumalt, sumalt$Clade == "Cercopithecoidea")

indimodel.cer<-brm(HI~Sex+Age+ (1|gr(Treename, cov = A )), data=sumcer, data2 = list(A = A))
summary(indimodel.cer)
#plot(indimodel.cer)
#pp_check(indimodel.cer, type = "ecdf_overlay")



indimodel.abs.cer<-brm(AbsHI~Sex+Age+ (1|gr(Treename, cov = A )), data=sumcer, data2 = list(A = A))
summary(indimodel.abs.cer)
#plot(indimodel.abs.cer)
#pp_check(indimodel.abs.cer, type = "ecdf_overlay")



sumhom<- subset(sumalt, sumalt$Clade == "Hominoidea")
sumhom<- subset(sumhom, sumhom$Species != "sapiens" & sumhom$Species != "troglodytes")

indimodel.hom<-brm(HI~Sex+Age+ (1|gr(Treename, cov = A )), data=sumhom, data2 = list(A = A),iter = 7000)
summary(indimodel.hom)
#plot(indimodel.hom)
#pp_check(indimodel.hom, type = "ecdf_overlay")



indimodel.abs.hom<-brm(AbsHI~Sex+Age+ (1|gr(Treename, cov = A )), data=sumhom, data2 = list(A = A),iter = 7000)
summary(indimodel.abs.hom)
#plot(indimodel.abs.hom)
#pp_check(indimodel.abs.hom, type = "ecdf_overlay")




