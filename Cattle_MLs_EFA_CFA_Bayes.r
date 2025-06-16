
#Set the folder for data analyses
setwd("/Users/____/Desktop/Data")

#MLs (AA, RF, ROC)
#######################################################################################################
#For AA (Figs.3a, S4, and S5) 
#memory clear
rm(list=ls(all=TRUE))
invisible(replicate(20, gc()))
library(dplyr)

#Set the folder for data analyses
setwd("/Users/____/Desktop/Data")

#Set the R library
if (!require("arules")) install.packages("arules")
if (!require("arulesViz")) install.packages("arulesViz")
require(arules) 
require(arulesViz)

dir.create("./Result_AA_RF_ROC", showWarnings = TRUE, recursive = FALSE, mode = "0777")

#Read the binarized file for AA (stored "Bibarized data for association analysis" in the sheet name "Fig.4a" of DataFCPCfinal241115.xlsx)
readfile="AA_binarized_raw_data.csv" # please check the sheetname "Fig. S4a and S5a" in the "Data_file.xlsx"

x=read.csv(readfile, header=T, colClasses="factor")
rulesAp1 <- apriori(x, parameter=list(support=0.3,confidence=0.5,maxlen=2))
rule_Ap1<-rulesAp1[quality(rulesAp1)$lift>1.3]

write(rule_Ap1,"./Result_AA_RF_XG/AA_calculated_data.csv",sep=",")

#For RF (Figs.3a, S4, and S5) 
#Read the raw data file for AA-selected components
library(randomForest)
FCPC.train <- read.csv("AA-selected_raw_data.csv") 
dim(FCPC.train) #[1] 91 X
set.seed(131)
train.x<- FCPC.train[,2:X] #dim(FCPC.train) #[1] 91 X
train.y<-as.factor(FCPC.train[,1])
model.rf<-tuneRF(train.x,train.y,doBest=T)
pred<-predict(model.rf,FCPC.test[,2:X]) #dim(FCPC.train) #[1] 91 X
table(FCPC.test[,1],pred)
print(model.rf$importance /sum(model.rf$importance))

write.csv(print(model.rf$importance /sum(model.rf$importance)),"randomForest_raw.csv")

rf_pred <- table(FCPC.test[,1],pred)
write.csv(rf_pred,"./Result_AA_RF_XG/randomForest_pred.csv")

FCPC.test_n <- cbind(y_1, train.x)

set.seed(22)
model = randomForest(y_1 ~ ., data = FCPC.test_n, importance = TRUE, proximity = TRUE)
print(model)
print(varImpPlot(model))

write.csv(print(varImpPlot(model)),"./Result_AA_RF_XG/randomForest_pred_importance_Gini.csv")

par(mar=c(100, 20, 30, 40))
rpp2 <- varImpPlot(model)
varFileName <- paste("./Result_AA_RF_XG/randomforest_tree_var.png",sep="") 
par(mar=c(100, 20, 30, 600)) 
png(file=varFileName, res=125, w=750, h=750)
rpp2 <- varImpPlot(model)
dev.off()


#ROC
library(Epi)

Data=read.csv("AA-selected_raw_data.csv") #
par("mar"=c(5,5,5,5))
ROC(test=Data$component_name, stat=Data$Category_name, plot="ROC") #

#######################################################################################################

#Confirm each sheet (Data_file.xlsx) 

#whole shapiro wilk test
input_r <- read.csv("file.csv") #col(column), each conditions; row, faecal bacteria or metabolites

dim(input_r) 
input=input[,-1] #deleted first line
library(MVN)
mvn(input, univariateTest="SW")
MVN_Cal=mvn(input, univariateTest="SW")

dir.create("./Result_whole_Statics", showWarnings = TRUE, recursive = FALSE, mode = "0777")
sink('./Result_whole_Statics/Shapiro_whole_MVN.txt', append = TRUE)
print(MVN_Cal)
sink()

#For CA (Fig.S10)
#Create the folder for Correlation network
dir.create("./CA_result", showWarnings = TRUE, recursive = FALSE, mode = "0777")

#Make the file with the components selected by ELA, AA, RF, and XGBoost
#File source: sheet name "Fig.3b" in "Data_file.xlsx"
d <- read.delim("Fig3b.txt", header=TRUE, row.name=1, sep="\t", fileEncoding="UTF-8",check.names=F)
d$name <- NULL
dat_t <- t(d)
no_zero<-subset(d,apply(d, 1, sum)>0)
dat_spearman<-cor(dat_t,dat_t,method="spearman")
write.csv(dat_spearman, "./CA_result/spearman_Pra.csv")

#For EFA (Fig. 4a)
dir.create("./EFA_result", showWarnings = TRUE, recursive = FALSE, mode = "0777")

#File source: sheet name "Figs.4a, S6 and S7" in "Data_file.xlsx"

library( psych )
library( GPArotation )
data=read.csv("Figs4a_S6_S7.csv") #
KMO(data)

KMO(data)

k_d=KMO(data)
sink('./EFA_result/KMO_Cattle_new.txt', append = TRUE)
print (k_d)
sink()

par("mar"=c(5,5,5,5))
fa.parallel(data, fa = "fa", use = "complete.obs")
abline(h = 0)
parallel=fa.parallel(data, fa = "fa", use = "complete.obs")
print(parallel)
abline(h = 0)
sink('./EFA_result/parallel.txt', append = TRUE)
print (parallel)
sink()

MAPminres <- vss(data, fm="minres")
print(MAPminres)

sink('./EFA_result/VSS_MAPminres.txt', append = TRUE)
print (MAPminres)
sink()

#Calculate as factor=7 based on the Kaiser-Guttman and vss
#Use rotate = "promax" based on the result using CA
result_7 = fa(data, nfactors = 7, fm = "minres", rotate = "promax", use = "complete.obs" )
print( result_7, digits = 3, sort = T)

sink('./EFA_result/nfactors_7.txt', append = TRUE)
print (result_7)
sink()

#factor analysis
result = fa( data, nfactors = 7, fm = "minres", rotate = "promax", use = "complete.obs" )
fa.diagram( result )

library(heatmaply)
heatmaply(fa(data, nfactors = 7, fm = "minress", rotate = "promax")$loadings,grid_gap = 1,subplot_widths = c(0.3, 0.2),subplot_heights = c(0.20, 0.70))


#For SEM (Fig.5a)

dir.create("./SEM_result", showWarnings = TRUE, recursive = FALSE, mode = "0777")

input=read.csv("Figs4a_S6_S7.csv") #

#Model
createModel10 <- function() {
  +     # variable ∼ dependent variable
    +     # + dependent variable
    return ("
    Hydroquinone ~ Caldibacillus
    AIB + Octadecanol ~ Caldibacillus + Muribaculum
    Glyceraldehyde  ~ AIB
            ")
}

#or

createModel10 <- function() {
  +     # variable ∼ dependent variable
    +     # + dependent variable
    return ("
    Caldibacillus ~ Butyrate
    AIB ~ Butyrate + Caldibacillus
    AIB  ~ Butyrate + Propionate
       ")
}

model.l10 <- createModel10()

library(lavaan)

# Select robust estimation (estimator = "MLR") based on the statistical result (p<0.05) using Shapiro-Wilk test

res.l10 <- lavaan(model.l10, data = input, estimator = "MLR", auto.var = TRUE)
summary(res.l10, fit.measure=TRUE,　standardized = TRUE)
fitMeasures(res.l10)

p0=model.l10
sink('./SEM_result/lavaanstaticsmodel.l10.txt', append = TRUE)
print (p0)
sink()

p1=summary(res.l10, fit.measure=TRUE,　standardized = TRUE)
sink('./SEM_result/lavaanstaticssummary.txt', append = TRUE)
print (p1)
sink()

p2=fitMeasures(res.l10)
sink('./SEM_result/lavaanstaticsfitMeasures.txt', append = TRUE)
print (p2)
sink()

#optinal evaluation:|SR|<2.58
sr <- residuals(res.l10, type = "standardized")
sink('./SEM_result/residuals.txt', append = TRUE)
print (sr)
sink()

# bootstrap based on ML

res.l12 <- lavaan(model.l10, data = input, auto.var = TRUE,se="bootstrap",bootstrap=1000)
P4=summary(res.l12, fit.measure=TRUE,　standardized = TRUE)
fitMeasures(res.l12)

sink('./SEM_result/lavaanstaticsmodel.lbootres12.txt', append = TRUE)
print (P4)
sink()

#Visualize the SEM
library(OpenMx)
library(psych)
library(semPlot)

semPaths(res.l10, what="stand", layout="tree", style="lisrel",
         shapeMan="rectangle", shapeLat="ellipse",
         sizeMan=6, residScale=5, posCol="darkgreen",
         negCol="violet", fade=FALSE, edge.label.cex=0.8) #Type1

semPaths(res.l10, what="stand", layout="tree", style="lisrel",
         shapeMan="rectangle", shapeLat="ellipse",
         sizeMan=8, residScale=5, posCol="darkgreen",
         negCol="violet", fade=FALSE, edge.label.cex=1.0) #Type2

semPaths(res.l10, what="stand", layout="tree", style="lisrel",
         shapeMan="rectangle", shapeLat="ellipse",
         sizeMan=10, residScale=5, posCol="darkgreen",
         negCol="violet", fade=FALSE, edge.label.cex=1.0) #gType3


#For BayesLiNGAM
dir.create("./Bayes_result", showWarnings = TRUE, recursive = FALSE, mode = "0777")

df <- read.csv("Figs4a_S6_S7.csv")

# Please install the commands for BayesLiNGAM based on the URL (https://www.cs.helsinki.fi/group/neuroinf/lingam/bayeslingam/)
setwd("/Users/_____/Bayeslingamtest")

library('fastICA')
source('bayeslingam/main/loud.R')

loud()

#Calculation
x <- df$Caldibacillus
y <- df$AIB
z <- df$Glyceraldehyde
w <- df$Hydroquinone
w1 <- df$Octadecanol
w2 <- df$Muribaculum
d <- data.frame(x1=x,x2=y,x3=z,x4=w,x5=w1,x6=w2)
result <- greedybayeslingam(d,model='GL')

#Data Export
sink('./Bayes_result/Bayes_Raw_data.txt', append = TRUE)
print (result)
sink()

#Visualize the plot
par("mar"=c(1,1,1,1))

#reference https://qiita.com/tomiyou/items/ca7032b1e0f1bf2b437b
library(igraph)
par(mfrow=c(2,3))
prob <- round(result$prob,digits=4) * 100
n.node <- max(result$components$node)
case <- nrow(result$DAGs)
node <- result$DAGs[,1:n.node]
res <- as.matrix(result$DAGs[,-c(1:n.node)],nrow=case)
name <-paste("X",1:n.node,sep="")

for(i in order(prob,decreasing=T)[1:6]){
  amat <- matrix(0,n.node,n.node)
  index <- node[i,]
  amat[lower.tri(amat,diag=FALSE)] <- res[i,]
  amat <- t(amat[index,index])
  g <- graph.adjacency(amat,mode="directed")
  E(g)$label <- NA
  pos <- layout.circle(g)
  rot <- matrix(c(0,1,-1,0),nrow=2,byrow=T)
  la <- pos %*% rot
  if(n.node == 2)la <- matrix(c(-1,0.5,1,0.5),,nrow=2,byrow=T)
  plot(g,layout=la,edge.arrow.size = 0.5,vertex.size=30,vertex.color = "lightgreen",vertex.label=name)     #Change the color dependent upon the conditions
  mtext(paste(prob[i],"%"),line=0)
}





