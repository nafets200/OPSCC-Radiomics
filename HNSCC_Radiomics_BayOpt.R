
# The purpose of this script is to define a UNIFIED BAYESIAN HYPERPARAMETER OPTIMIZATION FRAMEWORK 
# for machine learning classification problems. 


# *************** install R packages ************************

#install.packages("car")
#install.packages("caret")
#install.packages("caTools")
#install.packages("compareC")
#install.packages("e1071")
#install.packages("effsize")
#install.packages("flexsurv")
#install.packages("ggfortify")
#install.packages("ggplot2")
#install.packages("grid")
#install.packages("KMsurv")
#install.packages("MASS")
#install.packages("MASS")
#install.packages("My.stepwise")
#install.packages("naivebayes")
#install.packages("neuralnet")
#install.packages("pROC")
#install.packages("randomForest")
#install.packages("randomForestSRC")
#install.packages("ranger")
#install.packages("readxl")
#install.packages("rms")
#install.packages("rpart.plot")
#install.packages("rpart")
#install.packages("stepPlr")
#install.packages("stepPlr")
#install.packages("survival")
#install.packages("readxl")
#install.packages("mltools")
#install.packages("xgboost")
#install.packages("mRMRe")
#install.packages("praznik")
#install.packages("glmnet")
#install.packages("OneR")
#install.packages("kernlab")
#install.packages("rBayesianOptimization")
#install.packages("resample")
#install.packages("plsRglm")
#install.packages("gridExtra")
#install.packages("OptimalCutpoints")
#install.packages("robustbase")
#install.packages("data.table")
#install.packages("matrixStats")





# *****************************************************************
# **********************Feature Selection loop*********************
# *****************************************************************

AllSelectTrainD <- list()
AllSelectTestD <- list()

for (s in seeds){
    
set.seed(s)
    
Allradiom_shuffle <- Allradiom[sample(nrow(Allradiom)),]
Allradiom_shuffle <- data.frame(ID=Allradiom_shuffle$ID, 
                                HPV_status=Allradiom_shuffle$HPV_status,
                                Allradiom_shuffle[,sample(3:ncol(Allradiom_shuffle))]
)
    
if(anaMod == "primT"){ stratInd <- as.vector(as.numeric(as.character(Allradiom_shuffle$HPV_status))) 
caseInd_unique <- data.frame(ID = Allradiom_shuffle$ID)}
    
if(anaMod == "LN" | anaMod == "Consensus_LN" | anaMod == "Consensus_primT_LN"){
      
  caseInd <- data.frame(ID=Allradiom_shuffle$ID, HPV_status=Allradiom_shuffle$HPV_status)
  caseInd$ID <- substr(as.character(caseInd$ID), 1, nchar(as.character(caseInd$ID))-5)
      
  caseInd_unique <- caseInd[!duplicated(caseInd$ID),]
      
  stratInd <- as.vector(as.numeric(as.character(caseInd_unique$HPV_status)))
}    
    
    
folds <- folds(stratInd,
               nfolds=nfolds, 
               stratified=TRUE, 
               seed = s)
    
    
for (x in 1:nfolds){
      
print("FS, the test fold is:")
print(x)
      
if(anaMod == "primT"){ 
  testD=Allradiom_shuffle[which(folds==x),] 
  trainD=Allradiom_shuffle[which(!folds==x),] }
      
if(anaMod == "LN" | anaMod == "Consensus_LN" | anaMod == "Consensus_primT_LN"){
  testD=Allradiom_shuffle[ grep( paste(caseInd_unique$ID[folds==x], collapse = "|"), Allradiom_shuffle$ID, invert = FALSE), ] 
  trainD=Allradiom_shuffle[ grep( paste(caseInd_unique$ID[folds==x], collapse = "|"), Allradiom_shuffle$ID, invert = TRUE), ] 
}
      
      

for(j in (colnames(trainD)[3:length(colnames(trainD))]) ){
  mean <- mean(trainD[,j])
  SD <- sd(trainD[,j])
    
  trainD[,j] <- as.numeric((trainD[,j]-mean)/SD)
  testD[,j] <- as.numeric((testD[,j]-mean)/SD)
}

      
if(anaMod == "Consensus_LN" | anaMod == "Consensus_primT_LN"){
      
  Cohort_median <- colMedians(as.matrix( trainD[,3:ncol(trainD)] ), 
                              na.rm = FALSE, 
                              hasNA = FALSE, 
                              keep.names=TRUE)
  
  
  Consensus_trainD <- data.frame(ID = caseInd_unique$ID[!folds==x], 
                                 matrix(0, nrow = length(which(!folds==x)), ncol = ncol(trainD)-1 ) )
  colnames(Consensus_trainD) <- colnames(trainD)
  Consensus_trainD$HPV_status <- ordered(Consensus_trainD$HPV_status, levels=c(0, 1))
  
  Consensus_testD <- data.frame(ID = caseInd_unique$ID[folds==x],
                                matrix(0, nrow = length(which(folds==x)), ncol = ncol(testD)-1) )
  colnames(Consensus_testD) <- colnames(testD)
  Consensus_testD$HPV_status <- ordered(Consensus_testD$HPV_status, levels=c(0, 1))
  
  
  for(j in caseInd_unique$ID[!folds==x] ){
    
    All_ROI <- trainD[ trainD$ID %in% Allradiom_shuffle$ID[caseInd$ID == j], ]
    HPV_status <- All_ROI$HPV_status[1]
    All_ROI <- t(All_ROI[,3:ncol(All_ROI)])
    
    max_Delta <- abs( All_ROI - matrix(Cohort_median, nrow = length(Cohort_median), ncol = ncol(All_ROI)) )
    
    max_Delta_index <- (apply(max_Delta,1,which.max))
      
    Consensus_trainD[Consensus_trainD$ID == j, 2:ncol(Consensus_trainD)] <- data.frame(HPV_status = HPV_status,
                                                                                           t(All_ROI[ cbind(seq_along(max_Delta_index), max_Delta_index) ]) )
    }
        
  for(j in caseInd_unique$ID[folds==x] ){
          
    All_ROI <- testD[ testD$ID %in% Allradiom_shuffle$ID[caseInd$ID == j], ]
    HPV_status <- All_ROI$HPV_status[1]
    All_ROI <- t(All_ROI[,3:ncol(All_ROI)])
        
    max_Delta <- abs( All_ROI - matrix(Cohort_median, nrow = length(Cohort_median), ncol = ncol(All_ROI)) )
      
    max_Delta_index <- (apply(max_Delta,1,which.max))
          
    Consensus_testD[Consensus_testD$ID == j, 2:ncol(Consensus_testD)] <- data.frame(HPV_status= HPV_status, 
                                                                                          t(All_ROI[ cbind(seq_along(max_Delta_index), max_Delta_index) ]) )
  }
        
  trainD <- Consensus_trainD
  testD <- Consensus_testD
        
}
      
      
      
      
# ******************** Standardized Interface **********************
  
# This script will supply 2 data frames to the to feature selection section: 
# "trainD" and "testD" are both structured as follows: 
# 1st column = case ID, 2nd column = target variable (HPV_status), all other columns = features
# Every row represents a case. 

      
####################################################################
# FEATURE SELECTION
# Supply code for any feature selection method below
####################################################################
      
      
 

####################################################################

AllSelectTrainD[[paste(s, x, sep = "-")]] <- SelectTrainD
AllSelectTestD[[paste(s, x, sep = "-")]] <- SelectTestD
      

}
}









# *****************************************************************
# ********************Bayesian Optimization loop*******************
# *****************************************************************

capture_sd_AUC <- c()

passed_params <- data.frame()


BayOpt_cv <- function(params
                      ) {
  

passed_params <<- rbind(passed_params,
                        data.frame(nFeatures=nFeatures,
                                   mtry=mtry,
                                   cutoff=cutoff,
                                   maxnodes=maxnodes))

captureModelOutput <- data.frame()


for (s in seeds){
    
    
set.seed(s)
    
print("BayOpt, the seed is:")
print(s)
    
    
for (x in 1:nfolds){
  
print("BayOpt, the test fold is:")
print(x)


SelectTrainD <- as.data.frame(AllSelectTrainD[[ paste(s, x, sep = "-") ]])
SelectTestD <- as.data.frame(AllSelectTestD[[ paste(s, x, sep = "-") ]])

SelectTrainD <- SelectTrainD[, 1: (nFeatures+2) ]
SelectTestD <- SelectTestD[, 1: (nFeatures+2) ]




####################################################################
# MACHINE LEARNING CLASSIFIER
# Supply code for any ML classification method below
####################################################################



###################################################################   

AUC_test <- as.numeric(auc_roc(probtest_test, SelectTestD$HPV_status))
captureModelOutput <- rbind(captureModelOutput,
                            data.frame(FS_method = replicate(length(probtest_test), FS_identifier), 
                                       ML_method = replicate(length(probtest_test), ML_identifier), 
                                       ROI = replicate(length(probtest_test), anaMod), 
                                       Imaging_Modality = replicate(length(probtest_test), imageMod), 
                                       Nr_of_CV_folds = replicate(length(probtest_test), nfolds),
                                       ID=SelectTestD$ID, 
                                       HPV_status=SelectTestD$HPV_status,
                                       probtest_test=probtest_test, 
                                       AUC_test=replicate(length(probtest_test), AUC_test), 
                                       seed=replicate(length(probtest_test), s), 
                                       Test_fold=replicate(length(probtest_test), x))
                            )
                            

}
}

capture_AUC <- c()

for(i in unique(captureModelOutput$seed)){
  
  for(k in 1:nfolds){
      capture_AUC <- c(capture_AUC, 
                     captureModelOutput$AUC_test[ captureModelOutput$seed == i & captureModelOutput$Test_fold == k][1])
    }}

mean_AUC <- mean(capture_AUC)



sd_AUC <- sd(capture_AUC)
capture_sd_AUC <<- c(capture_sd_AUC, sd_AUC)


ifelse(anaMod == "LN", (index = Allradiom$ID), (index = caseInd_unique$ID))

meanProb <- c()

for(j in index){
  meanProb <- c(meanProb, 
                mean(captureModelOutput[ captureModelOutput$ID == j, "probtest_test"]))
}


BayOpt_result <- list(Score =  mean_AUC, Pred = meanProb)

ML_identifier <<- ML_identifier

return(BayOpt_result)


}


OptHyperP <- BayesianOptimization(

BayOpt_cv,
 bounds = list(parambounds
  ),
  
init_points = 20,
n_iter =  n_iter,
acq = "ucb", 
kappa = 2.576, 
eps = 0.9,
verbose = TRUE
)

print(OptHyperP$Best_Value)
print(OptHyperP$Best_Par)
print(OptHyperP$History)

OptHyperP_history <- data.frame(OptHyperP$History, SD_AUC = capture_sd_AUC)
OptHyperP_history <- rbind(OptHyperP_history, 
                           c(replicate(ncol(OptHyperP_history), NA)),
                           OptHyperP_history[order(OptHyperP_history$Value, decreasing = TRUE)[1:5],])

