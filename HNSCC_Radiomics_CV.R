
# The purpose of this script is to define a UNIFIED K-FOLD CROSS VALIDATION FRAMEWORK for 
# machine learning classification problems. 



# *************** install R packages ************************

install.packages("car")
install.packages("caret")
install.packages("caTools")
install.packages("compareC")
install.packages("e1071")
install.packages("effsize")
install.packages("flexsurv")
install.packages("ggfortify")
install.packages("ggplot2")
install.packages("KMsurv")
install.packages("MASS")
install.packages("MASS")
install.packages("My.stepwise")
install.packages("naivebayes")
install.packages("neuralnet")
install.packages("pROC")
install.packages("randomForest")
install.packages("randomForestSRC")
install.packages("ranger")
install.packages("readxl")
install.packages("rms")
install.packages("rpart.plot")
install.packages("rpart")
install.packages("stepPlr")
install.packages("stepPlr")
install.packages("survival")
install.packages("readxl")
install.packages("mltools")
install.packages("xgboost")
install.packages("mRMRe")
install.packages("praznik")
install.packages("glmnet")
install.packages("OneR")
install.packages("kernlab")
install.packages("resample")
install.packages("plsRglm")
install.packages("gridExtra")
install.packages("OptimalCutpoints")
install.packages("robustbase")
install.packages("data.table")
install.packages("matrixStats")



captureTuning <- data.frame()
capture_test_ClassProb <- data.frame()
capture_train_ClassProb <- data.frame()

captureAUC_test <- data.frame()
captureAUC_train <- data.frame()
captureOptCut <- data.frame()



for (s in c(seeds)){

  set.seed(s)
  
  print("the seed is:")
  print(s)

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

print("the test fold is:")
print(x)
    
if(anaMod == "primT"){ 
  testD=Allradiom_shuffle[which(folds==x),] 
  trainD=Allradiom_shuffle[which(!folds==x),] }

if(anaMod == "LN" | anaMod == "Consensus_LN" | anaMod == "Consensus_primT_LN"){
  testD=Allradiom_shuffle[ grep( paste(caseInd_unique$ID[folds==x], collapse = "|"), Allradiom_shuffle$ID, invert = FALSE), ] 
  trainD=Allradiom_shuffle[ grep( paste(caseInd_unique$ID[folds==x], collapse = "|"), Allradiom_shuffle$ID, invert = TRUE), ] 
 }
 
zscore <- preProcess(subset(trainD, select=-c(ID, HPV_status)), method = c("center", "scale"))
trainD[,3:ncol(trainD)] <- predict(zscore, trainD[,3:ncol(trainD)])
testD[,3:ncol(testD)] <- predict(zscore, testD[,3:ncol(testD)])
 

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
 
 



####################################################################
# FEATURE SELECTION
# Supply code for any feature selection method below
####################################################################
 


####################################################################
# MACHINE LEARNING CLASSIFIER
# Supply code for any ML classification method below
####################################################################
  

      
###################################################################   
      
AUC_test <- as.numeric(auc_roc(probtest_test, SelectTestD$HPV_status))
AUC_train <- as.numeric(auc_roc(probtest_train, SelectTrainD$HPV_status))
 
captureAUC_test[1,ncol(captureAUC_test)+1] <- as.numeric(AUC_test)
colnames(captureAUC_test)[ncol(captureAUC_test)] <- x
    
captureAUC_train[1,ncol(captureAUC_train)+1] <- AUC_train
colnames(captureAUC_train)[ncol(captureAUC_train)] <- x

capture_ConfMatrix <- data.frame()

for(i in seq(0, 1, 0.001)){
ConfMatrix <- confusionMatrix(data = factor(ifelse(probtest_test>=i  , 1, 0), levels = c(0,1)), 
                              reference = as.factor(SelectTestD$HPV_status), 
                              positive = "1",
                              mode = "everything"
                              )

capture_ConfMatrix <- rbind(capture_ConfMatrix, 
                            data.frame(cutoff=i, 
                                       Sensitivity=as.numeric(ConfMatrix$byClass[1]), 
                                       Specificity=as.numeric(ConfMatrix$byClass[2]), 
                                       Delta=as.numeric(ConfMatrix$byClass[1])-as.numeric(ConfMatrix$byClass[2]))
                            )}

if( any(  capture_ConfMatrix$Sensitivity == capture_ConfMatrix$Specificity) ){
    optimal_cutoff_balanced = median(
    capture_ConfMatrix$cutoff[capture_ConfMatrix$Sensitivity == capture_ConfMatrix$Specificity])
    } else {
          optimal_cutoff_balanced = mean(c(
          capture_ConfMatrix$cutoff[ which(capture_ConfMatrix$Delta<0)[1] ], 
          capture_ConfMatrix$cutoff[ (which(capture_ConfMatrix$Delta<0)[1])-1 ]
    ) )}



ConfMatrix_balanced <- confusionMatrix(data = factor(ifelse(probtest_test>=optimal_cutoff_balanced, 1, 0), levels = c(0,1)), 
                                      reference = as.factor(SelectTestD$HPV_status), 
                                      positive = "1",
                                      mode = "everything")

captureOptCut <- rbind(captureOptCut, 
                       data.frame(cutoff_balanced=optimal_cutoff_balanced, 
                                  Sensitivity_bal_cut=as.numeric(ConfMatrix_balanced$byClass[1]), 
                                  Specificity_bal_cut=as.numeric(ConfMatrix_balanced$byClass[2]), 
                                  Accuracy_bal_cut=as.numeric(ConfMatrix_balanced$overall[1]), 
                                  seed=s, 
                                  test_fold=x)  )
  
  
 
ClassProb_aggregated_test <- data.frame(FS_method = replicate(length(probtest_test), FS_identifier),
                                   ML_method = replicate(length(probtest_test), ML_identifier),
                                   ROI = replicate(length(probtest_test), anaMod),
                                   Imaging_Modality = replicate(length(probtest_test), imageMod),
                                   Nr_of_CV_folds = replicate(length(probtest_test), nfolds),
                                   ID=SelectTestD$ID, 
                                   HPV_status=SelectTestD$HPV_status,
                                   probtest_test=probtest_test, 
                                   AUC_test=replicate(length(probtest_test), AUC_test), 
                                   seed=replicate(length(probtest_test), s), 
                                   test_fold=replicate(length(probtest_test), x), 
                                   tuning_param_value=replicate(length(probtest_test), y))

capture_test_ClassProb <- rbind(capture_test_ClassProb, ClassProb_aggregated_test)


ClassProb_aggregated_train <- data.frame(FS_method = replicate(length(probtest_train), FS_identifier),
                                   ML_method = replicate(length(probtest_train), ML_identifier),
                                   ROI = replicate(length(probtest_train), anaMod),
                                   Imaging_Modality = replicate(length(probtest_train), imageMod),
                                   Nr_of_CV_folds = replicate(length(probtest_train), nfolds),
                                   ID=SelectTrainD$ID, 
                                   HPV_status=SelectTrainD$HPV_status,
                                   probtest_train=probtest_train, 
                                   AUC_train=replicate(length(probtest_train), AUC_train), 
                                   seed=replicate(length(probtest_train), s), 
                                   TEST_fold=replicate(length(probtest_train), x), 
                                   tuning_param_value=replicate(length(probtest_train), y))

capture_train_ClassProb <- rbind(capture_train_ClassProb, ClassProb_aggregated_train)



  }
}



mean_AUC_test <- rowMeans(captureAUC_test)
sd_AUC_test <- sd(as.numeric(captureAUC_test[which.max(mean_AUC_test),]))
mean_AUC_train <- rowMeans(captureAUC_train)
sd_AUC_train <- sd(as.numeric(captureAUC_train[which.max(mean_AUC_test),]))
mean_cutoff_balanced <- mean(captureOptCut$cutoff_balanced)


captureTuning[nrow(captureTuning)+1, "FS_method"] <- FS_identifier
captureTuning[nrow(captureTuning), "ML_method"] <- ML_identifier
captureTuning[nrow(captureTuning), "ROI"] <- anaMod
captureTuning[nrow(captureTuning), "Imaging_Modality"] <- imageMod
captureTuning[nrow(captureTuning), "Nr_of_CV_folds"] <- nfolds

captureTuning[nrow(captureTuning), "mean_AUC_test"] <- mean_AUC_test
captureTuning[nrow(captureTuning), "sd_AUC_test"] <- sd_AUC_test
captureTuning[nrow(captureTuning), "mean_AUC_train"] <- mean_AUC_train[as.numeric(which.max(mean_AUC_test))]
captureTuning[nrow(captureTuning), "sd_AUC_train"] <- sd_AUC_train

captureTuning[nrow(captureTuning), "mean_cutoff_balanced"] <- mean_cutoff_balanced


