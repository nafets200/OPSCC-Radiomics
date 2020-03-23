
# The purpose of this script is to perform final feature selection
# and final model training on a dedicated training cohort and final 
# model testing on two dedicated test cohorts (trainD_test & trainD_test_indep)


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




if(anaMod == "Consensus_primT_LN"){ 
  Allradiom_train$ID[grep("_LN_", Allradiom_train$ID, invert = TRUE)] <- paste(Allradiom_train$ID[grep("_LN_", Allradiom_train$ID, invert = TRUE)], "_____", sep="") 
  Allradiom_test$ID[grep("_LN_", Allradiom_test$ID, invert = TRUE)] <- paste(Allradiom_test$ID[grep("_LN_", Allradiom_test$ID, invert = TRUE)], "_____", sep="") 
  if(imageMod == "PET"){ Allradiom_test_indep$ID[grep("_LN_", Allradiom_test_indep$ID, invert = TRUE)] <- paste(Allradiom_test_indep$ID[grep("_LN_", Allradiom_test_indep$ID, invert = TRUE)], "_____", sep="") }
  }


zscore <- preProcess(subset(Allradiom_train, select=-c(ID, HPV_status)), method = c("center", "scale"))

Allradiom_train[,3:ncol(Allradiom_train)] <- predict(zscore, Allradiom_train[,3:ncol(Allradiom_train)])
trainD <- Allradiom_train

Allradiom_test[,3:ncol(Allradiom_test)] <- predict(zscore, Allradiom_test[,3:ncol(Allradiom_test)])
testD <- Allradiom_test

if(imageMod == "PET"){ 
  Allradiom_test_indep[,3:ncol(Allradiom_test_indep)] <- predict(zscore, Allradiom_test_indep[,3:ncol(Allradiom_test_indep)]) 
  testD_indep <- Allradiom_test_indep
  }


if(anaMod == "Consensus_LN" | anaMod == "Consensus_primT_LN"){
  
  caseInd_trainD <- trainD$ID
  caseInd_trainD <- substr(as.character(caseInd_trainD), 1, nchar(as.character(caseInd_trainD))-5)
  caseInd_unique_trainD <- caseInd_trainD[!duplicated(caseInd_trainD)]
  
  caseInd_testD <- testD$ID
  caseInd_testD <- substr(as.character(caseInd_testD), 1, nchar(as.character(caseInd_testD))-5)
  caseInd_unique_testD <- caseInd_testD[!duplicated(caseInd_testD)]
  
  if(imageMod == "PET"){
    caseInd_testD_indep <- testD_indep$ID
    caseInd_testD_indep <- substr(as.character(caseInd_testD_indep), 1, nchar(as.character(caseInd_testD_indep))-5)
    caseInd_unique_testD_indep <- caseInd_testD_indep[!duplicated(caseInd_testD_indep)]
    }
  
  
  Cohort_median <- colMedians(as.matrix( trainD[,3:ncol(trainD)] ), 
                              na.rm = FALSE, 
                              hasNA = FALSE, 
                              keep.names=TRUE)
  
  
  
  Consensus_trainD <- data.frame(ID = caseInd_unique_trainD, matrix(0, nrow = length(caseInd_unique_trainD), ncol = ncol(trainD)-1 ) )
  colnames(Consensus_trainD) <- colnames(trainD)
  Consensus_trainD$HPV_status <- ordered(Consensus_trainD$HPV_status, levels=c(0, 1))
  
  Consensus_testD <- data.frame(ID = caseInd_unique_testD, matrix(0, nrow = length(caseInd_unique_testD), ncol = ncol(testD)-1) )
  colnames(Consensus_testD) <- colnames(testD)
  Consensus_testD$HPV_status <- ordered(Consensus_testD$HPV_status, levels=c(0, 1))
  
  if(imageMod == "PET"){
    Consensus_testD_indep <- data.frame(ID = caseInd_unique_testD_indep, matrix(0, nrow = length(caseInd_unique_testD_indep), ncol = ncol(testD_indep)-1) )
    colnames(Consensus_testD_indep) <- colnames(testD_indep)
    Consensus_testD_indep$HPV_status <- ordered(Consensus_testD_indep$HPV_status, levels=c(0, 1))  
    }
  
  
  for(j in caseInd_unique_trainD ){
    
    All_ROI <- trainD[ caseInd_trainD %in% j, ]
    HPV_status <- All_ROI$HPV_status[1]
    All_ROI <- t(All_ROI[,3:ncol(All_ROI)])
    
    max_Delta <- abs( All_ROI - matrix(Cohort_median, nrow = length(Cohort_median), ncol = ncol(All_ROI)) )
    
    max_Delta_index <- (apply(max_Delta,1,which.max))
    
    Consensus_trainD[Consensus_trainD$ID == j, 2:ncol(Consensus_trainD)] <- data.frame(HPV_status = HPV_status,
                                                                                       t(All_ROI[ cbind(seq_along(max_Delta_index), max_Delta_index) ]) )
  }
  
  for(j in caseInd_unique_testD ){
    
    All_ROI <- testD[ caseInd_testD %in% j, ]
    HPV_status <- All_ROI$HPV_status[1]
    All_ROI <- t(All_ROI[,3:ncol(All_ROI)])
    
    max_Delta <- abs( All_ROI - matrix(Cohort_median, nrow = length(Cohort_median), ncol = ncol(All_ROI)) )
    
    max_Delta_index <- (apply(max_Delta,1,which.max))
    
    Consensus_testD[Consensus_testD$ID == j, 2:ncol(Consensus_testD)] <- data.frame(HPV_status= HPV_status, 
                                                                                    t(All_ROI[ cbind(seq_along(max_Delta_index), max_Delta_index) ]) )
  }
 
  if(imageMod == "PET"){
    
    for(j in caseInd_unique_testD_indep ){
      
      All_ROI <- testD_indep[ caseInd_testD_indep %in% j, ]
      HPV_status <- All_ROI$HPV_status[1]
      All_ROI <- t(All_ROI[,3:ncol(All_ROI)])
      
      max_Delta <- abs( All_ROI - matrix(Cohort_median, nrow = length(Cohort_median), ncol = ncol(All_ROI)) )
      
      max_Delta_index <- (apply(max_Delta,1,which.max))
      
      Consensus_testD_indep[Consensus_testD_indep$ID == j, 2:ncol(Consensus_testD_indep)] <- data.frame(HPV_status= HPV_status, 
                                                                                                        t(All_ROI[ cbind(seq_along(max_Delta_index), max_Delta_index) ]) )
      }
  }
  
  trainD <- Consensus_trainD
  testD <- Consensus_testD
  if(imageMod == "PET"){ testD_indep <- Consensus_testD_indep }
  
}


trainD <- trainD[sample(nrow(trainD)),]
testD <- testD[sample(nrow(testD)),]
if(imageMod == "PET"){ testD_indep <- testD_indep[sample(nrow(testD_indep)),] }

index <- sample(colnames(trainD)[3:length(colnames(trainD))])

trainD <- data.frame(ID = trainD$ID, 
                     HPV_status = trainD$HPV_status, 
                     trainD[, index ]  )

testD <- data.frame(ID = testD$ID, 
                     HPV_status = testD$HPV_status, 
                     testD[, index ]  )

if(imageMod == "PET"){
  testD_indep <- data.frame(ID = testD_indep$ID, 
                            HPV_status = testD_indep$HPV_status, 
                            testD_indep[, index ]  )
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
      
captureModelOutput <- data.frame(FS_method = replicate(length(probtest_train), FS_identifier),
                                 ML_method = replicate(length(probtest_train), ML_identifier),
                                 ROI = replicate(length(probtest_train), anaMod),
                                 Imaging_Modality = replicate(length(probtest_train), imageMod),
                                 data_set = replicate(length(probtest_train), "trainD"),
                                 ID = SelectTrainD$ID, 
                                 HPV_status = SelectTrainD$HPV_status,
                                 result_type = replicate(length(probtest_train), "probtest_train"), 
                                 result_value_numeric = probtest_train, 
                                 result_value_char = NA)

captureModelOutput <- rbind(captureModelOutput, 
                            data.frame(FS_method = replicate(length(probtest_test), FS_identifier),
                                 ML_method = replicate(length(probtest_test), ML_identifier),
                                 ROI = replicate(length(probtest_test), anaMod),
                                 Imaging_Modality = replicate(length(probtest_test), imageMod),
                                 data_set = replicate(length(probtest_test), "testD"),
                                 ID = SelectTestD$ID, 
                                 HPV_status = SelectTestD$HPV_status,
                                 result_type = replicate(length(probtest_test), "probtest_test"), 
                                 result_value_numeric = probtest_test,
                                 result_value_char = NA))

if(imageMod == "PET"){
captureModelOutput <- rbind(captureModelOutput, 
                            data.frame(FS_method = replicate(length(probtest_test_indep), FS_identifier),
                                 ML_method = replicate(length(probtest_test_indep), ML_identifier),
                                 ROI = replicate(length(probtest_test_indep), anaMod),
                                 Imaging_Modality = replicate(length(probtest_test_indep), imageMod),
                                 data_set = replicate(length(probtest_test_indep), "testD_indep"),
                                 ID = SelectTestD_indep$ID, 
                                 HPV_status = SelectTestD_indep$HPV_status,
                                 result_type = replicate(length(probtest_test_indep), "probtest_test_indep"), 
                                 result_value_numeric = probtest_test_indep, 
                                 result_value_char = NA))
}




captureModelOutput <- rbind(captureModelOutput, 
                            data.frame(FS_method = replicate((ncol(SelectTrainD)-2), FS_identifier), ML_method = replicate((length(colnames(SelectTrainD))-2), ML_identifier), ROI = replicate((length(colnames(SelectTrainD))-2), anaMod), Imaging_Modality = replicate((length(colnames(SelectTrainD))-2), imageMod),
                                 data_set = replicate((ncol(SelectTrainD)-2), NA), ID = NA, HPV_status = NA,
                                 result_type = replicate((ncol(SelectTrainD)-2), "selected_feature"), 
                                 result_value_numeric = NA, 
                                 result_value_char = colnames(SelectTrainD)[3:ncol(SelectTrainD)] ) )



AUC_train_mlt <- as.numeric(auc_roc(probtest_train, SelectTrainD$HPV_status))
pROC_train <- roc(response = SelectTrainD$HPV_status, predictor = probtest_train, levels= c(0, 1), percent=FALSE, na.rm=FALSE, direction= "<", quiet = FALSE, smooth=FALSE, auc=TRUE, ci=TRUE, ci.method="delong", conf.level= 0.95, plot=FALSE)
pROC_train_bootstrap <- roc(response = SelectTrainD$HPV_status, predictor = probtest_train, levels= c(0, 1), percent=FALSE, na.rm=FALSE, direction= "<", quiet = FALSE, smooth=FALSE, auc=TRUE, ci=TRUE, ci.method="bootstrap", boot.n = 2000, boot.stratified = TRUE, conf.level= 0.95, plot=FALSE)
AUC_train_pROC <- as.numeric(pROC_train$auc)
AUC_CI_train_Delong <- as.numeric(pROC_train$ci)
AUC_CI_train_bootstrap <- as.numeric(pROC_train_bootstrap$ci)

AUC_test_mlt <- as.numeric(auc_roc(probtest_test, SelectTestD$HPV_status))
pROC_test <- roc(response = SelectTestD$HPV_status, predictor = probtest_test,levels= c(0, 1), percent=FALSE, na.rm=FALSE,direction= "<", quiet = FALSE,smooth=FALSE, auc=TRUE, ci=TRUE, ci.method="delong",conf.level= 0.95,plot=FALSE)
pROC_test_bootstrap <- roc(response = SelectTestD$HPV_status,predictor = probtest_test,levels= c(0, 1), percent=FALSE, na.rm=FALSE,direction= "<", quiet = FALSE,smooth=FALSE, auc=TRUE, ci=TRUE, ci.method="bootstrap", boot.n = 2000, boot.stratified = TRUE, conf.level= 0.95,plot=FALSE)
AUC_test_pROC <- as.numeric(pROC_test$auc)
AUC_CI_test_Delong <- as.numeric(pROC_test$ci)
AUC_CI_test_bootstrap <- as.numeric(pROC_test_bootstrap$ci)

if(imageMod == "PET"){
AUC_test_indep_mlt <- as.numeric(auc_roc(probtest_test_indep, SelectTestD_indep$HPV_status))
pROC_test_indep <- roc(response = SelectTestD_indep$HPV_status, predictor = probtest_test_indep, levels= c(0, 1), percent=FALSE, na.rm=FALSE, direction= "<", quiet = FALSE, smooth=FALSE, auc=TRUE, ci=TRUE, ci.method="delong", conf.level= 0.95, plot=FALSE)
pROC_test_indep_bootstrap <- roc(response = SelectTestD_indep$HPV_status, predictor = probtest_test_indep, levels= c(0, 1), percent=FALSE, na.rm=FALSE, direction= "<", quiet = FALSE, smooth=FALSE, auc=TRUE, ci=TRUE, ci.method="bootstrap", boot.n = 2000, boot.stratified = TRUE, conf.level= 0.95, plot=FALSE)
AUC_test_indep_pROC <- as.numeric(pROC_test_indep$auc)
AUC_CI_test_indep_Delong <- as.numeric(pROC_test_indep$ci)
AUC_CI_test_indep_bootstrap <- as.numeric(pROC_test_indep_bootstrap$ci)
}

captureModelOutput <- rbind(captureModelOutput, 
                            data.frame(FS_method = replicate(6, FS_identifier), ML_method = replicate(6, ML_identifier), ROI = replicate(6, anaMod), Imaging_Modality = replicate(6, imageMod),
                                       data_set = replicate(6, "trainD"), ID = NA, HPV_status = NA,
                                       result_type = c("AUC_mlt", "AUC_pROC", "AUC_95%_CI_lower_bound_Delong", "AUC_95%_CI_upper_bound_Delong", "AUC_95%_CI_lower_bound_bootstrap", "AUC_95%_CI_upper_bound_bootstrap"), 
                                       result_value_numeric = c(AUC_train_mlt, AUC_train_pROC, AUC_CI_train_Delong[1], AUC_CI_train_Delong[3], AUC_CI_train_bootstrap[1], AUC_CI_train_bootstrap[3]), 
                                       result_value_char = NA))

captureModelOutput <- rbind(captureModelOutput, 
                            data.frame(FS_method = replicate(6, FS_identifier), ML_method = replicate(6, ML_identifier), ROI = replicate(6, anaMod), Imaging_Modality = replicate(6, imageMod),
                                       data_set = replicate(6, "testD"), ID = NA, HPV_status = NA,
                                       result_type = c("AUC_mlt", "AUC_pROC", "AUC_95%_CI_lower_bound_Delong", "AUC_95%_CI_upper_bound_Delong", "AUC_95%_CI_lower_bound_bootstrap", "AUC_95%_CI_upper_bound_bootstrap"), 
                                       result_value_numeric = c(AUC_test_mlt, AUC_test_pROC, AUC_CI_test_Delong[1], AUC_CI_test_Delong[3], AUC_CI_test_bootstrap[1], AUC_CI_test_bootstrap[3]), 
                                       result_value_char = NA))


if(imageMod == "PET"){
  captureModelOutput <- rbind(captureModelOutput, 
                              data.frame(FS_method = replicate(6, FS_identifier), ML_method = replicate(6, ML_identifier), ROI = replicate(6, anaMod), Imaging_Modality = replicate(6, imageMod),
                                         data_set = replicate(6, "testD_indep"), ID = NA, HPV_status = NA,
                                         result_type = c("AUC_mlt", "AUC_pROC", "AUC_95%_CI_lower_bound_Delong", "AUC_95%_CI_upper_bound_Delong", "AUC_95%_CI_lower_bound_bootstrap", "AUC_95%_CI_upper_bound_bootstrap"), 
                                         result_value_numeric = c(AUC_test_indep_mlt, AUC_test_indep_pROC, AUC_CI_test_indep_Delong[1], AUC_CI_test_indep_Delong[3], AUC_CI_test_indep_bootstrap[1], AUC_CI_test_indep_bootstrap[3]), 
                                         result_value_char = NA))
  }



ConfMatrix <- function(cutoff, data, aim){
  
  if(data=="trainD"){probtest = probtest_train 
                     reference = SelectTrainD$HPV_status
  }else if(data=="testD"){probtest = probtest_test
                          reference = SelectTestD$HPV_status
  }else if(data=="testD_indep"){probtest = probtest_test_indep
                                reference = SelectTestD_indep$HPV_status}
  
  ConfMatrix <- confusionMatrix(data = factor(ifelse(probtest>=cutoff, 1, 0), levels = c(0,1)), reference = reference, positive = "1", mode = "everything")
  
  return(rbind(
         data.frame(result_type = aim, 
                    result_value_numeric = cutoff), 
         data.frame(result_type = c(paste(aim, names(ConfMatrix$overall), sep="_"), paste(aim, names(ConfMatrix$byClass), sep = "_")),
                    result_value_numeric = c(as.numeric(ConfMatrix$overall), as.numeric(ConfMatrix$byClass))))
        )
}

captureModelOutput <- rbind(captureModelOutput, 
                            data.frame(FS_method = replicate(19, FS_identifier), ML_method = replicate(19, ML_identifier), ROI = replicate(19, anaMod), Imaging_Modality = replicate(19, imageMod),
                                       data_set = replicate(19, "trainD"), ID = NA, HPV_status = NA,
                                       rbind(ConfMatrix(optimal_cutoff_balanced, "trainD", "balanced_cutoff")),
                                       result_value_char = NA))

captureModelOutput <- rbind(captureModelOutput, 
                            data.frame(FS_method = replicate(19, FS_identifier), ML_method = replicate(19, ML_identifier), ROI = replicate(19, anaMod), Imaging_Modality = replicate(19, imageMod),
                                       data_set = replicate(19, "testD"), ID = NA, HPV_status = NA,
                                       rbind(ConfMatrix(optimal_cutoff_balanced, "testD", "balanced_cutoff")),
                                       result_value_char = NA))

if(imageMod == "PET"){
captureModelOutput <- rbind(captureModelOutput, 
                            data.frame(FS_method = replicate(19, FS_identifier), ML_method = replicate(19, ML_identifier), ROI = replicate(19, anaMod), Imaging_Modality = replicate(19, imageMod),
                                       data_set = replicate(19, "testD_indep"), ID = NA, HPV_status = NA,
                                       rbind(ConfMatrix(optimal_cutoff_balanced, "testD_indep", "balanced_cutoff")),
                                       result_value_char = NA))  }

