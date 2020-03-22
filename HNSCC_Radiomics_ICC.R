
# The purpose of this script is to calculate the INTRA- / INTERCLASS CORRELATION COEFFICIENT (ICC)
# for Radiomics features multiple delineation stability assessment. 


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
install.packages("grid")
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
install.packages("psych")
install.packages("rel")



captureICC <- data.frame()

for(i in colnames(Allradiom_OS11[3:ncol(Allradiom_OS11)])){


  
Feature_for_ICC <- cbind(Allradiom_OS11[,i], Allradiom_OS12[,i])

ICC_1 <- ICC(x = Feature_for_ICC,
             missing=TRUE,
             alpha=.05,
             lmer=FALSE,
             check.keys=FALSE
             )



Feature_for_ICC <- cbind(Allradiom_OS11[,i], Allradiom_OS2[,i])

ICC_2 <- ICC(x = Feature_for_ICC,
             missing=TRUE,
             alpha=.05,
             lmer=FALSE,
             check.keys=FALSE
            )


captureICC <- rbind(captureICC, data.frame(as.character(i), 
                                           
                                           as.numeric(ICC_1$results$`lower bound`[2]), 
                                           as.numeric(ICC_1$results$ICC[2]),
                                           as.numeric(ICC_1$results$`upper bound`[2]),
                                           
                                           as.numeric(ICC_2$results$`lower bound`[2]),
                                           as.numeric(ICC_2$results$ICC[2]),
                                           as.numeric(ICC_2$results$`upper bound`[2])  ) 
                    )

captureICC[,1] <- as.character(captureICC[,1])
captureICC[,2] <- as.numeric(captureICC[,2])
captureICC[,3] <- as.numeric(captureICC[,3])
captureICC[,4] <- as.numeric(captureICC[,4])
captureICC[,5] <- as.numeric(captureICC[,5])
captureICC[,6] <- as.numeric(captureICC[,6])
captureICC[,7] <- as.numeric(captureICC[,7])

}

colnames(captureICC) <- c("Feature_Name", 
                          "INTRArater_ICC_lower_bound", 
                          "INTRArater_ICC", 
                          "INTRArater_ICC_upper_bound", 
                          
                          "INTERrater_ICC_lower_bound",
                          "INTERrater_ICC",
                          "INTERrater_ICC_upper_bound")



