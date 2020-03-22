
# The purpose of this script is to store standardized code for FEATURE SELECTION TECHNIQUES ("FS")
# for use in machine learning classification problems. 


# *****************************************************************************************
# ************************ FS methods *****************************************************
# *****************************************************************************************


###########################################################################################


 
# *************** Minimum Redundancy Maximum Relevance - mRMR  ***************************

FS_identifier <- "MRMR"

trainD_formrmr <- mRMR.data(data=(subset(trainD, select=-c(ID))))

Features <- mRMR.classic("mRMRe.Filter",
                         feature_count = feature_count,
                         data = trainD_formrmr, 
                         target_indices=(which(trainD_formrmr@feature_names=="HPV_status")) 
)



###########################################################################################


# ******************** Mutual Information Maximization Filter ****************************

FS_identifier <- "MIM"


trainD_forPraznik <- trainD[,3:ncol(trainD)]

Features <- MIM(X = trainD_forPraznik, 
                Y = as.factor(trainD$HPV_status), 
                k = k
                )


###########################################################################################


# ***************** Logistic Regression with RIDGE regularization for FS *****************

FS_identifier <- "RIDGE"

RegReg <-  cv.glmnet(# cv.glmnet arguments   
  
                    x= as.matrix(trainD[3:ncol(trainD)])  ,
                    y=as.factor(trainD$HPV_status),
                    type.measure="auc",
                    nfolds=10,
                    family="binomial",
                    alpha=0,
                    standardize = FALSE
                    )



Coef <- predict(RegReg,
                s="lambda.min",
                type="coefficients"
                )  

###########################################################################################



# ******************* Hierarchical Clustering for FS ****************

FS_identifier <- "HClust"

dist <- dist(x = t(trainD[,3:ncol(trainD)]), 
             method = "euclidean", 
             diag = FALSE, 
             upper = FALSE)

hClust <- hclust(
                d = dist,
                method = "ward.D2", 
                members = NULL
                )

CuthClust <- cutree(hClust, 
                    k = k)

trainD_Meta <- data.frame(ID=trainD$ID, 
                          HPV_status=trainD$HPV_status)

testD_Meta <- data.frame(ID=testD$ID, 
                         HPV_status=testD$HPV_status)

if(ScriptApp == "final_train_test" & imageMod == "PET"){
  testD_indep_Meta <- data.frame(ID=testD_indep$ID, 
                                 HPV_status=testD_indep$HPV_status)  }


for(i in 1:k){
  
  
  Meta_Feature_trainD <- data.frame(
                                    rowMeans(
                                        data.frame( 
                                            trainD[, ( names(CuthClust[as.numeric(CuthClust)==i]) )    ]
                                    )))
  
  Meta_Feature_testD <- data.frame(
                                  rowMeans(
                                      data.frame( 
                                          testD[, ( names(CuthClust[as.numeric(CuthClust)==i]) )    ]
                                    )))
  
  if(ScriptApp == "final_train_test" & imageMod == "PET"){
    Meta_Feature_testD_indep <- data.frame(
                                        rowMeans(
                                            data.frame( 
                                                testD_indep[, ( names(CuthClust[as.numeric(CuthClust)==i]) )    ]
                                      )))  }
  
  
  colnames(Meta_Feature_trainD) <- paste(c("Cluster"), i, sep = "_", collapse = NULL)
  colnames(Meta_Feature_testD) <- paste(c("Cluster"), i, sep = "_", collapse = NULL)
  if(ScriptApp == "final_train_test" & imageMod == "PET"){   colnames(Meta_Feature_testD_indep) <- paste(c("Cluster"), i, sep = "_", collapse = NULL)   }
    
  
  trainD_Meta <- cbind(trainD_Meta, Meta_Feature_trainD)
  testD_Meta <- cbind(testD_Meta, Meta_Feature_testD)
  if(ScriptApp == "final_train_test" & imageMod == "PET"){   testD_indep_Meta <- cbind(testD_indep_Meta, Meta_Feature_testD_indep)   }
  
}

###########################################################################################



# ************* PCA-based feature selection proposed by Song et al. ***********************

FS_identifier <- "PCA"

PCA <- prcomp(x = trainD[,3:ncol(trainD)], 
              retx = TRUE, 
              center = FALSE, 
              scale. = FALSE,
              tol = NULL, 
              rank. = NULL )

PCAVar = PCA$sdev^2/sum(PCA$sdev^2)*100

CumSum_PCAVar <- cumsum(PCAVar)


ev <- PCA$rotation

Sev <- ev[,1:m]

C_features <- data.frame(rowSums(abs(Sev)))



###########################################################################################
 
    
# ******* Unsupervised Pearson Correlation redundancy reduction - fast *******************

PearCor <- cor(trainD[,3:ncol(trainD)], 
               method="pearson")

PearCor[upper.tri(PearCor)] <- 0
diag(PearCor) <- 0


p <- 0.9


trainD <- data.frame(ID=trainD$ID, 
                     HPV_status=trainD$HPV_status, 
                     trainD[,(which((!apply(PearCor,2,function(a) any(abs(a) > p)))==TRUE)+2)]
)                          


testD <- data.frame(ID=testD$ID, 
                    HPV_status=testD$HPV_status, 
                    testD[,(which((!apply(PearCor,2,function(a) any(abs(a) > p)))==TRUE)+2)]
) 

if(ScriptApp == "final_train_test" & imageMod == "PET"){
  testD_indep <- data.frame(ID=testD_indep$ID, 
                            HPV_status=testD_indep$HPV_status, 
                            testD_indep[,(which((!apply(PearCor,2,function(a) any(abs(a) > p)))==TRUE)+2)]
  ) }



###########################################################################################


# ******************************** No Feature Selection **********************************

FS_identifier <- "noFS"

###########################################################################################
