
# The purpose of this script is to store standardized  MACHINE LEARNING CLASSIFIERS ("MLC") code.


# *****************************************************************************************
# ************************ Machine Learning Classsifiers **********************************
# *****************************************************************************************


###########################################################################################


# ****************** Logistic Regression with Elastic Net Regularization ******************

ML_identifier <- "ElNet"

ElNet <-  cv.glmnet(
                    x= as.matrix(SelectTrainD[3:ncol(SelectTrainD)])  ,
                    y=as.factor(SelectTrainD$HPV_status),
                    type.measure="auc",
                    nfolds=10,
                    family="binomial",
                    alpha=alpha,
                    standardize = FALSE
                    )

probtest_train<- predict(object = ElNet, 
                         newx = as.matrix(SelectTrainD[,3:ncol(SelectTrainD)]), 
                         s= c("lambda.min"),
                         type = "response")

probtest_test <- predict(object = ElNet, 
                         newx = as.matrix(SelectTestD[,3:ncol(SelectTestD)]), 
                         s= c("lambda.min"),
                         type = "response")


###########################################################################################


# ************************* Random Forest ********************************

ML_identifier <- "RF"


rf = randomForest(x = SelectTrainD[,3:ncol(SelectTrainD)], 
                  y = SelectTrainD$HPV_status, 
                  ntree=ntree,
                  mtry=mtry,
                  maxnodes = maxnodes,
                  importance=TRUE, 
                  norm.votes=TRUE, 
                  do.trace=FALSE,
                  keep.forest=TRUE
)


probtest_train <- predict(object=rf, 
                          newdata=SelectTrainD[,3:ncol(SelectTrainD)], 
                          type="prob",
                          norm.votes=TRUE, 
                          predict.all=FALSE, 
                          proximity=FALSE, 
                          nodes=FALSE)


probtest_test <- predict(object=rf, 
                         newdata=SelectTestD[,3:ncol(SelectTestD)], 
                         type="prob",
                         norm.votes=TRUE, 
                         predict.all=FALSE, 
                         proximity=FALSE, 
                         nodes=FALSE)



###########################################################################################


# ************************* XGBoost ********************************

ML_identifier <- "XGB"

SelectTrainD_forXGB <- xgb.DMatrix(data=as.matrix(SelectTrainD[,3:ncol(SelectTrainD)]),
                                   label=as.matrix(SelectTrainD$HPV_status),
                                   missing=NA)

SelectTestD_forXGB <- xgb.DMatrix(data=as.matrix(SelectTestD[,3:ncol(SelectTestD)]),
                                  label=as.matrix(SelectTestD$HPV_status),
                                  missing=NA)

xgb <- xgb.train(
  
        params = list(
                booster="gbtree", 
                eta=eta,
                gamma=gamma,
                max_depth=max_depth,
                min_child_weight=min_child_weight, 
                subsample=subsample,
                colsample_bytree=colsample_bytree,
                lambda=lambda,
                objective="binary:logistic",
                eval_metric="auc",
                eval_metric="error"
                ), 
                data = SelectTrainD_forXGB,
                nrounds=nrounds,
                watchlist= list(test=SelectTestD_forXGB, train=SelectTrainD_forXGB),
                verbose=0,
                )



probtest_train <- predict(object=xgb, 
                          newdata=as.matrix(SelectTrainD[,3:ncol(SelectTrainD)]))


probtest_test <- predict(object=xgb, 
                         newdata=as.matrix(SelectTestD[,3:ncol(SelectTestD)]))



###########################################################################################


# ************************* Naive Bayes ****************************

ML_identifier <- "NBayes"
naiveB <- naive_bayes(x = SelectTrainD[,3:ncol(SelectTrainD)], 
                      y = SelectTrainD$HPV_status,
                      laplace = 0,
                      usekernel = FALSE
)

probtest_train <- predict(naiveB, 
                          SelectTrainD[,3:ncol(SelectTrainD)], 
                          type = "prob")

probtest_test <- predict(naiveB, 
                         SelectTestD[,3:ncol(SelectTestD)], 
                         type = "prob")


###########################################################################################



# ***************** Support Vector Machine - SVM - with radial kernel*********************

ML_identifier <- "SVM"

SVM <- svm(x = as.matrix(SelectTrainD[,3:ncol(SelectTrainD)]), 
           y = as.factor(SelectTrainD$HPV_status), 
           scale = FALSE, 
           type = "C-classification", 
           kernel = "radial",
           gamma = gamma,
           cost = cost , 
           class.weights = c("0" = 1-sum((as.factor(SelectTrainD$HPV_status)==0))/length(as.factor(SelectTrainD$HPV_status)), "1" = sum((as.factor(SelectTrainD$HPV_status)==0))/length(as.factor(SelectTrainD$HPV_status))),
           probability = TRUE,
           na.action = na.fail 
)

probtest_train <- predict(object = SVM, 
                          newdata = as.matrix(SelectTrainD[,3:ncol(SelectTrainD)]), 
                          decision.values = FALSE,
                          probability = TRUE,  
                          na.action = na.fail)


probtest_test <- predict(object = SVM, 
                         newdata = as.matrix(SelectTestD[,3:ncol(SelectTestD)]), 
                         decision.values = FALSE,
                         probability = TRUE,  
                         na.action = na.fail)

###########################################################################################

