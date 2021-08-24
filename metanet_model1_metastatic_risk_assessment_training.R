library(xgboost)
library(parallel)
library(dplyr)

rm(list=ls())

# xgboost functions
pvsm_cvparam <- function(x,y,foldid,logfilename){
  ns <- nrow(x)
  nf <- ncol(x)

  dtrain <- xgb.DMatrix(x, label = y, missing = NA)
  
  cv_folds <- list()
  fid <- as.integer(as.factor(foldid))
  nfd <- max(fid)
  for (i in 1:nfd){cv_folds[[i]] <- which(fid == i)}

  best_param <- list()
  best_seednumber <- 9999
  best_auc <- 0
  best_auc_index <- 0
  best_cvmod <- 0

  max_round <- 1000

  for (iter in 1:10000) {
    param <- list(objective = "binary:logistic",
                  max_depth = sample(2:6, 1),
                  eta = runif(1, .001, .050),  # Learning rate, default: 0.3
                  subsample = runif(1, .5, 1),
                  colsample_bytree = runif(1, .5, 1),
                  min_child_weight = sample(1:5, 1), # These two are important
                  max_delta_step = sample(1:5, 1),    # Can help to focus error into a small range.
                  scale_pos_weight = sum(y == 0)/sum(y == 1), 
                  lambda = runif(1, .001, .2), # L2
                  alpha = runif(1, .4, .6)) # L1

    seed.number <- sample.int(10000, 1) # set seed for the cv
    set.seed(seed.number)
    nco <- parallel::detectCores() - 1
    bstcv <- xgb.cv(params = param,
                    data = dtrain,
                    nrounds = max_round,
                    nfold = nfd,
                    folds = cv_folds,
                    metrics = 'auc',
                    early_stopping_rounds = 10,
                    nthread = nco,
                    verbose = T)

    max_auc_index <- bstcv$best_iteration
    max_auc <- bstcv$evaluation_log[max_auc_index]$test_auc_mean

    # Print to log.txt
    cat("iter=",iter,
        "ncores=",nco,
        "max_auc_index=",max_auc_index,
        "max_auc=",max_auc,
        "max_depth=",param$max_depth,
        "eta=",param$eta,
        "subsample=",param$subsample,
        "colsample_bytree=",param$colsample_bytree,
        "min_child_weight=",param$min_child_weight,
        "max_delta_step=",param$max_delta_step,
        "lambda=",param$lambda,
        "alpha=",param$alpha,
        "\n",
        sep = "\t",
        file = logfilename,
        append = TRUE)

    if (max_auc > best_auc) {
      best_auc <- max_auc
      best_auc_index <- max_auc_index
      best_seednumber <- seed.number
      best_param <- param
      best_cvmod <- bstcv
    }
  }

  res <- list(bestauc = best_auc, bestparam = best_param, bestcvmod = best_cvmod)
  return(res)
}

pvsm_train <- function(X,y1,params,nrounds){
  # Training all the training samples without cross validation:
  bst <- xgboost(xgb.DMatrix(X,label = y1,missing = NA),
                 params = params,
                 nrounds = nrounds,
                 eval_metric = "auc")
  return(bst)
}

pvsm_pred_shap <- function(bst,X){
  # Testing:
  ypred <- predict(bst, xgb.DMatrix(X, missing = NA))
  ypred <- as.data.frame(ypred)
  rownames(ypred) <- rownames(X)

  # calculate shap
  shap_contrib <- predict(bst,xgb.DMatrix(X, missing = NA),
                          predcontrib = TRUE,approxcontrib = FALSE)
  shap_contrib <- as.data.frame(shap_contrib)
  rownames(shap_contrib) <- rownames(X)

  res <- list(ypred = ypred, shap_contrib = shap_contrib)
  return(res)
}

pvsm_output <- function(ypred,shap,X,outdir,outfilemark){
  ypred$sampleorder <- match(rownames(ypred), rownames(X))
  ypred <- ypred[order(ypred$sampleorder),]
  write.table(data.frame("sampleid"=rownames(ypred),ypred),
              paste0(outdir,'Ypred1_',outfilemark,'.txt'),
              row.names=FALSE,sep = '\t',quote = F)

  shap$sampleorder <- match(rownames(shap), rownames(X))
  shap <- shap[order(shap$sampleorder),]
  write.table(data.frame("sampleid"=rownames(shap),shap),
              paste0(outdir,'SHAP1_',outfilemark,'.txt'),
              row.names=FALSE,sep = '\t',quote = F)
}

# Main function starts here:
# Reading...
indir <- 'xgbin/'
X <- as.matrix(read.delim(paste0(indir,'input_X.txt'),na.strings = c("NaN"),fill = TRUE,row.names = 1))
y1 <- as.matrix(read.delim(paste0(indir,'input_y1.txt'),na.strings = c("NaN"),fill = TRUE,row.names = 1))
foldid <- as.matrix(read.delim(paste0(indir,'input_foldid.txt'),na.strings = c("NaN"),fill = TRUE,row.names = 1))
Xval <- as.matrix(read.delim(paste0(indir,'input_Xval.txt'),na.strings = c("NaN"),fill = TRUE,row.names = 1))
testfold <- max(foldid)

# clinic:1-6, histology:1-9
NF <- list("C" = 6, "CH" = 9, "CHG" = ncol(X))
outmarker <- 'CH'
sf <- NF[[outmarker]]
X <- X[,1:sf]
Xval <- Xval[,1:sf]

outdir <- paste0('xgbout_',outmarker,'/')

# Start training...
# Search for the best parameters set

for (testfold in 1:5){
  logfilename <- paste0(outdir,'log_5cv_',outmarker,testfold,'.txt')
  # Parameter
  resxgb <- pvsm_cvparam(X[foldid != testfold,],
                         y1[foldid != testfold,],
                         foldid[foldid != testfold,],
                         logfilename)
  # Model
  bst <- pvsm_train(X[foldid != testfold,],
                    y1[foldid != testfold,],
                    resxgb$bestparam,
                    resxgb$bestcvmod$niter)
  xgb.save(bst, paste0(outdir,'xgb.pvm.',outmarker,testfold,'.model'))
  
  # Output
  res_test <- pvsm_pred_shap(bst,X[foldid == testfold,])
  if (testfold == 1){ytest <- res_test$ypred}else{ytest <- rbind(ytest, res_test$ypred)}
  if (testfold == 1){shap_test <- res_test$shap_contrib}else{shap_test <- rbind(shap_test, res_test$shap_contrib)}
}
pvsm_output(ytest,shap_test,X,outdir,outmarker)

# TCGA Application
logfilename <- paste0(outdir,'log_all_',outmarker,'.txt')
resxgb <- pvsm_cvparam(X,
                       y1,
                       foldid,
                       logfilename)
bst <- pvsm_train(X,
                  y1,
                  resxgb$bestparam,
                  resxgb$bestcvmod$niter)

xgb.save(bst, paste0(outdir,'xgb.pvm.',outmarker,'.all.model'))

res_train <- pvsm_pred_shap(bst,X)
pvsm_output(res_train$ypred,res_train$shap_contrib,X,outdir,paste0(outmarker,'_mmm'))

res_val <- pvsm_pred_shap(bst,Xval)
pvsm_output(res_val$ypred,res_val$shap_contrib,Xval,outdir,paste0(outmarker,'_TCGA'))
