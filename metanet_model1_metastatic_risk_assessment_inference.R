library(xgboost)

# Make sure the model file is in the working directory
xgbmodel <- xgb.load('xgb.model')

# Read a small example cohort of 10 samples from MSKCC
X <- as.matrix(read.delim('input_X.txt',na.strings = c("NaN"),fill = TRUE,row.names = 1))

# Predict metastatic risk
mrisk <- predict(xgbmodel, xgb.DMatrix(X, missing = NA))
mrisk <- as.data.frame(mrisk)
rownames(mrisk) <- rownames(X)
head(mrisk)

# Calculate SHAP value
shap_contrib <- predict(xgbmodel,xgb.DMatrix(X, missing = NA),
                        predcontrib = TRUE,approxcontrib = FALSE)
shap_contrib <- as.data.frame(shap_contrib)
rownames(shap_contrib) <- rownames(X)
head(shap_contrib)