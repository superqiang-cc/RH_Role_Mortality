# Project: Regional Variation in the Role of Humidity on City-level Heat-Related Mortality
# Author: Dr.GUO Qiang, The University of Tokyo
# Contact: qiang@rainbow.iis.u-tokyo.ac.jp
# This script is developed and tested in R 4.0.3 and Linux platform, and can be run also in Windows with
# proper environment settings.
# Description:
# This script is used to implement the Random Forest using R

# Note: be careful the index starts from 0 in Python while from 1 in R.

rm(list=ls())

# load libraries
library(reticulate)
library(caTools)
library(randomForest)


# load X (739, 13) and y (739) for random forest

X_org <- read.csv(file=paste0("RF_features.csv"))
X_org <- X_org[, -which((names(X_org) == 'y'))]
X_org <- X_org[, -which((names(X_org) == 'city_name'))]

# remove sh and use rh
X_org <- X_org[, -which((names(X_org) == 'sh_mean'))]
X_org <- X_org[, -which((names(X_org) == 'sh_std'))]
X_org <- X_org[, -which((names(X_org) == 'corr_t_sh'))]

# remove rh and use sh
# X_org <- X_org[, -which((names(X_org) == 'rh_mean'))]
# X_org <- X_org[, -which((names(X_org) == 'rh_std'))]
# X_org <- X_org[, -which((names(X_org) == 'corr_t_rh'))]

y_str_org <- read.csv(file=paste0("TableS4_Guo_etal_1008.csv"))[, 12] # BFI is in column 12

# Regarding Tw, Ts, Twbg, Tswbgt the heavy humid-heat
y_str_org[y_str_org == 'Tw'] <- 'hh'
y_str_org[y_str_org == 'Ts'] <- 'hh'
y_str_org[y_str_org == 'TWBG'] <- 'hh'
y_str_org[y_str_org == 'TsWBG'] <- 'hh'


# Only focus on classify tair and hh
y_str <- y_str_org
X <- X_org
y_str <- subset(y_str, (y_str_org=="Tair")|(y_str_org=="hh")) 
X <- subset(X, (y_str_org=="Tair")|(y_str_org=="hh")) 

y_str <- factor(y_str)


# Repeat the implementing the Random Forest 500 times to ensure the robustness
rp_n = 500

imp_n <- array(NA, dim=c(13, rp_n))
cf_mtx_n <- array(NA, dim=c(2, 2, rp_n))
accuracy_n <- array(NA, dim=c(rp_n))
precision_n <- array(NA, dim=c(2, rp_n))
recall_n <- array(NA, dim=c(2, rp_n))

for (sd in seq(rp_n)){
    set.seed(sd)  # setting seed
    # Splittig data into train and test
    split <- sample.split(y_str, SplitRatio=0.7)
    train_X <- subset(X, subset=split=="TRUE")
    test_X <- subset(X, subset=split=="FALSE")
    train_y_str <- subset(y_str, subset=split=="TRUE")
    test_y_str <- subset(y_str, subset=split=="FALSE")

    # Fitting Random Forest to the train dataset

    cls_RF <- randomForest(x=train_X,
                            y=train_y_str,
                            ntree=500,
                            mtry=4,
                            nodesize=7)

    # cls_RF

    # Predicting the Test set results
    y_pred_v = predict(cls_RF, newdata=test_X)
    # Confusion Matrix
    conf_mtx_v = table(test_y_str, y_pred_v)

    cf_mtx_n[, , sd] <- conf_mtx_v

    accuracy_n[sd] <- (conf_mtx_v[1, 1] + conf_mtx_v[2, 2]) / sum(conf_mtx_v) * 100
    precision_n[1, sd] <- conf_mtx_v[1, 1] / (conf_mtx_v[1, 1] + conf_mtx_v[2, 1]) * 100
    precision_n[2, sd] <- conf_mtx_v[2, 2] / (conf_mtx_v[2, 2] + conf_mtx_v[1, 2]) * 100
    recall_n[1, sd] <- conf_mtx_v[1, 1] / (conf_mtx_v[1, 1] + conf_mtx_v[1, 2]) * 100
    recall_n[2, sd] <- conf_mtx_v[2, 2] / (conf_mtx_v[2, 2] + conf_mtx_v[2, 1]) * 100

    # Plotting model
    # plot(cls_RF)

    # Importance
    imp <- importance(cls_RF)
    imp_ratio <- imp / sum(imp) * 100

    imp_n[, sd] <- imp_ratio

    # Variable importance plot
    # varImpPlot(cls_RF)

    print(paste0(sd, " seed finished."))
}

rownames(imp_n) <- rownames(imp)
print(mean(accuracy_n))
print(rowMeans(imp_n))

# save(imp_n, 
#      cf_mtx_n,
#      accuracy_n,
#      precision_n,
#      recall_n,
#      file=paste0(stats_under_dir, "RF_stats_local.Rdata"))

print('All Finished.')












