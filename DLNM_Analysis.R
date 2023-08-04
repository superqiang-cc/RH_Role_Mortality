# Project: Regional Variation in the Role of Humidity on City-level Heat-Related Mortality
# Author: Dr.GUO Qiang, The University of Tokyo
# Contact: qiang@rainbow.iis.u-tokyo.ac.jp
# This script is developed and tested in R 4.0.3 and Linux platform, and can be run also in Windows with
# proper environment settings.
# Description:
# This script is used to conduc the heat-motality association analysis using 
# Distributed Non-Linear Model

# Note: be careful the index starts from 0 in Python while from 1 in R.

rm(list=ls())

# load libraries
library(dlnm)
library(mixmeta)
library(splines)
library(reticulate)
library(lubridate)
library(MASS)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
# library(ggtext)

# some information
np <- import("numpy")
hsi_nm <- c("at", "tw", "ts", "wbgt", "swbgt", "hx",  "apt", "utci", "hi")
hsi_nm_weight <- c("at", "hi", "utci", "apt", "hx", "swbgt",  "wbgt", "ts", "tw")
max_mean <- c("max", "mean")

# Directories
mcc_dir <- "/"
mcc_info <- "/"
clm_dir <- "/"
main_dir <- "/"
stats_dir <- "/"

# load city and country information
city_info <- read.csv(file=paste0(mcc_info, "mcc_city_info_new.csv"))
country_info <- read.csv(file=paste0(mcc_info, "country_info_new.csv"))
mcc_city_sel <- city_info$city_sel_no + 1  # change python index to R index 
hsi_city_sel <- city_info$coord_sel_no + 1  # change python index to R index 
city_sel_name <- city_info$city_sel
city_full_name <- city_info$city_name_sel
print("cityinfo and countryinfo load finished.")


# load MCC mortality data (dlist_clean, city_otl_rt, city_otl_rm, city_avb_len)
load(paste0(stats_dir, "dlist_clean.Rdata"))  
print("MCC mortality data load finished.")


# load warm season statistics
hsi_warm_hfcsv <- np$load(paste0(stats_dir, "hsi_warm_hfcsv.npz"))$f[["hsi_warm_hfcsv"]]  # (764, 6)
indsummer <- np$load(paste0(stats_dir, "indsummer.npz"))$f[["indsummer"]]  # (764, 14610)
print("Warm season statistics finished.")

# load date matrix
clm_date <- np$load(paste0(stats_dir, "date_mx.npz"))$f[["date_mx"]]  # (14610, 3)

# remove the cities whose mortality has too much missing value or too short
# 764 - 24 = 740
mcc_city_sel <- subset(mcc_city_sel, !city_otl_rm)
hsi_city_sel <- subset(hsi_city_sel, !city_otl_rm)
city_sel_name <- subset(city_sel_name, !city_otl_rm)
city_full_name <- subset(city_full_name, !city_otl_rm)
hsi_warm_hfcsv <- subset(hsi_warm_hfcsv, !city_otl_rm)
indsummer <- subset(indsummer, !city_otl_rm)
dlist_clean <- subset(dlist_clean, !city_otl_rm)

###############################################################################################
# Set the parameters for the DLNM model
# A quasi-Poisson regression is fitted for heat-mortality association
# 4 factors are considered
# Heat stress:      A natural spline (quadratic B-spline) function with two internal knots 
#                   at the 50th and 90th percentile of the warm season 
#                   (consecutive hottest 4 months)
# Time lag:         A natural spline function with two internal knots at equally spaced 
#                   values in the log scale over 10 days of lag
# Seasonality:      A natural spline function with 4 degrees of freedom of day of the year
# Long-term trends: A natural spline function of time with aproximately one knot every 10 years 
# Day of week:      
# An interaction between seasonality and year to allow different seanal trends
###############################################################################################

## q-AIC Function
fqaic <- function(model) {
  loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
  phi <- summary(model)$dispersion
  qaic <- -2*loglik + 2*summary(model)$df[3]*phi
  return(qaic)
}

## q-BIC Function
fqbic <- function(model) {
  loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
  phi <- summary(model)$dispersion
  qbic <- -2*loglik + log(length(model$y))*summary(model)$df[3]*phi
  return(qbic)
}


# Specification of the exposure function
varfun <- "ns"
vardegree <- NULL
varper <- c(50, 90)

# Specification of the lag function
lag <- 10
lagnk <- 2

# Degree of freedom for seasonality
dfseas <- 4

# Degree of freedom for trend
dftrend <- 1

# Model formula
formula <- death ~ cb + dow + ns(doy, df=dfseas):factor(year) + 
           ns(date, df=round(length(unique(year)) / dftrend / 10))


ncoef <- length(varper) + ifelse(varfun == "bs", vardegree, 1)


# Create empty objects to store the qAIC and mis-fit flag
qaicall <- array(NA, dim=c(9, length(mcc_city_sel)))
qbicall <- array(NA, dim=c(9, length(mcc_city_sel)))


for (mm in seq(2)){

    # load HSI data
    clm_hsi <- np$load(paste0(clm_dir, "mcc_hsi_", max_mean[mm], 
                            ".npz"))$f[["mcc_hsi"]]  # (9, 764, 14610)
    clm_hsi <- clm_hsi[, !city_otl_rm, ]

    for (hsi_i in seq(9)){

        clm_hsi_loop <- clm_hsi[hsi_i, ,]
        
        for (ct in seq(length(mcc_city_sel))) {  # length(mcc_city_sel)

            # Match the date, combine mortality and HSI data
            datacity <- dlist_clean[[ct]]

            # Get the climate and mortality date index
            if (datacity$year[1] < 1980){
                clm_index_s <- 1
                mot_index_s <- which((datacity$year == 1980) & 
                                     (datacity$month == 1) & 
                                     (datacity$day == 1))

            } else {
                clm_index_s <- which((clm_date[, 1] == datacity$year[1]) & 
                                     (clm_date[, 2] == datacity$month[1]) & 
                                     (clm_date[, 3] == datacity$day[1]))
                mot_index_s <- 1
            }
            
            if (datacity$year[dim(datacity)[1]] > 2019){
                clm_index_e <- dim(clm_date)[1]
                mot_index_e <- which((datacity$year == 2019) & 
                                     (datacity$month == 12) & 
                                     (datacity$day == 31))
            } else {
                clm_index_e <- which((clm_date[, 1] == datacity$year[dim(datacity)[1]]) & 
                                     (clm_date[, 2] == datacity$month[dim(datacity)[1]]) & 
                                     (clm_date[, 3] == datacity$day[dim(datacity)[1]]))
                mot_index_e <- dim(datacity)[1]
            }
            

            # Cut the data based on date index
            datacity <- datacity[mot_index_s: mot_index_e, ]
            hsicity <- clm_hsi_loop[ct, clm_index_s: clm_index_e]
            slcity <- indsummer[ct, clm_index_s: clm_index_e]
            warmcity <- hsi_warm_hfcsv[ct, ]

            # Check whether the date length matches
            if (dim(datacity)[1] != length(hsicity)) {
                print(paste0(hsi_nm[hsi_i], " ", max_mean[mm], " ", ct, " ", 
                             city_sel_name[ct], ": the date length does not match."))
                next
            }

            datacity$hsi <- round(hsicity, 2)
            datacity$indsummer <- slcity

            # Subset to the summer period
            datasea <- subset(datacity, month %in% warmcity)

            # Define the crossbasis
            argvar <- list(
                fun = varfun,
                knots = quantile(datasea$hsi, varper / 100, na.rm = T),
                Bound = range(datasea$hsi, na.rm = T)
            )
            arglag <- list(knots = logknots(lag, lagnk))

            cb <- crossbasis(datasea$hsi, lag = lag, argvar = argvar, arglag = arglag,
                group = datasea$indsummer)

            # Run the model and obtain predictions
            model <- glm(formula, datasea, family = quasipoisson, na.action = "na.exclude")

            # Initial predict to get the minimum mortality value
            pre_pdt <- crosspred(cb, model, cen=mean(datasea$hsi, na.rm=T), by=0.01)
            cen <- pre_pdt$predvar[which.min(pre_pdt$allRRfit)]

            # Calculate the qAIC and get the mis-fit flag
            qaicall[hsi_i, ct] <- fqaic(model)
            qbicall[hsi_i, ct] <- fqbic(model)


            print(paste0(hsi_nm[hsi_i], " ", max_mean[mm], " ", ct, " ", 
                         city_sel_name[ct], " finished."))
        }


    }
    save(qaicall, qbicall, misfit,
    file=paste0(stats_dir, "dlnm_qaic_", max_mean[mm], "_local.Rdata"))
}


print('All Finished.')