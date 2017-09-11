#!/usr/bin/env Rscript

source('../library_SG_models_test.R')

inject_anomalies2 = function(X, bldg, from, duration, magnitude){
  # X: the slice of the data matrix in which anomalies will be injected
  # bldg: the buildings to inject anomalies into
  # from: inject anomalies starting at time point 'from'
  # duration: duration (in time slots of 2-mins) of anomalies
  # magnitude: magnitude of injected anomalies in terms of the bldg's standard deviation. Eg, magintude = 0.5 will inject
  #            anomalies of size 0.5 * sigma_bldg
  N = dim(X)[1]; k = dim(X)[2]
  stopifnot(duration < N)
  stopifnot(bldg <= k)
  X[(from):(from+duration-1), bldg] = X[(from):(from+duration-1), bldg] + magnitude * sd(diff(X[, bldg]))
  cat("Injecting ", magnitude * sd(diff(X[, bldg])), " Watts", "(actual series has sigma ", sd(X[, bldg]), "Watts) \n")
  return(X) # returns the data with injected anomalies
}

get_granger_based_cluster = function(bldg, granger_input_file, K=3){
  G = read.csv(granger_input_file)
  g = G[,bldg] 
  s = sort(g[which(g< 0.01)], index.return = TRUE)
  if (length(which(g< 0.01)) < K ){
    s = sort(g[which(g < 1)], index.return = TRUE)
    return(which(g<1)[s$ix[1]])
  }else{return(which(g<0.01)[s$ix[1:K]])}
}


get_ar_zscores = function(fit, Ya, ATTACK_END, ANOMALOUS_BLDG, DURATION, offset = 0, verbose = FALSE){
  from_day = 14; to_day = 15
  Ycc_test = (Ya[((offset+from_day)*720+1):((offset+to_day)*720), ])
  Ycc_test.diff = diff(Ycc_test)
  Ycc_test.diff[ATTACK_END, ANOMALOUS_BLDG] = Ycc_test.diff[ATTACK_END+1, ANOMALOUS_BLDG] # voiding the "2nd spike" of the injected attack that appears, as a result of differencing, at the end of the attack.
  capture.output(md_AR <- do_detection(Ycc_test.diff, bldg = ANOMALOUS_BLDG, fit), file='/dev/null')
  # EWMA-based monitoring:
  zscores.thres =  md_AR$zcores
  zscores.thres[which(abs(zscores.thres)>100)] = 100*sign(zscores.thres[which(abs(zscores.thres)>100)]) # Apply a threshold to avoid getting a (+/-)infty zscore
  return(zscores.thres)
}

ar_detection_eval = function(injected_anomalies, zscores, lam, L,  DURATION){
  alerts2 = ewma(zscores = zscores, lam = lam, L = L, two_in_a_row_rule = FALSE) +1  # add +1 to account for differencing
  R2 = single_attack_detection_accuracy(injected_anomalies, alerts2, duration = DURATION)
  FPR = R2$FP / (720-DURATION) # FPR = false positive rate
  TPR = (DURATION - R2$del) / (DURATION) # TPR = true positive rate
  F1 = 2 * R2$prec * R2$rec / (R2$prec + R2$rec)
  return(list("f1"= if (is.nan(F1)) 0 else F1, "FPR"=FPR, "TPR"=TPR,"alerts"=alerts2))
}

get_var_zscores = function(fit, Ya, CLUSTER, INDEX, ATTACK_END, ANOMALOUS_BLDG, DURATION, offset = 0, verbose = FALSE){
  from_day = 14; to_day = 15
  Ycc_test = (Ya[((offset+from_day)*720+1):((offset+to_day)*720), ])
  Ycc_test.diff = diff(Ycc_test)
  Ycc_test.diff[ATTACK_END, ANOMALOUS_BLDG] = Ycc_test.diff[ATTACK_END+1, ANOMALOUS_BLDG] # voiding the "2nd spike" of the injected attack that appears, as a result of differencing, at the end of the attack.
  capture.output(md_VAR <- do_VAR_detection(Ycc_test.diff, cluster = CLUSTER, fitted_model = fit), file='/dev/null')
  zscores.thres = md_VAR$res[,INDEX] / sqrt(m_VAR$Sigma[INDEX,INDEX])
  zscores.thres[which(abs(zscores.thres)>100)] = 100*sign(zscores.thres[which(abs(zscores.thres)>100)]) # Apply a threshold to avoid getting a (+/-)infty zscore
  return(zscores.thres)
}


var_detection_eval = function(injected_anomalies, zscores, lam, L, DURATION, fit){
  alerts1 = ewma(zscores = zscores, lam = lam, L = L, two_in_a_row_rule = FALSE) + 1 + fit$order # add +1 to account for differencing and +m_VAR$order (=p) to account for the fact that we moit the first 'p' data points
  R1 = single_attack_detection_accuracy(injected_anomalies, alerts1, duration = DURATION)
  FPR = R1$FP / (720-DURATION) # FPR = false positive rate
  TPR = (DURATION - R1$del) / (DURATION) # TPR = true positive rate
  F1 =  2 * R1$prec * R1$rec / (R1$prec + R1$rec)
  return(list("f1"= if (is.nan(F1)) 0 else F1, "FPR"=FPR, "TPR"=TPR, "alerts"=alerts1))
}

get_dfm_zscores = function(fit, Yc, Ya, ATTACK_END, INDEX, DURATION, offset = 0, verbose = FALSE){
  Ycc_test = Ya %*% diag(1/apply(Yc, 2, sd)) # Standardize data
  Ycc_test.diff = diff(Ycc_test)
  Ycc_test.diff[ATTACK_END, ANOMALOUS_BLDG] = Ycc_test.diff[ATTACK_END+1, ANOMALOUS_BLDG] # voiding the "2nd spike" of the injected attack that appears, as a result of differencing, at the end of the attack.
  capture.output(md_DFM <- dynamic_factors_detection(Ycc_test.diff, fit = fit, significance_level = 0.05), file='/dev/null')
  zscores.thres = md_DFM$zscoresPerBldg[,INDEX]
  zscores.thres[which(abs(zscores.thres)>100)] = 100*sign(zscores.thres[which(abs(zscores.thres)>100)]) # Apply a threshold to avoid getting a (+/-)infty zscore
  return(zscores.thres)
}

dfm_detection_eval = function(injected_anomalies, zscores, lam, L, DURATION, fit){
  alerts = ewma(zscores = zscores, lam = lam, L = L, two_in_a_row_rule = FALSE) + 1 + fit$p # add +1 to account for differencing and +DFM$p (=p) to account for the fact that we ignore the first 'p' data points
  R1 = single_attack_detection_accuracy(injected_anomalies, alerts, duration = DURATION, verbose = FALSE)
  FPR = R1$FP / (720-DURATION) # FPR = false positive rate
  TPR = (DURATION - R1$del) / (DURATION) # TPR = true positive rate
  F1 =  2 * R1$prec * R1$rec / (R1$prec + R1$rec)
  return(list("f1"= if (is.nan(F1)) 0 else F1, "FPR"=FPR, "TPR"=TPR, "alerts"=alerts))
}

url_KW='../data_granular.csv'

# Example of use:
X = get_data(url_KW)
Y = get_good_buildings(X, from=1, to=12*720*7, beta = 0.01)
sigma.Y = apply(Y, 2, sd)
Y = Y[,which(sigma.Y < 1e2)] # Remove buildings with "anomalously" large variance
n.bldgs = dim(Y)[2]

# Train the DFM model
INDEX=1; DURATION = 30;  OFFSET = 0;
Yc = Y[((OFFSET+10)*720):((OFFSET+14)*720), c(1:n.bldgs)] 
Yc.std = Yc %*% diag(1/apply(Yc, 2, sd)) # Standardize data
Ycc = diff((Yc.std)) # First differences of our data to remove trend
capture.output(m_DFM <- dynamic_factor_modeling(Ycc, r = 30, max.p = 3), file='/dev/null')
# Now we are ready to start the evaluation
ar_alerts = matrix(rep(0, n.bldgs*730), nrow=n.bldgs)
var_alerts = matrix(rep(0, n.bldgs*730), nrow=n.bldgs)
dfm_alerts = matrix(rep(0, n.bldgs*730), nrow=n.bldgs)
iter = 1
pairs = list(c(.53, 3.714))
cat("iter", "model ", "ewma_lam", "ewma_L", "ANOMALOUS_BLDG", "mag", "F1", "FPR", "TPR\n")
for (mag in c(30)){ # Injections are in KW
    for (ANOMALOUS_BLDG in 1:n.bldgs){
    ANOMALY_END = 530;
    # Experiment setting
    CLUSTER =  c(ANOMALOUS_BLDG, get_granger_based_cluster(ANOMALOUS_BLDG,'../granger_matrices/GrangerTest_matrix_from10_to14_granular_data.csv', K=3))
    Yc = Y[((OFFSET+10)*720):((OFFSET+14)*720), CLUSTER] 
    Ycc = diff((Yc)) # First differences of our data to remove trend
    # Train the AR model
    capture.output(m_AR <- arma_modeling(Ycc[, INDEX]), file='/dev/null')
    # Train the VAR model
    maxp = 25
    capture.output(m1 <- VARorder(Ycc, maxp = maxp), file='/dev/null') # order selection
    capture.output(m_VAR <- model_selection(Ycc, m1$bic, maxp,  log.file = paste(c('./demo-logs/var_model_for_bldg_', ANOMALOUS_BLDG,'.log'), collapse = '')), file='/dev/null')
    # Inject anomalies;  clusters to try: CLUSTER = c(11, 12), CLUSTER = c(1, 19)
    injected_anomalies = c((ANOMALY_END-DURATION+1):(ANOMALY_END)) # TODO: for multiple anomalies do: injected_anomalies = list (a1=c(), a2= c(), ... ak = ...) contains k sets of anomalies
    Ya = inject_anomalies2(Y, bldg = ANOMALOUS_BLDG, from = 14*720 + ANOMALY_END-DURATION+1, duration = DURATION, magnitude = mag) 

    ar_zscores = get_ar_zscores(m_AR, Ya, ANOMALY_END, ANOMALOUS_BLDG, DURATION, OFFSET)
    for (e in pairs){
      ewma_lam = e[1]; ewma_L = e[2] # EWMA-based monitoring
      ar_res = ar_detection_eval(injected_anomalies, ar_zscores, ewma_lam, ewma_L,  DURATION) 
      ar_alerts[ANOMALOUS_BLDG, ar_res$alerts] = 1
      cat(iter, "AR ", ewma_lam, ewma_L, ANOMALOUS_BLDG, mag, ar_res$f1, ar_res$FPR, ar_res$TPR, "\n")
    }
    var_zscores = get_var_zscores(m_VAR, Ya, CLUSTER, INDEX, ANOMALY_END, ANOMALOUS_BLDG, DURATION, OFFSET)
    for (e in pairs){
      ewma_lam = e[1]; ewma_L = e[2] # EWMA-based monitoring
      var_res = var_detection_eval(injected_anomalies, var_zscores, ewma_lam, ewma_L, DURATION, m_VAR)
      var_alerts[ANOMALOUS_BLDG, var_res$alerts] = 1
      cat(iter, "VAR ", ewma_lam, ewma_L, ANOMALOUS_BLDG, mag, var_res$f1, var_res$FPR, var_res$TPR, "\n")
    }
    Yc = Y[((OFFSET+10)*720):((OFFSET+14)*720), c(1:n.bldgs)] 
    Ycc_test = Ya[((OFFSET+14)*720+1):((OFFSET+15)*720), c(1:n.bldgs)]
    dfm_zscores = get_dfm_zscores(m_DFM, Yc, Ycc_test, ANOMALY_END, ANOMALOUS_BLDG, DURATION, offset = OFFSET)
    for (e in pairs){
      ewma_lam = e[1]; ewma_L = e[2] # EWMA-based monitoring
      dfm_res = dfm_detection_eval(injected_anomalies, dfm_zscores, ewma_lam, ewma_L, DURATION, m_DFM)
      dfm_alerts[ANOMALOUS_BLDG, dfm_res$alerts] = 1
      cat(iter, "DFM ", ewma_lam, ewma_L, ANOMALOUS_BLDG, mag, dfm_res$f1, dfm_res$FPR, dfm_res$TPR, "\n")
    } 
    iter = iter + 1
  }
}

write.csv(ar_alerts, 'ar_alerts_voiding.csv', row.names = FALSE)
write.csv(var_alerts, 'var_alerts_voiding.csv', row.names = FALSE)
write.csv(dfm_alerts, 'dfm_alerts_voiding.csv', row.names = FALSE)

