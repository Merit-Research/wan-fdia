#
# Prereqs
#
# install.packages("MTS")
library(MTS)

#
# General functionality methods
#

get_data=function(url){
  # Function which reads the data matrix
  # input: 
  #   url=url of the csv file with data
  # output:
  #   X=the final data matrix
  # example of use:
  #   X=get_data()
  KW=read.csv(url)
  n=dim(KW)[1]
  p=dim(KW)[2]
  X=matrix(as.numeric(as.matrix(KW[,1:p])),nrow=n)
  colnames(X)=colnames(KW)
  return (X)
}

get_good_buildings = function(X, from, to, beta = 0.05){
  # Exclude buildings with many duplicate values (ie, more than 100*beta%). Many consequent duplicate values is
  # an indication that the smart meter was not reporting values properly
  # example:
  # Y = get_good_buildings(X, 1, 720*30, beta = 0.05)
  X = X[from:to,] #TODO; add a stopifnot() check that asserts that our length N is a multiple of 720
  N <- dim(X)[1]; p <- dim(X)[2]
  bnames = colnames(X)
  good = c()
  for(b in 1:p){
    v = X[,b]
    vs = sort(v); # sort the values
    re = rle(vs) # get the run length encoding (similar to unix's command uniq -c) to help us find duplicates
    if ((max(re$lengths) < beta*length(v)) & (length(which(v<0)) == 0)){
      good = cbind(b, good)
    }
  }
  good = sort(good)
  Y = X[,good]
  colnames(Y) = bnames[good]
  return(Y)
}

view_buildings<-function(url_KW='um_buildings_watts_aggregates_with_names.csv',
                         url_Bnames='buildings_UM_analysis_code/Bnames.csv', days=7, bldgs=c(1:10), savefigs =  FALSE){
  #   Function that displays the consumption time series of some buildings
  #   url_KW=url of the kilowatt usage file
  #   url_Bnames=url of building name file with their corresponding numbers
  #   days=number of days to be used
  #   bldgs=which buildings to see
  # example of use:
  #   view_buildings(days=7, bldgs=c(1,2))
  
  d = 720 # this is the diurnal period
  
  # Load the data: our data starts on 1420174800 UTC aka 1/2/2015, 12:00:00 AM Eastern
  X=get_data(url_KW)
  
  for (b in bldgs){
    if (savefigs==TRUE){
      pdf(paste(c("building", b, ".pdf"), collapse = ""), width=6,height=4)
    }
    if (days==7){
      par(xaxt="n")
    }
    plot.ts(smooth(smooth(X[(10*d):((10+days)*d),b])), col="red", main = colnames(X)[b], ylab = 'Power Consumption (Watts)')
    if (days==7){
      # draw an axis on the bottom 
      par(xaxt="s")
      axis(1, at=seq(0, d*8-1, by = d)  ,labels=c("Mon, 01/12", "", "Wed, 01/14", "", "Fri, 01/16", "", "Sun, 01/18", ""), las=1)
    }
    if (savefigs==TRUE){
      dev.off()
    }
  }
  
}

high_pass_filter <- function(X, type=c("diff", "diff2", "spencer15")){
  # Function that applies a linear high-pass filter to the data in order to remove slow-moving trends
  # X: the Nxp data matrx, p is the number of buildings we have
  # type = {"diff", "spencer.15"}. The type of filters supported
  #        "diff": is just the first difference operator X_t - X_{t-1}
  #        "spencer15": the spencer.15 filter; see Brockwell 1987 (Time Series:Theory and Methods), $1.4, page 18
  
  type <- match.arg(type)
  
  print(type)
  
  if (type=="diff"){
    return(diff(X))
  }
  if (type=="diff2"){
    return((diff(X,difference=2)))
  }
  if (type=="spencer15"){
    require(signal)
    return(X - spencer(X))
  }
} 

explore_acf <- function(X, bldg=1, days=1, deseason = FALSE){
  
  d = 720 # period
  b = bldg
  
  par(mfrow=c(4,1))
  plot.ts(smooth(smooth(X[(10*d):((10+days)*d),b])), col="red", main = colnames(X)[b], ylab = 'Power Consumption (Watts)')
  Yd = high_pass_filter(X[(10*d):((10+days)*d),b],type="diff") 
  #Ys = high_pass_filter(X[(10*d):((10+days)*d),b],type="spencer15") 
  
  if (days>1 && deseason == TRUE){
    print(days)
    plot.ts(diff(Yd,lag=720), col="red", main = colnames(X)[b], ylab = 'Detrended & De-seasoned Power Consumption (Watts)')
  }else{
    plot.ts(Yd, col="red", main = colnames(X)[b], ylab = 'Detrended Power Consumption (Watts)')
  }
  acf(Yd, 720, ylim = c(-0.2, .4), main = "ACF of detrended (with diff) signal") 
  pacf(Yd[8:(length(Yd)-8)], 720, ylim = c(-0.2, .4), main = "Partial-ACF of detrended (with diff) signal") 
  #acf(Ys[8:(length(Ys)-8)], 720, ylim = c(-0.2, .4), main = "ACF of detrended (with Spencer15) signal") # add and subtract 8 to avoid the N values that Spencer15 filtering yields
  #pacf(Ys[8:(length(Ys)-8)], 720, ylim = c(-0.2, .4), main = "Partial-ACF of detrended (with Spencer15) signal") 
}

explore_and_safe_acf <- function(X, bldg=c(1:8), days=1){
  
  d = 720 # period
  
  for (b in bldg){
    
    pdf(paste(c("./figures/acf-building-", b, ".pdf"), collapse = ""))
    par(mfrow=c(4,1))
    #par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
    plot.ts(smooth(smooth(X[(10*d):((10+days)*d),b])), col="red", main = colnames(X)[b], ylab = 'Power Consumption (Watts)')
    Yd = high_pass_filter(X[(10*d):((10+days)*d),b],type="diff") 
    
    if (days>1){
      print(days)
      plot.ts(diff(Yd,lag=720), col="red", main = colnames(X)[b], ylab = 'Detrended & De-seasoned Power Consumption (Watts)')
    }else{
      plot.ts(Yd, col="red", main = colnames(X)[b], ylab = 'Detrended Power Consumption (Watts)')
    }
    acf(Yd, 720, ylim = c(-0.2, .4), main = "ACF of detrended (with diff) signal") 
    pacf(Yd[8:(length(Yd)-8)], 720, ylim = c(-0.2, .4), main = "Partial-ACF of detrended (with diff) signal") 
    dev.off()
  }
}

get_building_profiles = function(X, normalize = TRUE){
  N <- dim(X)[1]; p <- dim(X)[2]
  bnames = colnames(X)
  P = c(); # this will be our array of profiles for each building
  for(b in 1:p){
    tmp = matrix(X[,b], ncol = 720, byrow = TRUE)
    profile = colMeans(tmp)
    if (normalize == TRUE)
      profile = profile/sqrt(sum(profile^2)) # normalizing the profile of each building
    P = rbind(P, profile)
  }
  P = as.matrix(P)
  rownames(P) = bnames
  return(P)
}

cluster_buildings = function(P, k = 10){
  d <- dist(P) # find distance matrix
  hc <- hclust(d) # hierarchical clustering with 'complete' (default) linkage
  plot(hc)
  return(cutree(hc, k))
}

normalize = function(A, normalize=TRUE){
  # Nomralize columns of matrix A
  D = sqrt(colSums(A^2))
  if (normalize==TRUE)
    return(t(t(A)/D))
  if (normalize==FALSE)
    return(A)
}

ewma = function(zscores, lam = 0.84, L = 3.719, two_in_a_row_rule = FALSE){
  #
  # zscores <- an array of zscores
  # (lam, L) <- the EWMA pair; see Lucas and Saccucci
  #
  lambda = lam; sigma_l = sqrt(lambda/(2 - lambda)); THRESHOLD = sigma_l*L
  z = zscores
  if (two_in_a_row_rule == TRUE){
  #cat("I'm in the two-in-a-row rule block now. Greetings!\n")
  # Hypothesis testing
  Sn_1 = 0; Sn = 0; alerts = rep(0,length(z)); two_inrow_flag = 0
  for (t in c(1:length(z))){
    Zt = z[t]
    if (Zt <= 4 && two_inrow_flag == 0){
      # Update EWMA:
      Sn = (1-lambda)*Sn_1 + lambda*Zt
      if (abs(Sn) > THRESHOLD){
        alerts[t] = 1
      }
      two_inrow_flag = 0
    }else if (Zt > 4  &&  two_inrow_flag == 0){
      # Do NOT update EWMA now
      two_inrow_flag = 1
    }else if (Zt <= 4  &&  two_inrow_flag == 1){
      two_inrow_flag = 0
      # Update EWMA:
      Sn = (1-lambda)*Sn_1 + lambda*Zt
      if (abs(Sn) > THRESHOLD){
        alerts[t] = 1
      }
    }else{
      alerts[t] = 1
      two_inrow_flag = 1
      # Update EWMA:
      Sn = (1-lambda)*Sn_1 + lambda*Zt
    }
    Sn_1 = Sn
  }
  return(which(alerts==1)) # returns time indices that alerts were set to 1
  }
  # Otherwise, execure regular EWMA without the two-in-row rule
  S = rep(0,length(z))
  for (i in 2:length(S)) {
    S[i] = (1 - lambda)*S[i-1] + lambda*z[i]
  }
  return(which(abs(S) > (sigma_l*L))) # returns time indices that out-of-control data points are detected
}

inject_anomalies = function(X, bldg, from, duration, magnitude){
  # X: the slice of the data matrix in which anomalies will be injected
  # bldg: the buildings to inject anomalies into
  # from: inject anomalies starting at time point 'from'
  # duration: duration (in time slots of 2-mins) of anomalies
  # magnitude: magnitude of injected anomalies in terms of the bldg's standard deviation. Eg, magintude = 0.5 will inject
  #            anomalies of size 0.5 * sigma_bldg
  N = dim(X)[1]; k = dim(X)[2]
  stopifnot(duration < N)
  stopifnot(bldg <= k)
  X[(from):(from+duration-1), bldg] = X[(from):(from+duration-1), bldg] + magnitude * sd(X[, bldg])
  return(X) # returns the data with injected anomalies
}

inject_anomalies2 = function(X, bldg, from, duration, magnitude){
  # X: the slice of the data matrix in which anomalies will be injected
  # bldg: the buildings to inject anomalies into
  # from: inject anomalies starting at time point 'from'
  # duration: duration (in time slots of 2-mins) of anomalies
  # magnitude: magnitude of injected anomalies -- constant shift
  N = dim(X)[1]; k = dim(X)[2]
  stopifnot(duration < N)
  stopifnot(bldg <= k)
  X[(from):(from+duration-1), bldg] = X[(from):(from+duration-1), bldg] + magnitude 
  #cat("Injecting ", magnitude, " Watts", "(actual series has sigma ", sd(X[, bldg]), "Watts) \n")
  return(X) # returns the data with injected anomalies
}

inject_anomalies3 = function(X, bldg, from, duration, magnitude){
  # X: the slice of the data matrix in which anomalies will be injected
  # bldg: the buildings to inject anomalies into
  # from: inject anomalies starting at time point 'from'
  # duration: duration (in time slots of 2-mins) of anomalies
  # magnitude: magnitude of injected anomalies in terms of the bldg's ROBUST standard deviation (R function mad()). Eg, magintude = 0.5 will inject
  #            anomalies of size 0.5 * robust_sigma_bldg
  N = dim(X)[1]; k = dim(X)[2]
  stopifnot(duration < N)
  stopifnot(bldg <= k)
  X[(from):(from+duration-1), bldg] = X[(from):(from+duration-1), bldg] + magnitude * mad(X[, bldg])
  return(X) # returns the data with injected anomalies
}

single_attack_detection_accuracy = function(injected_anomalies, alerts, duration=30, verbose = FALSE){
  true_alerts = intersect(injected_anomalies, alerts)
  false_alerts = setdiff(alerts, true_alerts)
  if (verbose!=FALSE){
    cat("__Debug__All alerts: \n"); print(alerts)
    cat("__Debug__True alerts: \n"); print(true_alerts)
    cat("__Debug__False alerts: \n"); print(false_alerts)
  }
  detect_delay = duration
  if (length(true_alerts) > 0){
    #Tp = 1; Fn = 0 # true positives = 1; false negatives = 0
    detect_delay = true_alerts[1] - injected_anomalies[1]
    Tp = duration - detect_delay; Fn = duration - Tp # (= detect_delay)
  }else{
    #Tp = 0; Fn = 1  # true positives = 0; false negatives = 1
    Tp = 0; Fn = duration  
  }
  Fp = length(false_alerts)
  
  if (Tp==0){
    prec = 0
  }else{
    prec = Tp / (Tp + Fp)
  }
  
  return(list("prec" = prec, "rec" = Tp / (Tp + Fn), "del" = detect_delay, "FP" = Fp))
}

regularize_cov=function(X){ # covariance matrix estimation with regularization. See http://perso.ens-lyon.fr/patrick.flandrin/LedoitWolf_JMA2004.pdf
  n=dim(X)[1]
  p=dim(X)[2]
  S=cov(X)
  I=diag(p)
  u=frob_norm(S,I)
  d=frob_norm(S,u*I)
  val=NULL
  for (k in 1:n){
    mx=X[k,]%*%t(X[k,])-S
    val[k]=frob_norm(mx,mx)
  }
  b.tmp=sum(val)/(n*n)
  b=min(b.tmp,d)
  a=d-b;
  S1=(b*u/d)*I+(a/d)*S
}

#
# Univariate (AR) methods
#


goodness_of_fit <- function(X, bldg=c(1:8), days=1){
  require(forecast)
  d = 720 # period
  for (b in bldg){
    pdf(paste(c("./figures/tsdiag-", b, ".pdf"), collapse = ""))
    Yd = high_pass_filter(X[(10*d):((10+days)*d),b],type="diff")
    fit = auto.arima(Yd)
    tsdiag(fit)
    dev.off()
  }
}

ar_predictions <- function(X, bldg=1, train_days, test_days){
  d = 720
  days = train_days + test_days
  Yd = high_pass_filter(X[(10*d):((10+days)*d),bldg],type="diff") # remove trend
  Ydd = diff(Yd,lag=d) # remove diurnal pattern (periodicity)
  fit_ar = ar(Ydd[1:((train_days-1)*d)], aic=TRUE, order.max = 30)
  # perform the actual predictions on the test set
  newdata = Ydd[((train_days-1)*d):length(Ydd)]
  ts.plot(newdata)
  pfit = predict(fit_ar, newdata, n.ahead = test_days*d, se.fit = TRUE)
  return(list("D"=newdata, "fit"=fit_ar, "P"=pfit))
}

arma_modeling <- function (data, bldg = bldg, max.p = 20, plot.tsdiag = FALSE){
  Yd_train = data
  
  # Try several AR(p) models up to order max.p
  models <- vector("list", max.p)
  for (p in 1:max.p){
    fit <- arima(Yd_train, order = c(1,0,0), include.mean = FALSE) # This is the "backup" model that can be used if we ever encounter an error in the arima() call
    fit <- tryCatch({
      cat("Fitting model with p = ", p, " ...")
      fit <- arima(Yd_train, order = c(p,0,0), include.mean = FALSE) # since we have a differenced series, no need for intercept term -> inldude.mean = FALSE
    }, warning = function(war) {
      # warning handler picks up where error was generated
      print(paste("MY_WARNING:  ",war))
    }, error = function(err) {
      # error handler picks up where error was generated
      print(paste("MY_ERROR:  ",err))
      return(fit)
    }, finally = {
      cat("DONE!\n")
    }) # END tryCatch  
    models[[p]] = fit
  }
  stopifnot(length(models)==max.p)
  
  # Show AIC and BIC of above models in a Table
  models.eval = c()
  for (i in 1:length(models)){ 
    m = models[[i]]; 
    models.eval <- c(models.eval, BIC(m))
    cat("AR(", i, "): AIC = ", AIC(m), " BIC = ", BIC(m), "\n")
  }
  s = sort(models.eval, index.return = 1)
  
  # Elect the best model based on BIC -- the selected model should pass the Box.test criterion
  best = -1
  for (i in 1:length(models)){
    cat("Trying AR model that ranked ", i, " with BIC ", s$x[i], "...\n")
    best = i
    m = models[[s$ix[i]]]
    for (lag in 1:20){
      m.lb.test = Box.test(m$residuals, type = "Ljung-Box", lag=lag)
      if (m.lb.test$p.value < 0.05){
        cat("Goodness-of-fit for model that ranked  ", i, "fails at lag ", lag, "with p-value ", m.lb.test$p.value, "\n")
        best = -1
        break
      }
    }
    if (best != -1){
      break # best is positive, which means a model passed all goodness-of-fit tests
    }
  }
  
  if (best > 0){
    cat("We have a good model! Horay! This is model that ranked ", best, " and is an AR(",m$arma[1],") model.\n")
    best.model = m
  }else{
    cat("All models failed goodness-of-fits tests. Using the one with higher BIC anyways...\n")
    best = 1
    best.model = models[[s$ix[1]]]
  }
  
  # At this time, we should have an AR(best) model to perform predictions; predict and compare with test dataset
  if (plot.tsdiag==TRUE){
    pdf(paste(c("./tsdiag-bestfit.pdf"), collapse = ""))
    tsdiag(best.model)
    dev.off()
  }
  
  #Note: To forecast using the same parameters on different data, you might try "refitting" 
  #      the same model on new data but fix the parameters 
  #      (using the fixed argument to arima()) at the values you estimated on a different data set.
  return(list("fit"=best.model)) # we return the fitted model 
}

do_forecast <- function (X, fitted_model, offset = 0, from_day, to_day){
  # fitted_model <- list returned by arma_modeling() above
  # from_day <- do forecasts from this day
  # to_day  <- do forecasts until this day
  
  fit = fitted_model$fit
  bldg = fitted_model$b
  stopifnot(from_day < to_day)
  from = (offset + from_day) *720; to = (offset + to_day) * 720; 
  Y_test = (X[from:to, bldg]);  X_test = X[from:to, bldg];
  x0  = X_test[1] # we need x0 for doing predictions
  
  #Note: To forecast using the same parameters on different data, you might try "refitting" 
  #      the same model on new data but fix the parameters 
  #      (using the fixed argument to arima()) at the values you estimated on a different data set.
  forecast1 = arima(Y_test, order=c(fit$arma[1], 0, 0), fixed = fit$coef, include.mean = FALSE)
  
  # See Brockwell, Chapter 9, Subsection " Forecasting ARIMA models":
  # pred(X_n+1) = X0 + Y1 + ... + Yn + pred(Y_n+1). The following lines are just a fancy way of getting
  # all the pred(X_n+1) for all n
  Y_pred = Y_test - forecast1$residuals # we assume that R defines residual = actual_value - estimate
  Y_obs = as.vector(c(0, Y_test[1:(length(Y_test) -1)]))
  X_pred = (x0 + cumsum(Y_obs)) + Y_pred
  
  cat("RMSE for building  =", bldg, " (diff'ed series) is: ", sqrt(mean((forecast1$residuals)**2)), "\n")
  
  # do some plotting
  X_test =  X_test[2:length(X_test)]
  level = 0.05; quantile = qnorm(1 - level/2) # for 95% conf intervals, quantile is about 1.96
  outliers = X_test[which( (X_test < X_pred - quantile * sqrt(fit$sigma2)) | (X_test > X_pred + quantile * sqrt(fit$sigma2))   )] 
  cat("Outliers = ", length(outliers), "\n")
  plot(c(1:length(X_pred)), smooth(X_pred)/1e3, xlab = 'Time (2-min intervals)', ylab = 'Power Consumption (KWatts)', col = "red", type = "l", ylim=c(0.8*min(X_pred/1e3), 1.2*max(X_pred/1e3)))
  x = c(1:length(X_pred)); y2 = smooth(X_pred + quantile * sqrt(fit$sigma2)); y1 = X_pred - quantile * sqrt(fit$sigma2)
  polygon(c(x,rev(x)),c(y2/1e3,rev(y1/1e3)),col="grey", border = NA)
  lines(c(1:length(X_pred)), smooth(X_pred)/1e3, col = "black", lwd=1.5, type = "l", ylim=c(0.9*min(X_pred/1e3), 1.1*max(X_pred/1e3)))
  points(which((X_test < X_pred - quantile * sqrt(fit$sigma2)) | (X_test > X_pred + quantile * sqrt(fit$sigma2))),
         outliers/1e3, pch=20, col="darkred", cex=0.8)
  
  # calculate RMSE
  cat("RMSE for building =", bldg, "is: ", sqrt(mean((X_pred - X_test)**2)), "\n")
  return(X_pred)
}

do_detection <- function (X, bldg = bldg, fitted_model,significance_level=0.05){
  # fitted_model <- list returned by arma_modeling() above
  # from_day <- do forecasts from this day
  # to_day  <- do forecasts until this day
  # significance_level <- the level of significance; eg, 0.05. 
  
  fit = fitted_model$fit
  Y_test = X[, bldg]
  
  #Note: To forecast using the same parameters on different data, you might try "refitting" 
  #      the same model on new data but fix the parameters 
  #      (using the fixed argument to arima()) at the values you estimated on a different data set.
  forecast1 = arima(Y_test, order=c(fit$arma[1], 0, 0), fixed = fit$coef, include.mean = FALSE)
  
  p_vals = pnorm(abs((forecast1$residuals)/sqrt(fit$sigma2)), lower.tail = FALSE) # get the upper tail aka p-value
  outliers = p_vals[which(p_vals < (significance_level/2))] # significance_level is usually 0.05; need to divide by two because Gaussian is symmetric
  
  return(list("res" = forecast1$residuals, "p_vals"=p_vals, "outliers"=outliers, "zcores" = forecast1$residuals/sqrt(fit$sigma2)))
}

ords_hist = function(X, train_days = 2, max.p = 15){
  # Returns (and saves) a histogram of the selected AR(p) orders
  B = dim(X)[2]
  ords = c()
  for (i in c(1:B)){
    cat("Fitting building", i, "\n")
    Yd = diff(X[(3*720):(5*720), i])
    capture.output(fit <- arma_modeling(Yd, max.p = max.p))
    ords = rbind(ords, fit$fit$arma[1])
  }
  hist(ords, plot=TRUE, xlab = "Selected AR order") 
  return(ords)
}


#
# VAR methods
#

mq2 = function (x, lag = 24, adj = 0) 
  #
  # Adjusted from R. Tsay's MTS package. Added return values.
  #
  # Borrowing the MTS:mq() function to perform model checking. Our motification is simply need to return the p-values of the test
{
  if (!is.matrix(x)) 
    x = as.matrix(x)
  nr = nrow(x)
  nc = ncol(x)
  g0 = var(x)
  ginv = solve(g0)
  qm = 0
  QM = NULL
  df = 0
  for (i in 1:lag) {
    x1 = x[(i + 1):nr, ]
    x2 = x[1:(nr - i), ]
    g = cov(x1, x2)
    g = g * (nr - i - 1)/(nr - 1)
    h = t(g) %*% ginv %*% g %*% ginv
    qm = qm + nr * nr * sum(diag(h))/(nr - i)
    df = df + nc * nc
    dff = df - adj
    mindeg = nc^2 - 1
    pv = 1
    if (dff > mindeg) 
      pv = 1 - pchisq(qm, dff)
    QM = rbind(QM, c(i, qm, dff, pv))
  }
  pvs = QM[, 4]
  #dimnames(QM) = list(names(pvs), c("  m  ", "    Q(m) ", "   df  ",  " p-value"))
  #cat("Ljung-Box Statistics: ", "\n")
  #printCoefmat(QM, digits = 3)
  return(list(qm = QM[,2], pvs = QM[, 4]))
}


model_selection = function(data, criterion, maxp, log.file = "/dev/null", refine = TRUE){
  Ycc = data
  s = sort(criterion, index.return = 1)
  # Elect the best model based on BIC -- the selected model should pass the Box.test criterion
  best = -1; maxp= maxp;
  for (i in 1:maxp){
    cat("Trying VAR model that ranked ", i, " with BIC ", s$x[i], "...\n")
    best = i
    # NOTE: The function MTS:VAR() is estimating the parameters using a Least-Squeares method. Type 'VAR' and press enter to see the code!
    m_tmp = VAR(Ycc,  s$ix[i]-1, output = FALSE, include.mean =FALSE)  # we need to subtract -1 from s$ix[i] because the function VARorder also considers order 'p' with p=0. 
    if (refine==TRUE){
      capture.output(m <- refVAR(m_tmp,thres=1), file =  "/dev/null" ) # refine the model to remove insignificant Phi parameters. Output for this function always goes to /dev/null
    }else{
      m = m_tmp
    }
    nr = dim(Ycc)[1]; nc = dim(Ycc)[2]
    lag = floor(sqrt(nr));  m.mq.test = mq2(m$residuals, lag=lag, adj = m$order * nc * nc - length(which(m$Phi==0))) # adjusting the df by the number of non-zero AR parameters
    cat("Test fails at ", length(which(m.mq.test$pvs < 0.05)), " lags\n")
    if (length(which(m.mq.test$pvs < 0.05)) > 1){ # if it fails only at 1 lag, we still accept the model
      best = -1
    }
    if (best != -1){
      break # best is positive, which means a model passed all goodness-of-fit tests
    }
  }
  
  if (best > 0){
    cat("We have a good model! Horay! This is model that ranked ", best, " and is a VAR(",m$order,") model.\n")
    best.model = m
  }else{
    cat("All models failed goodness-of-fits tests. Using the BIC one that ranked top anyways...\n")
    sink(log.file) # logging this to a file too
    cat("All models failed goodness-of-fits tests. Using the BIC one that ranked top anyways...\n")
    sink()
    set.order = s$ix[1]-1
    if (set.order == 0){set.order = 1} # This guards against the case of trying to fit a VAR(0)
    best.model = VAR(Ycc,  set.order, include.mean =FALSE)  
  } 
  return(best.model)
}

do_VAR_detection <- function (Y, cluster, fitted_model, significance_level=0.05){
  # X: the input data (should be differenced/detrended)
  # cluster: the set of buildings in the VAR model
  # fitted_model <- selected VAR model, eg. via the model_selection() function above
  # significance_level <- the level of significance; eg, 0.05. 
  # Example usage:
  # Yc = Y[1:(5*720), c(1,2)]; Ycc = diff((Yc))
  # m1 = VARorder(Ycc, maxp = 20) # order selection. Per Tsay p.66 (remark), VARorder calculates the Maximum Likelihood (ML) function in order to get the AIC,BIC,HQ statistics
  # m = model_selection(Ycc, m1$bic, maxp) 
  # res = do_VAR_detection(diff(Y), cluster = c(1,2), fitted_model = m, from_day = 5, to_day = 6)
  
  d = 720 # we have 720 observations in a day
  m = fitted_model
  bldg = cluster
  Y_test = Y[, bldg];  N = dim(Y_test)[1]; k = dim(Y_test)[2] 
  
  #Note: Use the m$Phi matrix coefficients for predictions
  p = m$order
  forecast = matrix(0, k, N-p) # creating an empty matrix of dimension  k x (N-p) to hold the predictions
  Phi = array(m$Phi, dim = c(k,k,p)) # Recall that m$Phi is a k-by-kp matrix of VAR coefficients in the form Phi=[Phi1, Phi2, ..., Phip]
  for (i in 1:(N-p)){
    pred = 0 # becuase we set include.mean=False, there is no m$Ph0. Otherwise, set pred = m$Ph0
    Y_o = t(Y_test[rev((i):(i+p-1)),]) # Y_o = [Y_{t-1}, [Y_{t-2}, ..., [Y_{t-p}] is a (kxp) matrix of past observations
    for (j in 1:p){
      if (p == 1){ # we have a VAR(1) model
        pred = pred + Phi[,,j] %*% t(Y_o) # in this case, Y_o is just a column vector     
      }else{
      pred = pred + Phi[,,j] %*% Y_o[,j]   
      }
    }
    forecast[, i] = pred # i runs from 1 to N-p
  }
  
  forecast =  t(forecast) # transpose to change the shape to (N-p) x k
  residuals = Y_test[(p+1):N,] - forecast
  
  r2 = rep(0, times = N-p)
  p_vals = rep(0, times = N-p)
  for (i in c(1:(N-p))){
    r2[i] = residuals[i,] %*% solve(m$Sigma) %*% residuals[i,]
    p_vals[i] = pchisq(r2[i], df = k, lower.tail = FALSE)  # get the upper tail aka p-value
  }
  
  outliers = p_vals[which(p_vals < (significance_level))] # significance_level is usually 0.05
  
  return(list("res"=residuals, "p_vals"=p_vals, "outliers"=outliers, "zscores" = qnorm(p_vals, lower.tail = TRUE)))
}

#
# Below we have functions for Dynamic Factor Models 
#

synthetic = function(predictors=100, measurements = 30*720, factors = 2, sigma_v = 1, b = 0.2){
  # Create synthetic data. See Bai and Ng, CONFIDENCE INTERVALS FOR DIFFUSION INDEX FORECASTS AND INFERENCE FOR FACTOR-AUGMENTED REGRESSIONS
  r = factors; N = predictors; T = measurements
  L = c() # generate the loadings matrix Nxr
  # set.seed(0) # Uncomment for Testing!
  for (i in 1:N){
    lambda_i = 10*runif(r) 
    L  = rbind(L, lambda_i)
  }
  # set.seed(0) # Uncomment for Testing!
  v = matrix(rnorm(T*N, sd = sigma_v), nrow=T)
  if (b!=0){
    row1 = c(); for (j in 1:N){if (j <=40) {row1[j]=b^j} else row1[j]=0}
    Omega = toeplitz(row1)
    C = chol(Omega)
    e = v %*% C # T x N matrix
  }else{
    e = v * 0
  }
  # Generate the factors Ft
  Fact = matrix(rep(0, T * r), nrow = T)
  F0 = runif(r)
  rho  = c(); for (j in 1:r){rho = cbind(rho, 0.8^j)}
  for (t in 1:T){
    for (j in 1:r){
      if (t==1){Fact[t,j] = rho[j] * F0[j] + (1 - rho[j])^2 * rnorm(1, sd = 1*sigma_v)}
      else{Fact[t,j] = rho[j] * Fact[t-1,j] + (1 - rho[j])^2 * rnorm(1, sd = 1*sigma_v)}
    }
  }
  
  X = Fact %*% t(L) + e
  return(list("X" = as.matrix(X), "F" = as.matrix(Fact), "Lam" = as.matrix(L)))
}


synthetic2=function(n,r,u.val,rho,d,L,tau,T){
  # input:
  #  n = dim of Et (error process)
  #  r = dim of Ft (factor process)
  #  u.val = parameter of for uniform distribution [u.val,1-u.val]
  #  rho, L = parameters for the A matrix where A(L)Ft=ut
  #  d, L parameters for the D matrix where D(L)Et=vt
  #  tau = parameter for the covariance matrix of vt
  #  T = time span on which Xt is generated as Xt = Lam Ft + Et
  # ouput:
  #  X=T x n matrix of observations
  # example of use:
  #  X=synthetic2(n=100,r=10,u.val=0.2,rho=0.2,d=0.2,L=4,tau=0.7,T=1000)
  library(MASS)
  X=matrix(0,T,n);  Factors=matrix(0,T,r)
  Lam=matrix(rnorm(n*r,0,1),n,r)
  for (t in 1:T){
    beta=runif(n,u.val,(1-u.val))
    alpha=(beta/(1-beta))*rowSums(Lam^2)
    Tau=matrix(0,n,n)
    for (i in 1:n){
      for (j in 1:n){
        Tau[i,j]=sqrt(alpha[i]*alpha[j])*(tau^(abs(i-j)))*(1-d^2)
      }
    }
    A=(1-rho*L)*diag(r)
    D=(1-d*L)*diag(n)
    mu.u=rep(0,r)
    sig.u=(1-rho^2)*diag(r)
    mu.v=rep(0,n)
    sig.v=Tau
    u=mvrnorm(1,mu.u,sig.u)
    v=mvrnorm(1,mu.v,sig.v)
    F=solve(A)%*%u
    Factors[t,]=F
    E=solve(D)%*%v
    X[t,]=Lam%*%F+E
  }
  return(list("X" = as.matrix(X), "F" = as.matrix(Factors), "Lam" = as.matrix(Lam)))
}

factor_analysis = function(X, r = 2, max.p = 5){
  # X <- the data as an (TxN) matrix
  # r <- the number of factors
  # max.p <- the maximum order to be tried for fitting the VAR(p) model
  T = dim(X)[1]; N = dim(X)[2]
  #Xs = X %*% diag(1/apply(X, 2, sd)) # Standardize data
  S = 1/T*t(X)%*%X # Sample covariance matrix for our data
  s.v.d = svd(S) # perform PCA via SVD: S = U D U' (where V = U in this case)
  P_hat = s.v.d$u[,1:r]; D_hat = s.v.d$d[1:r]
  Lambda_hat = P_hat %*% diag(sqrt(D_hat)) # estimate of the factor loadings
  Psi_d = diag(S - Lambda_hat %*% t(Lambda_hat))
  G_hat = X %*% P_hat %*% diag(sqrt(1/D_hat)) # G_hat.T = diag(sqrt(1/D_hat)) %*%t(P_hat) %*% t(X)
  # Now we fit a VAR(p) model on the estimated G_hat to get the autoregressive parameters
  stopifnot(dim(G_hat)[1] == dim(X)[1])
  max.p = max.p
  m1 = VARorder(G_hat, maxp = max.p) # order selection
  #browser() #for debugging, uncomment 
  m = model_selection(G_hat, m1$bic, maxp = max.p, log.file = "/tmp/factors_var_model_selection.log", refine = TRUE) 
  return(list("Lambda_hat" = Lambda_hat, "VAR" = m, "Psi_d" = Psi_d ))
}

run_kalman = function(data, A, H, Q, R, p, horizon){
  # Runs the Kalman filter algorithm
  # Xk = AXk-1 + w_k with w_k ~ N(0, Q)
  # Zk = HXk + v_k with v_k ~ N(0, R)
  rp = dim(A)[1] # When we work with a VAR(p=1) then rp = r x 1 = r 
  # initial conditions:
  x0 = as.vector(rep(0,rp)); P0 = diag(rp)
  X = c(); Z = data
  T = horizon
  for (k in p:T){
    #Time Update
    if (k==p){x_k.prior = A %*% x0}
    else{x_k.prior = A %*% xk_1.est}
    if (k==p){P_k.prior = A %*% P0 %*% t(A) + Q}
    else{P_k.prior = A %*% Pk_1 %*% t(A) + Q}
    #Measurement Update
    K.k = P_k.prior %*% t(H) %*% solve(H %*% P_k.prior %*% t(H) + R)
    xk.est = x_k.prior + K.k %*% (Z[k,] - H %*% x_k.prior)
    Pk = (diag(rp) - K.k %*% H) %*% P_k.prior
    # Populate matrix X with the kalman-based estimates
    X = rbind(X, t(xk.est))
    # make prev iteration estimates the "new" ones
    xk_1.est =  xk.est; Pk_1 = Pk
  }
  return(list("X" = X, "Pk" = Pk)) #returns estimates of the latent state (the factors, in our case), and the covariance matrix Pk
}

dynamic_factor_modeling = function(data, var.expl = 0, r = 2, max.p = 5){
  # Input: data <- detrended data, sliced at a training window of size T
  # Input: var.expl <- select the top 'r' components that explain at least a 'var.expl' fraction of variance 
  # Input: r <- the number of factors if the user wants to set this value manually and avoid basing this on the percnt. of variance explained. var.expl should be set to 0.
  # Input: max.p <- the maximum order to be tried for fitting the VAR(p) model
  T = dim(data)[1]; N = dim(data)[2]
  S = 1/T * t(data) %*% data; D = svd(S)$d
  cum.D = cumsum(D)/sum(D)
  if (var.expl == 0){
    r = r # setting this manually, after extensive experimentation with screeplot (try: pca = princomp(S); screeplot(pca)) 
  }
  else{r = which(cum.D > var.expl)[1]}
  #browser()
  
  # Learn the factor parameters
  factor.params =  factor_analysis(data, r = r, max.p = max.p)
  A = factor.params$VAR$Phi; Q = factor.params$VAR$Sigma
  H = factor.params$Lambda_hat; R = factor.params$Psi_d
  # Below we build the block matrices needed for Kalman filter and the VAR(p) parameters of factors
  # See "A two-step estimator for large approximate dynamic factor models based on Kalman filtering" by Doz et al., page 10
  A.ncol = dim(A)[2]; p = A.ncol / r
  if (p > 1){
    I = diag(r*(p-1)) 
    Z0 = matrix(rep(0, r*r*(p-1)), ncol = r )
    A = rbind(A ,cbind(I, Z0))
    H = cbind(H, matrix(rep(0, N * r * (p-1)), ncol = r * (p-1)))
    tmp = matrix(rep(0, p*p), ncol = p); tmp[1,1] = 1 
    Q = tmp %x% Q # doing a Kronecker product trick
  }
  # Run the Kalman filter for the "training" period
  kalman.state = run_kalman(data, A, H, Q, diag(R), p, horizon = T)
  # return the fitted model
  return(list("fp"=factor.params, "kp" = kalman.state, "p" = p))
}

dynamic_factors_detection = function(data, fit, significance_level){
  # Inputs:
  # data <- detrended multivariate series (matrix T x N)
  # fit <- fitted DFM model (the output of dynamic_factor_modeling)
  # significance_level <- significance level for declaring outliers/anomalies
  T = dim(data)[1]; N = dim(data)[2]
  A = fit$fp$VAR$Phi; Q = fit$fp$VAR$Sigma
  H = fit$fp$Lambda_hat; R = diag(fit$fp$Psi_d)
  r = dim(A)[1]; 
  # Below we build the block matrices needed for Kalman filter and the VAR(p) parameters of factors
  # See "A two-step estimator for large approximate dynamic factor models based on Kalman filtering" by Doz et al., page 10
  A.ncol = dim(A)[2]; p = A.ncol / r
  if (p > 1){
    I = diag(r*(p-1)) 
    Z0 = matrix(rep(0, r*r*(p-1)), ncol = r )
    A = rbind(A ,cbind(I, Z0))
    H = cbind(H, matrix(rep(0, N * r * (p-1)), ncol = r * (p-1)))
    tmp = matrix(rep(0, p*p), ncol = p); tmp[1,1] = 1 
    Q = tmp %x% Q # doing a Kronecker product trick
  }
  rp = dim(A)[1] # When we work with a VAR(p=1) then rp = r x 1 = r 
  # initial conditions:
  #x0 = fit$kp$X[dim(fit$kp$X)[1],]; P0 = fit$kp$Pk #use the final state of the training period as initial conditions
  x0 = as.vector(rep(0,rp)); P0 = diag(rp)
  X = c(); Z = data; Z_hat = c()
  r2 = rep(0, times = T - (p-1)); p_vals = rep(0, times = T - (p-1)); 
  fcast = matrix(0, nrow = T - (p-1), ncol = N); res = matrix(0, nrow = T - (p-1), ncol = N); zscoresPerBldg = matrix(rep(0, times = (T-(p-1))*N), ncol=N);
  for (k in p:T){
    #Time Update
    if (k==p){x_k.prior = A %*% x0}
    else{x_k.prior = A %*% xk_1.est}
    if (k==p){P_k.prior = A %*% P0 %*% t(A) + Q}
    else{P_k.prior = A %*% Pk_1 %*% t(A) + Q}
    #Measurement Update
    Sigma_z = H %*% P_k.prior %*% t(H) + R
    K.k = P_k.prior %*% t(H) %*% solve(Sigma_z) # K.k = P_k.prior %*% t(H) %*% solve(H %*% P_k.prior %*% t(H) + R)
    xk.est = x_k.prior + K.k %*% (Z[k,] - H %*% x_k.prior)
    Pk = (diag(rp) - K.k %*% H) %*% P_k.prior
    zk.est = H %*% x_k.prior # Measurements estimates zk.est are: zk.est = H * x_k.prior
    Z_hat = rbind(Z_hat, t(zk.est))
    # Populate matrix X with the kalman-based estimates
    X = rbind(X, t(xk.est))
    # make prev iteration estimates the "new" ones
    xk_1.est =  xk.est; Pk_1 = Pk
    # Residuals housekeeping
    residual = as.vector(Z[k-(p-1),] - zk.est) # NOTE: using an offset '(p-1)' to account for the fact that we start at k=p
    fcast[k-(p-1),] =  zk.est
    res[k-(p-1),] = residual  
    r2[k-(p-1)] = residual %*% solve(Sigma_z) %*% residual
    p_vals[k-(p-1)] = pchisq(r2[k-(p-1)], df = N, lower.tail = FALSE)  # get the upper tail aka p-value
    zscoresPerBldg[k-(p-1),] = residual / sqrt(diag(Sigma_z))
  }
  # Calculate outliers and return
  outliers = p_vals[which(p_vals < (significance_level))] # significance_level is usually 0.05
  return(list("fcast" = fcast, "p_vals"=p_vals, "outliers"=outliers, "res" = res, "zscores" = qnorm(p_vals, lower.tail = TRUE), "zscoresPerBldg" = zscoresPerBldg))
}

# Function for GrangerTest
granger_compute=function(Y,BLDG,from_day_tr,to_day_tr,n_period=720, max.p = 10){
  Ycc=diff(Y[(from_day_tr*n_period+1):(to_day_tr*n_period),])
  p=dim(Ycc)[2]
  G=rep(1,p)
  for (i in setdiff(1:p,BLDG)){
    maxp = max.p
    capture.output(m1 <- VARorder(Ycc[,c(BLDG,i)], maxp = maxp), file='NUL') # order selection
    capture.output(m_VAR <- model_selection(Ycc[,c(BLDG,i)], m1$bic, maxp), file='NUL')
    granger_outp=capture.output(GrangerTest(Ycc[,c(BLDG,i)],m_VAR$order,include.mean = F))
    granger_str=unlist(strsplit(granger_outp[[2]]," "))
    G[i]=as.numeric(granger_str[length(granger_str)])
  }
  return (G)
}


