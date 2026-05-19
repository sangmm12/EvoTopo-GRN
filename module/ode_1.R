
get_interaction <- function(data, col, reduction = FALSE ){
  
  if (nrow(data)==2) {
    return_obj = list(ind.name = rownames(data)[col],
                      dep.name = rownames(data)[-col],
                      coefficient = cor(t(data))[1,2])
    
  } else{
    data <- t(data); name <- colnames(data)
    y = as.matrix(data[,col])
    x = as.matrix(data[,-col])
    
    n <- ncol(x)
    if (reduction == TRUE) {
      vec <- abs(apply(x, 2, cor, y))
      if (all(is.na(vec))) {
        return_obj = list(ind.name = name[col],
                          dep.name = NA,
                          coefficient = 0)
      } else{
        x = x[,order(vec, decreasing = T)[1:(n/log(n))]]
      }
    }
    
    if ( all(y==0) |  all(y==1) ) {
      return_obj = list(ind.name = name[col],
                        dep.name = NA,
                        coefficient = 0)
    } else{
      ridge_cv <- try(cv.glmnet(x = x, y = y,alpha = 0))
      if ('try-error' %in% class(ridge_cv)) {
        return_obj = list(ind.name = name[col],
                          dep.name = NA,
                          coefficient = 0)
        
      } else{
        ridge_cv <- cv.glmnet(x = x, y = y, type.measure = "mse", nfolds = 10, alpha = 0)
        best_ridge_coef <- abs(as.numeric(coef(ridge_cv, s = ridge_cv$lambda.1se))[-1])
        
        fit <- cv.glmnet(x = x, y = y, alpha = 1, family = "gaussian", type.measure = "mse",
                         penalty.factor = 1/best_ridge_coef,
                         nfolds = 10, keep = TRUE, thresh=1e-10, maxit=1e6)
        lasso_coef <- coef(fit, s = fit$lambda.1se)
        return_obj = list(ind.name = name[col],
                          dep.name = lasso_coef@Dimnames[[1]][lasso_coef@i + 1][-1],
                          coefficient = lasso_coef@x[-1])
        if ( length(return_obj$dep.name)==0 ) {
          tmp = cor(x,y)
          return_obj$dep.name = rownames(tmp)[which.max(abs(tmp))]
          return_obj$coefficient = tmp[which.max(abs(tmp))]*1/3
        }
        
      }
      
    }
    
    
  }
  
  
  return(return_obj)
}




qdODEmod <- function(Time, State, Pars, power_par) {
  nn = length(Pars)
  ind_effect = paste0("alpha","*",names(State)[1])
  dep_effect = sapply(2:nn, function(c) paste0(paste0("beta",c-1),"*",names(State)[c]))
  dep_effect = paste0(dep_effect, collapse = "+")
  all_effect = paste0(ind_effect, "+", dep_effect)
  expr = parse(text = all_effect)
  
  with(as.list(c(State, Pars)), {
    dx = eval(expr)
    dy <- power_par[,1]*power_par[,2]*Time^(power_par[,2]-1)
    dind = alpha*x
    for(i in c(1:(nn-1))){
      tmp = paste0(paste0("beta",i),"*",paste0("y",i))
      expr2 = parse(text = tmp)
      assign(paste0("ddep",i),eval(expr2))
    }
    return(list(c(dx, dy, dind, mget(paste0("ddep",1:(nn-1))))))
  })
}

#' @title least-square fit for qdODE model
#' @importFrom deSolve ode
#' @param pars vector for unknown ODE parameters
#' @param data data contain independent effect as first row and dependent effect
#' @param Time vector of time point
#' @param power_par matrix of power equation parameters for dependent effect
#' @return mean-square error
qdODE_ls <- function(pars, data, Time, power_par){
  n = length(pars)
  power_par = as.matrix(power_par)
  if (n==2) {
    Pars = c(alpha = pars[1], beta1 = pars[2:n])
    power_par = t(power_par)
    State = c(x=data[1,1],y1 = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  } else{
    Pars = c(alpha = pars[1], beta = pars[2:n])
    State = c(x=data[1,1],y = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  }
  out = as.data.frame(ode(func = qdODEmod, y = State, parms = Pars,
                          times = Time, power_par = power_par))
  X = as.numeric(data[1,])
  fit = as.numeric(out[,2])
  ind = as.numeric(out[,(n+2)])
  
  if ( any(fit < 0) | any(ind < 0) ) {
    penalty = sum((fit[fit < 0])^2) + sum((ind[ind<0])^2)
  }
  else{
    penalty = 0
  }
  
  sse = sum(crossprod(X - fit)) + penalty 
  
  #sse = sum(crossprod(X-fit),sum((ind[ind<0])^2))
  
  return(sse)
}


#' @title legendre polynomials fit to qdODE model
#' @importFrom deSolve ode
#' @param pars vector of qdODE parameters
#' @param data dataframe of observed data
#' @param Time vector of time point
#' @param power_par matrix of power equation parameters for dependent effect
#' @param LOP_order scalar of LOP order
#' @param new_time vector produce new defined time point
#' @param n_expand scalar for how many interpolation needed
#' @return list contain legendre polynomials parameters, qdODE values and LOP fitted values
qdODE_fit <- function(pars, data, Time, power_par, LOP_order = 6, new_time = NULL, n_expand = 100){
  n = length(pars)
  if (n==2) {
    Pars = c(alpha = pars[1], beta1 = pars[2:n])
    power_par = t(power_par)
    State = c(x=data[1,1],y1 = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  } else{
    Pars = c(alpha = pars[1], beta = pars[2:n])
    State = c(x=data[1,1],y = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  }
  out = as.data.frame(ode(func = qdODEmod, y = State, parms = Pars,
                          times = Time, power_par = power_par))
  out2 = data.frame(x = out[,1], y = data[1,], y.fit = out[,2],
                    ind = out[,(n+2)], dep = out[,(n+3):(ncol(out))])
  colnames(out2)[4:ncol(out2)] = c(rownames(data)[1], rownames(data)[2:n])
  rownames(out2) = NULL
  
  all_LOP_par = sapply(2:ncol(out2),function(c)get_legendre_par(out2[,c], LOP_order, out2$x))
  
  if (is.null(new_time)) {
    time2 = seq(min(Time), max(Time), length = n_expand)
    out3 = apply(all_LOP_par, 2, legendre_fit, x = time2)
    out3 = cbind(time2, out3)
  } else{
    out3 = apply(all_LOP_par, 2, legendre_fit, x = new_time)
    out3 = cbind(new_time, out3)
  }
  colnames(out3) = colnames(out2)
  result = list(fit = out2,
                predict = data.frame(out3),
                LOP_par = all_LOP_par)
  return(result)
}



#' @title wrapper for qdODE model
#' @importFrom stats optim
#' @param result result from power_equation_fit
#' @param relationship list contain variable selection results
#' @param i scalar for which id used for qdODE solving, must <= nrow
#' @param init_pars scalar for initial parameters
#' @param LOP_order scalar of LOP order
#' @param method scalar of qdODE solving methodm, cuurent only support least square
#' @param new_time vector produce new defined time point
#' @param n_expand scalar for how many interpolation needed
#' @param maxit scalar of Optim iteration setting
#' @return list contain variable selection results and LOP parameters for every row
#' @export
qdODE_all <- function(result, relationship, i, init_pars = 1, LOP_order = 6, method = "ls",
                      new_time = NULL, n_expand = 100, maxit = 1e3){
  Time = as.numeric(colnames(result$power_fit))
  variable = c(relationship[[i]]$ind.name, relationship[[i]]$dep.name)
  data = result$power_fit[variable,]
  
  if (length(variable)<=1) {
    qdODE.est = NA
    result = NA
    return.obj <- append(result, list(ODE.value = NA,
                                      parameters = NA))
  } else{
    power_par = result$power_par[variable,][-1,]
    n = nrow(data)
    pars_int = c(init_pars,relationship[[i]]$coefficient)
    if (method == "ls") {
      qdODE.est <- optim(pars_int, qdODE_ls, data = data, Time = Time, power_par = power_par,
                         method = "L-BFGS-B",
                         lower = c(rep(-10,(length(pars_int)))),
                         upper = c(rep(10,(length(pars_int)))),
                         control = list(trace = TRUE, maxit = maxit))
      
      result <- qdODE_fit(pars = qdODE.est$par,
                          data = data,
                          power_par = power_par,
                          Time = Time)
      return.obj <- append(result, list(ODE.value = qdODE.est$value,
                                        parameters = qdODE.est$par))
    } else{
      qdODE.est <- optim(pars_int, qdODE_ls, data = data, Time = Time, power_par = power_par,
                         method = "L-BFGS-B",
                         lower = c(rep(-10,(length(pars_int)))),
                         #lower = c(0, rep(-10,(length(pars_int))-1)),
                         upper = c(rep(10,(length(pars_int)))),
                         control = list(trace = TRUE, maxit = maxit))
      
      result <- qdODE_fit(pars = qdODE.est$par,
                          data = data,
                          power_par = power_par,
                          Time = Time)
      return.obj <- append(result, list(ODE.value = qdODE.est$value,
                                        parameters = qdODE.est$par))
    }
  }
  return(return.obj)
}


#' @title wrapper for qdODE_all in parallel version
#' @param result result from power_equation_fit
#' @param reduction use n/log(n) dimension reduction
#' @param thread scales for how many threads used
#' @param maxit scalar of Optim iteration setting
#' @return list contain variable selection results and LOP parameters for every row
#' @export
qdODE_parallel <- function(result, reduction = FALSE, thread = 2, maxit = 1e3){
  data = result$original_data
  relationship = lapply(1:nrow(data),function(c)get_interaction(data, c, reduction = reduction))
  core.number <- thread
  cl <- makeCluster(getOption("cl.cores", core.number))
  #cl <- parallel::makeCluster(1, setup_strategy = "sequential")
  
  clusterEvalQ(cl, {require(orthopolynom)})
  clusterEvalQ(cl, {require(deSolve)})
  clusterExport(cl, c("qdODEmod", "qdODE_ls", "qdODE_fit", "qdODE_all","get_legendre_matrix",
                      "get_legendre_par","legendre_fit"), envir=environment())
  result = parLapply(1:nrow(data),function(c) qdODE_all(result = result,
                                                        relationship = relationship,
                                                        i = c,
                                                        maxit = maxit
  ), cl = cl)
  stopCluster(cl)
  names(result) = rownames(data)
  names(relationship) = rownames(data)
  return_obj <- list(ode_result = result,
                     relationship = relationship)
  return(return_obj)
}

#' @title convert qdODE results to plot data
#' @importFrom reshape2 melt
#' @param result list of qdODE all
#' @importFrom stats na.omit
qdODEplot_convert <- function(result){
  data = result$predict
  n = ncol(data)
  colnames(data)[4:n] = c(paste0("ind.",colnames(data)[4]),
                          paste0("dep.",colnames(data)[5:n]))
  
  plot.df = melt(data, id.vars = c("x"))
  
  name = levels(plot.df[,2])
  
  ind.name = name[grep("ind", name)]
  ind.name2 = strsplit(ind.name,split = "\\.")[[1]][2]
  ind.df <- subset(plot.df, plot.df[,2] == ind.name)
  ind.df$type = "ind"
  ind.df$variable = ind.name2
  
  depname = levels(plot.df[,2])[grep("dep",name )]
  dep.df <- subset(plot.df, plot.df[,2] %in% depname)
  dep.df$type = "dep"
  dep.df$variable = sapply(strsplit(as.character(dep.df$variable),"\\."),"[",2)
  
  
  original.df = subset(plot.df, plot.df[,2] == "y")
  original.df$type = "original"
  
  fit.df = subset(plot.df, plot.df[,2] == "y.fit")
  fit.df$type = "fit"
  
  plot.df2 = rbind(ind.df, dep.df,fit.df)
  
  name.df = subset(plot.df2, plot.df[,1] == max(plot.df2[,1]))
  name.df = name.df[-nrow(name.df),]
  name.df[,2][name.df[,2] == "y.fit"] = ind.name2
  
  name.df = name.df[-which(name.df[,4] == "fit"),]
  
  name.df[,1] = name.df[,1]*1.002
  return_obj = list(plot.df2 = plot.df2,
                    name.df = name.df,
                    ind.name2 = ind.name2)
  return(return_obj)
}



