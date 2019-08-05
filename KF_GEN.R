#SMAC
#Lionel Voirol
#July 2019
#Kalman-filter for the sum of any combinations of AR processes, RW, WN and DR

#Load libraries and define parameters ####
library(simts)
library(tidyverse)

#define parameteres used in the function for debugging purposes
n = 500
estimate_model = F
model = AR(phi= 0.82345009, sigma2 = 0.06499131) + AR(phi= -0.22345009, sigma2 = 0.06700869) + RW(5.85e-2)
y = gen_lts(n = n, model = model)
model_to_estimate = NULL


#define kalman_filter general function
kalman_filter = function(model, y, estimate_model = F, model_to_estimate = NULL, method = 'mle'){
  "
  Compute a kalman filter on a state space model which
  can be composed of the sum of AR(), RW(), DR() and WN() processes.
  Note that there can be any ammount of AR() processes while RW(), DR() and WN()
  processes are limited to 1. The function can both estimate the model parameter
  from a given model and apply the kalman filter or directly apply the kalman 
  filter with specified parameters. The function takes as input the
  model which is the defined model and y which is the observed time serie.
  If estimate_model is set to False (by default), the user need to provide a model
  and its parameters (a ts.model class object). If estimate_model is set to True, the user just need to provide
  the selected model to estimate and the estimation method (see function estimate)
  "
  
  #return error if estimate_model == T and model to estimate is empty
  if(estimate_model == T & is.null(model_to_estimate)) stop("No defined model to estimate")

  #extract info from model and build matrices
  #define specific order
  my_order = data.frame('sel_order' = c(1,2,3,4,4), 'process' = c('WN','DR','RW','AR','SIGMA2'))
  
  #If the user provides the model form and parameters
  if(estimate_model == F){
    df = data.frame('process'= model$process.desc, 'val' = model$theta)
  }
  #if user want to estimate a given model
  else if (estimate_model == T){
    #estimate model
    estimated_model = estimate(model = model_to_estimate, Xt = y , method = method)
    #store in dataframe
    df = data.frame('process.desc' = rownames(estimated_model$mod$estimate),
                    'theta' = estimated_model$mod$estimate)
    #rename colnames
    colnames(df) = c('process', 'val')
  }
  
  #join possible processes and model processes, had to do a little trick to avoid warnings
  combined = sort(union(levels(df$process), levels(my_order$process)))
  join_df <- right_join(mutate(my_order, process=factor(process, levels=combined)),
                 mutate(df, process=factor(process, levels=combined)), by = 'process')

  #isolate all processes and values
  
  #wn
  wn_process = join_df %>% arrange(sel_order) %>% filter(sel_order == 1)
  wn_val = wn_process$val
  
  #dr
  dr_process = join_df %>% arrange(sel_order) %>% filter(sel_order == 2)
  dr_process_val = dr_process$val
  
  #rw
  rw_process = join_df %>% arrange(sel_order) %>% filter(sel_order == 3)
  rw_process_val = rw_process$val
  
  #ar
  ar_processes = join_df %>% arrange(sel_order) %>% filter(sel_order == 4)
  dim_ar_processes = dim(ar_processes)[1]
  phi_vec_index = seq(1, dim_ar_processes, 2) #identify the phi parameters in the AR processes
  phi_vec = ar_processes[phi_vec_index,]$val  #extract the phi parameters in the AR processes
  
  #Define trans mat 
  dr_process_length = dim(dr_process)[1]
  dr_vec = rep(1, dr_process_length)
  rw_process_length = dim(rw_process)[1]
  rw_vec = rep(1, rw_process_length)
  trans_mat = diag(x = c(phi_vec, rw_vec, dr_vec))
  
  #state transition
  #trans_mat %* X_t + T * U where the T * U stands for the matrix multiplication for the DR() process
  
  #define T mat
  t_phi_rw = rep(0, length(c(phi_vec, rw_vec)))
  t_dr = rep(1, length(dr_vec))
  T_mat = diag(c(t_phi_rw, t_dr))
  
  #define U vec
  U_vec = c(rep(0, length(c(phi_vec, rw_vec))), dr_process_val)
  
  #initial measurment
  x_0 = y[1]
  
  #number of obs
  n = length(y) 
  
  #initial process_noise_cov_mat defined as the matrix multiplication of the process noise vector
  #only ar processes and rw processes have a wn terms added, dr will have a wn term of 0
  p_length = length(c(phi_vec, rw_vec, dr_vec))
  phi_wn_vec_index = seq(2, dim_ar_processes, 2)
  phi_wn_vec = ar_processes[phi_wn_vec_index,]$val
  dr_wn_vec = t_dr = rep(0, length(dr_vec))
  p_0 = matrix(nrow = p_length, data = c(phi_wn_vec, rw_process_val, dr_wn_vec)) %*% t(matrix(nrow = p_length, data = c(phi_wn_vec, rw_process_val, dr_wn_vec)))
  process_noise_cov_mat = p_0
  
  #define the matrix H 
  total_states_length = length(c(phi_vec, rw_vec, dr_vec))
  H = matrix(rep(1, total_states_length), ncol = total_states_length)
  
  #Define measurment error and if null set to 0
  measurment_error = wn_val
  if(length(measurment_error) == 0){
    measurment_error = 0
  }
  #creation of empty list of matrices and 2 vectors 
  X_t = list() 
  X_h = list()
  P_t = list()
  P_h = list()
  K   = list() #Kalman gain

  #Initialization
  #define initial X_t
  X_t_length = dim(trans_mat)[1]
  x0 = y[1]
  X_t[[1]] = X_h[[1]] = matrix(c(rep(x0, X_t_length)), nrow = X_t_length)
  P_t[[1]] = P_h[[1]] = process_noise_cov_mat
  
  for (k in seq(n-1)){
    # estimate of next state and error given observation untill t
    X_t[[k+1]] = trans_mat %*% X_h[[k]] + T_mat %*% U_vec
    P_t[[k+1]] = trans_mat %*% P_h[[k]] %*% t(trans_mat) + process_noise_cov_mat
    
    #update estimate given measurment t
    #Note that there is no measurment error in this model
    K[[k+1]]   = P_t[[k+1]] %*% t(H) %*% solve(H %*% P_t[[k+1]] %*% t(H) + measurment_error)  #update Kalman gain
    X_h[[k+1]] = X_t[[k+1]] + K[[k+1]] * as.numeric((y[k+1] - H %*% X_t[[k+1]])) #update state prediction
    P_h[[k+1]] = (diag(X_t_length)-K[[k+1]] %*% H) %*% P_t[[k+1]] #update process_noise_cov_mat
    
  }
  
  out = list("X_h" = X_h, "X_t" = X_t, "P_h" = P_h, "K" = K)
  class(out) = "KF"
  return(out)
}

#apply kalman_filter on savingrt

#load data
data("savingrt")

#kalman filter on savingrt with the estimates provided by R. Molinari
my_mod = AR(phi= 0.82345009, sigma2 = 0.06499131) + AR(phi= -0.22345009, sigma2 = 0.06700869) + RW(5.85e-2)

#apply kalman filter on savingrt with all three processes with estimates R. Molinari
my_res = kalman_filter(model = my_mod, y = savingrt)
estimated = sapply(my_res$X_h, FUN = sum) #using sapply to return a vector instead of a list w/ apply
plot(gts(savingrt))
lines(estimated, col = 'darkorange')
my_res$P_h

#estimate parameter for savingrt with gmwm with all three processes
res_gmwm = kalman_filter(estimate_model = T, method = 'gmwm', model_to_estimate =(AR(1) + AR(1) + RW()), y = savingrt)
estimated_gmwm = sapply(res_gmwm$X_h, FUN = sum)
plot(gts(savingrt))
lines(estimated_gmwm, col = 'darkorange')

#estimate parameter for savingrt with rgmwm with all three processes
res_rgmwm = kalman_filter(estimate_model = T, method = 'rgmwm', model_to_estimate =(AR(1) + AR(1) + RW()), y = savingrt)
estimated_rgmwm = sapply(res_rgmwm$X_h, FUN = sum)
plot(gts(savingrt))
lines(estimated_rgmwm, col = 'darkorange')

#only with the 2 AR processes with estimates given by R. Molinari
my_mod_only_ar = AR(phi= 0.82345009, sigma2 = 0.06499131) + AR(phi= -0.22345009, sigma2 = 0.06700869)
res_only_ar = kalman_filter(model = my_mod_only_ar, y = savingrt)
estimated_only_ar = sapply(res_only_ar$X_h, FUN = sum)
plot(gts(savingrt))
lines(estimated_only_ar, col = 'darkorange')


################################## Code JunkYard ##################################

#plot confidence interval
diagsum = function(x){
  sum(diag(x))
}
sum_cov = function(x){
  sum(x[lower.tri(x)])
}

#plot lines and estimates
plot(gts(savingrt))
lines(estimated, col = 'darkorange')

#adding confidence interval
var_vec = sapply(my_res$P_h, FUN = diagsum) + 2 * sapply(my_res$P_h, FUN = sum_cov)
fit_sd = sqrt(abs(var_vec))
n =length(my_res$P_h)
alpha = .05
polygon(x = c(1:n, rev(1:n)),    
        y = c(estimated + qnorm(1-alpha/2)*-fit_sd,
              rev(estimated + qnorm(1-alpha/2)*fit_sd)))    
plot(x=seq(100), y=rep(5,100), type = 'l')

ci_couleur = "#F8766D4D",
state_couleur = "#F8766DFF"
polygon()