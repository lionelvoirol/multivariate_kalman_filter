#SMAC
#Lionel Voirol
#Summer 2019
#Kalman-filter on savingrt

# A few example using data savingrt
#load data
data("savingrt")

#kalman filter on savingrt with the estimates provided by R. Molinari
my_mod = AR(phi= 0.82345009, sigma2 = 0.06499131) + AR(phi= -0.22345009, sigma2 = 0.06700869) + RW(5.85e-2)

#apply kalman filter on savingrt with all three processes with estimates R. Molinari
my_res = kalman_filter(model = my_mod, y = savingrt)
plot(my_res)

#estimate parameter for savingrt with gmwm with all three processes
res_gmwm = kalman_filter(estimate_model = T, method = 'gmwm', model_to_estimate =(AR(1) + AR(1) + RW()), y = savingrt)
plot(res_gmwm)

#estimate parameter for savingrt with rgmwm with all three processes
res_rgmwm = kalman_filter(estimate_model = T, method = 'rgmwm', model_to_estimate =(AR(1) + AR(1) + RW()), y = savingrt)
plot(res_rgmwm)

#only with the 2 AR processes with estimates given by R. Molinari
my_mod_only_ar = AR(phi= 0.82345009, sigma2 = 0.06499131) + AR(phi= -0.22345009, sigma2 = 0.06700869)
res_only_ar = kalman_filter(model = my_mod_only_ar, y = savingrt)
plot(res_only_ar)
