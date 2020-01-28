
nc <- nc_open('data/SeaFlow_SizeDist_regrid-15-5.nc')

PAR         <- ncvar_get(nc,'PAR')
w_obs       <- ncvar_get(nc,'w_obs')
m           <- ncvar_get(nc,'m')
delta_v_inv <- ncvar_get(nc,'delta_v_inv')
v_min       <- ncvar_get(nc,'v_min')
time        <- ncvar_get(nc,'time')

delta_v <- 1/delta_v_inv
v       <- v_min*2^(0:14*delta_v)

data$w_obs       <- ncvar_get(nc,'w_obs')
data$PAR         <- ncvar_get(nc,'PAR')
data$m           <- ncvar_get(nc,'m')
data$delta_v_inv <- ncvar_get(nc,'delta_v_inv')
data$v_min       <- ncvar_get(nc,'v_min')
data$time        <- ncvar_get(nc,'time')

data$nt     <- ndays*1440/data$dt
data$nt_obs <- length(data$time)

t                 <- (0:data$nt)*data$dt
data$E            <- approx(data$time,data$PAR,xout=t)$y
data$return_prior <- 0
data$obs          <- t(data$w_obs)
data$t_obs        <- data$time
data$w_ini        <- data$w_obs[1,]

ind_obs    <- data$t_obs > 3
data$t_obs <- data$t_obs[ind_obs]
data$obs   <- data$obs[,ind_obs]

if(limit_to_numdays > 0){
    thresh     <- limit_to_numdays*1440
    ind_obs    <- data$t_obs < thresh
    data$t_obs <- data$t_obs[ind_obs]
    data$obs   <- data$obs[,ind_obs]

    data$nt    <- thresh/data$dt
    data$E     <- data$E[t < thresh]
}

if(stride_t_obs > 0){
    data$t_obs <- data$t_obs[seq(1,length(data$t_obs),stride_t_obs)]
    data$obs   <- data$obs[,seq(1,ncol(data$obs),stride_t_obs)]
}

data$nt_obs <- dim(data$obs)[2]
