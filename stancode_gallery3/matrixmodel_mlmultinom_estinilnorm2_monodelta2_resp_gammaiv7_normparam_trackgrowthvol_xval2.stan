data {
    // size variables
    int<lower=0> m;         // number of size classes
    int<lower=0> nt;        // number of timesteps
    int<lower=0> nt_obs;    // number of timesteps with observations
    // model parameters and input
    int<lower=0> dt;        // delta t in minutes
    real<lower=0> E[nt];    // vector of incident radiation values
    real<lower=0> v_min;    // size in smallest size class in um^-3
    int<lower=0> delta_v_inv;  // inverse of delta_v 
    // observations
    real<lower=0,upper=nt*dt>  t_obs[nt_obs]; // the time of each observation
    real<lower=0> obs[m,nt_obs]; // observations
    int<lower=0> obs_count[m,nt_obs]; // count observations
    // for cross-validation
    int<lower=0, upper=1> i_test[nt_obs];
    int prior_only;
}
transformed data {
    int j;
    real v_range;
    real<lower=0> delta_v;
    real<lower=0> dt_days;  // dt in units of days
    real<lower=0> dt_norm;  // dt in units of days and doublings
    real<lower=0> v[m+1];   // vector of (minimum) sizes for each size class
    row_vector[m] v_mid;    // vector of sizes for each size class
    real<lower=0> v_diff[m-1];// vector of size-differences for first m-1 size classes
    int<lower=0> t[nt];     // vector of times in minutes since start 
    int<lower=1, upper=nt> it_obs[nt_obs]; // the time index of each observation
    int n_test = sum(i_test);
    real xi_max = 10.0;

    j = 1 + delta_v_inv; 
    delta_v = 1.0/delta_v_inv;
    dt_days = dt/1440.0;
    dt_norm = dt/(1440.0 * (2^delta_v - 1.0));
    for (i in 1:m+1){
        v[i] = v_min*2^((i-1)*delta_v);
    }
    for (i in 1:m){
        // using geometric mean
        v_mid[i] = sqrt(v[i] * v[i+1]);
    }
    v_range = v_mid[m]-v_mid[1];
    for (i in 1:m-1){
        // difference between the centers for each class
        v_diff[i] = 0.5*(v_mid[i+1]-v_mid[i]);
    }
    // populate time vector
    t[1] = 0;
    for (i in 2:nt){
        t[i] = (t[i-1] + dt);
    }
    // populate index vector it_obs
    for (k in 1:nt_obs){
        for (i in 1:nt){
            if (t_obs[k]>=t[i] && t_obs[k]<t[i]+dt){
                it_obs[k] = i;
                break;
            }
        }
    }
}
parameters {
    simplex[m-j+1] delta_incr;
    real<lower=0, upper=1.0/dt_days> delta_max;
    real<lower=0,upper=1.0/dt_norm> gamma_max;
    real<lower=0,upper=1.0/dt_norm> rho_max; 
    real<lower=0, upper=5000> E_star; 
    real<lower=1e-10> sigma; 
    real<lower=-xi_max,upper=xi_max> xi;
    simplex[m] theta[nt_obs];
    simplex[m] w_ini;  // initial conditions
}
transformed parameters {
    real divrate;
    real<lower=0, upper=1.0/dt_days> delta[m-j+1];
    matrix<lower=0>[m,nt_obs] mod_obspos;
    real<lower=0> resp_vol_loss[nt];    // record volume loss due to respiration
    real<lower=0> growth_vol_gain[nt];  // record volume gain due to cell growth
    real<lower=0> total_vol[nt];        // record total volume
    real<lower=0> cell_count[nt];       // record relative cell count for each time step 
    {
        // helper variables
        vector[m] w_curr; 
        vector[m] w_next;
        real delta_i = 0.0;
        real gamma;
        real a;
        real rho;
        real sizelim_gamma[m];
        real x;
        int ito = 1;
      
        // populate delta using delta_incr
        delta[1] = delta_incr[1] * delta_max;
        for (i in 1:m-j){
            delta[i+1] = delta[i] + delta_incr[i+1] * delta_max;
        }
        
        // pre-compute size-limitations
        if (xi > 0){
            for (i in 1:m){ // size-class loop
                sizelim_gamma[i] = exp(xi*(v_mid[i]-v_mid[m])/v_range);
            }
        } else {
            for (i in 1:m){ // size-class loop
                sizelim_gamma[i] = exp(xi*(v_mid[i]-v_mid[1])/v_range);
            }
        }

        w_curr = w_ini;

        for (it in 1:nt){ // time-stepping loop
            // record current solution 
            // here is another place where normalization could be performed
            if (it == it_obs[ito]){
                mod_obspos[,ito] = w_curr;
                ito += 1;
                if (ito > nt_obs){
                    // just needs to be a valid index
                    // cannot use break because resp_vol_loss needs to be recorded
                    ito = 1;
                }
            }
            
            w_next = rep_vector(0.0, m);
            resp_vol_loss[it] = 0.0;
            growth_vol_gain[it] = 0.0;
            total_vol[it] = v_mid * w_curr;
            cell_count[it] = sum(w_curr);
            for (i in 1:m){ // size-class loop
                // compute delta_i
                if (i >= j){
                    delta_i = delta[i-j+1] * dt_days;
                }
                // compute gamma_i
                gamma = dt_norm * sizelim_gamma[i] * gamma_max * (1.0 - exp(-E[it]/E_star));
                // compute rho_i
                rho = dt_norm * rho_max;
                
                // fill superdiagonal (respiration)
                if (i >= j){
                    //A[i-1,i] = rho * (1.0-delta_i);
                    a = rho * (1.0-delta_i);
                    w_next[i-1] += a * w_curr[i];
                    resp_vol_loss[it] += a * w_curr[i] * v_diff[i-1];
                } else if (i > 1){
                    //A[i-1,i] = rho;
                    a = rho;
                    w_next[i-1] += a * w_curr[i];
                    resp_vol_loss[it] += a * w_curr[i] * v_diff[i-1];
                }
                // fill subdiagonal (growth)
                if (i == 1){
                    //A[i+1,i] = gamma;
                    a = gamma;
                    w_next[i+1] += a * w_curr[i];
                    growth_vol_gain[it] += a * w_curr[i] * v_diff[i];
                } else if (i < j){
                    //A[i+1,i] = gamma * (1.0-rho);
                    a = gamma * (1.0-rho);
                    w_next[i+1] += a * w_curr[i];
                    growth_vol_gain[it] += a * w_curr[i] * v_diff[i];
                } else if (i < m){
                    //A[i+1,i] = gamma * (1.0-delta_i) * (1.0-rho);
                    a = gamma * (1.0-delta_i) * (1.0-rho);
                    w_next[i+1] += a * w_curr[i];
                    growth_vol_gain[it] += a * w_curr[i] * v_diff[i];
                }
                // fill (j-1)th superdiagonal (division)
                if (i >= j){
                    //A[i+1-j,i] = 2.0*delta_i;
                    a = 2.0*delta_i;
                    w_next[i+1-j] += a * w_curr[i];
                }
                // fill diagonal (stasis)
                if (i == 1){
                    //A[i,i] = (1.0-gamma);
                    a = (1.0-gamma);
                    w_next[i] += a * w_curr[i];
                } else if (i < j){
                    //A[i,i] = (1.0-gamma) * (1.0-rho);
                    a = (1.0-gamma) * (1.0-rho);
                    w_next[i] += a * w_curr[i];
                } else if (i == m){
                    //A[i,i] = (1.0-delta_i) * (1.0-rho);
                    a = (1.0-delta_i) * (1.0-rho);
                    w_next[i] += a * w_curr[i];
                } else {
                    //A[i,i] = (1.0-gamma) * (1.0-delta_i) * (1.0-rho);
                    a = (1.0-gamma) * (1.0-delta_i) * (1.0-rho);
                    w_next[i] += a * w_curr[i];
                }
            }
            /*
            // use this check to test model consistency when division is set to "a = delta_i;" (no doubling)
            if (fabs(sum(w_next)-1.0)>1e-12){
                reject("it = ",it,", sum(w_next) = ", sum(w_next), ", rho =", rho)
            }
            */
            // do not normalize population here
            w_curr = w_next;
        }
        divrate = log(sum(w_curr))*60*24/(nt*dt); // daily division rate
    }
}
model {
    vector[m] alpha;
    
    // priors
    
    gamma_max ~ normal(10.0, 10.0) T[0,1.0/dt_norm];
    rho_max ~ normal(3.0, 10.0) T[0, 1.0/dt_norm];
    E_star ~ normal(1000.0,1000.0) T[0,];
    sigma ~ lognormal(1000.0, 1000.0) T[1,];
    xi ~ normal(0.0, 0.1);

    // fitting observations
    if (prior_only == 0){
        for (it in 1:nt_obs){
            if(i_test[it] == 0){
                alpha = mod_obspos[:,it]/sum(mod_obspos[:,it]) * sigma + 1;
                theta[it] ~ dirichlet(alpha);
                obs_count[:,it] ~ multinomial(theta[it]);
            }
        }
    }
}
