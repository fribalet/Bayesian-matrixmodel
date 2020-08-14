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
    // initial conditions
    vector[m] w_ini;
}
transformed data {
    int j;
    real v_max;
    real<lower=0> delta_v;
    real<lower=0> dt_days;  // dt in units of days
    real<lower=0> dt_norm;  // dt in units of days and doublings
    real<lower=0> v[m+1];   // vector of (minimum) sizes for each size class
    row_vector[m] v_mid;    // vector of sizes for each size class
    real<lower=0> v_diff[m-1];// vector of size-differences for first m-1 size classes
    int<lower=0> t[nt];     // vector of times in minutes since start 
    int<lower=1, upper=nt> it_obs[nt_obs]; // the time index of each observation
    int n_test = sum(i_test);
    int mitosis_maxlen_min = 60*5; // maximum length of mitosis in units of minutes
    int nt_mitosis = mitosis_maxlen_min/dt;

    j = 1 + delta_v_inv; 
    delta_v = 1.0/delta_v_inv;
    dt_days = dt/1440.0;
    dt_norm = dt/(1440.0 * (2^delta_v - 1.0));
    for (i in 1:m+1){
        v[i] = v_min*2^((i-1)*delta_v);
    }
    v_max = v[m];
    for (i in 1:m){
        v_mid[i] = 0.5*(v[i]+v[i+1]);
    }
    for (i in 1:m-1){
        // difference between the centers for each class
        v_diff[i] = 0.5*(v[i+2]-v[i]); 
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
    real<lower=0> delta_mu; 
    real<lower=0> delta_sigma; 
    real<lower=0,upper=1.0/dt_days> delta[m-j+1]; 
    real<lower=0,upper=1.0/dt_norm> gamma_max;
    real<lower=0,upper=1.0/dt_norm> rho_max; 
    real<lower=0, upper=5000> E_star; 
    real<lower=dt,upper=dt*nt_mitosis> zeta_offset;
    real<lower=1e-10> sigma; 
    simplex[m] theta[nt_obs];
}
transformed parameters {
    real divrate;
    matrix<lower=0>[m,nt_obs] mod_obspos;
    real<lower=0> resp_vol_loss[nt];    // record volume loss due to respiration
    real<lower=0> growth_vol_gain[nt];  // record volume gain due to cell growth
    real<lower=0> total_vol[nt];        // record total volume
    real<lower=0> cell_count[nt];       // record relative cell count for each time step 
    {
        // helper variables
        matrix[m,nt_mitosis] w_curr; 
        matrix[m,nt_mitosis] w_next;
        real zeta = 0.0;
        int imit_next;
        real delta_i = 0.0;
        real gamma;
        real a;
        real rho;
        real x;
        int ito = 1;
      
        // for now, assuming they all start not in mitosis
        w_curr = rep_matrix(0.0, m, nt_mitosis);
        for (i in 1:m){
            w_curr[i,1] = w_ini[i];
        }

        for (it in 1:nt){ // time-stepping loop
            // record current solution 
            // here is another place where normalization could be performed
            if (it == it_obs[ito]){
                for (i in 1:m){
                    mod_obspos[i,ito] = 0.0;
                    for (imit in 1:nt_mitosis){
                        mod_obspos[i,ito] += w_curr[i,imit];
                    }
                }
                ito += 1;
                if (ito > nt_obs){
                    // just needs to be a valid index
                    // cannot use break because resp_vol_loss needs to be recorded
                    ito = 1;
                }
            }
            
            // compute gamma and rho
            gamma = gamma_max * dt_norm * (1.0 - exp(-E[it]/E_star)) - rho_max * dt_norm;
            if (gamma > 0){
                rho = 0.0;
            } else {
                rho = -gamma;
                gamma = 0.0;
            }

            resp_vol_loss[it] = 0.0;
            growth_vol_gain[it] = 0.0;
            total_vol[it] = 0.0;
            for (imit in 1:nt_mitosis){
                total_vol[it] += v_mid * w_curr[:,imit];
            }
            cell_count[it] = sum(w_curr);
            w_next = rep_matrix(0.0, m, nt_mitosis);
            
            //
            // non-mitotic cells
            //
            
            for (i in 1:m){ // size-class loop
                // compute delta_i
                if (i >= j){
                    delta_i = delta[i-j+1] * dt_days;
                }
                
                // respiration
                if (i >= j){
                    //A[i-1,i] = rho * (1.0-delta_i);
                    a = rho * (1.0-delta_i);
                    w_next[i-1,1] += a * w_curr[i,1];
                    resp_vol_loss[it] += a * w_curr[i,1] * v_diff[i-1];
                } else if (i > 1){
                    //A[i-1,i] = rho;
                    a = rho;
                    w_next[i-1,1] += a * w_curr[i,1];
                    resp_vol_loss[it] += a * w_curr[i,1] * v_diff[i-1];
                }
                // growth
                if (i == 1){
                    //A[i+1,i] = gamma;
                    a = gamma;
                    w_next[i+1,1] += a * w_curr[i,1];
                    growth_vol_gain[it] += a * w_curr[i,1] * v_diff[i];
                } else if (i < j){
                    //A[i+1,i] = gamma * (1.0-rho);
                    a = gamma * (1.0-rho);
                    w_next[i+1,1] += a * w_curr[i,1];
                    growth_vol_gain[it] += a * w_curr[i,1] * v_diff[i];
                } else if (i < m){
                    //A[i+1,i] = gamma * (1.0-delta_i) * (1.0-rho);
                    a = gamma * (1.0-delta_i) * (1.0-rho);
                    w_next[i+1,1] += a * w_curr[i,1];
                    growth_vol_gain[it] += a * w_curr[i,1] * v_diff[i];
                }
                // enter mitosis
                if (i >= j){
                    a = delta_i;
                    w_next[i,2] += a * w_curr[i,1];
                }
                // stasis
                if (i == 1){
                    //A[i,i] = (1.0-gamma);
                    a = (1.0-gamma);
                    w_next[i,1] += a * w_curr[i,1];
                } else if (i < j){
                    //A[i,i] = (1.0-gamma) * (1.0-rho);
                    a = (1.0-gamma) * (1.0-rho);
                    w_next[i,1] += a * w_curr[i,1];
                } else if (i == m){
                    //A[i,i] = (1.0-delta_i) * (1.0-rho);
                    a = (1.0-delta_i) * (1.0-rho);
                    w_next[i,1] += a * w_curr[i,1];
                } else {
                    //A[i,i] = (1.0-gamma) * (1.0-delta_i) * (1.0-rho);
                    a = (1.0-gamma) * (1.0-delta_i) * (1.0-rho);
                    w_next[i,1] += a * w_curr[i,1];
                }
            }

            //
            // mitotic cells
            //

            for (imit in 2:nt_mitosis-1){ // mitosis stage loop
                imit_next = imit + 1;
                // compute zeta
                if ((imit-1)*dt >= zeta_offset){
                    zeta = 1.0;
                } else {
                    zeta = 0.0;
                }
                for (i in j:m){ // size-class loop (note the starting index j)
                    // divide and leave mitosis
                    a = 2.0 * zeta;
                    //a = zeta; // for consistency check
                    w_next[i+1-j,1] += a * w_curr[i,imit];
                    
                    // respiration (not permitted to shrink below size class j)
                    if (i > j){
                        //A[i-1,i] = rho * (1.0-zeta);
                        a = rho * (1.0-zeta);
                        w_next[i-1,imit_next] += a * w_curr[i,imit];
                        resp_vol_loss[it] += a * w_curr[i,imit] * v_diff[i-1];
                    } 
                    // growth
                    if (i == j){
                        //A[i+1,i] = gamma * (1.0-zeta);
                        a = gamma * (1.0-zeta);
                        w_next[i+1,imit_next] += a * w_curr[i,imit];
                        growth_vol_gain[it] += a * w_curr[i,imit] * v_diff[i];
                    } else if (i < m){
                        //A[i+1,i] = gamma * (1.0-zeta) * (1.0-rho);
                        a = gamma * (1.0-zeta) * (1.0-rho);
                        w_next[i+1,imit_next] += a * w_curr[i,imit];
                        growth_vol_gain[it] += a * w_curr[i,imit] * v_diff[i];
                    }
                    // stasis
                    if (i == j){
                        //A[i,i] = (1.0-gamma) * (1.0-zeta);
                        a = (1.0-gamma) * (1.0-zeta);
                        w_next[i,imit_next] += a * w_curr[i,imit];
                    } else if (i == m){
                        //A[i,i] = (1.0-zeta) * (1.0-rho);
                        a = (1.0-zeta) * (1.0-rho);
                        w_next[i,imit_next] += a * w_curr[i,imit];
                    } else {
                        //A[i,i] = (1.0-gamma) * (1.0-zeta) * (1.0-rho);
                        a = (1.0-gamma) * (1.0-zeta) * (1.0-rho);
                        w_next[i,imit_next] += a * w_curr[i,imit];
                    }
                }
            }
            // forced division at maximum length of mitosis
            for (i in j:m){ // size-class loop (note the starting index j)
                w_next[i+1-j,1] += 2.0 * w_curr[i,nt_mitosis];
                //w_next[i+1-j,1] += w_curr[i,nt_mitosis]; // for consistency check
            }
            // use this check to test model consistency when division is set to "a = delta_i;" (no doubling)
            /*
            if (fabs(sum(w_next)-1.0)>1e-12){
                for (imit in 1:nt_mitosis){
                    print("in stage ",imit,": ",sum(w_next[:,imit])," before: ",sum(w_curr[:,imit]))
                }
                reject("it = ",it,", sum(w_next) = ", sum(w_next), ", zeta_offset =", zeta_offset)
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
    
    delta_mu ~ normal(3.0, 1.0);
    delta_sigma ~ exponential(1.0);
    delta ~ normal(delta_mu, delta_sigma);
    gamma_max ~ normal(10.0, 10.0) T[0,1.0/dt_norm];
    rho_max ~ normal(1.5, 0.06) T[0, 1.0/dt_norm]; // copied from exp_zs_20200701_g2_ext results
    E_star ~ normal(86, 22) T[0,]; // copied from exp_zs_20200701_g2_ext results
    sigma ~ lognormal(1000.0, 1000.0) T[1,];
    zeta_offset ~ normal(4*60,20);

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
