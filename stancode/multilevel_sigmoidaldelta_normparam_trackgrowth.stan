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
    vector<lower=0>[m] w_ini;  // initial conditions 
    // observations
    real<lower=0,upper=nt*dt>  t_obs[nt_obs]; // the time of each observation
    real<lower=0> obs[m,nt_obs]; // observations
}
transformed data {
    int j;
    real<lower=0> delta_v;
    real<lower=0> dt_days;  // dt in units of days
    real<lower=0> dt_norm;  // dt in units of days and doublings
    real<lower=0> v[m];     // vector of (minimum) sizes for each size class
    int<lower=0> t[nt];     // vector of times in minutes since start 
    int<lower=1, upper=nt> it_obs[nt_obs]; // the time index of each observation
    int ndays = 0;

    j = 1 + delta_v_inv; 
    delta_v = 1.0/delta_v_inv;
    dt_days = dt/1440.0;
    dt_norm = dt/(1440.0 * (2^delta_v - 1.0));
    for (i in 1:m){
        v[i] = v_min*2^((i-1)*delta_v);
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
    while (ndays*1440 < t[nt]){
        ndays += 1;
    }
}
parameters {
    real<lower=0> delta_max_mu;
    real<lower=0> delta_max_sigma;
    vector[ndays] delta_max_raw; 
    real<lower=0> gamma_max_mu;
    real<lower=0> gamma_max_sigma;
    vector[ndays] gamma_max_raw;
    real<lower=v[1],upper=v[m]> sig_offset_mu;
    real<lower=0> sig_offset_sigma;
    vector[ndays] sig_offset_raw;
    real<lower=0> sig_steepness_mu;
    real<lower=0> sig_steepness_sigma;
    vector[ndays] sig_steepness_raw;
    real<lower=5> E_star_mu; // introducing a lower bound of 5 here
    real<lower=0> E_star_sigma; 
    vector[ndays] E_star_raw; 
    real<lower=1e-10> sigma; 
}
transformed parameters {
    matrix[m,nt_obs] mod_obspos;
   
    // it appears like bounds are unnecessary here due to the truncation imposed in the model block
    // in fact, imposing bounds here will lead to error messages
    vector[ndays] delta_max = delta_max_mu + delta_max_sigma * delta_max_raw;
    vector[ndays] gamma_max = gamma_max_mu + gamma_max_sigma * gamma_max_raw;
    vector[ndays] sig_offset = sig_offset_mu + sig_offset_sigma * sig_offset_raw;
    vector[ndays] sig_steepness = sig_steepness_mu + sig_steepness_sigma * sig_steepness_raw;
    vector[ndays] E_star = E_star_mu + E_star_sigma * E_star_raw; 
    {
        // helper variables
        vector[m] w_curr; 
        vector[m] w_next;
        real delta_i = 0.0;
        real gamma;
        real a;
        real tmp;
        int ito = 1;
        int iday = 1;
        
        w_curr = w_ini;

        for (it in 1:nt){ // time-stepping loop
            // record current solution 
            // here is another place where normalization could be performed
            if (it == it_obs[ito]){
                mod_obspos[,ito] = w_curr;
                ito += 1;
                if (ito > nt_obs){
                    break;
                }
            }
            // increment iday
            if (t[it] > iday*1440){
                iday += 1;
            }
            // compute gamma
            gamma = gamma_max[iday] * dt_norm * (1.0 - exp(-E[it]/E_star[iday]));

            w_next = rep_vector(0.0, m);
            for (i in 1:m){ // size-class loop
                // compute delta_i
                if (i >= j){
                    delta_i = delta_max[iday]/(1.0+exp(-sig_steepness[iday]*(v[i]-sig_offset[iday]))) * dt_days;
                }
                
                // fill subdiagonal (growth)
                if (i < j){
                    //A[i+1,i] = gamma;
                    a = gamma;
                    w_next[i+1] += a * w_curr[i];
                } else if (i < m){
                    //A[i+1,i] = gamma * (1.0-delta_i);
                    a = gamma * (1.0-delta_i);
                    w_next[i+1] += a * w_curr[i];
                }
                // fill (j-1)th superdiagonal (division)
                if (i >= j){
                    //A[i+1-j,i] = 2.0*delta_i;
                    a = 2.0*delta_i;
                    w_next[i+1-j] += a * w_curr[i];
                }
                // fill diagonal (stasis)
                if (i < j){
                    //A[i,i] = (1.0-gamma);
                    a = (1.0-gamma);
                    w_next[i] += a * w_curr[i];
                } else if (i == m){
                    //A[i,i] = (1.0-delta_i);
                    a = (1.0-delta_i);
                    w_next[i] += a * w_curr[i];
                } else {
                    //A[i,i] = (1.0-gamma) * (1.0-delta_i);
                    a = (1.0-gamma) * (1.0-delta_i);
                    w_next[i] += a * w_curr[i];
                }
            }
            // do not normalize population here
            w_curr = w_next;
        }
    }
}
model {
    real diff;
    real popsum;
    int iday;
    
    // priors
    delta_max_mu ~ normal(3.5, 2.0);
    delta_max_sigma ~ exponential(1.0);

    sig_offset_mu ~ uniform(v[1],v[m]);
    sig_offset_sigma ~ exponential(1.0);
    
    sig_steepness_mu ~ lognormal(5.0,1.0);
    sig_steepness_sigma ~ exponential(1.0);

    gamma_max_mu ~ uniform(0.0,1440.0/dt);
    gamma_max_sigma ~ exponential(1.0);

    E_star_mu ~ normal(700.0,300.0);
    E_star_sigma ~ exponential(0.1);

    sigma ~ exponential(1000.0);
    
    // Outcomes in truncated distributions must be univariate, hence the loop.
    for (i in 1:ndays){
        delta_max_raw[i] ~ normal(0.0,1.0) T[-delta_max_mu/delta_max_sigma,];
        gamma_max_raw[i] ~ normal(0.0,1.0) T[-gamma_max_mu/gamma_max_sigma,];
        //sig_offset_raw[i] ~ normal(0.0,1.0) T[(v[1]-sig_offset_mu)/sig_offset_sigma,(v[m]-sig_offset_mu)/sig_offset_sigma];
        sig_offset_raw[i] ~ normal(0.0,1.0);
        sig_steepness_raw[i] ~ normal(0.0,1.0) T[-sig_steepness_mu/sig_steepness_sigma,];
        E_star_raw[i] ~ normal(0.0,1.0) T[(5.0-E_star_mu)/E_star_sigma,];
    }

    // fitting observations
    for (it in 1:nt_obs){
        diff = 0.0;
        // normalization is now happening here but could also be moved to transformed parameters block
        popsum = sum(mod_obspos[,it]);
        for (iv in 1:m){
            diff += fabs((mod_obspos[iv,it]/popsum) - obs[iv,it]);
        }
        diff = diff/sigma;
        diff ~ normal(0.0, 1.0) T[0,];
    }
}
