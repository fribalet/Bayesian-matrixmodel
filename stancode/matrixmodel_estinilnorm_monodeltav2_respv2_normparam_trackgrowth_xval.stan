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
    // for cross-validation
    int<lower=0, upper=1> i_test[nt_obs];
}
transformed data {
    int j;
    real v_max;
    real<lower=0> delta_v;
    real<lower=0> dt_days;  // dt in units of days
    real<lower=0> dt_norm;  // dt in units of days and doublings
    real<lower=0> v[m+1];   // vector of (minimum) sizes for each size class
    int<lower=0> t[nt];     // vector of times in minutes since start 
    int<lower=1, upper=nt> it_obs[nt_obs]; // the time index of each observation
    int n_test = sum(i_test);

    j = 1 + delta_v_inv; 
    delta_v = 1.0/delta_v_inv;
    dt_days = dt/1440.0;
    dt_norm = dt/(1440.0 * (2^delta_v - 1.0));
    for (i in 1:m+1){
        v[i] = v_min*2^((i-1)*delta_v);
    }
    v_max = v[m];
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
    simplex[m-j+1] delta_max_frac; 
    real<lower=0> delta_max_m;

    real<lower=0> gamma_max;
    real<lower=0> respiration; 
    real<lower=0, upper=5000> E_star; 
    real w_ini_mu;
    real<lower=0> w_ini_sigma;
    real<lower=1e-10> sigma; 
}
transformed parameters {
    real divrate;
    real delta_max[m-j+1]; 
    matrix[m,nt_obs] mod_obspos;
    vector<lower=0>[m] w_ini;  // initial conditions 
    {
        // helper variables
        vector[m] w_curr; 
        vector[m] w_next;
        real delta_i = 0.0;
        real gamma;
        real a;
        real rho;
        real x;
        int ito = 1;
      
        // populate delta_max using delta_max_incr
        delta_max[1] = delta_max_frac[1] * delta_max_m;
        for (i in 1:m-j){
            delta_max[i+1] = delta_max[i] + delta_max_frac[i] * delta_max_m;
        }
        
        for (i in 1:m){
            x = 0.5*(v[i]+v[i+1]) - v[1];
            w_ini[i] = exp(-(log(x)-w_ini_mu)^2/(2*w_ini_sigma^2 * pi()))/(x*w_ini_sigma*sqrt(2*pi()));
        }
        w_ini = w_ini/sum(w_ini);

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
            
            // compute gamma and rho
            gamma = gamma_max * dt_norm * (1.0 - exp(-E[it]/E_star)) - respiration * dt_norm;
            if (gamma > 0){
                rho = 0.0;
            } else {
                rho = -gamma;
                gamma = 0.0;
            }

            w_next = rep_vector(0.0, m);
            for (i in 1:m){ // size-class loop
                // compute delta_i
                if (i >= j){
                    delta_i = delta_max[i-j+1] * dt_days;
                }
                
                // fill superdiagonal (respiration)
                if (i >= j){
                    //A[i-1,i] = rho * (1.0-delta_i);
                    a = rho * (1.0-delta_i);
                    w_next[i-1] += a * w_curr[i];
                } else if (i > 1){
                    //A[i-1,i] = rho;
                    a = rho;
                    w_next[i-1] += a * w_curr[i];
                }
                // fill subdiagonal (growth)
                if (i == 1){
                    //A[i+1,i] = gamma;
                    a = gamma;
                    w_next[i+1] += a * w_curr[i];
                } else if (i < j){
                    //A[i+1,i] = gamma * (1.0-rho);
                    a = gamma * (1.0-rho);
                    w_next[i+1] += a * w_curr[i];
                } else if (i < m){
                    //A[i+1,i] = gamma * (1.0-delta_i) * (1.0-rho);
                    a = gamma * (1.0-delta_i) * (1.0-rho);
                    w_next[i+1] += a * w_curr[i];
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
    real diff;
    real popsum;
    
    // priors
    
    // delta_max for (last) size class m 
    delta_max_m ~ normal(3.0,3.0);
    delta_max_frac ~ dirichlet(rep_vector(10.0,m-j+1));

    gamma_max ~ uniform(0.0,1440.0/dt);
    //respiration ~ uniform(0.0,10.0);
    respiration ~ uniform(0.0,10.0);
    E_star ~ normal(1000.0,1000.0);
    sigma ~ exponential(1000.0);

    w_ini_mu ~ normal(-3.0, 1.0);
    w_ini_sigma ~ uniform(0.03, 3.0);

    // fitting observations
    for (it in 1:nt_obs){
        if (i_test[it] == 0){
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
}
generated quantities{
    real log_like_test = 0;
    {
        real diff;
        real popsum;

        for(it in 1:nt_obs){
            if(i_test[it] == 1){
                diff = 0.0;
                popsum = sum(mod_obspos[,it]);
                for(iv in 1:m){
                    diff += fabs(mod_obspos[iv,it]/popsum - obs[iv,it]);
                }
                diff = diff/sigma;
                log_like_test += normal_lpdf(diff | 0.0, 1.0);
            }
        }
        log_like_test = log_like_test/n_test;
    }
}
