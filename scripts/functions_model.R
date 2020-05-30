#' Delta function. Fraction of cells that divide between t and t + dt
#' 
#' @param volbins volbins
#' @param dmax dmax
#' @param b b
#' @return delta
#' @export 
delta <- function(volbins, dmax, b){
	
	# NOTE: most values of volbins need to be < 1 for compatibility issue (HYNES et al. 2015)
	volbins <- as.numeric(volbins)
	dmax <- as.numeric(dmax)
	b <- as.numeric(b)

	#find the size that represent twice the size of the smallest size class
	j <- findInterval(2 * volbins[1], volbins)
	
	# scale volbins to allow smooth transition over volbins
	v <- volbins - volbins[j-1] 
	
	# rescale volbins so values < 1, for compatibility with delta function (Hynes et al. 2015)
	v.norm <- v/max(v)
	
	# calculate probability of division
    del <-   dmax * v.norm^(b*10) / (1 + v.norm^(b*10))
	d <- matrix(del, 1, length(volbins))
	
	# Small phytoplankton (i=1,..., j-1) are less than twice as big as the smallest size class, and so are prohibited to divide
	d[1:(j-1)] <- 0
	
	return(d)
}


#' Gamma function. Fraction of cells that grow into next size class between t and t + dt
#' 
#' @param Edata Edata
#' @param gmax gmax
#' @param E_star E_star
#' @return gamma
#' @export 
gamma_t <- function(Edata, gmax, E_star){
	
	Edata <- as.numeric(Edata)
	gmax <- as.numeric(gmax)
	E_star <- as.numeric(E_star)

	# y <- (1-exp(-Edata/E_star)) * gmax # original 2003 model
	y <- (gmax/(1000*E_star)) * Edata 
	y[which(Edata >= 1000*E_star)] <- gmax
	return(y)
}


#' Rho function. Fraction of cells that shrink between t and t + dt
#' Assumptions:
#' 1) rate of respiration is a function of photosynthetic rate (gamma function)
#' 2) rate of respiration constant (h-1) constant over the nighttime period 
#'	(polysaccharide is drawn down linearly over the nighttime period)
#' @param y Gamma values
#' @return rho
#' @export 
rho_t <- function(y){

	y <- as.numeric(y)
	
	# 30% total fixed carbon is respired over a 24h-period
	r <-  0 * 0.3 * mean(y[which(y > 0)]) - y 

	# probablity to decrease size class = 0 when growth > respiration
	r[which(r < 0)] <- 0 

	return(r)
}


#' Construct matrix A(t) for each time step within an hour based on delta and gamma at each 10 minute time #' intervals then construct B(t) which is A(t)'s multiplied for the given hour
#' multiply B(t)*w(t) to get the projection to next hour dist if desired
#'
#' @param hr blabla
#' @param Edata blabla
#' @param volbins blabla
#' @param gmax blabla
#' @param dmax blabla
#' @param b blabla
#' @param E_star blabla
#' @param resol blabla
#' @param breaks blabla
#' @return Matrix B(t)
#' @examples
#' \dontrun{
#' }
#' @export
matrix_conct_fast <- function(hr, Edata, volbins, gmax, dmax, b, E_star, resol){

		# time interval 
		dt <- resol/60

		# no division during the first 't.nodiv' hours after dawn (only for Synechococcus model)
		t.nodiv <- 0 
		
		#find the sie that represent twice the size of the smallest size class
		j <- findInterval(2 * volbins[1], volbins)

 		# dimensions of the squared matrix
		m <- length(volbins)

		## growth
		y <- gamma_t(Edata=Edata, gmax=gmax, E_star=E_star)

		## respiration
		r <- rho_t(y=y)

		## division
		del <- delta(volbins, dmax=dmax, b=b)
			
			if(hr < t.nodiv){div  <- matrix(data=0, 1, m)
			}else{div <- matrix(del, 1, m)}
    

  
		## CONSTRUCTION SPARSE MATRIX 
    	# Diagonal (0) stasis 
		stasis_ind <- seq(1,m^2,by=m+1) 
		# Subdiagonal (-1) growth 
		growth_ind <- seq(2,m^2,by=m+1) 
		# Superdiagonal (+1) respiration 
		resp_ind <- seq(m+1, m^2, by=m+1) 
		# Superdiagonal (j-1) division
		div_ind <- seq((((j-1)*m)+1), m^2, by=m+1) 

		for(t in 1:(1/dt)){

			A <- matrix(data=0,nrow=m, ncol=m)
			growth <- y[t+hr/dt]
			resp <- r[t+hr/dt]


			# Stasis (main diagonal)
			A[stasis_ind] <- (1-div)*(1-growth)*(1-resp) # the hr/dt part in the indexing is because each hour is broken up into dt segments for the irradiance spline
			A[m,m] <- (1-div[m])*(1-resp)
			A[1,1] <- (1-growth)

			# Cell growth (subdiagonal)
			A[growth_ind] <- growth*(1-div[1:(m-1)])

			# Division (first row and superdiagonal j-1)
			A[div_ind] <- 2 * div[j:m] # The cell division terms for large (i > = j) phytoplankton

			# Respiration (superdiagonal)
			A[resp_ind] <- resp*(1-div[2:m])*(1-growth)

					if(t == 1){B <- A}else{B <- A %*% B}
			}

		return(B)
}


#' This function calculates the sum of squares of the of the differences between the hourly observations and 
#' the model given the specified parameters
#' This function returns a column vector - called by "determine.opt.para" for the optimization.
#'
#' @param params blabla
#' @param Edata blabla
#' @param distribution blabla
#' @param resol blabla
#' @return goodness of fit
#' @examples
#' \dontrun{
#' }
#' @export
sigma_lsq <- function(params=params, Edata=Edata, distribution=distribution, resol=resol){
               	# gmax <- 0.2 
                # dmax <- 0.2
                # b <- 0.4
                # E_star <- 0.4
				time.interval <- median(diff(distribution$time))
				res <- which(diff(distribution$time) == time.interval)# select time that have at least 2 consecutive time points, required for comparing the projection to the next time point
                PSD <- t(distribution[,-c(1)])
				dim <- dim(PSD)
				sigma <- matrix(NA, dim[1], dim[2]-1) # preallocate sigma
				TotN <- as.matrix(colSums(PSD))
                volbins <- as.numeric(rownames(PSD))
        
                gmax <- as.numeric(params[1]) 
                dmax <- as.numeric(params[2])
                b <- as.numeric(params[3]) 
                E_star <- as.numeric(params[4])

			for(hr in res){
					B <- matrix_conct_fast(hr=hr-1, Edata=Edata, volbins=volbins, gmax=gmax, dmax=dmax, b=b, E_star=E_star, resol=resol)
					wt <- B %*% PSD[,hr]/sum(PSD[,hr], na.rm=T) # calculate the projected size-frequency distribution
					wt.norm <- wt/sum(wt, na.rm=T) # normalize distribution
					sigma[,hr] <-  (PSD[,hr+1] - round(TotN[hr+1]*wt.norm))^2 # ABSOLUTE observed value - fitted value
					}
			sigma <- colSums(sigma)/colSums(PSD[,-1])
			sigma <- mean(sigma, na.rm=T)
			return(sigma)

}



#' This function calculates the sum of squares of the of the differences between the hourly observations 
#' and the model given the specified parameters using Hubert Loss (1954) approach.
#'
#' @param params blabla
#' @param Edata blabla
#' @param distribution blabla
#' @param resol blabla
#' @return goodness of fit
#' @examples
#' \dontrun{
#' }
#' @export
sigma_hl <- function(params=params, Edata=Edata, distribution=distribution, resol=resol){

			    time.interval <- median(diff(distribution$time))
				res <- which(diff(distribution$time) == time.interval)# select time that have at least 2 consecutive time points, required for comparing the projection to the next time point
                PSD <- t(distribution[,-c(1)])
				dim <- dim(PSD)
				sigma <- matrix(NA, 1, dim[2]-1) # preallocate sigma
				TotN <- as.matrix(colSums(PSD))
                volbins <- as.numeric(rownames(PSD))
        
                gmax <- as.numeric(params[1]) 
                dmax <- as.numeric(params[2]) 
                b <- as.numeric(params[3])
                E_star <- as.numeric(params[4])

        d <- 1.345
			for(hr in res){
					B <- matrix_conct_fast(hr=hr-1, Edata=Edata, volbins=volbins, gmax=gmax, dmax=dmax, b=b, E_star=E_star, resol=resol)
					wt <- B %*% PSD[,hr]/sum(PSD[,hr], na.rm=T) # calculate the projected size-frequency distribution
					wt.norm <- wt/sum(wt, na.rm=T) # normalize distribution
                  # Huber loss calculation
                    a <- PSD[,hr+1] - round(TotN[hr+1]*wt.norm)
                    loss <- ifelse(abs(a) <= d,
                                       0.5 * a^2,
                                       d * (abs(a) - 0.5 * d))
                    sigma[,hr] <- mean(loss, na.rm=T)
                     }
            sigma <- sum(sigma, na.rm=T)/100 ## HUBER loss

			return(sigma)

}


#' model optimization usin Diff Evolution (require DEoptim package)
#' 
#' @param distribution blabla
#' @param Edata blabla
#' @param resol blabla
#' @return optimization
#' @examples
#' \dontrun{
#' }
#' @export
determine_opt_para <- function(distribution=distribution,Edata=Edata,resol=resol){

		require(DEoptim)
		#require(cmaes)

		distributuion <- as.matrix(distribution)
		resol <- as.numeric(resol)
		Edata <- as.vector(Edata)

		##################
		## OPTIMIZATION ##
		##################
		print("Optimizing gmax, dmax, b and E_star")
		start <- Sys.time()

		f <- function(params) sigma_lsq(params, Edata, distribution, resol)
		#f <- function(params) sigma_hl(params, Edata, distribution, resol)

		opt <- DEoptim::DEoptim(f, lower=c(1e-6,1e-6,1e-6,1e-6), upper=c(1,1,1,1), 
								control=DEoptim.control(itermax=1000, 
														reltol=1e-6, 
														trace=10, 
														steptol=100, 
														strategy=2, 
														parallelType=0))

		print(Sys.time()-start)

		params <- opt$optim$bestmem
        gmax <- as.numeric(params[1]) 
        dmax <- as.numeric(params[2])
        b <- as.numeric(params[3]) 
        E_star <- as.numeric(params[4])
		resnorm <- opt$optim$bestval
		
		# opt <- cma_es(par=c(0.5,0.5,0.5,0.5),f, lower=c(0,0,0,0), upper=c(1,1,1,1))
		# params <- opt$par
        # gmax <- as.numeric(params[1]) 
        # dmax <- as.numeric(params[2])
        # b <- as.numeric(params[3]) 
        # E_star <- as.numeric(params[4])
		# resnorm <- opt$value

		####################################################
		## Calculate projections from best fit parameters ##
		####################################################
		print(params)
		time.interval <- median(diff(distribution$time))
		res <- which(diff(distribution$time) == time.interval)# select time that have at least 2 consecutive time points, required for comparing the projection to the next time point
        Nproj <- t(distribution[,-c(1)])
        volbins <- as.numeric(rownames(Nproj))
        
		for(hr in res){
					B <- matrix_conct_fast(hr=hr-1, Edata, volbins, gmax, dmax, b, E_star, resol)
					Nproj[,hr+1] <- round(B %*% Nproj[,hr]) # calculate numbers of individuals
					}
		
		PSD <- data.frame(t(Nproj))
		PSD <- add_column(PSD, time=distribution$time, .before=T)
		colnames(PSD) <- colnames(distribution)

		#############################
		## Growth rate calculation ##
		#############################
		mu_N <- diff(log(rowSums(PSD[,-c(1)], na.rm=T))) / as.numeric(diff(PSD$time))
		print("hourly growth rates:")
		print(round(mu_N[-c(25)],3)) # last two values are NA , ?? due to NA from Edata

		parameters <- data.frame(cbind(gmax,dmax,b,E_star,resnorm), row.names=NULL)

		output <- list(parameters, PSD)
		names(output) <- c("parameters","PSD")
		return(output)
}
