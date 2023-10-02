#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#		Calibrate model
#		
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
source("stocking_mod.R")

###################
## Calibration ####
###################
# tuning catchability (match desired Depletion (SSB_y/SSB_0))
print(input_list$q)

# calibrate PWU functions (3 parameters)
# 1 - max PWU (1.5 for all)
# 2 - steepness (larger = steeper)
# 3 - inflection point (typically mean value)

## values for PWU curve  
Ut_list <- list()
Ut_list[[1]] <- input_list$var.cpue <-  c(1.5, 4, 1.979228)
Ut_list[[2]] <- input_list$var.harv <- c(1.5, 10, 1.025765)
Ut_list[[3]] <- input_list$var.size <- c(1.5, 0.002, 509.4874)
Ut_list[[4]] <- input_list$var.dist <- c(1.5, -0.03, 25)
Ut_list[[5]] <- input_list$var.crowd <- c(1.5, -0.0008, 1475.775)


#################
## Run model ####
#################
years <- input_list$years
nsites <- input_list$nsites
H_null <- matrix(1, nrow=years, ncol=nsites)
stock_null <- matrix(rep(0,nsites), nrow=years, ncol=nsites, byrow=T)
test <- stocking_mod(input_list,
                    sig1e = 1,
                    Ls = 100,
                    H = H_null,
                    stock = stock_null,
                    isopleth_plots = FALSE,
                    DD_flag = TRUE,
                    M = 0.1,
                    discard = 0.1,
                    recK = 15,
                    qt = 0.000823,
                    DD_sd = 0.01)
mean(test$WSB[input_list$years - input_list$years_init,])
# paste to values above (inflection point)
mean(test$cpue)
mean(test$hpue)
mean(test$avg_size)
mean(test$et)


################
## Plotting ####
################
# check shape of utility functions
# cpue and hpue - logistic increase
# average size - linear increase
# distance and crowding (i.e., effort) - linear decrease (disutilities)
sig1e <- c(0.5,1,200)
vect <- list()

Tot_init <- test$U.init
vect[[1]] <- seq(0, Ut_list[[1]][3]*2, length.out=500)
vect[[2]] <- seq(0, Ut_list[[2]][3]*2, length.out=500)
vect[[3]] <- seq(0.01, Ut_list[[3]][3]*2, length.out=500)
vect[[4]] <- seq(0.01, Ut_list[[4]][3]*2, length.out=500)
vect[[5]] <- seq(0.01, Ut_list[[5]][3]*2, length.out=500)
Ut_tot_sum <- seq(0.01, 2*Tot_init, length.out = 500)

pmax_eff <- list()
for(i in 1:length(sig1e)) pmax_eff[[i]] <- 1/(1+exp(-(Ut_tot_sum-Tot_init)/(sig1e[i]*Tot_init)))

## Plotting
ylab_list <- c("PWU Catch per unit effort", "PWU Harvest per unit effort", "PWU Average size", "PWU Distance", "PWU Crowding", "Proportion of maximum effort")
xlab_list <- c("Catch per unit effort", "Harvest per unit effort", "Average size (cm)", "Distance", "Total effort", "Total utility")
lty_type <- 2:4
col_type <- c("blue","black","darkgreen")

export <- FALSE # export plot
if(export) jpeg(filename="utility_functions.jpeg", height = 5, width = 8, units = "in", res = 500)
  par(mfrow = c(2,3), mar = c(4.3,4.3,2,2), oma =c(0.1,0.1,1,0.1))
    for(i in 1:5){
      Ut <- Ut_list[[i]][1]/(1+exp(Ut_list[[i]][2]*(Ut_list[[i]][3]-vect[[i]])))
      plot(vect[[i]], Ut, ylab = ylab_list[i], xlab = xlab_list[i], xlim = c(0, vect[[i]][length(vect[[i]])]), type = "l")
      abline(v = Ut_list[[i]][3], col = "red", lty = 3)
      if(i == 2) mtext("Part-worth utility (PWU) functions and effort allocation", side = 3, line = 1.6, cex = 0.8)
    }
  plot(Ut_tot_sum, seq(0.01,max(unlist(pmax_eff))*1.1,length.out=500), ylab = ylab_list[6], xlab = xlab_list[6], type = "n")
    for(i in 1:3){
      lines(Ut_tot_sum, pmax_eff[[i]], col = col_type[i], lty = lty_type[i])
    }
    abline(v = Tot_init, col = "red", lty = 3)
    legend("bottomright", legend = c("high", "moderate", "low"), col = col_type, lty = 2:4, bty = "n", cex = 0.8)
if(export) dev.off()
