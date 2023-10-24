#' @title stocking_mod
#'
#' @description function to run a stocking model (for largemouth bass)
#'
#' INPUT LIST PARAMETERS (set here)
#' @param years how long to run model (include years_init)
#' @param years_init burnin period
#' @param nsites number of sites/lakes
#' @param Amax maximum age
#' @param Linf_mu mean asymptotic length (across lakes)
#' @param vbk Brody growth coefficient
#' @param t0 time at length = 0
#' @param M natural mortality
#' @param lorenzc Lorenzen mortality allometric shape parameter
#' @param recK compensation ratio
#' @param TLr reference length for survival function
#' @param alfa weight length conversion a parameter
#' @param bet weight length conversion b parameter
#' @param R0_sing unfished recruitment at single lake
#' @param city_pen "penalty" for being clsoer to city
#' @param dist relative travel distances assuming single human population center near site 1
#' @param qt catchability
#' @param reg_eff approximate effort (can be informed by data)
#' @param cost cost function for gravity model ("null", "multi", "exp")
#' @param cstar default = 1; logical parameter for whether larval are assumed to be capable of locating recruitment habitat (1) or not (0)
#' @param persis default = 1; stickyness parameter for effort dynamics (seen in van Poorten and Camp 2019)
#' @param kill rate of killed (on purpose)
#' @param discard discard rate (regulatory or voluntary)
#' @param Vcat_low minimum size for vulnerability to capture
#' @param Vcat_high maximum size for vulnerability to capture
#' @param Vcat_sd standard deviation of capture vulnerability
#' @param Vret_low minimum size for retention (retention)
#' @param Vret_sd standard deviation of retention curve
#' @param So.92 back-scaled survival of fish size
#' first parm is the max pwu, second determines steepness (larger numbers = steeper), third is inflection point
#' @param var.cpue utility parameters for CPUE
#' @param var.harv utility parameters for HPUE
#' @param var.size utility parameters for average size of fish
#' @param var.dist dis-utility parameters for distance
#' @param var.crowd dis-utility parameters for crowding (function of effort)
#' @param fit_hat hypothetical fitness of hatchery fish (defined as wild born offspring of hatchery originated fish)
#' @param fit_st hypothetical fitness of stocked fish
#' @param hert_hat heritability of life history traits
#' @param hert_st heritability of life history traits
#' @param ism instaneous stocking mortality--the proportion of stocked fish that die immediately upon stocking from shock/transport stress
#' @param L0 length at beginning of DD process, assumed to be settlement, from Stunz et al. 2002
#' @param M1 natural mortality/year of 10 mm fish, roughly 15
#' @param ts time to grow from L0 to Lr
#' @param Lr Length at recruitment, ballpark guess
#'

## Largemouth bass life history parameter
input_list <- list(
  seed = 24,
  years = 150,
  years_init = 50,
  nsites = 5,
  Amax = 12,
  Linf_mu = 643,
  vbk = 0.26,
  t0 = -0.024,
  lorenzc = 1,
  TLr = 450,
  alfa = 0.00000339,
  bet = 3.24,
  R0_sing = 4e4,
  city_pen = c(0.5, 0.6, 0.7, 0.8, 0.9),
  dist = c(10, 20, 30, 40, 50),
  tot_eff = 10255 * (1 + (1 - 0.27)), # 10255 hours accounted for 27% of total effort (FWC 2021 report)
  cost_type = "null",
  cstar = 1,
  persis = 0.5,
  kill = 0.37,
  Vcat_low = 250,
  Vcat_high = 850,
  Vcat_sd = 0.1,
  Vret_low = 350,
  Vret_sd = 0.1,
  So.92 = 0.89,
  var.cpue = c(1.5, 4, 2.373142),
  var.harv = c(1.5, 10, 1.2314),
  var.size = c(1.5, 0.002, 499.1246),
  var.dist = c(1.5, -0.03, 25),
  var.crowd = c(1.5, -0.0008, 1505.835),
  grav_pow = 1,
  fit_hat = 1,
  fit_st = 0.85,
  hert_hat = 0.2,
  hert_st = 0,
  ism = 0.15,
  L0 = 25,
  M1 = 15,
  ts = 0.7,
  Lr = 128,
  trophy_age = 8
)

################
## FUNCTION ####
################
#' parameters that will change between scenarios
#' @param input_list list of life history values (not changing)
#' @param sig1e parameter that changes response of dynamic effort (low numbers = high response, high numbers = low response)
#' @param Ls length at stocking
#' @param H # habitat matrix
#' @param stock stocking scenario
#' @param isopleth_plots if TRUE, export results for isopleths (WSB, total utility, CTB, and effort summed across sites)
#' @param DD_flag density dependence flag; TRUE is on, FALSE is off
#' @param M natural mortality
#' @param discard discard mortality rate
#' @param recK compensation ratio
#' @param qt catchability
#' @param DD_sd standard deviation for Linf at various lakes (assuming normal distribution, mean = Linf mean)

stocking_mod <- function(input_list,
                         sig1e,
                         Ls,
                         H,
                         stock,
                         isopleth_plots = FALSE,
                         DD_flag = TRUE,
                         M = 0.1,
                         discard = 0.1,
                         recK = 15,
                         qt = 0.00115,
                         DD_sd = 0.01) {
  with(input_list, {
    ################
    ## Settings ####
    ################
    set.seed(seed)
    ## Life history
    Age <- 1:Amax
    Linf_mat <- matrix(NA, nrow = years, ncol = nsites) # matrix of asymptotic lengths (different across lakes and years)
    if (DD_flag) Linf_mat[1, ] <- rnorm(nsites, Linf_mu, Linf_mu * DD_sd) # if DD_flag == FALSE, this should be NAs
    TL_global <- Linf_mu * (1 - exp(-vbk * Age - t0)) # global length at age
    TL_bar <- array(NA, c(years, Amax, nsites)) # length at age for each lake and year (modeling density dependent growth)
    for (a in 1:Amax) TL_bar[1, a, ] <- Linf_mat[1, ] * (1 - exp(-vbk * Age[a])) # for just first year (estimated in time dynamics loop)
    Wt <- alfa * TL_global^bet # weight at age
    wmat <- alfa * TLr^bet # refers to maturity (used to calculate fecundity)
    Fec <- ifelse(Wt < wmat, 0, Wt - wmat)
    Mat <- Fec / max(Fec) # maturity at age (logistic function)
    R0_site <- rep(R0_sing / nsites, nsites) # R0 (unfished recruitment) across lakes
    So <- exp(-(M * TLr / TL_global))^lorenzc # survival
    So_st <- So
    So_hat <- So
    S <- c(1, cumprod(So[1:(Amax - 1)])) # survivorship (used for recruitment calculations)
    S_mat <- matrix(rep(S, each = nsites), nrow = Amax, ncol = nsites, byrow = TRUE)

    ## Selectivities
    # vulnerability to capture
    Vulcat_low <- 1 / (1 + exp(-(TL_global - Vcat_low) / (Vcat_sd * Vcat_low)))
    Vulcat_high <- 1 / (1 + exp(-(TL_global - Vcat_high) / (Vcat_sd * Vcat_high)))
    Vulcat <- Vulcat_low - Vulcat_high
    Vulcat <- Vulcat / max(Vulcat)
    Vulcat_mat <- matrix(rep(Vulcat, each = nsites), nrow = Amax, ncol = nsites, byrow = TRUE)
    # retention
    ret_sigma <- 0.1 * Vcat_low # sd from CV
    Vulret <- pnorm((TL_global - Vcat_low) / ret_sigma)
    Vulret <- Vulret / max(Vulret)
    # vulnerability to discard/release
    Vuldis <- Vulcat * (1 - Vulret)
    Vuldis <- Vuldis / max(Vuldis)
    Vuldis_mat <- matrix(rep(Vuldis, each = nsites), nrow = Amax, ncol = nsites, byrow = TRUE)

    ## Socioeconomic and fishery parameters
    # costs
    cost <- switch(cost_type,
      null = 50 + dist * 0, # low differences in costs between sites, so high movement of anglers
      mult = 50 + dist * 3,
      exp = 5 + dist^1.5 # great differences in costs between sites, so low movement of anglers
    )



    ###################
    ## Recruitment ####
    ###################
    EPR0 <- sum(Fec * S) # eggs per recruit unfished
    E0 <- (R0_sing / nsites) * EPR0 # total eggs unfished
    # h = CR/(CR+4), CR = 4h/(1-h)
    h <- recK / (recK + 4) # steepness conversion
    a <- (4 * h) / (EPR0 * (1 - h)) # Beverton-Holt a parameter
    b <- (5 * h - 1) / ((1 - h) * R0_site * EPR0) # Beverton-Holt b parameter

    ## habitat specific parameters
    ahab <- a * H^cstar
    bhab <- t(t(H^cstar) * b) / H
    R0_mat <- (ahab * EPR0 - 1) / t(t(bhab) * EPR0)

    ## unpacking
    aa <- ahab # local Beverton-Holt parameters for unpacking
    bb <- ahab / bhab # local Beverton-Holt parameters for unpacking
    Ls1 <- rep(Ls, nsites) # length at stocking across lakes
    DD <- (Ls1 - L0) / (Lr - L0) #
    v <- (Lr - L0) / ts # linear growth rate through recruitment period
    S1 <- (L0 / Ls1)^(M1 / v) # base survival for phase 1 if stage 2
    S2 <- (Ls1 / Lr)^(M1 / v) # base survival for phase 2 of stage 2
    S1S2 <- S1 * S2 # cumulative base survival for recruitment period
    SR <- t(t(aa) / S1S2) # survival rate of larvae
    a1 <- S1 # survival for phase 1 of stage 2
    b1 <- t(t(aa) * DD) / bb # DD component of survival for phase 1 of stage 2
    a2 <- S2 # survival for phase 2 of stage 2
    b2 <- (aa / bb - b1) / (t(t(SR) * a1)) ## DD component of survival for phase 2 of stage 2

    ## hatchery and stocking parameters
    a1_hat <- a1 * fit_hat
    a2_hat <- a2 * fit_hat
    a2_st <- a2 * fit_st



    ######################
    ## Initialization ####
    ######################
    et <- matrix(0, nrow = years, ncol = nsites) # effort
    hr <- matrix(0, nrow = years, ncol = nsites) # harvest rate
    FM_tot <- array(0, c(years, Amax, nsites)) # total fishing mortality (direct+discard)
    FM_direct <- array(0, c(years, Amax, nsites)) # direct fishing mortality (from fishery)
    FM_discard <- array(0, c(years, Amax, nsites)) # discard (regulatory and release) fishing mortality
    st <- matrix(0, nrow = years, ncol = nsites) # stocked fish
    larv <- matrix(NA, nrow = years, ncol = nsites) # wild larvae
    larv_hat <- matrix(NA, nrow = years, ncol = nsites) # hatchery larvae
    larv_tot <- matrix(NA, nrow = years, ncol = nsites) # total larvae (wild+hatchery)
    eggs <- matrix(NA, nrow = years, ncol = nsites) # wild eggies
    eggs_hat <- matrix(NA, nrow = years, ncol = nsites) # hatchery eggies
    R <- matrix(NA, nrow = years, ncol = nsites) # wild recruits
    R_hat <- matrix(NA, nrow = years, ncol = nsites) # hatchery recruits
    R_st <- matrix(NA, nrow = years, ncol = nsites) # stocked recruits
    eff_dens <- matrix(NA, nrow = years, ncol = nsites) # effective density (for DD growth)
    N1 <- matrix(NA, nrow = years, ncol = nsites) # stage 1 total recruits
    N1_hat <- matrix(NA, nrow = years, ncol = nsites) # stage 1 hatchery recruits
    N1_w <- matrix(NA, nrow = years, ncol = nsites) # stage 1 wild recruits
    N2_tot <- matrix(NA, nrow = years, ncol = nsites) # stage 2 total recruits
    nage <- array(0, c(years, Amax, nsites)) # numbers at age of wild
    nage_hat <- array(0, c(years, Amax, nsites)) # numbers at age of hatchery
    nage_st <- array(0, c(years, Amax, nsites)) # numbers at age of stocked
    grav_wt <- matrix(NA, nrow = years, ncol = nsites) # gravity weights (calculated as function of total utility)
    pmax_eff <- NULL # max proportion of effort allocated
    persis_pmax_eff <- NULL # stickyness of effort
    catch_num <- matrix(NA, nrow = years, ncol = nsites) # catch in numbers (wild, hatchery, stocked)
    catch_numage <- array(NA, c(years, Amax, nsites)) # catch at age in numbers
    catch_numage_hat <- array(NA, c(years, Amax, nsites)) # catch at age in numbers
    catch_numage_st <- array(NA, c(years, Amax, nsites)) # catch at age in numbers
    harvest_num <- matrix(NA, nrow = years, ncol = nsites) # harvest in numbers (includes discarded)
    hpue <- matrix(NA, nrow = years, ncol = nsites) # harvest per unit of effort
    avg_size <- matrix(NA, nrow = years, ncol = nsites) # average size
    cpue <- matrix(NA, nrow = years, ncol = nsites) # catch per unit of effort
    bt <- matrix(NA, nrow = years, ncol = nsites) # DD growth parameter
    Linft <- matrix(NA, nrow = years, ncol = nsites) # DD asymptotic length
    U.cpue <- matrix(NA, nrow = years, ncol = nsites) # utility for CPUE
    U.hpue <- matrix(NA, nrow = years, ncol = nsites) # utility for HPUE
    U.size <- matrix(NA, nrow = years, ncol = nsites) # utility for size
    U.dist <- matrix(NA, nrow = years, ncol = nsites) # utility for distance
    U.crd <- matrix(NA, nrow = years, ncol = nsites) # utility for crowding (as a function of effort)
    U.st <- matrix(0.25, nrow = years, ncol = nsites) # utility for stocking
    U.tot_add <- matrix(NA, nrow = years, ncol = nsites) # total utility (adding)
    U.tot <- matrix(NA, nrow = years, ncol = nsites) # total utility
    U.tot_sum <- NULL # total utility across sites
    val <- matrix(NA, nrow = years, ncol = nsites) # value (utility*effort)

    ## First year - initialize
    # dynamic effort
    pmax_eff[1] <- sum(et[1, ]) / tot_eff # proportion of effort allocated in year
    persis_pmax_eff[1] <- pmax_eff[1] # no stickyness first year

    # fishing mortality
    FM_direct[1, , ] <- t(t(Vulcat_mat * Vuldis) * et[1, ]) * qt
    FM_discard[1, , ] <- Vuldis_mat * discard
    FM_tot[1, , ] <- FM_direct[1, , ] + FM_discard[1, , ]

    # numbers at age, eggs, recruitment
    nage[1, , ] <- t(R0_mat[1, ] * t(S_mat)) # unfished recruitment and survivorship
    eggs[1, ] <- sum(Fec * nage[1, , ])
    # larv[1,] <- eggs[1,]
    nage_hat[1, , ] <- 0
    eggs_hat[1, ] <- 0
    larv_hat[1, ] <- 0
    R_hat[1, ] <- 0
    nage_st[1, , ] <- 0
    R_st[1, ] <- 0

    # catch and harvest
    catch_num[1, ] <- 0 # no fishing first year
    catch_numage[1, , ] <- 0

    # utility
    cpue[1, ] <- colSums(nage[1, , ] * Vulcat)
    hpue[1, ] <- colSums(nage[1, , ] * Vulcat * Vulret)
    ## avg_size - DD_flag
    avg_size[1, ] <- colSums(nage[1, , ] * Vulcat * TL_global) / colSums(nage[1, , ] * Vulcat)
    if (DD_flag) avg_size[1, ] <- colSums(nage[1, , ] * Vulcat * TL_bar[1, , ]) / colSums(nage[1, , ] * Vulcat)
    U.cpue[1, ] <- var.cpue[1] / (1 + exp(var.cpue[2] * (var.cpue[3] - cpue[1, ])))
    U.hpue[1, ] <- var.harv[1] / (1 + exp(var.harv[2] * (var.harv[3] - hpue[1, ])))
    U.size[1, ] <- var.size[1] / (1 + exp(var.size[2] * (var.size[3] - avg_size[1, ])))
    U.dist[1, ] <- var.dist[1] / (1 + exp(var.dist[2] * (var.dist[3] - dist)))
    U.crd[1, ] <- var.crowd[1] / (1 + exp(var.crowd[2] * (var.crowd[3] - et[1, ])))
    U.tot_add[1, ] <- U.cpue[1, ] + U.hpue[1, ] + U.size[1, ] + U.dist[1, ] + U.crd[1, ] + U.st[1, ]
    U.tot[1, ] <- U.tot_add[1, ]
    U.tot_sum[1] <- sum(U.tot[1, ])
    U.init <- 1 * U.tot_sum[1]

    # DD growth
    eff_dens[1, ] <- colSums(nage[1, , ] * TL_global^2)
    if (DD_flag) eff_dens[1, ] <- colSums(nage[1, , ] * TL_bar[1, , ]^2) # effective density
    gam1 <- 1.25 * Linf_mat[1, ]
    gam2 <- ((gam1 / Linf_mat[1, ]) - 1) / eff_dens[1, ]



    #####################
    ## Time Dynamics ####
    #####################
    for (y in 2:years) {
      for (k in 1:nsites) grav_wt[y, k] <- (U.tot[y - 1, k] / cost[k])^grav_pow

      # proportion of max effort that goes out in each year
      pmax_eff[y] <- 1 / (1 + exp(-(U.tot_sum[y - 1] - U.init) / (sig1e * U.init)))
      # accounts for "stickyness" in effort in each year
      persis_pmax_eff[y] <- pmax_eff[y] * (1 - persis) + pmax_eff[y - 1] * persis

      ## stocking
      st[y, ] <- 0
      # stocking occurs after burnin and 30 years, then reduced by instantaneous mortality
      if (y >= 30 + years_init) st[y, ] <- stock[y, ] * (1 - ism)

      for (k in 1:nsites) {
        if (st[y, k] > 0) U.st[y, k] <- U.st[y, k] * (0.25 * 1 + stock[y, k] / max(stock[y, ]))
        # Effort
        et[y, k] <- (tot_eff * persis_pmax_eff[y] * grav_wt[y, k]) / sum(grav_wt[y, ])
        hr[y, k] <- 1 - exp(-et[y, k] * qt)
      }

      ## fishing mortalities (direct, discard, and total)
      FM_direct[y, , ] <- qt * t(t(Vulcat_mat * Vuldis) * et[y, ])
      FM_discard[y, , ] <- Vuldis_mat * discard
      FM_tot[y, , ] <- FM_direct[y, , ] + FM_discard[y, , ]

      ## DD growth
      Linf_mat[y, ] <- gam1 / (1 + gam2 * eff_dens[y - 1, ])
      # calculate length at age per year and site
      if (DD_flag) TL_bar[y, 1, ] <- Linf_mat[y, ] * (1 - exp(-vbk))

      ## Recruitment
      # no dispersal (no movement)
      larv[y, ] <- eggs[y - 1, ]
      larv_hat[y, ] <- eggs_hat[y - 1, ]
      larv_tot[y, ] <- larv[y, ] + larv_hat[y, ]
      # first stage of DD
      N1_hat[y, ] <- (larv_hat[y, ] * (1 - hert_hat)) * SR[y, ] * a1_hat / (1 + b1[y, ] * larv_tot[y, ])
      N1_w[y, ] <- (larv[y, ] + (hert_hat * larv_hat[y, ])) * SR[y, ] * a1_hat / (1 + b1[y, ] * larv_tot[y, ])
      # second stage of DD
      N2_tot[y, ] <- N1_hat[y, ] + N1_w[y, ] + st[y, ]
      R_hat[y, ] <- N1_hat[y, ] * a2_hat / (1 + b2[y, ] * N2_tot[y, ])
      R_st[y, ] <- st[y, ] * a2_st / (1 + b2[y, ] * N2_tot[y, ])
      R[y, ] <- N1_w[y, ] * a2 / (1 + b2[y, ] * N2_tot[y, ])
      # mortality on recruits before age 1
      nage[y, 1, ] <- R[y, ] * So.92
      nage_hat[y, 1, ] <- R_hat[y, ] * So.92
      nage_st[y, 1, ] <- R_st[y, ] * So.92

      ## Continuous fishing
      for (a in 2:Amax) {
        # wild fish
        nage[y, a, ] <- nage[y - 1, a - 1, ] * exp(-(M + FM_tot[y - 1, a - 1, ]))
        # hatchery fish
        nage_hat[y, a, ] <- nage_hat[y - 1, a - 1, ] * exp(-(M + FM_tot[y - 1, a - 1, ]))
        # stocked fish
        nage_st[y, a, ] <- nage_st[y - 1, a - 1, ] * exp(-(M + FM_tot[y - 1, a - 1, ]))
        # DD growth
        if (DD_flag) TL_bar[y, a, ] <- Linf_mat[y, ] + (TL_bar[y - 1, a - 1, ] - Linf_mat[y, ]) * exp(-vbk)
      }

      ## end of season calculations
      eff_dens[y, ] <- colSums(nage[y, , ] * TL_global^2)
      if (DD_flag) eff_dens[y, ] <- colSums(nage[y, , ] * TL_bar[y, , ]^2)
      catch_numage[y, , ] <- (FM_tot[y, , ] / (FM_tot[y, , ] + M)) * nage[y, , ] * (1 - exp(-(M + FM_tot[y, , ])))
      catch_numage_hat[y, , ] <- (FM_tot[y, , ] / (FM_tot[y, , ] + M)) * nage_hat[y, , ] * (1 - exp(-(M + FM_tot[y, , ])))
      catch_numage_st[y, , ] <- (FM_tot[y, , ] / (FM_tot[y, , ] + M)) * nage_st[y, , ] * (1 - exp(-(M + FM_tot[y, , ])))
      catch_num[y, ] <- colSums(catch_numage[y, , ] + catch_numage_hat[y, , ] + catch_numage_st[y, , ])
      harvest_num[y, ] <- colSums(catch_numage[y, , ] * Vulret) + colSums(catch_numage_hat[y, , ] * Vulret) + colSums(catch_numage_st[y, , ] * Vulret)
      cpue[y, ] <- catch_num[y, ] / et[y, ]
      cpue[is.na(cpue)] <- 0
      hpue[y, ] <- harvest_num[y, ] / et[y, ]
      hpue[is.na(hpue)] <- 0
      avg_size[y, ] <- (colSums(nage[y, , ] * TL_global * Vulcat) + colSums(nage_hat[y, , ] * TL_global * Vulcat) + colSums(nage_st[y, , ] * TL_global * Vulcat)) /
        (colSums(nage[y, , ] * Vulcat) + colSums(nage_hat[y, , ] * Vulcat) + colSums(nage_st[y, , ] * Vulcat))
      if (DD_flag) {
        avg_size[y, ] <- (colSums(nage[y, , ] * TL_bar[y, , ] * Vulcat) + colSums(nage_hat[y, , ] * TL_bar[y, , ] * Vulcat) + colSums(nage_st[y, , ] * TL_bar[y, , ] * Vulcat)) /
          (colSums(nage[y, , ] * Vulcat) + colSums(nage_hat[y, , ] * Vulcat) + colSums(nage_st[y, , ] * Vulcat))
      }
      avg_size[is.na(avg_size)] <- 0
      # realized utility
      U.cpue[y, ] <- var.cpue[1] / (1 + exp(var.cpue[2] * (var.cpue[3] - cpue[y, ])))
      U.hpue[y, ] <- var.harv[1] / (1 + exp(var.harv[2] * (var.harv[3] - hpue[y, ])))
      U.size[y, ] <- var.size[1] / (1 + exp(var.size[2] * (var.size[3] - avg_size[y, ])))
      U.dist[y, ] <- var.dist[1] / (1 + exp(var.dist[2] * (var.dist[3] - dist)))
      U.crd[y, ] <- var.crowd[1] / (1 + exp(var.crowd[2] * (var.crowd[3] - et[y, ])))
      U.tot_add[y, ] <- U.cpue[y, ] + U.hpue[y, ] + U.size[y, ] + U.dist[y, ] + U.crd[y, ] + U.st[y, ]
      U.tot[y, ] <- U.tot_add[y, ]
      U.tot_sum[y] <- sum(U.tot[y, ])
      # eggs
      eggs[y, ] <- colSums(Fec * nage[y, , ])
      eggs_hat[y, ] <- colSums(Fec * nage_hat[y, , ]) + colSums(Fec * nage_st[y, , ])
    } # end of time dynamics loop [y]



    ################################
    ## Out of loop calculations ####
    ################################
    val <- U.tot * et
    # caught trophy numbers
    CTB <- sapply(1:nsites, function(x) rowSums(catch_numage[, trophy_age:Amax, x]))
    # caught trophy numbers (alive)
    CTB_alive <- sapply(1:nsites, function(x) rowSums(nage[, trophy_age:Amax, x]))
    # trophy biomass from population
    # CTB <- sapply(1:nsites, function(x) rowSums(nage[,trophy_age:Amax,x]*Wt))
    WSB <- t(t(eggs) / E0)

    export_years <- (years_init + 1):years
    res <- list(
      et = et, hr = hr, WSB = WSB, st = st, R = R, R_hat = R_hat, R_st = R_st,
      eggs = eggs, catch = catch_num, cpue = cpue,
      harvest = harvest_num, hpue = hpue, avg_size = avg_size,
      U.tot = U.tot, U.dist = U.dist, U.crd = U.crd, U.cpue = U.cpue, U.hpue = U.hpue, U.size = U.size, U.st = U.st,
      val = val, Linf_mat = Linf_mat, CTB = CTB, CTB_alive = CTB_alive
    )
    for (i in 1:length(res)) res[[i]] <- t(sapply(export_years, function(x) res[[i]][x, ]))
    res$TL_bar <- TL_bar[export_years, , ]
    res$nage <- nage[export_years, , ]
    res$nage_hat <- nage_hat[export_years, , ]
    res$nage_st <- nage_st[export_years, , ]
    res$U.init <- U.init

    if (isopleth_plots) res <- c(mean(res$CTB), mean(res$U.tot), mean(res$cpue), mean(res$et), mean(res$WSB))

    return(res)
  }) # end of with function
} # end of function
