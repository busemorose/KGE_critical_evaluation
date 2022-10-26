## ---------------------------
##
## Author: Guillaume Cinkus
##
## Date Created: 2021-09-09
##
## Email: guillaume.cinkus@umontpellier.fr
##
## ---------------------------
##
## Notes:
##   
## This routine represents HBV Snow Routine and is written by Pianosi et al. (2015) 
## and modified by Chen et al. (2018).
## The temperature index approach by Hock (1999) is applied for snow melt calculation
## with consideration of potential clear-sky solar radiation.
## The routine is suitable for the calculation at daily and hourly time step.
## For the daily time step calculation, no data of potential clear-sky solar radiation
## need to be given.
## 
## Input:
## temp   = time series of temperature                            - vector (T,1)
## prec   = time series of precipitation                          - vector (T,1)
## srad   = time series of potential clear-sky solar radiation    - vector (T,1)
## param  = model parameters                                      - vector (T,5)
##            1. Ts    = threshold temperature [C]
##            2. MF    = melt factor [mm/C]
##            3. CFR   = refreezing factor [-]
##            4. CWH   = Water holding capacity of snow [-]   
##            5. alfa  = radiation coefficient [-]
##
## Output:
##      P = time series of simulated water leaving the            - vector (T,1)
##          routine [mm/timestep]        
## STATES = time series of simulated storages (all in mm)         - matrix (T,2)
##          1st column: water content of snowpack (snow component)
##          2nd column: water content of snowpack (liquid component)
##
## ---------------------------

snow_routine <- function(temp, prec, srad, param) {

#---------------------
# Preparing parameters
#---------------------

Ts <- param[1] # threshold temperature (?C)
MF <- param[2] # melt factor (mm/C)
CFR <- param[3] # refreezing factor (-)
CWH <- param[4] # water holding capacity of snow (-)
alfa <- param[5] # radiation coefficient (-)

N <- 1:length(prec) # length of the time series

#----------------
# Running routine
#----------------

P <- rep(0, length(N)) # water leaving the routine/recharge to the soil (mm/timestep)
rain <- ifelse(temp < Ts, 0, prec) # (mm/timestep)
snow <- ifelse(temp >= Ts, 0, prec) # (mm/timestep)
Ta <- ifelse(temp < Ts, 0, temp - Ts) # active temperature for snowmelt
Tn <- ifelse(temp >= Ts, 0, Ts - temp) # active temperature for refreezing
m <- rep(0, length(N)) # snowmelt (mm/timestep)
rfz <- rep(0, length(N)) # refreezing (mm/timestep)
v <- rep(0, length(N) + 1) # snowpack depth (mm): solid component
vl <- rep(0, length(N) + 1) # snowpack depth (mm): liquid component

for (t in rep(N)) {
  
  m[t] <- min((MF * Ta[t] + alfa * srad[t] * Ta[t]), v[t])
  rfz[t] <- min(CFR * MF * Tn[t], vl[t])
  
  # snowpack dynamics: solid component
  
  v[t + 1] <- v[t] - m[t] + snow[t] + rfz[t]
  
  # snowpack dynamics: liquid component
  
  vl[t + 1] <- vl[t] + m[t] + rain[t] - rfz[t]
  
  if (vl[t + 1] > CWH * v[t + 1]) { # if the liquid component exceeds the snowpack
  
   # holding capacity
    
    P[t] <- vl[t + 1] - CWH * v[t + 1]
    vl[t + 1] <- CWH * v[t + 1]
    
  } else {
    
    P[t] = 0
    
  }
}

output_snow_routine <- list("P" = P, "v" = v, "vl" = vl, "snow" = snow)
return(output_snow_routine)

}