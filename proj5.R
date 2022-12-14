################################################################################
###################### 5: UK excess deaths 2020-22 #############################
################################################################################

################################ Contributor ################################## 
##---------------------------  Kartik (s2407270)------------------------------##

##------------------------------ Description ---------------------------------##
## In this practical, we would like to estimate the excess deaths from 2020-22
## with relative to data for some previous period.  Given these data we can 
## simply apply the dying and ageing process to update the population in each 
## age group and obtain the expected number of deaths each week. Using this 
## approach, we have created a model in which we are predicting the deaths 
## each week from the start of 2020 to 2022. By obtaining these predictions, we
## can get excess Deaths each week by subtracting them from the actual deaths
## happened in that time. We can then using jags, build a simple time series
## model and collected 10,0000 samples. The model is useful for highlighting the 
## parts of the series that are most unusual such as weeks around Christmas and 
## new year because many deaths are reported late. By calculating the 
## posterior expected values for 10,000 samples we can get pretty accurate
## predictions for deaths happening each week.



##-------------------------------- CODE --------------------------------------##

## Loading ggplot
library(ggplot2)    

## load rjags (and coda)
library(rjags)

estimate_death <- function (d, mf, mm, Nf, Nm){
  ## This function predicts the number of deaths that's are happening in a 
  ## specific week and updating the population according to deaths happened in
  ## that week. The function takes arguments as d(death rate), mf and mm 
  ## instantaneous per ca-pita death rate per year of females and males 
  ## respectively and Nf & Nm represents populations of females and males 
  ## respectively.
  
  Nm_star_init <- Nm[1]                   ## Initial value for the N0* (males) 
  
  Nf_star_init <- Nf[1]                   ## Initial value for the N0* (females)
  
  qm <- 1 - exp(-mm/52)                   ## The expected proportion of males
                                          ## dying in a week  
  
  qf <- 1 - exp(-mf/52)                   ## The expected proportion of females
                                          ## dying in a week 
  
  pred_D <- rep(0, length.out= length(d)) ## Creating a vector to store death
                                          ## predictions
  
  for (j in c(1: length(d))){
  ## This loop is iterating through forward length(d) weeks 
    
      ## For Male Population
      Dm <-  0.9885 * d[j] * qm  * Nm     ## Death predictions for males
      Nm_star <- Nm - Dm                  ## Subtracting the died males from
                                          ## the population
      
      ## Creating a new vector with N_0 as the initial value
      Nm_star0 <- c(Nm_star_init, Nm_star)
      Nm_star0 <- Nm_star0[1:length(Nm_star0)-1]
      
      ## Calculating aging male population
      Nm_plus <- Nm_star * 51/52 + Nm_star0 / 52
      Nm <- Nm_plus                      ## Updating the male population
     
      
      ## For Female Population
      Df <-  0.9885 * d[j] * qf  * Nf    ## Death predictions for females
      Nf_star <- Nf - Df                 ## Subtracting the died females from
                                         ## the population
      
      ## Creating a new vector with N_0 as the initial value
      Nf_star0 <- c(Nf_star_init, Nf_star)
      Nf_star0 <- Nf_star0[1:length(Nf_star0)-1]
      
      ## Calculating the aging female population
      Nf_plus <- Nf_star * 51/52 + Nf_star0 / 52 
      
      
      pred_D[j] <- sum(Dm + Df)          ## getting sum of predicted death for 
                                         ## male & female in that week
      
      
      Nf <- Nf_plus                      ## Updating the female population
    
    

    
  }
    return(pred_D)
  }
  

## Reading Data from the files
Death <- read.table("death1722uk.dat", header=TRUE)
it <- read.table("lt1720uk.dat")

## extracting the mortality rate
d <- Death[, "d"]    

##getting the instantaneous per capita death rate(females and males respectively)
mf <- it[ ,"mf"]                            
mm <- it[,"mm"]

## Extracting the female and male population respectively
Nf <- it[,"fpop20"]
Nm <- it[,"mpop20"]

## Calling the function to get predictions for deaths
pred_deaths <- estimate_death(d[157:305], mf, mm, Nf, Nm) 
                                            

## Calculating the total excess Deaths 
excess_D <-  sum(Death[157:nrow(Death), 1] - pred_deaths)

## Computing Excess Death for 2020
excess_2020D <- sum(Death[157:208, 1] - pred_deaths[1:52])
                         


## Calculating the excess deaths from year 2020 to end of data
x <-  Death[157:nrow(Death), 1] - pred_deaths
    

## Converting them to Data Frame to make plotting easier
week <- c(157:nrow(Death))                      
prediction <- data.frame(week, pred_deaths) ## adding the weeks and predictions
                                            ## into DataFrame




## Creating a plot using ggplot library using geom_point to make scatter plot
## and geom_line to make continuous line(of predictions)  
plot1 <- ggplot(Death[157:305,], aes(x=week , y=deaths)) + 
  geom_point() +geom_line(color='red', data = prediction, 
  aes(x= week, y = pred_deaths)) + ylim(0,max(Death))

plot1 <- plot1 + labs(title = "Actual & Predicted Deaths 2020-2022",
         subtitle = paste("Excess Deaths in 2020 are:", floor(excess_2020D), 
         "& Overall Deaths are:", floor(excess_D)))  
                                          
print(plot1)                                ## Explicitly Printing the plot

CummulativeDeaths <- cumsum(x)              ## Compute the cummulative deaths
                                            ## with cumsum function

cummul_D <- data.frame(CummulativeDeaths)   ## Converting data to data Frame

## Plotting using ggplot with the converted Data frame                       
plot2 <- ggplot(cummul_D, aes(x= week, y = CummulativeDeaths))+geom_line()
                                            
plot2 <- plot2 + labs(title="Cummulative Deaths 2020-2022")                    
print(plot2)                                ## Explicitly printing the plot                          

xs <- x                                     ## x is vector of excess Deaths
idx <- c(51, 52, 53, 105, 106)              ## setting Christmas and new year 
xs[idx] <- NA                               ## predicted deaths to NA


## Calling the jags model file to call the model
mod <- jags.model("model.jags", data = list(x=xs, N=length(xs))) 

## Extracting the samples using coda available in rjags
sam.coda <- coda.samples(mod, variable.names=c("mu", "rho", "k"), n.iter = 10000)



## Plotting Trace Plot for rho
traceplot(sam.coda[[1]][,"rho"], main = "Trace of rho" )

## Plotting histogram for rho
hist(sam.coda[[1]][,"rho"] , main = "Histogram for rho", xlab = "rho")

sam_mat <- as.matrix(sam.coda)               ## Creating a matrix of samples
                                             ## extracted from coda i.e. mu's

col_idx <- c(-c(1), -c(151), -c(152))        ## getting the column indices to
                                             ## remove from the matrix 

sam_mat <- sam_mat[ ,col_idx]                ## Removing the indices 

E_mu <- colMeans(sam_mat)                    ## Calculating the posterior 
                                             ## expected value which is simply
                                             ## means of columns
            

## Getting every 50th sampled mu vector using the seq function
sampled_mu <- sam_mat[seq(1, nrow(sam_mat), 50),] 

## created a new vector to make the all indices zero other than the special ones
x_new <- x                                
x_new[-idx] <- NA                              

## (5th plot) Using matplot to plot the matrix  
matplot(week, t(sampled_mu), ylab="Sampled(mu)",type = "l", col="grey")
     
lines(week, E_mu, col="blue")                 ## Drawing the posterior expected
                                              ## values against week in same

points(week, x, pch=20)                       ## x(Excess Deaths) values vs week

points(week, x_new, col='red', pch=8)         ## modified x_new vs week as point


legend(x="topright",legend=c("50th Sample of E(mu)", "Expectation for mu",
  "Excess Deaths", "Excess Deaths(weeks not used) ") , 
   col=c("grey","blue","black","red"), lty=1:2, pch = c(20,8) ,cex=0.8,bty='n')          


## Calculating the residuals
residual <- x - E_mu                          

residual <- data.frame(residual)              ## Converting them to DataFrame

## (6th plot) plotting residuals vs weeks (using ggplot)
plot6 <- ggplot(residual, aes(x = week , y = residual)) + geom_point()
plot6 <- plot6 + labs(title="Residual vs weeks")
print(plot6)
